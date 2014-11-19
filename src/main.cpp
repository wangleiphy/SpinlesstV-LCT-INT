#include "interaction_expansion.hpp"

#include <alps/ngs.hpp>
#include <alps/parseargs.hpp>
#include <alps/stop_callback.hpp> 
//#include <alps/check_schedule.hpp> 
#include "mycheck_schedule.hpp"
#include <boost/chrono.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <functional>  

//#include <alps/hdf5.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>

class MpiSimulation : public InteractionExpansion
{

    public:
        typedef InteractionExpansion Base;
        typedef check_schedule ScheduleChecker;
        
        MpiSimulation(
              parameters_type& parm
            , boost::mpi::communicator const & comm
            , ScheduleChecker const & check = ScheduleChecker()
        )
            : Base(parm, comm.rank())
            , communicator(comm)
            , schedule_checker(check)
            , binnumber(parm["BINNUMBER"] | 128)
        {}

        bool run(boost::function<bool ()> const & stop_callback) 
        {

            bool done = false, stopped = false;
            while( !done ) 
            {

                update();
                measure();
                
                if(stopped || schedule_checker.pending())
                {
                    stopped = stop_callback(); 
                    double fraction = stopped ? 1. : Base::fraction_completed()*communicator.size();

                    std::pair<double, double> timing = schedule_checker.update(fraction);
                    done = (fraction >= 1.);

                    if( communicator.rank() == 0 )
                        std::cout << "Completed " << 100*fraction << "%. Next check in " << timing.first << "s. " 
                                  << "Finish in " << timing.second << "s. " 
                                  << "Time used " << schedule_checker.timespend() << "s. " 
                                  << "Sweeps " << progress()  << " "
                                  << "PertOrder " << pertorder()  << " "
                                  << "Block " << block()  << " "
                                  << "Cycles " << cycle()  << std::endl; 

                }
            }
            std::cout << "Rank " << communicator.rank() << " stopping after doing " << 100*Base::fraction_completed() << "% of the work in " << schedule_checker.timespend() << "s." << std::endl;

            measurements["Walltime"] << schedule_checker.timespend(); //wall time in unit of seconds 

            return !stopped;
        }

        results_type collect_results() const 
        {
            results_type results = Base::collect_results();
            for(results_type::iterator it = results.begin(); it != results.end(); ++it)
                it->second = it->second.reduce(communicator, binnumber);
            return results;
        }

    private:
        boost::mpi::communicator communicator;
        ScheduleChecker schedule_checker;
        unsigned binnumber;
};


int main(int argc, char** argv){

    try {


        // Init MPI environment
        boost::mpi::environment env(argc, argv);
        boost::mpi::communicator comm;

        // Parse command line options
        alps::parseargs options(argc, argv);


        //params file from input 
        alps::params params; 
        if (comm.rank()== 0){
            params = alps::params(options.input_file); 
            //std::cout << "lattice name = " << params["LATTICE"] << std::endl;
            //std::cout << "L,W= " << params["L"] << ","<< params["W"] << std::endl;
        }
        broadcast(comm, params);

        std::string filename = boost::lexical_cast<std::string>(params["filename"]);  
        std::string h5output_file = filename.substr(0, filename.find_last_of('.')) + ".out.h5"; // hdf5 output file 
        std::string checkpoint_path = filename.substr(0, filename.find_last_of('.'))  + ".chkp"; 
        std::string checkpoint_file = checkpoint_path +  "/clone"+ boost::lexical_cast<std::string>(comm.rank()) + ".h5";

        if (!boost::filesystem::exists(checkpoint_path) && comm.rank() ==0)
            boost::filesystem::create_directory(checkpoint_path);

      // Run simulation
      MpiSimulation sim(params, comm, check_schedule(options.tmin, options.tmax));
    
      if (options.resume && boost::filesystem::exists(checkpoint_file)){
        sim.load(checkpoint_file);

        if (comm.rank()== 0)
            std::cout << "loaded configuration from " << checkpoint_path  << std::endl;
      }else{
        sim.initialize_tvlist(); 

        if (comm.rank()== 0)
            std::cout << "will write configuration to " << checkpoint_path  << std::endl;
      }

      sim.run(alps::stop_callback(options.timelimit));  

      sim.save(checkpoint_file);

      //Collect results  
      MpiSimulation::results_type results = alps::collect_results(sim);
 

     
      if (comm.rank() ==0) 
      {
               sim.evaluate(results);

               alps::hdf5::archive ar(h5output_file, "w");
               ar["/parameters"] << params;
               ar["/simulation/results"] << results;
               ar.close(); 

               std::cout << results << std::endl; 

               /*
               std::ofstream outfile;
               outfile.open(filename.c_str());
               outfile.precision(12); 

               MpiSimulation::result_names_type obslist = sim.result_names(); //get a list of observables names
               //header in the .dat file  
               outfile << "#V, "; 
               for (unsigned i =0; i<obslist.size(); ++i){
                    std::string obsname = obslist[i]; 
                    if (results[obsname].is_type<double>()){
                        std::string obsname = obslist[i]; 
                        outfile << obsname << ", "; 
                    }
               }

               //the actual content    
               outfile <<  boost::lexical_cast<std::string>(params["V"]); 
               for (unsigned i =0; i<obslist.size(); ++i){
                    std::string obsname = obslist[i]; 
                    if (results[obsname].is_type<double>()){
                         double obs_val = results[obsname].mean<double>();
                         double obs_err = results[obsname].error<double>();
                         std::cout << obsname << " " << obs_val
                                              << " +- "<< obs_err << std::endl; 

                         outfile << " " << obs_val
                                 << " " << obs_err; 
                     
                    }
               }
               outfile << std::endl; 
               outfile.close();
               */
      }

    } catch (std::exception const & e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
