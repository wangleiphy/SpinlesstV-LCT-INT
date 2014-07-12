#ifndef MYCHECK_SCHEDULE_HPP
#define MYCHECK_SCHEDULE_HPP

#include <boost/chrono.hpp>

    class check_schedule
    {
    public:
        typedef boost::chrono::high_resolution_clock clock;
        typedef boost::chrono::duration<double> duration;
        typedef boost::chrono::time_point<clock, duration> time_point;

        check_schedule(const double tmin=300., const double tmax=600.)
        :   min_check_(tmin)
        ,   max_check_(tmax)
        ,   check_duration_(min_check_)
        ,   time_start_(clock::now())
        ,   last_check_time_(clock::now())
        ,   old_fraction_(0.)
        {
        }

        bool pending() const
        {
            time_point now = clock::now();
            return now > (last_check_time_ + check_duration_);
        }

        double timespend() const 
        {
            duration timespend = clock::now() - time_start_; 
            return timespend.count(); 
        }
    

        std::pair<double, double> update(const double fraction)
        {
            time_point now = clock::now();
            duration remaining; 
            
            if( fraction > old_fraction_ )
            {
                // estimate remaining time; propose to run 1/4 of that time
                remaining = (1. - fraction) * (now-last_check_time_)/ (fraction - old_fraction_);
                old_fraction_ = fraction;  
                check_duration_ = 0.25 * remaining;  
                if( check_duration_ < min_check_ ) check_duration_ = min_check_;
                if( check_duration_ > max_check_ ) check_duration_ = max_check_;
            }
            else
                check_duration_ = min_check_;
            
            last_check_time_ = now;

            return std::make_pair(check_duration_.count(), remaining.count()); 
        }

    private:
        duration min_check_;
        duration max_check_;
        duration check_duration_;

        time_point time_start_; 
        time_point last_check_time_;

        double old_fraction_;
    
    };

#endif // !defined MYCHECK_SCHEDULE_HPP
