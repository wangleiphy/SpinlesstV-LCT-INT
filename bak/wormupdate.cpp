#include "interaction_expansion.hpp"

double InteractionExpansion::wormcreate_impl(const double tau, const std::vector<site_t>& sites, const bool compute_only_weight)
{//create a worm at (time, sites)

  assert(M.matrix().rows() == M.matrix().cols()); 
  assert(M.matrix().rows() == M.creators().size()); 

  unsigned int Msize = M.matrix().rows();
  unsigned int Morder = sites.size(); 

  Eigen::MatrixXd Stilde = Eigen::MatrixXd::Zero(Morder, Morder); // diagonal term is always zero 
  Eigen::MatrixXd RM(Morder, Msize), R(Morder, Msize);
  Eigen::MatrixXd Q(Msize, Morder), MQ(Msize, Morder);

  double parity = 1.; 
  for (unsigned int i=0; i<Morder; ++i){
      parity *= lattice.parity(sites[i]); 
  }

  for(unsigned int i=0; i<Morder; ++i){
     for(unsigned int j=i+1; j<Morder; ++j){
         Stilde(i,j) = bare_green_itime(0, sites[i], sites[j]); 
         Stilde(j,i) = -Stilde(i,j)* lattice.parity(sites[i])* lattice.parity(sites[j]);   
     }
  }

  for(unsigned int i=0; i<Morder; ++i){
   for(unsigned int j=0; j< Msize; ++j){
        Q(j,i) = green0_spline(M.creators()[j].t()-tau, M.creators()[j].s(), sites[i]);
        R(i,j) = -lattice.parity(sites[i]) * M.creators()[j].parity()* Q(j,i);//anti-symmetrization 
    }
  }

  if(Msize>0){
    MQ.noalias() = M.matrix() * Q; 
    Stilde.noalias() -=  R * MQ; 
  }

  //return weight if we have nothing else to do
  if(compute_only_weight){
    return Stilde.determinant()*parity;// we have not yet perform the inverse, so it is actually 1./det(Stilde)
  }
  
  Stilde = Stilde.inverse().eval(); 

  if (Msize>0)
      RM.noalias() = R*M.matrix(); 

  M.matrix().conservativeResize(Msize+Morder, Msize+Morder);//conservativeResize keep the content 

  //perform the actual update  
  if(Msize>0){
    M.matrix().topRightCorner(Msize, Morder).noalias()  = - MQ * Stilde;
    M.matrix().topLeftCorner(Msize, Msize).noalias()    -=   M.matrix().topRightCorner(Msize, Morder) * RM; 
    M.matrix().bottomLeftCorner(Morder, Msize).noalias() = -Stilde * RM;
  }
  M.matrix().bottomRightCorner(Morder, Morder) = Stilde; 
 
  //we append the worm to the vertex list 
  for (unsigned int i=0; i< Morder; ++i)
    M.creators().push_back(creator(sites[i], lattice.parity(sites[i]), tau)); 
  
  wormlength += Morder;  
  return 1./Stilde.determinant()*parity;
}

double InteractionExpansion::wormdestroy_impl(const bool compute_only_weight, const unsigned int Morder)
{// this will only get call when there is a worm 

  assert(wormlength>0); 
  assert(M.matrix().rows() == M.matrix().cols()); 
  assert(M.matrix().rows() == M.creators().size()); 
  assert(M.matrix().rows() >= wormlength);

  unsigned int Msize = M.matrix().rows()-Morder;//always remove the last worm 

  double parity = 1.; 
  for (unsigned int i=0; i<Morder; ++i){
      parity *= M.creators()[Msize+i].parity(); 
  }

  Eigen::MatrixXd Stilde = M.matrix().bottomRightCorner(Morder, Morder); // the block we want to remove 

  if(compute_only_weight){
    return Stilde.determinant()*parity;
  }

  //now perform fastupdate of M
  Eigen::MatrixXd Qtilde(Msize, Morder);
  Eigen::MatrixXd Rtilde(Morder, Msize);

  if(Msize>0){
    Qtilde = M.matrix().topRightCorner(Msize, Morder) ; //block(i,j,rows,cols)
    Rtilde = M.matrix().bottomLeftCorner(Morder, Msize) ; 
  }

  M.matrix().conservativeResize(Msize, Msize); 

  if(Msize>0)
    M.matrix().noalias() -= Qtilde * Stilde.inverse() * Rtilde; 

  //get rid of the worm
  for (unsigned int i=0; i<Morder; ++i)
    M.creators().pop_back();

  wormlength -= Morder; 
  return Stilde.determinant()*parity;
}


//add a new vertext at (tau, sites) in the presence of worm 
double InteractionExpansion::wormadd_impl(const std::vector<double>& taus, const std::vector<site_t>& sites, const bool compute_only_weight)
{// this is *almost* identical to add in the absence of wroms
 // except after add vertex, we swap the new vertex with worm, thus keep the worm alway at end of the matrix 

  assert(M.matrix().rows() == M.matrix().cols()); 
  assert(M.matrix().rows() == M.creators().size()); 
  assert(taus.size() == sites.size());
  assert(wormlength>0);

  unsigned int Msize = M.matrix().rows();
  unsigned int n = taus.size()/2; //number of vertices to add  

  assert(n>=1);  
 
  //std::cout << "wormadd 1" << std::endl; 
  Eigen::MatrixXd Stilde = Eigen::MatrixXd::Zero(2*n, 2*n); // diagonal term is always zero 
  Eigen::MatrixXd RM(2*n, Msize), R(2*n, Msize);
  Eigen::MatrixXd Q(Msize, 2*n), MQ(Msize, 2*n);

   for(unsigned int i=0; i<2*n; ++i){
     //Stilde(i,i) = lattice.parity(sites[i])*delta; //staggered on site potential  
     for(unsigned int j=i+1; j<2*n; ++j){
         Stilde(i,j) = green0_spline(taus[i]-taus[j], sites[i], sites[j]);
         Stilde(j,i) = -Stilde(i,j)* lattice.parity(sites[i])* lattice.parity(sites[j]);   
     }   
  }

  //std::cout << "wormadd 2" << std::endl; 

  for(unsigned int i=0; i<2*n; ++i){
   for(unsigned int j=0; j< Msize; ++j){
        Q(j,i) = green0_spline(M.creators()[j].t()-taus[i], M.creators()[j].s(), sites[i]);
        R(i,j) = -lattice.parity(sites[i]) * M.creators()[j].parity()* Q(j,i);//anti-symmetrization 
    }
  }

  if(Msize>0){
    MQ.noalias() = M.matrix() * Q; 
    Stilde.noalias() -=  R * MQ; 
  }

  //std::cout << "wormadd 3"<< Msize << " " << n << std::endl; 

  //return weight if we have nothing else to do
  if(compute_only_weight){
    return Stilde.determinant();// we have not yet perform the inverse, so it is actually 1./det(Stilde)
  }
  
  Stilde = Stilde.inverse().eval(); 

  if (Msize>0)
      RM.noalias() = R*M.matrix(); 
 
  M.matrix().conservativeResize(Msize+2*n, Msize+2*n); //conservativeResize keep the content 

  //std::cout << "wormadd 4"<< Msize << " " << n << std::endl; 

  //perform the actual update  
  if(Msize>0){
    M.matrix().topRightCorner(Msize, 2*n).noalias()  = - MQ * Stilde ;      
    M.matrix().topLeftCorner(Msize, Msize).noalias() -= M.matrix().topRightCorner(Msize, 2*n) * RM;
    M.matrix().bottomLeftCorner(2*n, Msize).noalias() = -Stilde * RM;  
  }

  M.matrix().bottomRightCorner(2*n, 2*n) = Stilde; 
 
  //add new vertex to the end 
 for(unsigned int i=0; i<2*n; ++i)
    M.creators().push_back(creator(sites[i], lattice.parity(sites[i]), taus[i])); 

  for (int i=1; i<=wormlength/2; ++i){//swap with worm, from back to front 
      unsigned int pos = Msize-2*i;
      unsigned int pos_j = Msize+2*n -2*i;

      if (pos > pos_j){
        std::cout << "i,wormlength,Msize " << i << " " << wormlength << " " << Msize << std::endl; 
        std::cout << "pos, pos_j " << pos << " " << pos_j << std::endl; 
        abort();  
      }

      assert(pos <= pos_j); 

      M.matrix().row(pos).swap(M.matrix().row(pos_j));
      M.matrix().row(pos+1).swap(M.matrix().row(pos_j+1));

      M.matrix().col(pos).swap(M.matrix().col(pos_j));
      M.matrix().col(pos+1).swap(M.matrix().col(pos_j+1));
 
      std::swap(M.creators()[pos], M.creators()[pos_j]);
      std::swap(M.creators()[pos+1], M.creators()[pos_j+1]);
  }
  
  //std::cout << "wormadd 5"<< Msize << " " << n << std::endl; 

  M.num_vertices() += n;    
  return 1./Stilde.determinant();
}


double InteractionExpansion::wormremove_impl(const std::vector<unsigned int>& vertices, const bool compute_only_weight)
{// this *almost* identical to the remove in the absence of worm 
 // except to remove the vetex , we swap it with the worm, remove it from the end, then swap the worm back to the end 

  unsigned int Msize = M.matrix().rows();
  unsigned int n = vertices.size(); //number of vertices need to remove 
  
  Eigen::MatrixXd Stilde(2*n, 2*n); // the block we want to remove 
  //std::cout << "wormremove 1" << std::endl; 
  for (unsigned int i=0; i< n; ++i){
    for (unsigned int j=0; j< n; ++j){
        Stilde.block<2,2>(2*i,2*j) = M.matrix().block<2,2>(2*vertices[i], 2*vertices[j]);
   }   
  }

  if(compute_only_weight){
    return Stilde.determinant();
  }

  //std::cout << "wormremove 1" << std::endl; 

  //move vertices to end 
 std::vector<unsigned int> wormpos; //keep track of the worm position
 for (unsigned int i =0; i< wormlength/2; ++i)
    wormpos.push_back(Msize-wormlength+ 2*i);  

  //std::cout << "wormremove 2" << std::endl; 
  //std::cout << "vertices to remove" << std::endl; 
  //std::ostream_iterator<unsigned int> out_it (std::cout ,", ");
  //copy(vertices.begin(), vertices.end(), out_it );
  //std::cout << std::endl; 

  //std::cout << "wormpos" << std::endl; 
  //copy(wormpos.begin(), wormpos.end(), out_it );
  //std::cout << std::endl; 

  for (int i=1; i<=n; ++i){//swap, from back to front 
      unsigned int pos = 2*vertices[n-i]; 
      unsigned int pos_j = Msize-2*i; 

      assert(pos <= pos_j); 
      assert(std::find(wormpos.begin(), wormpos.end(), pos) == wormpos.end() ); 
     
      std::replace(wormpos.begin(), wormpos.end(), pos_j, pos); 

      M.matrix().row(pos).swap(M.matrix().row(pos_j));
      M.matrix().row(pos+1).swap(M.matrix().row(pos_j+1));

      M.matrix().col(pos).swap(M.matrix().col(pos_j));
      M.matrix().col(pos+1).swap(M.matrix().col(pos_j+1));

      std::swap(M.creators()[pos], M.creators()[pos_j]);
      std::swap(M.creators()[pos+1], M.creators()[pos_j+1]);
  }
  //std::cout << "wormremove 3" << std::endl; 
  //std::cout << "wormpos" << std::endl; 
  //copy(wormpos.begin(), wormpos.end(), out_it );
  //std::cout << std::endl; 

  //now perform fastupdate of M
  Msize -= 2*n; 
  Eigen::MatrixXd Qtilde(Msize, 2*n), Rtilde(2*n, Msize);

  //std::cout << "wormremove 4" << std::endl; 

  if(Msize>0){
    Qtilde = M.matrix().topRightCorner(Msize, 2*n) ; //block(i,j,rows,cols)
    Rtilde = M.matrix().bottomLeftCorner(2*n, Msize) ; 
  }

  M.matrix().conservativeResize(Msize, Msize); 

  if(Msize>0)
    M.matrix().noalias() -= Qtilde * Stilde.inverse() * Rtilde; 

  //std::cout << "wormremove 5" << std::endl; 

  //get rid of the operators
  for (unsigned int i=0; i< n; ++i){
    M.creators().pop_back();
    M.creators().pop_back();
  }
  

  //now swap the worms back to the end
  for (unsigned int i =0; i< wormlength/2; ++i){
    unsigned int pos = wormpos[i]; 
    unsigned int pos_j = Msize-wormlength+2*i; 
   
    //std::cout << "Msize, wormlength,i " <<  Msize <<  " " << wormlength << " " << i << std::endl; 
    //std::cout << "pos, pos_j " <<  pos <<  " " << pos_j << std::endl; 

    //if the last one is also a worm, we keep track its position 
    std::replace(wormpos.begin(), wormpos.end(), pos_j, pos); 

    M.matrix().row(pos).swap(M.matrix().row(pos_j));
    M.matrix().row(pos+1).swap(M.matrix().row(pos_j+1));
    
    M.matrix().col(pos).swap(M.matrix().col(pos_j));
    M.matrix().col(pos+1).swap(M.matrix().col(pos_j+1));
 
    std::swap(M.creators()[pos],   M.creators()[pos_j]);
    std::swap(M.creators()[pos+1], M.creators()[pos_j+1]);
  }
  

  //std::cout << "wormremove 6" << std::endl; 
  M.num_vertices() -= n; 
  return Stilde.determinant();
}


double InteractionExpansion::wormshift_impl(const double tau, const std::vector<site_t>& sites, const bool compute_only_weight)
{//shift a worm to a new (tau, sites)
 //the way to do it is we (temperatory) remove the worm, then add it back 

  assert(M.matrix().rows() == M.matrix().cols()); 
  assert(wormlength>0); 
  assert(M.matrix().rows()>=wormlength); 

  unsigned int Msize = M.matrix().rows()-wormlength;// the worm position 

  // the old worm 
  double oldparity = 1.; 
  for (unsigned int i=0; i<wormlength; ++i){
      oldparity *= M.creators()[Msize+i].parity() ;
  }

  Eigen::MatrixXd oldStilde = M.matrix().bottomRightCorner(wormlength, wormlength); //worm is always at the very end 
  Eigen::MatrixXd oldM(Msize, Msize), Qtilde(Msize, wormlength), Rtilde(Msize, wormlength);// oldM is a unfortunature name but it just means the M matrix if we remove worm 

  if(Msize>0){
    Qtilde = M.matrix().topRightCorner(Msize, wormlength) ; //block(i,j,rows,cols)
    Rtilde = M.matrix().bottomLeftCorner(wormlength, Msize) ; 
    oldM = M.matrix().topLeftCorner(Msize, Msize) ; //block(i,j,rows,cols)
    oldM.noalias() -= Qtilde * oldStilde.inverse() * Rtilde; 
  }

  // the new worm 
  Eigen::MatrixXd Stilde = Eigen::MatrixXd::Zero(wormlength, wormlength); // diagonal term is always zero 
  Eigen::MatrixXd RM(wormlength, Msize), R(wormlength, Msize);
  Eigen::MatrixXd Q(Msize, wormlength), MQ(Msize, wormlength);

  double parity = 1.; 
  for (unsigned int i=0; i<wormlength; ++i){
      parity *= lattice.parity(sites[i]); 
  }

  for(unsigned int i=0; i<wormlength; ++i){
     for(unsigned int j=i+1; j<wormlength; ++j){
         Stilde(i,j) = bare_green_itime(0, sites[i], sites[j]); 
         Stilde(j,i) = -Stilde(i,j)* lattice.parity(sites[i])* lattice.parity(sites[j]);   
     }
  }

  for(unsigned int i=0; i<wormlength; ++i){
   for(unsigned int j=0; j< Msize; ++j){
        Q(j,i) = green0_spline(M.creators()[j].t()-tau, M.creators()[j].s(), sites[i]);
        R(i,j) = -lattice.parity(sites[i]) * M.creators()[j].parity()* Q(j,i);//anti-symmetrization 
    }
  }

  if(Msize>0){
    MQ.noalias() = oldM * Q; 
    Stilde.noalias() -=  R * MQ; 
  }

  if(compute_only_weight){
    return Stilde.determinant() * oldStilde.determinant()*parity*oldparity;
  }

  Stilde = Stilde.inverse().eval(); 

  if (Msize>0)
      RM.noalias() = R*oldM; 

  //perform the actual update  
  if(Msize>0){
    M.matrix().topRightCorner(Msize, wormlength).noalias()  = - MQ * Stilde ;   
    M.matrix().topLeftCorner(Msize, Msize).noalias() = oldM - M.matrix().topRightCorner(Msize, wormlength) * RM; 
    M.matrix().bottomLeftCorner(wormlength, Msize).noalias()  = -Stilde * RM;
  }
  M.matrix().bottomRightCorner(wormlength, wormlength) = Stilde; 
  
  //get rid of oldworm and replace it by the new one  
  for(unsigned int i=0; i<wormlength; ++i){
    M.creators().pop_back();
  }

  for(unsigned int i=0; i<wormlength; ++i){
        M.creators().push_back(creator(sites[i], lattice.parity(sites[i]), tau)); 
  }
  
  return 1./Stilde.determinant() * oldStilde.determinant()*parity*oldparity;

}

