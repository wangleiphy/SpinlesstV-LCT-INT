#include "interaction_expansion.hpp"


/*implement add vertex*/
double InteractionExpansion::add_impl(const double tau, const std::vector<site_t>& sites, const bool compute_only_weight)
{

  assert(M.matrix().rows() == M.matrix().cols()); 
  assert(M.matrix().rows() == M.creators().size()); 

  unsigned int Msize = M.matrix().rows();
  unsigned int n = taus.size()/2; //number of vertices to add  

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

  //std::cout << "Stilde:\n"<< Stilde << std::endl; 
     
  for(unsigned int i=0; i<2*n; ++i){
   for(unsigned int j=0; j< Msize; ++j){
        Q(j,i) = green0_spline(M.creators()[j].t()-taus[i], M.creators()[j].s(), sites[i]);
        R(i,j) = -lattice.parity(sites[i]) * M.creators()[j].parity()* Q(j,i);//anti-symmetrization 
    }
  }

  if(Msize>0){
    MQ.noalias() = M.matrix() * Q; 
    Stilde.noalias() -= R * MQ; 
  }

  //return weight if we have nothing else to do
  if(compute_only_weight){
    return Stilde.determinant();// we have not yet perform the inverse, so it is actually 1./det(Stilde)
  }
  
  Stilde = Stilde.inverse().eval(); 

  if (Msize>0)
      RM.noalias() = R*M.matrix(); 

  M.matrix().conservativeResize(Msize+2*n, Msize+2*n);//conservativeResize keep the content 

  //perform the actual update  
  if(Msize>0){
    M.matrix().topRightCorner(Msize, 2*n).noalias()   = -MQ * Stilde ;       
    M.matrix().topLeftCorner(Msize, Msize).noalias() -=  M.matrix().topRightCorner(Msize, 2*n) * RM; 
    M.matrix().bottomLeftCorner(2*n, Msize).noalias() = -Stilde * RM; 
  }

  M.matrix().bottomRightCorner(2*n, 2*n) = Stilde; 

  for(unsigned int i=0; i<2*n; ++i)
     M.creators().push_back(creator(sites[i], lattice.parity(sites[i]), taus[i])); 

  M.num_vertices() += n; 

  return 1./Stilde.determinant();
}


/*implement remove vertex in partition funciton sector*/
double InteractionExpansion::remove_impl(const std::vector<unsigned int>& vertices, const bool compute_only_weight)
{// this will only get call when there is a vertex 
 // vertices contains from small to large indices 

  assert(M.matrix().rows() == M.matrix().cols()); 
  assert(M.matrix().rows() == M.creators().size()); 
  if (vertices.size()==2)
    assert(vertices[0]< vertices[1]); 
 
  unsigned int Msize = M.matrix().rows();
  unsigned int n = vertices.size(); //number of vertices need to remove  

  //std::cout << "Msize,n: " <<Msize << " "<<  n << std::endl; 
  //std::cout <<  M.matrix().block<2,2>(0, 0) << std::endl;
  //std::ostream_iterator<unsigned int> out_it (std::cout," ");
  //std::copy(vertices.begin(), vertices.end(), out_it );

  Eigen::MatrixXd Stilde(2*n, 2*n); // the block we want to remove 
  for (unsigned int i=0; i< n; ++i){
    for (unsigned int j=0; j< n; ++j){
        //std::cout << "i,j,M"<< i << " " << j << std::endl; 
        //std::cout << M.matrix().block<2,2>(2*vertices[i], 2*vertices[j]);

        Stilde.block<2,2>(2*i,2*j) = M.matrix().block<2,2>(2*vertices[i], 2*vertices[j]);
   }
  }

  if(compute_only_weight){
    return Stilde.determinant();
  }

  for (int i=n-1; i>=0; --i){//swap, from back to front 
      unsigned pos = 2*vertices[i]; 
      unsigned pos_j = Msize-2*(n-i); 

      assert(pos <= pos_j); 

      M.matrix().row(pos).swap(M.matrix().row(pos_j));
      M.matrix().row(pos+1).swap(M.matrix().row(pos_j+1));

      M.matrix().col(pos).swap(M.matrix().col(pos_j));
      M.matrix().col(pos+1).swap(M.matrix().col(pos_j+1));

      std::swap(M.creators()[pos], M.creators()[pos_j]);
      std::swap(M.creators()[pos+1], M.creators()[pos_j+1]);
  }

  //now perform fastupdate of M
  Msize -= 2*n; 
  Eigen::MatrixXd Qtilde(Msize, 2*n), Rtilde(2*n, Msize);

  if(Msize>0){
    Qtilde = M.matrix().topRightCorner(Msize, 2*n) ; //block(i,j,rows,cols)
    Rtilde = M.matrix().bottomLeftCorner(2*n, Msize) ; 
  }

  //std::cout << "remove:5" << std::endl; 

  M.matrix().conservativeResize(Msize, Msize); 

  if(Msize>0)
    M.matrix().noalias() -= Qtilde * Stilde.inverse() * Rtilde; 

  //std::cout << "remove:6" << std::endl; 
  //get rid of operators
  for (unsigned int i=0; i< n; ++i){
    M.creators().pop_back();
    M.creators().pop_back();
  }
  M.num_vertices() -= n; 

  return Stilde.determinant();
}
