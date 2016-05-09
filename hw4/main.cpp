#define NDEBUG
#include "md_multi.hpp"
#include "param.hpp"
#include <cmath>
#define LJMIN 1.122462


int main() {
  // Parameters. For usage see "md_multi.hpp"
  MD::NThread         =12;
  MD::BoxSize[0]      =10.0;
  MD::BoxSize[1]      =20.0;
  MD::Delta           =1.0e-12;
  MD::TimeStep        =1e-3;
  MD::TotalSteps      =500000;
  MD::OutputInterval  =100;
  MD::BoundaryConditions[0]=MD::boundary_condition_t::Periodic;
  MD::BoundaryConditions[1]=MD::boundary_condition_t::Reflective;
  MD::Cooling         = 1.0 - 2.0e-3;
  MD::Mass            = {1,1,1};
  MD::Fixed           = {1,0,0};
  for (int i=0; i<2; ++i) {
    for (int j=0; j<2; ++j) {
      MD::LJEplison[i][j] = EP1;
      MD::LJSigma[i][j]   = SI1;
    }
  }
  for (int i=0; i<2; ++i) {
    MD::LJEplison[i][2] = MD::LJEplison[2][i] = EP2;
    MD::LJSigma[i][2]   = MD::LJSigma[2][i]   = SI2;
  }
  MD::LJEplison[2][2] = EP3;
  MD::LJSigma[2][2]   = SI3;
  
  std::default_random_engine rng;
  std::normal_distribution<double> gasdev;
  vector<coordinate_vector> my_pos_mat;
  vector<coordinate_vector> my_vel_mat;
  vector<atom_type_t>       my_type_vec;
  
  int n[2] = {5,10};
  for (int i=-n[0];i<n[0];++i) {
    for (int j=0;j<n[1];++j) {
      coordinate_vector p;
      p[0] = i*MD::BoxSize[0]/(n[0]+0.5);
      p[1] = j*MD::BoxSize[1]/(n[1]+0.5);
      my_pos_mat.push_back(p);
      
      coordinate_vector v;
      v[0]=gasdev(rng)*LAMBDA;
      v[1]=gasdev(rng)*LAMBDA;
      my_vel_mat.push_back(v);
      
      my_type_vec.push_back(2);
    }
  }
  
  n[0] = MD::BoxSize[0]/LJMIN/MD::LJSigma[0][0];
  for (int i=-n[0];i<n[0];i++) {
          coordinate_vector p;
          p[0] = i*MD::LJSigma[0][0]*LJMIN;
          p[1] = -MD::BoxSize[1]+1*MD::LJSigma[0][0];
          my_pos_mat.push_back(p);
          p[1] = -MD::BoxSize[1]+3*MD::LJSigma[0][0];
          my_pos_mat.push_back(p);
          p[0] = (i+0.5)*MD::LJSigma[0][0]*LJMIN;
          p[1] = -MD::BoxSize[1]+2*MD::LJSigma[0][0];
          my_pos_mat.push_back(p);
          p[1] = -MD::BoxSize[1]+4*MD::LJSigma[0][0];
          my_pos_mat.push_back(p);
          
          coordinate_vector v;
          v[0]=0;
          v[1]=0;
          my_vel_mat.push_back(v);
          my_vel_mat.push_back(v);
          my_vel_mat.push_back(v);
          my_vel_mat.push_back(v);
          
          my_type_vec.push_back(0);
          my_type_vec.push_back(1);
          my_type_vec.push_back(0);
          my_type_vec.push_back(1);
  }

  
  // Run simulation
  MD::simulate(my_pos_mat, my_vel_mat, my_type_vec);
}
