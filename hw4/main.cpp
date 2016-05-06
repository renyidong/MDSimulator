#define NDEBUG
#include "md_multi.hpp"

int main() {
  // Parameters. For usage see "md_multi.hpp"
  MD::NThread         =4;
  MD::LJEplison[0][0] =1.0e-1;
  MD::LJSigma[0][0]   =1.0;
  MD::BoxSize         =10.0;
  MD::Delta           =1.0e-12;
  MD::TimeStep        =3e-4;
  MD::TotalSteps      =1000;
  MD::OutputInterval  =10;
  MD::BoundaryConditions.fill(MD::boundary_condition_t::Periodic);
  MD::Mass[0]         =1.0;
  
  std::default_random_engine rng;
  std::normal_distribution<double> gasdev;
  vector<coordinate_vector> my_pos_mat;
  vector<coordinate_vector> my_vel_mat;
  vector<atom_type_t>       my_type_vec;
  
  double Lambda = 20;
  int n = 4;
  for (int i=-n;i<n;++i) {
    for (int j=-n;j<n;++j) {
      coordinate_vector p;
      p[0] = i*MD::BoxSize/(n+0.5);
      p[1] = j*MD::BoxSize/(n+0.5);
      my_pos_mat.push_back(p);
      
      coordinate_vector v;
      v[0]=gasdev(rng)*Lambda;
      v[1]=gasdev(rng)*Lambda;
      my_vel_mat.push_back(v);
      
      my_type_vec.push_back(0);
    }
  }
  
  // Run simulation
  MD::simulate(my_pos_mat, my_vel_mat, my_type_vec);
}
