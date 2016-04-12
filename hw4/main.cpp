#define NDEBUG
#include "md_multi.hpp"

int main() {
  // Parameters. For usage see "md_multi.hpp"
  MD::NThread         =23;
  MD::LJEplison       =1.0e-1;
  MD::LJSigma         =1.0;
  MD::BoxSize         =15.0;
  MD::Delta           =1.0e-12;
  MD::TimeStep        =1e-5;
  MD::TotalSteps      =5000;
  MD::OutputInterval  =10;
  
  
  std::default_random_engine rng;
  std::normal_distribution<double> gasdev;
  vector<coordinate_vector> my_pos;
  vector<coordinate_vector> my_vel;
  vector<double> my_mass;
  
  double mass   = 1.0;
  double Lambda = 20;
  int n = 11;
  for (int i=-n;i<n+1;++i) {
    for (int j=-n;j<n+1;++j) {
      coordinate_vector p;
      p.push_back(i*MD::BoxSize/(n+0.5));
      p.push_back(j*MD::BoxSize/(n+0.5));
      my_pos.push_back(p);
      
      coordinate_vector v;
      v.push_back(gasdev(rng)*Lambda);
      v.push_back(gasdev(rng)*Lambda);
      my_vel.push_back(v);
      
      my_mass.push_back(mass);
    }
  }
  
  // Run simulation
  MD::simulate(my_pos, my_vel, my_mass);
}
