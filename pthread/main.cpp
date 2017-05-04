#define NDEBUG
#include "md_multi.hpp"

int main() {
  // Parameters. For usage see "md_multi.hpp"
  MD::NThread         =2;
  MD::LJEplison       =1.0e-1;
  MD::LJSigma         =1.0;
  MD::BoxSize         =15.0;
  MD::Delta           =1.0e-12;
  MD::TimeStep        =1e-5;
  MD::TotalSteps      =100;
  MD::OutputInterval  =10;
  
  
  std::default_random_engine rng;
  std::normal_distribution<double> gasdev;
  vector<coordinate_vector> my_pos;
  vector<coordinate_vector> my_vel;
  vector<double> my_mass;
  
  double mass   = 1.0;
  double Lambda = 20;
  int n = 6;
  for (int i=-n;i<n;++i) {
    for (int j=-n;j<n;++j) {
      coordinate_vector p;
      p[0] = i*MD::BoxSize/(n+0.5);
      p[1] = j*MD::BoxSize/(n+0.5);
      my_pos.push_back(p);
      
      coordinate_vector v;
      v[0]=gasdev(rng)*Lambda;
      v[1]=gasdev(rng)*Lambda;
      my_vel.push_back(v);
      
      my_mass.push_back(mass);
    }
  }
  
  // Run simulation
  MD::simulate(my_pos, my_vel, my_mass);
}
