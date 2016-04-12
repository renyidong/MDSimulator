#define NThread         5       // Number of threads to spawn, sugguest n_CPU + 1~2
#define LJEplison       1.0e-1  // Eplison in LJ potential equation
#define LJSigma         1.0     // Sigma in LJ potential equation, equilibrium distance
#define Lambda          10      // atom velocities factor
#define Mass            1.0     // Mass of each atom
#define BoxSize         15.0    // Half size of simulated box, simulation space is [-BoxSize,BoxSize]
#define Delta           1.0e-4  // Delta to calculate derivative of potential, force
#define TimeStep        5e-4    // TimeStep of each calculation
#define TotalSteps      2500    // Number of TimeSteps to be simulated
#define OutputInterval  10      // Print output each nth step

#define NDEBUG
#include "md_multi.hpp"

int main() {
    vector<coordinate_vector> my_pos;
    vector<coordinate_vector> my_vel;
    vector<double> my_mass;
    std::default_random_engine rng;
    std::normal_distribution<double> gasdev;
    for (int i=-11;i<12;++i) {
        for (int j=-11;j<12;++j) {
            coordinate_vector p;
            p.push_back(i*LJSigma);
            p.push_back(j*LJSigma);
            my_pos.push_back(p);
            
            coordinate_vector v;
            v.push_back(gasdev(rng)*Lambda);
            v.push_back(gasdev(rng)*Lambda);
            my_vel.push_back(v);
            
            my_mass.push_back(Mass);
        }
    }
    
    print_pos(my_pos, energy(my_pos, my_vel, my_mass), 0);
    for (unsigned int i=1; i<=TotalSteps; ++i) {
        velocity_verlet(my_pos, my_vel, my_mass);
        if (i%OutputInterval==0) {
            print_pos(my_pos, energy(my_pos, my_vel, my_mass), i);
        }
    }
}
