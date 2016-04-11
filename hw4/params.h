#define NDEBUG

#define NThread         6       // Number of threads to spawn, sugguest 1.2*n_CPU or n_CPU+2
#define LJEplison       1.0     // Eplison in LJ potential equation
#define LJSigma         1.0     // Sigma in LJ potential equation, equilibrium distance
#define Lambda          1.0     // atom velocities factor, sqrt( 3 (N-1) k_b T / sum(m v^2) )
#define Mass            1.0     // Mass of each atom
#define BoxSize         12.0    // Half size of simulated box, simulation space is [-BoxSize,BoxSize]
#define Delta           1e-4    // Delta to calculate derivative of potential, force
#define TimeStep        1e-2    // TimeStep of each calculation
#define TotalSteps      10000   // Number of TimeSteps to be simulated
#define OutputInterval  10      // Print output each nth step
