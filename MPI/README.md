Usage: mpirun -n NRank ./a.out Steps NRow NThread

Steps is number of timesteps to simulate.

NRow is number of atoms in a row. There will be NRow^3 atoms in total.
The running time is expected to be O(NRow^6)
Warning:
1. NRow must be divisible by 2
2. NRow^3 must be divisible by NRank * NThread. 

NThread is number of threads computing
