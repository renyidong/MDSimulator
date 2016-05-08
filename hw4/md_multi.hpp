#include <cstdlib>
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <pthread.h>
#include <array>


#define DIMENSION 2
#define MAX_TYPES 5

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::array;
typedef array<double, DIMENSION> coordinate_vector;

typedef unsigned int atom_type_t;


namespace MD {

double  Delta,                          // Delta to calculate derivative of potential, force, suggest: max_distance/Delta ~ 10^14
        Cooling,
        TimeStep;                       // TimeStep of each calculation
unsigned int  NThread,                  // Number of threads to spawn, sugguest n_CPU + 1~2
              TotalSteps,               // Number of TimeSteps to be simulated
              OutputInterval;           // Print output each nth step

enum class boundary_condition_t {Periodic, Reflective};
array<boundary_condition_t, DIMENSION> BoundaryConditions;
array<double, DIMENSION> BoxSize;       // Half size of simulated box, simulation space is [-BoxSize,BoxSize]

array<array<double, MAX_TYPES>, MAX_TYPES> LJEplison;
array<array<double, MAX_TYPES>, MAX_TYPES> LJSigma;
array<double, MAX_TYPES> Mass;
array<bool, MAX_TYPES> Fixed;

double potentials(const vector<coordinate_vector> &pos_mat, const vector<atom_type_t> &type_vec) {
  double pot = 0;
  for (size_t i=0; i < pos_mat.size(); ++i) {
    for (size_t j=i+1; j < pos_mat.size(); ++j) {
      double dist2 = 0.0;
      for (unsigned int k=0; k<DIMENSION; ++k) {
        double d = pos_mat[i][k]-pos_mat[j][k];
        if (BoundaryConditions[k] == boundary_condition_t::Periodic ) {
          if      (d >  BoxSize[k]) d -= 2*BoxSize[k];
          else if (d < -BoxSize[k]) d += 2*BoxSize[k];
        }
        dist2 += d * d;
      }
      atom_type_t i_type = type_vec[i],
                  j_type = type_vec[j];
      double ljs = LJSigma[i_type][j_type];
      pot += 4*LJEplison[i_type][j_type] * ( pow(ljs*ljs/dist2, 6) - pow(ljs*ljs/dist2, 3) );
    }
  }
  return pot;
}

struct accel_thread_arg_t {
  const vector<coordinate_vector> *pos_mat_p;
  const vector<atom_type_t> *type_vec_p;
  double orig_pot;
  size_t range_begin, range_end;
  vector<coordinate_vector> *result_vec_p;
};

void* accel_thread(void *args_p) {
  struct accel_thread_arg_t * args = (struct accel_thread_arg_t *)args_p;
  vector<coordinate_vector> pos_mat(*(args->pos_mat_p));
  double orig_pot(args->orig_pot);
  vector<coordinate_vector> &result_vec = *(args->result_vec_p);
  const vector<atom_type_t> &type_vec = *(args->type_vec_p);

  for (size_t i=args->range_begin; i<args->range_end; ++i) {
    for (unsigned int j=0; j<DIMENSION; ++j) {
      double orig_pos = pos_mat[i][j];
      pos_mat[i][j] += Delta;
      double new_pot = potentials(pos_mat, type_vec);
      pos_mat[i][j]  = orig_pos;
      result_vec[i][j] = (orig_pot-new_pot)/Delta/Mass[type_vec[i]];
    }
  }
  return NULL;
}

vector<coordinate_vector> acceleration(const vector<coordinate_vector> &pos_mat,
                                       const vector<atom_type_t> &type_vec) {
  double orig_pot = potentials(pos_mat, type_vec);

  // Allocate result space
  vector<coordinate_vector> result_vec(pos_mat.size());
  
  // Prepare arguments for each thread and spawn accel threads
  pthread_t children[NThread];
  vector<accel_thread_arg_t> accel_thread_args(NThread);
  for (unsigned int i=0; i<NThread; ++i) {
    accel_thread_args[i].pos_mat_p  = &pos_mat;
    accel_thread_args[i].orig_pot = orig_pot;
    accel_thread_args[i].range_begin = pos_mat.size()*i/NThread;
    if (i != NThread-1) accel_thread_args[i].range_end = pos_mat.size()*(i+1)/NThread;
    else                accel_thread_args[i].range_end = pos_mat.size();
    accel_thread_args[i].result_vec_p = &result_vec;
    accel_thread_args[i].type_vec_p   = &type_vec;
    
    pthread_create(&children[i], NULL, *accel_thread, &accel_thread_args[i]);
  }
  
  // Join threads and calculate acceleration
  for (unsigned int i=0; i<NThread; ++i) {
    pthread_join(children[i], NULL);
  }
  return result_vec;
}

void velocity_verlet(vector<coordinate_vector> &pos_mat,
                     vector<coordinate_vector> &vel_mat,
                     const vector<atom_type_t> &type_vec,
                     const vector<coordinate_vector> &raw_pos_mat) {
  //Step 1: pos_mat += vel_mat*TimeStep + accel_mat*TimeStep*TimeStep/2
  vector<coordinate_vector> accel_mat = acceleration(pos_mat, type_vec);
  for (size_t i=0; i<pos_mat.size(); ++i) {
    for (unsigned int j=0; j<DIMENSION; ++j) {
//       if (Fixed[type_vec[i]]) {
//         double d_pos = pos_mat[i][j]-raw_pos_mat[i][j];
//         if (d_pos >  BoxSize[j]) d_pos -= 2*BoxSize[j];
//         if (d_pos < -BoxSize[j]) d_pos += 2*BoxSize[j];
//         accel_mat[i][j] -= Restore * d_pos;
//       }
      pos_mat[i][j] += vel_mat[i][j]*TimeStep + accel_mat[i][j]*TimeStep*TimeStep/2;
      if (BoundaryConditions[j] == boundary_condition_t::Periodic ) {
        if      (pos_mat[i][j] >  BoxSize[j] ) pos_mat[i][j] -= 2*BoxSize[j];
        else if (pos_mat[i][j] < -BoxSize[j] ) pos_mat[i][j] += 2*BoxSize[j];
      }
    }
  }
  
  //Step 3: new_vel_mat = vel_mat + (accel_mat+new_accel_mat)*TimeStep/2
  vector<coordinate_vector> new_accel_mat = acceleration(pos_mat, type_vec);
  for (size_t i=0; i<vel_mat.size(); ++i) {
    for (unsigned int j=0; j<DIMENSION; ++j) {
      if (BoundaryConditions[j] == boundary_condition_t::Reflective ) {
        if (abs(pos_mat[i][j]) >  BoxSize[j] ) vel_mat[i][j] = -vel_mat[i][j];
      }
      vel_mat[i][j] += (accel_mat[i][j]+new_accel_mat[i][j])*TimeStep/2;
      if (Fixed[type_vec[i]]) vel_mat[i][j] *= Cooling;
    } 
  }
}

double energy(const vector<coordinate_vector> &pos_mat,
              const vector<coordinate_vector> &vel_mat,
              const vector<atom_type_t>       &type_vec) {
  double ener = 0;
  for (size_t i=0; i<vel_mat.size(); ++i) {
    double v2 = 0;
    for (unsigned int j=0; j<DIMENSION; ++j) {
      v2 += vel_mat[i][j]*vel_mat[i][j];
    }
    ener += Mass[type_vec[i]]*v2;
  }
  ener /= 2;
  ener += potentials(pos_mat, type_vec);
  return ener;
}

void print_pos(const vector<coordinate_vector> &pos_mat,
               const vector<atom_type_t> &type_vec,
               double energy,
               size_t step) {
  cout <<pos_mat.size() <<'\n'
       <<"Time " <<step*TimeStep <<"\tEnerygy " <<energy <<'\n';
  cerr <<step <<endl;
  for (size_t i=0; i<pos_mat.size(); ++i) {
    cout <<(char)('A'+type_vec[i]);
    for (unsigned int j=0; j<DIMENSION; ++j) {
      cout <<'\t' <<pos_mat[i][j];
    }
    cout <<endl;
  }
}

void simulate(vector<coordinate_vector> &pos_mat,
              vector<coordinate_vector> &vel_mat,
              const vector<atom_type_t> &type_vec) {
  print_pos(pos_mat, type_vec, energy(pos_mat, vel_mat, type_vec), 0);
  vector<coordinate_vector> raw_pos_mat(pos_mat);
  for (size_t i=1; i<=TotalSteps; ++i) {
    velocity_verlet(pos_mat, vel_mat, type_vec, raw_pos_mat);
    if (i%OutputInterval==0) {
      print_pos(pos_mat, type_vec, energy(pos_mat, vel_mat, type_vec), i);
    }
  }      
}

}
