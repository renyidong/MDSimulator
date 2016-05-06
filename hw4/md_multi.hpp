#include <cstdlib>
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <pthread.h>
#include <array>


#define DIMENSION 2

using std::vector;
using std::cout;
using std::array;
typedef array<double, DIMENSION> coordinate_vector;


namespace MD{

double  LJEplison,                      // Eplison in LJ potential equation
        LJSigma,                        // Sigma in LJ potential equation, equilibrium distance
        Lambda,                         // atom velocities factor
        Mass,                           // Mass of each atom
        BoxSize,                        // Half size of simulated box, simulation space is [-BoxSize,BoxSize]
        Delta,                          // Delta to calculate derivative of potential, force, suggest: max_distance/Delta ~ 10^14
        TimeStep;                       // TimeStep of each calculation
unsigned int  NThread,                  // Number of threads to spawn, sugguest n_CPU + 1~2
              TotalSteps,               // Number of TimeSteps to be simulated
              OutputInterval;           // Print output each nth step


double potentials(const vector<coordinate_vector> &pos_mat) {
  double pot = 0;
  for (size_t i=0; i < pos_mat.size(); ++i) {
    for (size_t j=i+1; j < pos_mat.size(); ++j) {
      double dist2 = 0.0;
      for (unsigned int k=0; k<DIMENSION; ++k) {
        double d = pos_mat[i][k]-pos_mat[j][k];
        if      (d >  BoxSize) d -= 2*BoxSize;
        else if (d < -BoxSize) d += 2*BoxSize;
        dist2 += d * d;
      }
      pot += 4*LJEplison * ( pow(LJSigma*LJSigma/dist2, 6) - pow(LJSigma*LJSigma/dist2, 3) );
    }
  }
  return pot;
}

struct force_thread_arg_t {
  vector<coordinate_vector> pos_mat;
  double orig_pot;
  size_t range_begin, range_end;
  vector<coordinate_vector> *result_vec_p;
};

void* force_thread(void *args_p) {
  struct force_thread_arg_t * args = (struct force_thread_arg_t *)args_p;
  vector<coordinate_vector> &pos_mat = args->pos_mat;
  double orig_pot(args->orig_pot);
  vector<coordinate_vector> &result_vec = *(args->result_vec_p);

  for (size_t i=args->range_begin; i<args->range_end; ++i) {
    for (unsigned int j=0; j<DIMENSION; ++j) {
      double orig_pos = pos_mat[i][j];
      pos_mat[i][j] += Delta;
      double new_pot = potentials(pos_mat);
      pos_mat[i][j]  = orig_pos;
      result_vec[i][j] = (orig_pot-new_pot)/Delta;
    }
  }
  return NULL;
}

vector<coordinate_vector> acceleration(const vector<coordinate_vector> &pos_mat,
                                       const vector<double> &mass_mat) {
  double orig_pot = potentials(pos_mat);

  // Allocate result space
  vector<coordinate_vector> result_vec(pos_mat.size());
  
  // Prepare arguments for each thread and spawn force threads 
  pthread_t children[NThread];
  size_t BatchSize = pos_mat.size()/NThread;
  vector<force_thread_arg_t> force_thread_args(NThread);
  for (unsigned int i=0; i<NThread; ++i) {
    force_thread_args[i].pos_mat  = pos_mat;
    force_thread_args[i].orig_pot = orig_pot;
    force_thread_args[i].range_begin = i*BatchSize;
    if (i != NThread-1) force_thread_args[i].range_end = (i+1)*BatchSize;
    else                force_thread_args[i].range_end = pos_mat.size();
    force_thread_args[i].result_vec_p = &result_vec;
    pthread_create(&children[i], NULL, *force_thread, &force_thread_args[i]);
  }
  
  // Join threads and calculate acceleration
  for (unsigned int i=0; i<NThread; ++i) {
    pthread_join(children[i], NULL);
    for (size_t j=force_thread_args[i].range_begin; j<force_thread_args[i].range_end; ++j) {
      for (unsigned int k=0; k<DIMENSION; ++k) {
        result_vec[j][k] /= mass_mat[j];
      }
    }
  }
  return result_vec;
}

void velocity_verlet(vector<coordinate_vector> &pos_mat,
                     vector<coordinate_vector> &vel_mat,
                     const vector<double>      &mass_mat) {
  //Step 1: pos_mat += vel_mat*TimeStep + accel_mat*TimeStep*TimeStep/2
  vector<coordinate_vector> accel_mat = acceleration(pos_mat, mass_mat);
  for (size_t i=0; i<pos_mat.size(); ++i) {
    for (unsigned int j=0; j<DIMENSION; ++j) {
      pos_mat[i][j] += vel_mat[i][j]*TimeStep + accel_mat[i][j]*TimeStep*TimeStep/2;
      if      (pos_mat[i][j] >  BoxSize ) pos_mat[i][j] -= 2*BoxSize;
      else if (pos_mat[i][j] < -BoxSize ) pos_mat[i][j] += 2*BoxSize;
    }
  }

  //Step 2: new_half_vel = vel_mat + accel_mat*TimeStep/2
  vector<coordinate_vector> new_half_vel(vel_mat.size());
  for (size_t i=0; i<new_half_vel.size(); ++i) {
    for (unsigned int j=0; j<DIMENSION; ++j) {
      new_half_vel[i][j] = vel_mat[i][j] + accel_mat[i][j]*TimeStep/2;
    }
  }
  
  //Step 3: new_vel_mat = vel_mat + (accel_mat+new_accel_mat)*TimeStep/2
  vector<coordinate_vector> new_accel_mat = acceleration(pos_mat, mass_mat);
  for (size_t i=0; i<vel_mat.size(); ++i) {
    for (unsigned int j=0; j<DIMENSION; ++j) {
      vel_mat[i][j] += (accel_mat[i][j]+new_accel_mat[i][j])*TimeStep/2;
    } 
  }
}

double energy(const vector<coordinate_vector> &pos_mat,
              const vector<coordinate_vector> &vel_mat,
              const vector<double>            &mass_mat) {
  double ener = 0;
  for (size_t i=0; i<vel_mat.size(); ++i) {
    double v2 = 0;
    for (unsigned int j=0; j<DIMENSION; ++j) {
      v2 += vel_mat[i][j]*vel_mat[i][j];
    }
    ener += mass_mat[i]*v2;
  }
  ener /= 2;
  ener += potentials(pos_mat);
  return ener;
}

void print_pos(const vector<coordinate_vector> &pos_mat,
               double energy,
               size_t step) {
  cout <<pos_mat.size() <<'\n'
       <<"Time " <<step*TimeStep <<"\tEnerygy " <<energy <<'\n';
  for (size_t i=0; i<pos_mat.size(); ++i) {
    cout <<'A';
    for (unsigned int j=0; j<DIMENSION; ++j) {
      cout <<'\t' <<pos_mat[i][j];
    }
    cout <<'\n' <<std::flush;
  }
}

void simulate(vector<coordinate_vector> &pos_mat,
              vector<coordinate_vector> &vel_mat,
              vector<double>            &mass_mat) {
  print_pos(pos_mat, energy(pos_mat, vel_mat, mass_mat), 0);
  for (size_t i=1; i<=TotalSteps; ++i) {
    velocity_verlet(pos_mat, vel_mat, mass_mat);
    if (i%OutputInterval==0) {
      print_pos(pos_mat, energy(pos_mat, vel_mat, mass_mat), i);
    }
  }      
}

}
