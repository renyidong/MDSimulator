#include <cstdlib>
#include <cmath>
#include <vector>
#include <list>
#include <cassert>
#include <iostream>
#include <random>
#include <pthread.h>

using std::vector;
using std::list;
using std::cout;

#define LJEplison 1
#define LJSigma 1
#define Delta 0.0001
#define NThread 6


typedef vector<double> coordinate_vector;

double potentials(const vector<coordinate_vector> &pos_mat) {
  double pot = 0;
  for (unsigned int i=0; i < pos_mat.size(); ++i) {
    for (unsigned int j=i+1; j < pos_mat.size(); ++j) {
      assert(pos_mat[i].size()==pos_mat[j].size());
      double dist2 = 0.0;
      for (unsigned int k=0; k<pos_mat[i].size(); ++k) {
        dist2 += pow(pos_mat[i][k]-pos_mat[j][k], 2);
      }
      pot += 4*LJEplison * ( pow(LJSigma, 12)/pow(dist2, 6) - pow(LJSigma, 6)/pow(dist2, 3) );
    }
  }
  return pot;
}

struct force_thread_arg_t {
  vector<coordinate_vector> pos_mat;
  double orig_pot;
  unsigned int range_begin, range_end;
  vector<coordinate_vector> *result_vec_p;
};

void* force_thread(void *args_p) {
  struct force_thread_arg_t * args = (struct force_thread_arg_t *)args_p;
  vector<coordinate_vector> &pos_mat = args->pos_mat;
  double orig_pot(args->orig_pot);
  vector<coordinate_vector> &result_vec = *(args->result_vec_p);

  for (unsigned int i=args->range_begin; i<args->range_end; ++i) {
    for (unsigned int j=0; j<pos_mat[i].size(); ++j) {
      double orig_pos = pos_mat[i][j];
      pos_mat[i][j] += Delta;
      double new_pot = potentials(pos_mat);
      pos_mat[i][j]  = orig_pos;
      result_vec[i][j] = (orig_pot-new_pot)/Delta;
    }
  }
}

vector<coordinate_vector> acceleration(const vector<coordinate_vector> &pos_mat,
                                       const vector<double> &mass_mat) {
  double orig_pot = potentials(pos_mat);

  // Allocate result space
  vector<coordinate_vector> result_vec(pos_mat.size());
  for (unsigned int i=0; i<pos_mat.size(); ++i) 
    result_vec[i].resize(pos_mat[i].size());
  
  // Prepare arguments for each thread
  unsigned int BatchSize = pos_mat.size()/NThread;
  vector<force_thread_arg_t> force_thread_args(NThread);
  for (unsigned int i=0; i<force_thread_args.size(); ++i) {
    force_thread_args[i].pos_mat  = pos_mat;
    force_thread_args[i].orig_pot = orig_pot;
    force_thread_args[i].range_begin = i*BatchSize;
    force_thread_args[i].range_end = (i+1)*BatchSize;
    force_thread_args[i].result_vec_p = &result_vec;
  }

  // Add all leftover atoms to the last thread
  force_thread_args.back().range_end = pos_mat.size();

  // Spawn threads
  pthread_t children[NThread];
  for (unsigned int i=0; i<NThread; ++i) {
    pthread_create(&children[i], NULL, *force_thread, &force_thread_args[i]);
  }

  // Join threads
  for (unsigned int i=0; i<NThread; ++i) {
    pthread_join(children[i], NULL);
  }
  
  for (unsigned int i=0; i<result_vec.size(); ++i) {
    for (unsigned int j=0; j<result_vec[i].size(); ++j) {
      result_vec[i][j] /= mass_mat[i];
    }
  }
  return result_vec;
}

void velocity_verlet(vector<coordinate_vector> &pos_mat,
                     vector<coordinate_vector> &vel_mat,
                     const vector<double>      &mass_mat,
                     double                    steptime,
                     int                       box_size) {
  //Step 1: pos_mat += vel_mat*TimeStep + accel_mat*TimeStep*TimeStep/2
  vector<coordinate_vector> accel_mat = acceleration(pos_mat, mass_mat);
  for (unsigned int i=0; i<pos_mat.size(); ++i) {
    for (unsigned int j=0; j<pos_mat[i].size(); ++j) {
      pos_mat[i][j] += vel_mat[i][j]*steptime + accel_mat[i][j]*steptime*steptime/2;
      if (abs(pos_mat[i][j]) > box_size) {
        vel_mat[i][j] = -vel_mat[i][j];
      }
    }
  }

  //Step 2: new_half_vel = vel_mat + accel_mat*TimeStep/2
  vector<coordinate_vector> new_half_vel(vel_mat.size());
  for (unsigned int i=0; i<new_half_vel.size(); ++i) {
    new_half_vel[i].reserve(vel_mat[i].size());
    for (unsigned int j=0; j<new_half_vel[i].size(); ++j) {
      new_half_vel[i].push_back(vel_mat[i][j] + accel_mat[i][j]*steptime/2);
    }
  }

  //Step 3: new_vel_mat = vel_mat + (accel_mat+new_accel_mat)*TimeStep/2
  vector<coordinate_vector> new_accel_mat = acceleration(pos_mat, mass_mat);
  for (unsigned int i=0; i<vel_mat.size(); ++i) {
    for (unsigned int j=0; j<vel_mat[i].size(); ++j) {
      vel_mat[i][j] += (accel_mat[i][j]+new_accel_mat[i][j])*steptime/2;
    } 
  }
}

double energy(const vector<coordinate_vector> &pos_mat,
              const vector<coordinate_vector> &vel_mat,
              const vector<double>            &mass_mat) {
  double ener = 0;
  for (unsigned int i=0; i<vel_mat.size(); ++i) {
    double v2 = 0;
    for (unsigned int j=0; j<vel_mat[i].size(); ++j) {
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
               unsigned int step) {
  cout <<pos_mat.size() <<'\n'
       <<"Step: " <<step <<" Energy: " <<energy <<'\n';
  for (unsigned int i=0; i<pos_mat.size(); ++i) {
    cout <<'A';
    for (unsigned int j=0; j<pos_mat[i].size(); ++j) {
      cout <<'\t' <<pos_mat[i][j];
    }
    cout <<'\n' <<std::flush;
  }
}

int main() {
  vector<coordinate_vector> my_pos;
  vector<coordinate_vector> my_vel;
  vector<double> my_mass;
  std::default_random_engine rng;
  std::normal_distribution<double> gasdev;
  for (int i=-11;i<11;++i) {
    for (int j=-11;j<11;++j) {
      coordinate_vector p;
      p.push_back(i);
      p.push_back(j);
      my_pos.push_back(p);

      coordinate_vector v;
      v.push_back(gasdev(rng));
      v.push_back(gasdev(rng));
      my_vel.push_back(v);

      my_mass.push_back(1);
    }
  }

  double my_box=15;
  unsigned int TotalStep=50;
  unsigned int OutputInterval=10;
  double steptime = 5e-3;

  print_pos(my_pos, energy(my_pos, my_vel, my_mass), 0);
  for (unsigned int i=1; i<=TotalStep; ++i) {
    velocity_verlet(my_pos, my_vel, my_mass, steptime, my_box);
    if (i%OutputInterval==0) {
      print_pos(my_pos, energy(my_pos, my_vel, my_mass), i);
    }
  }
}
