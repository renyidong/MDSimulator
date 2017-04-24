#include <stdlib.h> 
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <pthread.h>
#include <mpi.h>
#ifdef __bgq__
#include <hwi/include/bqc/A2_inlines.h>
#else
#warning Not on BG/Q, use MPI_Wtime to mimic GetTimeBase
uint64_t GetTimeBase() {
  return MPI_Wtime()*(1L<<30);
}
#endif


// ============================= Constants ====================================
#define DIMENSION 3
const double BoxSize = 100.0;
const double TimeStep = 1E-3;
const double LJEplison=1.0, LJSigma=1.0, Mass=1.0;
const double Delta = 1.0e-12;


// =========================== Global Variables ===============================
int mpi_commsize, mpi_my_rank;

double *global_position;
double *local_position;
double *velocity, *acceleration, *old_acceleration;
double old_potential;
size_t global_size, rank_begin, rank_end, local_size;

int NThread;
int NRow;
unsigned long TotalSteps;

// ======================== Other Functions ===================================
double dist_square(const double *x1, const double *x2) {
  double dist2 = 0.0;
  for (int i=0; i<DIMENSION; ++i) {
    double d = x1[i] - x2[i];
    if      (d >  BoxSize) d -= 2*BoxSize;
    else if (d < -BoxSize) d += 2*BoxSize;
    dist2 += d*d;
  }
  return dist2;
}

double potentials(size_t overlay_index, const double *overlay_position) {
  double pot=0;
  for (size_t i=0; i<global_size; ++i) {
    for (size_t j=i+1; j<global_size; ++j) {
      double dist2;
      if (i==overlay_index) dist2 = dist_square(overlay_position, global_position + DIMENSION*j);
      else if (j==overlay_index) dist2 = dist_square(global_position + DIMENSION*i, overlay_position);
      else dist2 = dist_square(global_position + DIMENSION*i, global_position + DIMENSION*j);
      
      // repeat multiplication is faster then pow()
      double sigma6 = LJSigma*LJSigma*LJSigma*LJSigma*LJSigma*LJSigma;
      double dist6 = dist2*dist2*dist2;
      double sigdist6 = sigma6/dist6;
      pot += 4 * LJEplison * (sigdist6*sigdist6 - sigdist6);
    }
  }
  return pot;
}
void update_acceleration_one(size_t index) {
  double overlay[DIMENSION];
  for (int j=0; j<DIMENSION; ++j) overlay[j] = global_position[DIMENSION*(rank_begin+index) + j];
  for (int j=0; j<DIMENSION; ++j) {
    overlay[j] += Delta;
    double new_pot = potentials(index, overlay);
    acceleration[index*DIMENSION+j] = (old_potential-new_pot)/Delta/Mass;
    overlay[j] = global_position[DIMENSION*(rank_begin+index) + j];
  }
}

typedef struct {
  size_t begin;
  size_t end;
  pthread_barrier_t *barrier;
} potential_thread_arg_t;

void* potential_thread(void * vargs) {
  potential_thread_arg_t *args = (potential_thread_arg_t*)vargs;
  size_t begin = args->begin;
  size_t end = args->end;
  pthread_barrier_t *barrier = args->barrier;
  while(1) {
    pthread_barrier_wait(barrier);
    for (size_t i=begin; i<end; ++i) {
      update_acceleration_one(i);
    }
    pthread_barrier_wait(barrier);
  }
  return NULL;
} 

void update_acceleration() {
  old_potential = potentials(0, global_position);
  {
    double *tmp=old_acceleration;
    old_acceleration = acceleration;
    acceleration = tmp;
  }
  
  size_t thread_size = local_size/NThread;
  static int thread_ready=0;
  static pthread_barrier_t *barrier=NULL;
  if (!thread_ready) {
    thread_ready = 1;
    barrier = malloc(sizeof(pthread_barrier_t));
    pthread_barrier_init(barrier, NULL, NThread);
    assert(local_size%NThread == 0);
    for (int i=1; i<NThread; ++i) {
      potential_thread_arg_t *arg = malloc(sizeof(potential_thread_arg_t));
      arg->begin = i*thread_size;
      arg->end = (i+1)*thread_size;
      arg->barrier = barrier;
      pthread_t thread;
      pthread_create(&thread, NULL, &potential_thread, arg);
      pthread_detach(thread);
    }
  }
  
  pthread_barrier_wait(barrier);
  for (size_t i=0; i<thread_size; ++i) {
    update_acceleration_one(i);
  }
  pthread_barrier_wait(barrier);
}

void step() {
  // Step 1: generate new positions
  for (size_t i=0; i<local_size*DIMENSION; ++i) {
    local_position[i] = global_position[DIMENSION*rank_begin+i] + TimeStep*velocity[i] + TimeStep*TimeStep/2*acceleration[i];
    if      (local_position[i] >  BoxSize) local_position[i] -= 2*BoxSize;
    else if (local_position[i] < -BoxSize) local_position[i] += 2*BoxSize;
  }
  
  // Exchange positions
  MPI_Allgather(local_position, local_size*DIMENSION, MPI_DOUBLE, global_position, local_size*DIMENSION, MPI_DOUBLE, MPI_COMM_WORLD);
  
  // Step 2: generate new acceleration
  update_acceleration();
  
  // Step 3: generate new velocity
  for (size_t i=0; i<local_size*DIMENSION; ++i) {
    velocity[i] += TimeStep/2 * (acceleration[i] + old_acceleration[i]);
  }
}

int main(int argc, char* argv[]) {
  MPI_Init(0, NULL);
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
  MPI_Comm_rank( MPI_COMM_WORLD, &mpi_my_rank);
  
  if (argc != 4) {
    printf("Usage: mpirun -n NRank %s Steps NRow NThread\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  TotalSteps = atol(argv[1]);
  NRow = atoi(argv[2]);
  NThread = atoi(argv[3]);
  
  // Allocate
  if (mpi_my_rank==0) {
    global_size = NRow*NRow*NRow;
  }
  
  MPI_Bcast(&global_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  local_size = global_size/mpi_commsize;
  rank_begin = local_size*mpi_my_rank;
  assert(global_size%mpi_commsize == 0);
  
  global_position = malloc(sizeof(double)*DIMENSION*global_size);
  local_position  = malloc(sizeof(double)*DIMENSION*local_size);
  velocity        = malloc(sizeof(double)*DIMENSION*local_size);
  acceleration    = malloc(sizeof(double)*DIMENSION*local_size);
  old_acceleration= malloc(sizeof(double)*DIMENSION*local_size);
  
  // Initialize positions
  size_t inited = 0;
  if (mpi_my_rank == 0) {
    for (int x=-NRow/2; x<NRow/2; ++x)
      for (int y=-NRow/2; y<NRow/2; ++y)
        for (int z=-NRow/2; z<NRow/2; ++z) {
          double *p = global_position + DIMENSION*inited;
          p[0] = x+0.5;
          p[1] = y+0.5;
          p[2] = z+0.5;
          inited++;
        }
  }
  MPI_Bcast(global_position, DIMENSION*global_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  // Initialize acceleration
  update_acceleration();
  
  // Initialize velocity
  srand(GetTimeBase());
  for (size_t i=0; i<local_size*DIMENSION; ++i) {
    velocity[i] = (float)rand()/(float)RAND_MAX;
  }
  
  // Run
  MPI_Barrier(MPI_COMM_WORLD);
  double time_begin = MPI_Wtime();
  for (unsigned int t=0; t<TotalSteps; ++t) {
    step();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  double time_end = MPI_Wtime();
  
  // End
  if (mpi_my_rank==0) {
    printf("%f\n", time_end-time_begin);
  }
  
  MPI_Finalize();
  return 0;
}
