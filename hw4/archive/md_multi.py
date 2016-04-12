#!/bin/python2
from __future__ import print_function
from math import sqrt
from random import random
from copy import deepcopy
from multiprocessing import cpu_count
from random import gauss as gasdev
from multiprocessing import Pool as ProcessPool
import signal

import sys
# ================ input ====================
LJSigma = 1        # Sigma in LJ
LJEplison = 1   # Eplison in LJ
Delta = 1e-8       # Delta x for calculate potential differencial
TimeStep = 1e-2
TotalTime = 99990
OutputInterval = 5
NProc = 4
BatchSize = 9     # Assert: BatchSize % NAtoms == 0
my_pos = []
for i in range(-4,5):
    for j in range(-4,5):
            my_pos.append([i*1.5,
                           j*1.5])
my_vel = [ [gasdev(0,1),
            gasdev(0,1)] for _ in my_pos]
my_mass = [1 for _ in my_pos]
my_box = 6
# ================ input ====================

def potentials(pos_mat):
    pot = 0.0
    for i in range(len(pos_mat)):
        for j in range(i+1,len(pos_mat)):
            assert(len(pos_mat[i])==len(pos_mat[j]))
            dist2 = 0.0
            for k in range(len(pos_mat[i])):
                dist2 += (pos_mat[i][k]-pos_mat[j][k])**2
            pot += 4*LJEplison * ( (LJSigma**2/dist2)**6 - (LJSigma**2/dist2)**3 )
    return pot
    
def force_single(args):
    i_pos_mat, orig_pot, i = args
    pos_mat = deepcopy(i_pos_mat)
    batch_result = []
    for i in range(BatchSize*i, BatchSize*(i+1)):
        force_line = [0.0 for direction in pos_mat[i]]
        for j in range(len(pos_mat[i])):
            orig_pos = pos_mat[i][j] 
            pos_mat[i][j] += Delta
            new_pot = potentials(pos_mat)
            pos_mat[i][j] = orig_pos
            force_line[j] = (orig_pot-new_pot)/Delta
        batch_result.append(force_line)
    return batch_result

POOL=ProcessPool(NProc)
def forces(pos_mat):
    orig_pot = potentials(pos_mat)
    assert(len(pos_mat)%BatchSize == 0)
    args = [(pos_mat, orig_pot, i) for i in range(len(pos_mat)/BatchSize)]
    force_mat = []
    for r in POOL.map(force_single, args):
        force_mat.extend(r)
    return force_mat

def acceleration(force_mat, mass_mat):
    assert(len(force_mat)==len(mass_mat))
    accel_mat =[[] for _ in force_mat]
    for i in range(len(force_mat)):
        accel_mat[i] = [ f/mass_mat[i] for f in force_mat[i] ]
    return accel_mat

def velocity_verlet(pos_mat, vel_mat, mass_mat, box_size):
    # Step 1
    accel_mat = acceleration(forces(pos_mat), mass_mat)
    #new_pos = pos_mat + vel_mat*TimeStep + accel_mat*TimeStep*TimeStep/2
    new_pos = [[0.0 for j in i] for i in pos_mat]
    for i in range(len(pos_mat)):
        for j in range(len(pos_mat[i])):
            new_pos[i][j] = pos_mat[i][j] + vel_mat[i][j]*TimeStep + accel_mat[i][j]*TimeStep*TimeStep/2
            if abs(new_pos[i][j]) > box_size:
                # reverse velocity / bounce back
                vel_mat[i][j] = -vel_mat[i][j]
    
    # Step 2
    #new_half_vel = vel_mat + accel_mat*TimeStep/2
    new_half_vel = [[0.0 for j in i] for i in vel_mat]
    for i in range(len(vel_mat)):
        for j in range(len(vel_mat[i])):
            new_half_vel[i][j] = vel_mat[i][j] + accel_mat[i][j]*TimeStep/2
    
    # Step 3
    new_accel_mat = acceleration(forces(new_pos), mass_mat)
    #new_vel_mat = vel_mat + (accel_mat+new_accel_mat)*TimeStep/2
    new_vel_mat = [[0.0 for j in i] for i in vel_mat]
    for i in range(len(vel_mat)):
        for j in range(len(vel_mat[i])):
            new_vel_mat[i][j] = vel_mat[i][j] + (accel_mat[i][j]+new_accel_mat[i][j])*TimeStep/2
    return (new_pos, new_vel_mat)

def energy(pos_mat, vel_mat, mass_mat):
    ener = 0.0
    for i in range(len(vel_mat)):
        ener += mass_mat[i] * sum([v**2 for v in vel_mat[i]])
    ener /= 2
    ener += potentials(pos_mat)
    return ener
    
def print_pos(pos_mat,vel_mat, mass_mat, step):
    print(len(pos_mat))
    print('Atoms. Step: ', step, ' Energy:', energy(pos_mat, vel_mat, mass_mat))
    for p in pos_mat:
        print('A',end='')
        for i in p:
            print('\t', '%.6f'%i, end='')
        print('\n', end='')
    sys.stdout.flush()

print_pos(my_pos, my_vel, my_mass, 0)
for i in range(1,TotalTime):
    my_pos, my_vel = velocity_verlet(my_pos, my_vel, my_mass, my_box)
    if i%OutputInterval == 0: 
        print_pos(my_pos, my_vel, my_mass, i)
