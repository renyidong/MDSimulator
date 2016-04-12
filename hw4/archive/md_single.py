from math import sqrt
from random import random

import sys
# ================ input ====================
LJSigma = 0.3     # Sigma in LJ
LJEplison = 0.661  # Eplison in LJ
Delta = 1e-8       # Delta x for calculate potential differencial
TimeStep = 1e-8
my_pos = []
for i in range(-1,2):
    for j in range(-1, 2):
        for k in range(-1, 2):
            my_pos.append([i*0.38 + 0.2*random(),
                           j*0.38 + 0.2*random(),
                           k*0.38 + 0.2*random()])
my_vel = [ [0,0,0] for _ in my_pos]
my_mass = [.0000004 for _ in my_pos]
# ================ input ====================

def distances(pos_mat):
    dist_mat = [[0.0 for i in pos_mat] for j in pos_mat ]
    for i in range(len(pos_mat)):
        for j in range(i+1,len(pos_mat)):
            dist2 = 0.0
            assert(len(pos_mat[i])==len(pos_mat[j]))
            for k in range(len(pos_mat[i])):
                dist2 += (pos_mat[i][k]-pos_mat[j][k])**2
            dist_mat[i][j]=dist_mat[j][i]=sqrt(dist2)
    return dist_mat

def potentials(dist_mat):
    pot = 0.0
    for i in range(len(dist_mat)):
        for j in range(i+1,len(dist_mat)):
            pot += 4*LJEplison * ( (LJSigma/dist_mat[i][j])**12 - (LJSigma/dist_mat[i][j])**6 )
    return pot
    

def forces(pos_mat):
    force_mat = [[0.0 for direction in pos] for pos in pos_mat]
    orig_pot = potentials(distances(pos_mat))
    for i in range(len(pos_mat)):
        for j in range(len(pos_mat[i])):
            orig_pos = pos_mat[i][j] 
            pos_mat[i][j] += Delta
            new_pot = potentials(distances(pos_mat))
            pos_mat[i][j] = orig_pos
            force_mat[i][j] = (orig_pot-new_pot)/Delta
    return force_mat

def acceleration(force_mat, mass_mat):
    assert(len(force_mat)==len(mass_mat))
    accel_mat =[[] for _ in force_mat]
    for i in range(len(force_mat)):
        accel_mat[i] = [ f/mass_mat[i] for f in force_mat[i] ]
    return accel_mat

def velocity_verlet(pos_mat, vel_mat, mass_mat):
    # Step 1
    accel_mat = acceleration(forces(pos_mat), mass_mat)
    #new_pos = pos_mat + vel_mat*TimeStep + accel_mat*TimeStep*TimeStep/2
    new_pos = [[0.0 for j in i] for i in pos_mat]
    for i in range(len(pos_mat)):
        for j in range(len(pos_mat[i])):
            new_pos[i][j] = pos_mat[i][j] + vel_mat[i][j]*TimeStep + accel_mat[i][j]*TimeStep*TimeStep/2
    
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

def print_pos(pos_mat,time):
    print(len(pos_mat))
    print('Atoms. Timestep:', int(time))
    for p in pos_mat:
        print('He',end='')
        for i in p:
            print('\t', '%.6f'%(i*5), end='')
        print('\n', end='')
    sys.stdout.flush()

print_pos(my_pos, 0)
for i in range(1,99990):
    my_pos, my_vel = velocity_verlet(my_pos, my_vel, my_mass)
    if i%10 == 0: 
        print_pos(my_pos, i/10)
