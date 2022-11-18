#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 17:47:41 2019

@author: guangchen
"""

import numpy as np
import sys

#########################################
############ geometry ###################
#########################################
xbound = np.array([0.0, 20.0])
ybound = np.array([-0.5, 0.5])
Nptxy = np.array([41, 3], dtype=int)
nGxy = np.array([2,2], dtype=int)
nGBC = 2
dof = 2

npts = Nptxy[0]*Nptxy[1]
Nint = nGxy[0]*nGxy[1]              #integration points
Nen = (Nptxy[0]-1)*(Nptxy[1]-1)     #elements numbers
ElNode = 4                          #element nodes: quad-element
ELGauss = nGxy[0]*nGxy[1]           #element Gauss points: quad-element

#########################################
############## problem type #############
#########################################
# 1 = linear elasticity, small strain;
# 2 = nonlinear elasticity, small strain;
# 3 = linear viscoelasticity, small strain; 
# 4 = hard coding for linear viscoelasticity 

problem_type = 3


##########################################
########### plane stress/strain ##########
##########################################
#1 = plane stress
#2 = plane strain

i_stresstrain = 1


##########################################
########### plane stress/strain ##########
##########################################
#1 = Hughes-Simo integral form
#2 = proposed differential form

i_stressupdate = 1

#########################################
########### boundary condition ##########
#########################################
# nodal index, dof
# 1 = all EBC;
# 2 = left EBC;
# 3 = user defined; 

EBCtype = 2

# if = 3, then give the following array
#        nodeID, dof
# EBC = [[1,0],
#        [1,1],
#        [4,0]]


#########################################
########### EBC prescribed form #########
#########################################
# nodal index, dof
# 1 = patch test: linear (u = a + bx + cy), or [u,v]^T = M*[x,y]
# 2 = patch test: quadratic
# 3 = fixed: zeros (u = 0)

# EBCDisp_type = 1
# LinMat = np.array([[0, 0.2, 0.3],[0, 0.1, 0.4]])

# EBCDisp_type = 2
# QuadMat = np.array([[0, 0.12, 0.14, 0.16, 0.18, 0.20],[0, 0.11, 0.13, 0.15, 0.17, 0.19]])

EBCDisp_type = 3


#########################################
############# Body force form ###########
#########################################
# 1 = no body force
# 2 = constant body force, e.g. g=[0,-9.81]
# 3 = user defined function

BF_type = 1

# if BF_type = 2, give the following value
# bf = np.array([ , ])


#########################################
############# Load force form ###########
#########################################
# 1 = no traction force
# 2 = tip load: give the load value and the load applied node index
# 3 = user defined load: load added on a side

#load_type = 1

load_type = 2
P0 = 1.0
load = np.array([0,-P0])
loadInx = npts-1

#load_type = 3
#P0 = 1.0
#load = np.array([P0,0]) # tension
#loadSide = 'right'


#########################################
############# material   type ###########
#########################################
# 1 = elastic
# 2 = viscoelastic, Maxwell element
# 3 = viscoelastic, Generalized Maxwell model


material_type = 3


#for Generalized-Maxwell Model: E_infty, E_1, tau_1
# problem_type = 3

if (material_type == 1):
    # problem_type = 1
    Emod = 3E7      #elastic modulus
    nu = 0.34        #nu: poisson ratio
    
elif (material_type == 2):
    #for Maxwell Model: E_infty, E_1, tau_1
    # problem_type = 3
    E_infty = 0
    E_1 = 10 # eta=5
    tau_1 = 1.3

elif (material_type == 3):
    # for GMM   
    nu = 0.3
    E_infty = 0.10E+06
    E_1 = 0.40E+06
    tau_1 = 1.0
    
#    # for GMM   
#    nu = 0.3
#    E_infty = 0 #0.10E+06
#    E_1 = 0.10E+06
#    tau_1 = 2.6 #1.0
    
    i_bulk = 1 # if i_bulk = 0, NOT consider bulk relaxation effect; if = 1, consider it
    
else:
    sys.exit('Error type of material. Exiting...')


#########################################
############# time parameters ###########
#########################################

dt = 0.1            #time step
DuraTime = 10.0     #duration for tip load
TolTime = 40.0      #total simulation time
dt_out = 1.0        #time length to output into vtk files













