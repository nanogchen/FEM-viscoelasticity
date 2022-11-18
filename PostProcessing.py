#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 10:47:02 2019

@author: guangchen
"""

def SaveTxt(Filename, dispU):
    import numpy as np
#    import ControlDat as CD
    
#    strMesh = str(CD.Nptxy[0]) + '-' + str(CD.Nptxy[1]) + '_' + str(CD.nGxy[0]) + '-' + str(CD.nGxy[1]) 
#    Filename = Filename + strMesh + '.txt'
    Filename = Filename + '.txt'
    
    np.savetxt(Filename, dispU, delimiter = ' ', newline = '\n')
    
    return 

########################################################

def GetDispPlot():
    import numpy as np
    import matplotlib.pyplot as plt
    import ControlDat as CD
    
    npx = CD.Nptxy[0]
    npy = CD.Nptxy[1]
    
    #read numerical results
    solDisp = np.loadtxt('./data/CalDisp.txt')
    solU = solDisp[::2]
    solV = solDisp[1::2]
    
    solU = solU.reshape((npx,npy))
    solV = solV.reshape((npx,npy))
    
    #read exact results
    exactDisp = np.loadtxt('./data/ExactDisp.txt')
    exactU = exactDisp[::2]
    exactV = exactDisp[1::2]
    
    exactU = exactU.reshape((npx,npy))
    exactV = exactV.reshape((npx,npy))
    
    xpts = np.linspace(CD.xbound[0], CD.xbound[1], CD.Nptxy[0], endpoint=True)
    ypts = np.linspace(CD.ybound[0], CD.ybound[1], CD.Nptxy[1], endpoint=True)
    
    XX,YY = np.meshgrid(xpts,ypts)
    
    #
    plt.figure(1,figsize=(18,16),dpi=80)
    fig1, axs = plt.subplots(2,1)
    ax = axs[0,0]
    CS = ax.contourf(XX,YY,solU)
    ax.clabel(CS, inline=True)
    ax.axis('equal')
    ax.set_title('Approx value')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.legend('u')
    
    ax = axs[0,1]
    CS = ax.contourf(XX,YY,exactU)
    ax.clabel(CS, inline=True)
    ax.axis('equal')
    ax.set_title('Exact value')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.legend('u')
    
    fig1.show()
    fig1.savefig('./data/dispU.png')
    
    #save the results into .csv format
    
    return

########################################################
    
def ViscoElasPlot(dispU):
    import numpy as np
    import matplotlib.pyplot as plt
    import ControlDat as CD
    
    t0 = 0
    t1 = CD.TolTime
    dt = CD.dt_out
    
    t = np.linspace(t0, t1, int(t1/dt)+1, endpoint=True)
        
    plt.figure(1)
    plt.plot(t,dispU)
    plt.show()
    
    return 

########################################################
    
def vtkStepWriter(stepOut,np,npGauss,xynodes,xyGaussPts, Disp,strain, CauchyStress,HStress,SStress, Fint, Fext):
    
    NodalPtsWriter(stepOut,np,xynodes,Disp, Fint, Fext)
    GaussPtsWriter(stepOut,npGauss,xyGaussPts,strain, CauchyStress,HStress,SStress)
    
    return

########################################################
    
def NodalPtsWriter(stepOut,np,xynodes,Disp, Fint, Fext):
    
    
    
    
    
    return

########################################################
    
def GaussPtsWriter(stepOut,npGauss,xyGaussPts,strain, CauchyStress,HStress,SStress):
    
    
    
    return

########################################################
    
def vtkHeader(f, npts):
    
    f.writelines("{:26s}".format('# vtk DataFile Version 4.0') + '\n')
    f.writelines("{:12s}".format('MyRK 2D file') + '\n')
    f.writelines("{:5s}".format('ASCII') + '\n')
    f.writelines("{:25s}".format('DATASET UNSTRUCTURED_GRID') + '\n')
    f.writelines("{:7s}{:10d}{:15s}".format('POINTS ',npts,' double ') + '\n')
       
    return

########################################################
















