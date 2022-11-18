#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 11:20:00 2019

@author: guangchen
"""

# this is the pre-processing unit of the program

###############################################################################

#           3 -------------- 4
#           |                |
#           |       |        |
#           |       |_ _     |
#           |                |
#           |                |
#           |                |
#           1 -------------- 2

def ShapeFunc(LocalPtx,Quadxy):
    
    import numpy as np
    
    Nemat=np.zeros((2,8))
    Bemat=np.zeros((3,8))
    Jacob=np.zeros((2,2))
    
    Ne,Be = IsoMapping(LocalPtx)
    
    Jacob = np.matmul(Be,Quadxy)
    Jacobian = np.linalg.det(Jacob)
    Bg = np.matmul(np.linalg.inv(Jacob),Be)
    
    for J in range(8):
        if (J%2 == 0):        #2 is the dof
            index = int(J/2)
            Nemat[0,J] = Ne[index]
            
            Bemat[0,J] = Bg[0,index]
            Bemat[2,J] = Bg[1,index]
        else:
            Nemat[1,J] = Ne[index]
    
            Bemat[1,J] = Bg[1,index]
            Bemat[2,J] = Bg[0,index]    
    
    return Nemat, Bemat, Jacobian

##########################################################3

def IsoMapping(LocalPtx):
    import numpy as np
    
    xi = LocalPtx[0]
    eta = LocalPtx[1]
    
    Ne = np.zeros((4))
    Be = np.zeros((2,4))

    Ne[0] = 0.25*(1 - xi)*(1 - eta)
    Ne[1] = 0.25*(1 + xi)*(1 - eta)
    Ne[2] = 0.25*(1 - xi)*(1 + eta)
    Ne[3] = 0.25*(1 + xi)*(1 + eta)
    
    Be[0,0] = 0.25*(-1)*(1 - eta)
    Be[1,0] = 0.25*(1 - xi)*(-1)

    Be[0,1] = 0.25*(+1)*(1 - eta)
    Be[1,1] = 0.25*(1 + xi)*(-1)

    Be[0,2] = 0.25*(-1)*(1 + eta)
    Be[1,2] = 0.25*(1 - xi)*(+1)

    Be[0,3] = 0.25*(1)*(1 + eta)
    Be[1,3] = 0.25*(1 + xi)*(+1)
    
    return Ne, Be

###############################################################################

def GetTopoInfo(Nptxy, dof, Nen, ElNode):
    import numpy as np
    # this is to get the topology information for local-global stiffness matrix transformation
    # ID(dof, npts)
    # IEN(EleNo, Nen)
    # LEM
    
    npx = Nptxy[0]
    npy = Nptxy[1]
    npts = npx*npy
    
    ID = np.zeros((dof,npts), dtype=int)
    IEN = np.zeros((ElNode,Nen), dtype=int)
    LM = np.zeros((dof*ElNode, Nen), dtype=int)
    
    #get ID
    for i in range(dof):
        for j in range(npts):
            ID[i,j] = int(j*2 + i)
    
    #get IEN
    for I in range(npx-1):
        for J in range(npy-1):
            iel = J*(npx-1) + I   # the index of the element
            inode = J*npx + I     # the index of the left-bottom node 
            
            IEN[:,iel] = [inode, inode+1, inode+npx, inode+npx+1]
            
    #get LM
    for ii in range(Nen):
        for jj in range(ElNode):
            for kk in range(dof):
                irow = dof*jj + kk
                inode = IEN[jj,ii]
                LM[irow, ii] = ID[kk, inode]                
    
    return IEN, LM

###############################################################################

def GetNodeList():
    import numpy as np
    import sys
    import ControlDat as CD
    
    # 1=all fixedEBC;
    # 2=left fixed EBC;
    # 3=user defined; 
    
    EBCtype = CD.EBCtype
      
    npx = CD.Nptxy[0]
    npy = CD.Nptxy[1]
    dof = CD.dof
    npts = CD.npts
    
    #get MeshID
    MeshID = np.zeros((npy,npx), dtype=int)
    for II in range(npx):
        for JJ in range(npy):
            MeshID[JJ,II] = JJ*npx + II
    
    MeshID = np.flipud(MeshID)
    
    # boundary List
    BList = MeshID[-1,:]
    RList = MeshID[:,-1]
    UList = MeshID[0,:]
    LList = MeshID[:,0]
    
    RList = RList[::-1] #flip 
    LList = LList[::-1] #flip
    
    if (EBCtype == 1):
        Len = 2*(npx+npy) - 4
        EBC = np.zeros((dof*Len,2), dtype=int)
        
        EBCNodeList_temp = np.concatenate((BList,RList[1:-1], LList[1:-1]))
        EBCNodeList = np.concatenate((EBCNodeList_temp, UList))
        EBCNodeList = np.sort(EBCNodeList)
        
        for i in range(Len):
            for j in range(dof):
                EBC[dof*i,:] = EBCNodeList[i],0
                EBC[dof*i+1,:] = EBCNodeList[i],1
    
    elif(EBCtype == 2):
        Len = npy
        EBC = np.zeros((dof*Len,2), dtype=int)
        EBCNodeList = LList;
        
        for i in range(Len):
            for j in range(CD.dof):
                EBC[dof*i,:] = EBCNodeList[i],0
                EBC[dof*i+1,:] = EBCNodeList[i],1
    
    elif(EBCtype == 3):
        EBC = CD.EBC
        
    else:
        sys.exit('Error in the EBC type')
    
    # to create the NodeList
    TolDof = dof*npts
    
    NpEBC = EBC.shape[0]
    EBCList = np.zeros((NpEBC), dtype=int)
    NodeList = np.arange((TolDof), dtype=int)
    ElseList = np.zeros((TolDof-NpEBC), dtype=int)
    
    # obtain the EBC List
    for I in range(NpEBC):
        EBCList[I] = EBC[I,0]*CD.dof + EBC[I,1]
    
    # the total index
    NodeList = np.arange((TolDof), dtype=int)
    
    # obtain ElseList
    for i in range(TolDof):
        for j in range(NpEBC):
            if (NodeList[i] == EBCList[j]):
                NodeList[i] = 0
    idx = 0    
    for i in range(TolDof):
        if (NodeList[i] != 0):
            ElseList[idx] = NodeList[i]
            idx = idx + 1
            
    # ordered Node List
    NodeList = np.concatenate((ElseList, EBCList))
        
    return NodeList, NpEBC

###############################################################################

def GenerateMesh(xbound, ybound, Nptxy):
    import numpy as np    
    
    npx = Nptxy[0]
    npy = Nptxy[1]
    npts = npx*npy
    xynodes = np.zeros((npts,2))    
    
    xpts = np.linspace(xbound[0],xbound[1],npx,endpoint=True)
    ypts = np.linspace(ybound[0],ybound[1],npy,endpoint=True)
    
    for i in range(npy):
        for j in range(npx):
            k = i*npx + j
            xynodes[k,0] = xpts[j]
            xynodes[k,1] = ypts[i]        
    
    return xynodes

###############################################################################

def Generate2DGaussPts(xynGaussPts):
# this is to generate the 2D elemental Gauss points
    import numpy as np

    xGPts,xGWs = Generate1DGaussPts(xynGaussPts[0])
    yGPts,yGWs = Generate1DGaussPts(xynGaussPts[1])
    
    npGauss = xynGaussPts[0]*xynGaussPts[1]
    GaussPoints = np.zeros((npGauss,2))
    GaussWeights = np.zeros((npGauss))
    
    for i in range(xynGaussPts[1]):
        for j in range(xynGaussPts[0]):
                k = i*xynGaussPts[0]+j;
                GaussPoints[k,0] = xGPts[j];
                GaussPoints[k,1] = yGPts[i];        
                
                GaussWeights[k] = xGWs[j]*yGWs[i];


    return GaussPoints, GaussWeights

###############################################################################

def Generate1DGaussPts(npGaussPts):
    import sys
    import numpy as np

    if(npGaussPts == 1):
        GaussPts = np.array([0.0])
        GaussWs = np.array([2.0])
        
    elif(npGaussPts == 2):
        GaussPts = np.array([-0.577350269, 0.577350269])
        GaussWs = np.array([1,1])
        
    elif(npGaussPts == 3):
        GaussPts = np.array([-0.774596669, 0, 0.774596669])
        GaussWs = np.array([5/9, 8/9, 5/9])
        
    elif(npGaussPts == 4):
        GaussPts = np.array([-0.861136312, -0.339981044, 0.339981044, 0.861136312])
        GaussWs = np.array([0.347854845, 0.652145155,0.652145155, 0.347854845])
        
    elif(npGaussPts == 5):
        GaussPts = np.array([-0.906179846, -0.538469310, 0, 0.538469310, 0.906179846])
        GaussWs = np.array([0.236926885, 0.478628670, 0.568888889, 0.478628670, 0.236926885])
        
    else:
        sys.exit('Gauss Points NOT included yet')


    return GaussPts, GaussWs