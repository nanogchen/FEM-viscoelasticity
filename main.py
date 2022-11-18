#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:38:53 2019

@author: guangchen
"""

if __name__ == '__main__':
    import sys
    import ControlDat as CD
    import Solvers
    import preprocessing as pre
    import processing as pro
    import postprocessing as post

    # generate the mesh and local Gauss points
    xynodes = pre.GenerateMesh(CD.xbound, CD.ybound, CD.Nptxy)    
    GaussPts,GaussWts = pre.Generate2DGaussPts(CD.nGxy)
    
    # generate the topology information for local and global elements/nodes/dofs assemblage
    IEN, LM = pre.GetTopoInfo(CD.Nptxy, CD.dof, CD.Nen, CD.ElNode)
    
    # get the NodeList and No. of EBC dofs
    NodeList, NpEBC = pre.GetNodeList()
    
    # solve the given small strain linear elastic problem
    if (CD.problem_type == 1):
        # get the stiffness matrix
        KGmat, FGext = pro.GlobalAssembler(CD.Nen,IEN,LM,GaussPts,GaussWts,xynodes)    
        
        disp = Solvers.LinElasSolver(NpEBC, NodeList, xynodes, KGmat, FGext)
        post.SaveTxt('./data/Bending', disp)
        print(disp[-1])
        
    elif (CD.problem_type == 2):
        pass
    
    # viscoelastic problem
    elif (CD.problem_type == 3):
        # get the stiffness matrix
        KGmat, FGext = pro.GlobalAssembler(CD.Nen,IEN,LM,GaussPts,GaussWts,xynodes)    
        Solvers.Viscoelastic(NpEBC, NodeList, IEN, LM, GaussPts, GaussWts, xynodes, KGmat, FGext)
        
    # hard code viscoelastic problem
    elif (CD.problem_type == 4):
        strain,dev_strain,stress,HStress,SStress = Solvers.ViscoelasticHardCode()
         #find the end value
        print(HStress[0,-1])

    else:
        sys.exit('Error in setting the problem type. Exiting...')
    

    
    #post-processing
    #post.GetDispPlot()
    
    
    
 
    
    
    
    
    
    
    
    
    
    
    

    