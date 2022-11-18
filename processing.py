#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 11:31:07 2019

@author: guangchen
"""

# this is the processing unit of the program

######################################################################################

def GetGStrain(IEN,LM,GaussPts,xynodes,disp):
    import numpy as np
    import preprocessing as pre
    import ControlDat as CD
    
    nrow = 3
    Gstrain = np.zeros((nrow*CD.Nint, CD.Nen))
    
    for i in range(CD.Nen):
        for j in range(CD.Nint):
            quadxy = xynodes[IEN[:,i],:]
            disp_e = disp[LM[:,i]]
            
            Nemat, Bemat, Jacobian = pre.ShapeFunc(GaussPts[j,:],quadxy)            
            Gstrain[j*nrow:j*nrow+nrow,i] = np.matmul(Bemat, disp_e) #the strain
    
    return Gstrain

######################################################################################
    
def GetGStress_GC_noBulkRelax(strain, inc_strain, CauchyStress,HStress):
    import numpy as np
    import math
    import sys
    import ControlDat as CD

    nrow = 3
    CauchyStress_new = np.zeros((nrow*CD.Nint, CD.Nen))    
    HStress_new = np.zeros((4*CD.Nint, CD.Nen))
    
    dt=CD.dt
    tau_1=CD.tau_1
    
    Cmat, mu_1, mu, K_mod, gamma_inf, gamma_1, lamd_star, mu_star,K_star = GetGMM_GC()
    beta_1 = dt/tau_1/(1+dt/tau_1)
    
    for i in range(CD.Nen):
        for j in range(CD.Nint):
            strain_e = strain[nrow*j:nrow*j+nrow,i] 
            
# The HS scheme
            HS = HStress[4*j:4*j+4, i]

            if (CD.i_stresstrain == 1 ):  #plane stress
                strain_33 = (gamma_1*math.exp(-dt/(2*tau_1))*SS[2] - gamma_1*math.exp(-dt/tau_1)*HS[2] - lamd_star*(strain_e[0] + strain_e[1]))/(K_mod + 4./3.*mu_star)
                strain_ext = np.array([strain_e[0], strain_e[1], strain_33, strain_e[2]*0.5])                    
            elif(CD.i_stresstrain == 2): #plane strain            
                strain_ext = np.array([strain_e[0], strain_e[1], 0, strain_e[2]*0.5])
                
            else:
                sys.exit('Error in i_stresstrain. Exiting...')
                
            # the split of strain
            dev_strain, sph_strain = GetStrainSplit(strain_ext)
            inc_strain_e, inc_strain_b = GetStrainSplit(strain_ext)
            
            # will use e_n and e_n+1 for the increment dev strain
            HS_tr = HS + 2*mu_1*inc_strain_e # for the Maxwell element
            CS_new = K_mod*sph_strain*3 + 2*mu*dev_strain - beta_1*HS_tr
            HS_new = (1-beta_1)*HS_tr
            
            HStress_new[4*j:4*j+4, i] = HS_new 
            CauchyStress_new[nrow*j:nrow*j+nrow, i] = CS_new[[0,1,3]]           

    return CauchyStress_new, HStress_new

######################################################################################
    
def GetGStress_HS_noBulkRelax(strain, CauchyStress,HStress,SStress):
    import numpy as np
    import math
    import sys
    import ControlDat as CD

    nrow = 3
    CauchyStress_new = np.zeros((nrow*CD.Nint, CD.Nen))    
    HStress_new = np.zeros((4*CD.Nint, CD.Nen))
    SStress_new = np.zeros((4*CD.Nint, CD.Nen))
#    HStress_new = np.zeros((nrow*CD.Nint, CD.Nen))
#    SStress_new = np.zeros((nrow*CD.Nint, CD.Nen))
    
    dt=CD.dt
    tau_1=CD.tau_1
    
    Cmat, Cmat_e, mu, K_mod, gamma_inf, gamma_1, lamd_star, mu_star,K_star = GetGMM()
    
    for i in range(CD.Nen):
        for j in range(CD.Nint):
            strain_e = strain[nrow*j:nrow*j+nrow,i] 
            
# The HS scheme
            HS = HStress[4*j:4*j+4, i]
            SS = SStress[4*j:4*j+4, i]

            if (CD.i_stresstrain == 1 ):  #plane stress
                strain_33 = (gamma_1*math.exp(-dt/(2*tau_1))*SS[2] - gamma_1*math.exp(-dt/tau_1)*HS[2] - lamd_star*(strain_e[0] + strain_e[1]))/(K_mod + 4./3.*mu_star)
                strain_ext = np.array([strain_e[0], strain_e[1], strain_33, strain_e[2]*0.5])                    
            elif(CD.i_stresstrain == 2): #plane strain            
                strain_ext = np.array([strain_e[0], strain_e[1], 0, strain_e[2]*0.5])
                
            else:
                sys.exit('Error in i_stresstrain. Exiting...')
                
            # the split of strain
            dev_strain, sph_strain = GetStrainSplit(strain_ext)
            dev_SStress = 2*mu*dev_strain

            SS_new, vol_temp = GetStressSplit(dev_SStress)
            HS_new = HS*math.exp(-dt/tau_1)+ math.exp(-dt/(2.*tau_1))*(SS_new - SS)
            CS_new = K_mod*sph_strain*3 + gamma_inf*SS_new + gamma_1*HS_new
            
            SStress_new[4*j:4*j+4, i] = SS_new
            HStress_new[4*j:4*j+4, i] = HS_new 
            
#            CauchyStress_new[nrow*j:nrow*j+nrow, i] = CS_new  
            CauchyStress_new[nrow*j:nrow*j+nrow, i] = CS_new[[0,1,3]]           

    return CauchyStress_new, HStress_new, SStress_new


######################################################################################
    
def GetGStress_HS_bulkRelax(strain, CauchyStress,HStress,SStress, MStress, PStress):
    import numpy as np
    import math
    import sys
    import ControlDat as CD

    nrow = 3
    CauchyStress_new = np.zeros((nrow*CD.Nint, CD.Nen))    
    HStress_new = np.zeros((4*CD.Nint, CD.Nen))
    SStress_new = np.zeros((4*CD.Nint, CD.Nen))
    PStress_new = np.zeros((4*CD.Nint, CD.Nen))
    MStress_new = np.zeros((4*CD.Nint, CD.Nen))
#    HStress_new = np.zeros((nrow*CD.Nint, CD.Nen))
#    SStress_new = np.zeros((nrow*CD.Nint, CD.Nen))
    
    dt=CD.dt
    tau_1=CD.tau_1
    
    Cmat, Cmat_e, mu, K_mod, gamma_inf, gamma_1, lamd_star, mu_star,K_star = GetGMM()
    
    for i in range(CD.Nen):
        for j in range(CD.Nint):
            strain_e = strain[nrow*j:nrow*j+nrow,i] 
            
            
            # The HS scheme
            HS = HStress[4*j:4*j+4, i]
            SS = SStress[4*j:4*j+4, i]
            PS = PStress[4*j:4*j+4, i]
            MS = MStress[4*j:4*j+4, i]
            
            if (CD.i_stresstrain == 1 ):  #plane stress
                strain_33 = -((K_star - 2/3*mu_star)*(strain_e[0] + strain_e[1]))/(K_star + 4./3.*mu_star) - (gamma_1*math.exp(-dt/tau_1)*HS[2] - gamma_1*math.exp(-dt/(2*tau_1))*SS[2])/(K_star + 4./3.*mu_star) - (gamma_1*math.exp(-dt/tau_1)*MS[2] - gamma_1*math.exp(-dt/(2*tau_1))*PS[2])/(K_star + 4./3.*mu_star)
                strain_ext = np.array([strain_e[0], strain_e[1], strain_33, strain_e[2]*0.5])                    
            elif(CD.i_stresstrain == 2): #plane strain            
                strain_ext = np.array([strain_e[0], strain_e[1], 0, strain_e[2]*0.5])
                
            else:
                sys.exit('Error in i_stresstrain. Exiting...')
                
            # the split of strain
            dev_strain, sph_strain = GetStrainSplit(strain_ext) # 4*1 vector
            
            # update the dev stress
            dev_SStress = 2*mu*dev_strain            
            SS_new, vol_temp = GetStressSplit(dev_SStress)
            HS_new = HS*math.exp(-dt/tau_1)+ math.exp(-dt/(2.*tau_1))*(SS_new - SS)
            
            # update the bulk stress
            PS_new = K_mod*sph_strain*3
            MS_new = MS*math.exp(-dt/tau_1) + math.exp(-dt/(2.*tau_1))*(PS_new - PS)
            
            # update the local cauchy stress
            CS_new = (gamma_inf*PS_new + gamma_1*MS_new) + (gamma_inf*SS_new + gamma_1*HS_new)
            
#            print(CS_new)
            
            #update the global stress
            SStress_new[4*j:4*j+4, i] = SS_new
            HStress_new[4*j:4*j+4, i] = HS_new
            PStress_new[4*j:4*j+4, i] = PS_new
            MStress_new[4*j:4*j+4, i] = MS_new
            

            CauchyStress_new[nrow*j:nrow*j+nrow, i] = CS_new[[0,1,3]]
#            print(CS_new[[2]])
# The HS scheme              

    return CauchyStress_new, HStress_new, SStress_new, MStress_new, PStress_new

######################################################################################

def GetStrainSplit(strain):
    import numpy as np
    
    dev_strain = np.zeros((4,1))
    sph_strain = np.zeros((4,1))
    
    t_strain = strain[0] + strain[1] + strain[2]
    sph_strain = np.array([t_strain/3., t_strain/3., t_strain/3., 0])
    dev_strain = strain - sph_strain    
    
    return dev_strain, sph_strain

######################################################################################
    
def GetStressSplit(stress):
    import numpy as np
    
    dev_stress = np.zeros((4,1))
    sph_stress = np.zeros((4,1))

    t_stress = stress[0] + stress[1] + stress[2]
    sph_stress = np.array([t_stress/3., t_stress/3., t_stress/3., 0.])
    dev_stress = stress - sph_stress  
    
    return dev_stress, sph_stress

######################################################################################

def Solve(npEBC, NodeList, Kmat, fvec, disp_in):
    import numpy as np
    
    Len = Kmat.shape[0]
    dispU = np.zeros((Len))
    
    K1,f1,disp1,transMatRight = rearrange(NodeList, Kmat, fvec, disp_in)
    
    Knn = K1[0:Len-npEBC, 0:Len-npEBC]
    Kne = K1[0:Len-npEBC, Len-npEBC:Len]
    fn = f1[0:Len-npEBC]
    dispE = disp1[Len-npEBC:Len]
    
    invKnn = np.linalg.inv(Knn)
    disp_Else = np.matmul(invKnn, fn-np.matmul(Kne, dispE))
    
    dispU[0:Len-npEBC] = disp_Else
    dispU[Len-npEBC:Len] = dispE
    
    dispU = np.matmul(transMatRight,dispU)
    
    return dispU


######################################################################################
    
def rearrange(NodeList, K0, f0, disp0):
    import numpy as np
    
    Ntotdof = NodeList.shape[0]
    transMatLeft = np.zeros((Ntotdof,Ntotdof),dtype=int)      

    for i in range(Ntotdof):
        for j in range(Ntotdof):
            if (j == NodeList[i]):
                transMatLeft[i,j] = 1
            else:
                transMatLeft[i,j] = 0
                
    transMatRight = transMatLeft.T
    
    K = np.matmul(np.matmul(transMatLeft,K0),transMatRight)
    f = np.matmul(transMatLeft,f0)
    disp = np.matmul(transMatLeft,disp0)    
    
    return K,f,disp, transMatRight

########################################################

def CalExactDisp(Ptx):
    import numpy as np
    import sys
    import ControlDat as CD
    
    disp = np.zeros((CD.dof))
    
    x = Ptx[0]
    y = Ptx[1]
    
    if (CD.EBCDisp_type == 1):        # for linear essential boundary
        Coord = np.array([1, x, y])
        disp = np.matmul(CD.LinMat, Coord)
        
    elif (CD.EBCDisp_type == 2):
        Coord = np.array([1, x, y, x**2, x*y, y**2])
        disp = np.matmul(CD.QuadMat, Coord)
    
    elif (CD.EBCDisp_type == 3):
        pass
    
    else:
        sys.exit('Error in prescribing the displacement')
    
    return disp


########################################################
    
def CalExactDomainDisp(xynodes):
    import numpy as np
    import ControlDat as CD
    
    npts = xynodes.shape[0]    
    dispDomain = np.zeros((CD.dof*npts))
    
    for i in range(npts):
        disp = CalExactDisp(xynodes[i])
        dispDomain[CD.dof*i] = disp[0]
        dispDomain[CD.dof*i+1] = disp[1]        
        
    return dispDomain

######################################################################################
    
def GetFGint(IEN,LM,GaussPts,GaussWts,xynodes, CauchyStress):
    import numpy as np
    import ControlDat as CD
    
    npts = xynodes.shape[0]
    dof = CD.dof
    Nen = CD.Nen
    FGint = np.zeros((dof*npts))
    
    for i in range(Nen):
        quadxy = xynodes[IEN[:,i],:]
             
        FEint = GetEleFint(dof,CD.ElNode,GaussPts,GaussWts,quadxy, CauchyStress[:,i])
        for j in range (dof*CD.ElNode):
            II = LM[j,i] #global index
            
            FGint[II] = FGint[II] + FEint[j]
    
    return FGint

######################################################################################

def GetEleFint(dof,ElNode,GaussPts,GaussWts,quadxy, CauchyStress_e):
    import numpy as np
    import ControlDat as CD
    import preprocessing as pre
    
    FEint = np.zeros((dof*ElNode))
    nrow = 3

    for i in range(CD.Nint):
        Nemat, Bemat, Jacobian = pre.ShapeFunc(GaussPts[i,:],quadxy)
        
        FEint = FEint + np.matmul(Bemat.T, CauchyStress_e[nrow*i:nrow*i+nrow])*(Jacobian)*GaussWts[i]
    
    return FEint

######################################################################################

def GetGMM():
    import numpy as np
    import sys
    import math
    import ControlDat as CD

    E_inf = CD.E_infty
    E_1 = CD.E_1
    nu = CD.nu
    dt = CD.dt
    tau_1 = CD.tau_1

    K_inf = E_inf/(3.*(1.-2.*nu))
    K_1 = E_1/(3.*(1.-2.*nu))
    mu_inf = E_inf/(2.*(1.+nu))
    mu_1 = E_1/(2.*(1.+nu))
    
    gamma_inf = mu_inf/(mu_inf+mu_1)
    gamma_1 = mu_1/(mu_inf+mu_1)
    
    gstar = gamma_inf + gamma_1*math.exp(-dt/(2.*tau_1))
#    gstar = gamma_inf + gamma_1*(1-math.exp(-dt/tau_1))/(dt/tau_1)
    # both above equations are tested, the first one is more accurate
    
    K_mod = K_inf + K_1
    mu = mu_inf + mu_1
    
    lamd = K_mod - 2.*mu/3.
    lamd_star = K_mod*gstar - 2.*mu/3.*gstar
    mu_star = mu*gstar
    K_star = K_mod*gstar
    
    if (CD.i_stresstrain == 1): # plane stress
        lamd_bar = 2.*lamd_star*mu_star/(lamd_star + 2.*mu_star)  #viscoelastic
#        lamd_e_bar = 2.*lamd*mu/(lamd + 2.*mu)                    #just elastic
        
        Cmat_e = np.array([[lamd+2*mu, lamd, lamd, 0.], [lamd, lamd+2*mu,lamd,0.], [lamd,lamd, lamd+2*mu,0.], [0.,0.,0., mu]])
        Cmat = np.array([[lamd_bar+2*mu_star, lamd_bar, 0.], [lamd_bar, lamd_bar+2*mu_star,0.], [0.,0.,mu_star]])
        
        # The value is bigger than the ref paper, thus is stiffer with smaller disp.
        
        # from the ref paper Wiechert table
#        Cmat = np.array([[5.0762 ,1.5229, 0], [1.5229, 5.0762, 0], [0, 0, 1.7767]])
#        Cmat= Cmat*1E5

    elif (CD.i_stresstrain == 2): # plane strain
        
        Cmat_e = np.array([[lamd+2*mu, lamd, 0.], [lamd, lamd+2*mu,0.], [0.,0.,mu]])
        Cmat = np.array([[lamd_star+2*mu_star, lamd_star, 0.], [lamd_star, lamd_star+2*mu_star,0.], [0.,0.,mu_star]])
    
    else:
        sys.exit("Something wrong happen in choosing the plane stress/strain model.")  
    
    return Cmat, Cmat_e, mu, K_mod, gamma_inf, gamma_1, lamd_star, mu_star,K_star

######################################################################################

def GlobalAssembler(Nen,IEN,LM,GaussPts,GaussWts,xynodes):
    import numpy as np
    import ControlDat as CD
    
    npts = xynodes.shape[0]
    dof = CD.dof
    
    KGmat = np.zeros((dof*npts, dof*npts))
    FGext = np.zeros((dof*npts))    
    
    for i in range(Nen):
        quadxy = xynodes[IEN[:,i],:]
#        disp_e = dispU[LM[:,i]]
        
        KEmat, FEext = GetEleKF(dof,CD.ElNode,GaussPts,GaussWts,quadxy)
        for j in range (dof*CD.ElNode):
            II = LM[j,i] #global index
            
            for k in range(dof*CD.ElNode):
                JJ = LM[k,i] #global index
                KGmat[II, JJ] = KGmat[II,JJ] + KEmat[j,k]
            
            FGext[II] = FGext[II] + FEext[j]
    
    return KGmat, FGext

######################################################################################

def GetEleKF(dof,ElNode,GaussPts,GaussWts,quadxy):
    import numpy as np
    import sys
    import preprocessing as pre
    import ControlDat as CD
    
    KEmat = np.zeros((dof*ElNode,dof*ElNode))
    FEext = np.zeros((dof*ElNode))
    
    if (CD.problem_type == 1):
        Dmat = GetElasticMat()
        
    elif (CD.problem_type == 2):
        pass
    
    elif (CD.problem_type == 3):
        Dmat, Cmat_e, mu, K_mod, gamma_inf, gamma_1, lamdStar, muStar,K_star = GetGMM()
        
    else:
        sys.exit('Error in problem type. Exiting...')
    
    for i in range(CD.Nint):
        Nemat, Bemat, Jacobian = pre.ShapeFunc(GaussPts[i,:],quadxy)
        GPtx = np.matmul(Nemat, quadxy.reshape(dof*ElNode,1))
        ForceGauss = GetBodyForce(GPtx)
        
        KEmat = KEmat + np.matmul(np.matmul(Bemat.T, Dmat),Bemat)*(Jacobian)*GaussWts[i]
        FEext = FEext + np.matmul(Nemat.T, ForceGauss)*(Jacobian)*GaussWts[i]
        
    # traction force
    #TractionGauss = GetTraction()
    
    return KEmat, FEext

######################################################################################

def GetBodyForce(Ptx):
    # this is to obtain the body force at the given Gauss point
    import ControlDat as CD
    import numpy as np
    import sys
    
    force = np.zeros((2))
    
    if (CD.BF_type == 1):
        force = np.zeros((2)) #no body force
        
    elif (CD.BF_type == 2):
        force = CD.bf
        
    else:
        sys.exit('Error in the body force type.')
    # Equation for body force given here
    
    return force

######################################################################################

def GetTraction(Ptx):
    # this is to obtain the traction force at the given Gauss point
    import numpy as np
    
    force = np.zeros((2))
    
    # Equation given here
    
    return force
    
######################################################################################

def GetElasticMat():
    import numpy as np
    import ControlDat as CD
    import sys
    
    Emod = CD.Emod
    Pois = CD.nu
    
    if(CD.i_stresstrain == 1): #plane stress
        Dmat = Emod/(1.-Pois**2)*np.array([[1.0,Pois,0.0],[Pois, 1.0, 0.0],[0., 0., 0.5*(1-Pois)]])
        
    elif(CD.i_stresstrain == 2):
        Dmat = Emod*(1-Pois)/((1+Pois)*(1-2*Pois))*np.array([[1., Pois/(1-Pois), 0.],[Pois/(1-Pois), 1., 0.],  [0., 0., (1-2*Pois)/(2*(1-Pois))]])
    
    else:
        sys.exit('Error input of plane hypothesis')
        
    return Dmat
    
######################################################################################
    
def AddTipLoad(idx, load, FGext):
    import ControlDat as CD
    import numpy as np
    
    FGext[CD.dof*idx: (CD.dof*idx + CD.dof)] = FGext[CD.dof*idx: (CD.dof*idx + CD.dof)] + np.array(load)
    
    return FGext
    