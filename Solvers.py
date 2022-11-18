
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 13:45:41 2019

@author: guangchen
"""

###############################################################################

# the solver for linear elasticity
def LinElasSolver(NpEBC, NodeList, xynodes, KGmat, FGext):
    import sys
    import processing as pro
    import ControlDat as CD
    
    #Add tip load
    if (CD.load_type == 1):
        pass
    
    elif (CD.load_type == 2):        
        FGext = pro.AddTipLoad(CD.loadInx, CD.load, FGext)
        
    elif (CD.load_type == 3):
        sys.exit('Not yet implemented yet.')
        
    else:
        sys.exit('Error in the load type.')    
    
    dispU0 = pro.CalExactDomainDisp(xynodes)
    
    dispU = pro.Solve(NpEBC, NodeList, KGmat, FGext, dispU0)    
    
    return dispU

###############################################################################
    
# the solver for linear elasticity
def Viscoelastic(NpEBC, NodeList, IEN, LM, GaussPts, GaussWts, xynodes, KGmat, FGext):
    from datetime import datetime
    import numpy as np
    import math
    import time
    import sys
    import processing as pro
    import postprocessing as post
    import ControlDat as CD
    
    # the variables
    Res = np.zeros((CD.dof*CD.npts))
    FGint = np.zeros((CD.dof*CD.npts))
    CauchyStress = np.zeros((3*CD.Nint, CD.Nen))        
#    HStress = np.zeros((3*CD.Nint, CD.Nen))
#    SStress = np.zeros((3*CD.Nint, CD.Nen))
    HStress = np.zeros((4*CD.Nint, CD.Nen))
    SStress = np.zeros((4*CD.Nint, CD.Nen))
    PStress = np.zeros((4*CD.Nint, CD.Nen))
    MStress = np.zeros((4*CD.Nint, CD.Nen))
    disp = np.zeros((CD.dof*CD.npts))
    disp_inc = np.zeros((CD.dof*CD.npts))
    disp_inc0 = np.zeros((CD.dof*CD.npts))
    
    # the wanted displacement
    DISP = np.zeros((int(CD.TolTime/CD.dt_out)+1))
    
    timer = 0.
    output_timer = 0.0
    step = 0
    output_step = 0    
    
    #Add tip load
    if (CD.load_type == 1):
        pass # no need to treat
        
    elif (CD.load_type == 2):
        load = CD.load
        idx = CD.loadInx
        FGext = pro.AddTipLoad(idx, load, FGext)
        
    elif (CD.load_type == 3):
        if (CD.loadSide == 'right'):
            load = CD.load/(CD.Nptxy[1])
            for i in range(CD.Nptxy[1]):
                idx = (i+1)*CD.Nptxy[0] - 1
                FGext = pro.AddTipLoad(idx, load, FGext)
                
    else:
        sys.exit('Error in the load type. Exiting...')
        
    
    # begin the loop
    print('Calculation begins at:', datetime.now().time()) #current time
    time_begin = time.time()

    print("{:10s}{:10s}{:20s}{:20s}{:29s}".format('Outstep', 'Step', 'Calculated time', 'Time elapsed', 'Real time'))
    
    # get inc_disp
    Res = FGext - FGint
    disp_inc = pro.Solve(NpEBC, NodeList, KGmat, Res, disp_inc0)
    # update total disp
    disp = disp + disp_inc
    DISP[0] = disp[-1]
    print(disp[-1])
    name = format(output_step, '04d')
    post.SaveTxt('./data/tension_'+name, disp)
    
    while True:
        
        # get strain
        Gstrain = pro.GetGStrain(IEN,LM,GaussPts,xynodes,disp)
        Ginc_strain = pro.GetGStrain(IEN,LM,GaussPts,xynodes,disp_inc)
    
        # get stress
        if (CD.i_bulk == 0): # with no bullk relaxation
            
            if (CD.i_stressupdate == 1): # Hughes-Simo
                CauchyStress_new, HStress_new, SStress_new = pro.GetGStress_HS_noBulkRelax(Gstrain,CauchyStress,HStress,SStress)
            elif (CD.i_stressupdate == 2):
                CauchyStress_new, HStress_new, SStress_new = pro.GetGStress_GC_noBulkRelax(Gstrain,Ginc_strain, CauchyStress,HStress,SStress)
            else:
                sys.exit('Error in the stress update scheme. Exiting...')
            
             # update the stress
            CauchyStress = CauchyStress_new
            HStress = HStress_new
            SStress = SStress_new
            
        elif (CD.i_bulk == 1): # with bulk relaxation
            
            if (CD.i_stressupdate == 1): # Hughes-Simo
                CauchyStress_new, HStress_new, SStress_new, MStress_new, PStress_new = pro.GetGStress_HS_bulkRelax(Gstrain,CauchyStress,HStress,SStress, MStress, PStress)
            elif (CD.i_stressupdate == 2):
                CauchyStress_new, HStress_new, SStress_new = pro.GetGStress_GC_bulkRelax(Gstrain,CauchyStress,HStress,SStress)
            else:
                sys.exit('Error in the stress update scheme. Exiting...')
            
            # update the stress
            CauchyStress = CauchyStress_new
            HStress = HStress_new
            SStress = SStress_new
            MStress = MStress_new
            PStress = PStress_new
        
        else:
            sys.exit('Error in i_bulk type. Exiting...')
    
        # get Fint
        FGint = pro.GetFGint(IEN,LM,GaussPts,GaussWts,xynodes, CauchyStress)
        
        # get inc_disp
        Res = FGext - FGint
        disp_inc = pro.Solve(NpEBC, NodeList, KGmat, Res, disp_inc0)        
        # update total disp
        disp = disp + disp_inc
    
        # update the recorder
        timer = timer + CD.dt
        output_timer = output_timer + CD.dt
        step = step + 1
        
        if (math.isclose(timer, CD.DuraTime, rel_tol=1e-5)):
            FGext = np.zeros((CD.dof*CD.npts))        
        
        if (math.isclose(output_timer, CD.dt_out, rel_tol=1e-5)):
            output_step = output_step + 1
            output_timer = 0.
            
            # vtkStepWriter            
            DISP[output_step] = disp[-1]
            name = format(output_step, '04d') #or: "{0:04d}".format(output_step)
            post.SaveTxt('./data/tension_'+name, disp)
            
            time_end = time.time()
            print("{:7d}{:7d}{:21.4f}{:20.4f}".format(output_step, step, timer, time_end-time_begin), '    ', datetime.now().time())
            print(disp[-1])
            
        if (timer > CD.TolTime):
            break
            
    time_end = time.time()
    print('Total computation time is:', time_end - time_begin)    
    print('Calculation ends at:', datetime.now().time())

    post.SaveTxt('./data/DISP', DISP)
#    plot the result
    post.ViscoElasPlot(DISP)
    
    return

###############################################################################
    
def ViscoelasticHardCode():
    import numpy as np
    import math
    import processing as pro
    import ControlDat as CD
    import matplotlib.pyplot as plt
    
    nu = CD.nu
    Ttotal=1.0
    nStep=100+1
    dt = Ttotal/nStep;
    relax_t = CD.tau_1
    Time = dt*np.arange(nStep)
    strain = np.zeros((4,nStep))
    dev_strain = np.zeros((4,nStep))
    stress = np.zeros((4,nStep))
    SStress = np.zeros((4,nStep))
    HStress = np.zeros((4,nStep))
    strain_inc = 1./1000.*np.array([1.0, -nu, -nu, 0.0])
    unit_vector = np.array([1,1,1,0])
    
    Cmat, Cmat_e, mu, K_mod, gamma_inf, gamma_1 = pro.GetGMM()
    
    # prescribed the strain
    for i in range (1,nStep):
        strain[:,i] = strain_inc*Time[i]
    
    # do the HS scheme
    for i in range (1,nStep):
        dev_strain[:,i], sph_strain = pro.GetStrainSplit(strain[:,i])

        SStress[:,i], sph_stress = pro.GetStressSplit(2*mu*dev_strain[:,i])
        HStress[:,i] = math.exp(-dt/relax_t)*HStress[:,i-1] + math.exp(-dt/(2*relax_t))*(SStress[:,i] - SStress[:,i-1])
        stress[:,i] = K_mod*sph_strain*unit_vector + HStress[:,i]
        
    # plot eStress_11 with time
    plt.subplot(2,1,1)
    plt.plot(Time,dev_strain[0,:])
    plt.title('dev_strain')
    plt.subplot(2,1,2)
    plt.plot(Time,HStress[0,:])
    plt.title('dev_stress')
    
    return strain,dev_strain,stress,HStress,SStress

    return

###############################################################################
    














