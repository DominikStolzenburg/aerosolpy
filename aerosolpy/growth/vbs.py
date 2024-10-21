# -*- coding: utf-8 -*-

import numpy as np
import bisect
from scipy.integrate import solve_ivp
#import necessary base package functions and classes 
from aerosolpy.growth.h2so4 import SulfuricAcid


class VbsModel(SulfuricAcid):
    """
    a class for simulating monodisperse aerosol growth from volatility basis
    set (VBS), sulfuric acid and particle phase reactions

    Parameters
    ----------
    time : array_like
        time steps at which vapor concentrations are measured in [min],
        must be given as array even if vbs_traces and sa_traces are 1d/0d
    vbs_traces : array_like
        1d or 2d array of vapor concentrations in [cm-3] per VBS bin
    vbs_mass : array_like
        mean molecular mass of each VBS bin
    vbs_logC : array_like
        logarithm of saturation mass concentration [ug m-3] of each VBS bin
    sa_trace : string or array_like, optional
        whether to include sulfuric acid or not into the calculations
        default is 'None', if array_like 1d array of sulfurci acid 
        concentrations in [cm-3]
    activity : string or 3-tuple, optional
        whether to solve for particle phase activity coefficients or not
        default is 'unity', if tuple, it contains the O-C non-linearity 
        coefficient to be used, nC and nO for each volatility bin
    particle_phase : string or float or 3-tuple
        whether to calculate particle phase reactions or not.
        default is 'None'. 
    temp_kelvin : float, optional
        temperature in [K], default 300 K
    rh : float, optional
        relative humidity as fraction of 1, default 0.5
        only important if sa_traces is not 'None'

        
    Notes
    ----------
    ordinary differential equation solver using Eulerian forward integration
    solver operates with time in the scale of minutes, solver time step can
    be set by user in _solve function
    
    The system is solved in log-space to avoid negative growth rates.
    """
    def __init__(self,time,vbs_traces,vbs_mass,vbs_logC,
                 sa_trace='None',
                 activity='unity',
                 particle_phase='None',
                 temp_kelvin=300,rh=50,
                 **kwargs):
        
        # load time axis, needs to be equal for SA and VBS
        self.time = np.array(time)
        
        # fixed parameters, should be changed by experienced user only
        self.n_p = 1e3 
        self.dp_seed = 0.9

        if  isinstance(sa_trace, str):
            if sa_trace=='None':
                self.rho_sa = 1800 
                self.cs_sa_init = (np.pi/6.*(self.dp_seed*1e-9)**3
                                   *(self.n_p*1e6)*self.rho_sa*1e9)
                self.interact_sa = 0.75 
                self.cv_sa = np.zeros(len(self.time))
            else:
                raise ValueError(sa_trace,"needs to be 'None' or np.array")
        else:
            ## properties of SA 
            ## acting as a surrogate species for H2SO4+H2O+NH3 mixture      
            self.rho_sa = 1600      # H2SO4+H2O(+NH3) density [kg m-3]
            self.m_sa = 134.        # vapor cluster mass in amu 
            self.dp_sa = (((6./np.pi)*self.m_sa*1.6605e-27
                           /(self.rho_sa))**(1/3.)
                          *1e9) # size of SA vapor in nm
            # seed is made of initial condensed phase SA
            self.cs_sa_init = (np.pi/6.*(self.dp_seed*1e-9)**3
                               *(self.n_p*1e6)*self.rho_sa*1e9)
            self.interact_sa = 1
        
            # load H2SO4 input data 
            sa_trace = np.array(sa_trace)
            # extend input data if input is 0d
            if (sa_trace.ndim==0):
                sa_trace = np.repeat(sa_trace, len(self.time), axis=0)
            # convert SA concentration to saturation mass concentration
            self.cv_sa = sa_trace[:]*1e6/6.022e23*self.m_sa*1e6

        ## properties of organics
        ## simulated within the VBS scheme
        self.dk = 4.8*(300./temp_kelvin) # Kelvin-diameter in nm
        self.rho_org = 1400 # equal for all VBS bins
        # load organics input data
        self.log_c0 = np.array(vbs_logC)
        self.c0 = 10**self.log_c0
        self.m_org = np.array(vbs_mass)
        self.dp_org = (((6/np.pi)*(self.m_org*1.6605e-27/self.rho_org))**(1/3.)
                       *1e9) # in nm
        vbs_traces = np.array(vbs_traces)
        # extend gas-phase data if input is 0d
        if ((vbs_traces.ndim==1) or (vbs_traces.shape[0]==1)):
            vbs_traces = np.repeat(vbs_traces.flatten()[:,np.newaxis],
                                   len(self.time),axis=1).transpose()    
        # convert organic concentration to saturation mass concentration
        self.cv_org = vbs_traces[:,:]*1e6/6.022e23*self.m_org[:]*1e6
        self.n_bins = self.cv_org.shape[1]
        
        #activity coefficients
        if activity=='unity':
            self.solve_activity = False
        elif len(activity)==3:
            self.solve_activity = True
            self.b_co = activity[0]
            self.n_c = np.array(activity[1])
            self.n_o = np.array(activity[2])
        else:
            raise ValueError(activity,
                             ("if not 'unity' must be 3 element tuple/list"))
        
        self.particle_phase = particle_phase

        super().__init__(temp_kelvin=temp_kelvin, rh=rh, **kwargs)
    
    def _dynVBS(self,t,Ccomb):
        """
        set of differential equations for VBS dynamics
        
        Parameters
        ----------
        t : np.array
            time vector
        Ccomb : np.array
            combined vector of initial concentrations Cs,Cv
        
        Returns
        ----------
        np.array
            dCcombdt, differentials of vapor and condensed phase mass
            concentrations
        
        Notes
        ----------
        does not allow for changes in n_p and Cseed. 
        """

        # disentengle input
        CvSa = Ccomb[:1] # gaseous sulfuric acid
        CvOrg = Ccomb[1:self.numBins+1]
        # solution of the problem is achieved in the logCs space 
        # this avoids negative values. However, init values cannot be 0 
        # but must take a very small number if necessary
        CsSa = np.exp(Ccomb[1+self.numBins:self.numBins+2]) 
        CsOrg = np.exp(Ccomb[2+self.numBins:2*self.numBins+2])

        # raoult term
        COA = np.sum(CsOrg)
        Cpart = COA + self.interact_sa*CsSa
        
        #calculate mass-weighted organic activity coefficients
        if self.solve_activity==True:
            fC_s = np.ones(self.numBins)
            for i in range(self.numBins):
                nC_Cs_weighted = 0
                nO_Cs_weighted = 0
                for j in range(self.numBins):
                    if j!=i:
                        nC_Cs_weighted = nC_Cs_weighted+self.nC[j]*CsOrg[j]
                        nO_Cs_weighted = nO_Cs_weighted+self.nO[j]*CsOrg[j]
                fC_s[i] = 1./(1+nO_Cs_weighted/nC_Cs_weighted)
            #calculate solute carbon fractions
            fC_i = 1./(1+self.nO[:]/self.nC[:])
            # calculate activity coefficients
            gamma_i = np.exp(-2*(self.nC[:]+self.nO[:])
                             *((fC_i[:])**2+(fC_s[:])**2-2*fC_i[:]*fC_s[:])
                             *(self.bCO)*(690/self.T)
                             )
            gamma_i = np.nan_to_num(gamma_i,nan=1.0)
        else:
            gamma_i = np.ones(self.numBins)
        
        
        # particle volume,diameter and mass from condensed mass concentration
        Vp = ((CsSa/self.rho_sa + COA/self.rho_org)
              /(self.n_p*1e6*1e9)
              )
        Dp = ((6/np.pi)*Vp)**(1/3.) * 1e9
        

        beta_i_p_org = self.coll_kernel_vp(self.dpOrg[:], Dp,
                                           self.rho_org, self.rho_org)

        beta_i_p_sa = self.coll_kernel_vp_vdw(self.dp_sa, Dp,
                                              self.rho_sa, self.rho_sa,
                                              hamaker=5.2e-20,
                                              diff_coeff_v=self.diff_coeff_h2so4(),
                                              method='sceats') 

        #kinetic H2SO4+H2O condensation
        dlogCsSa_dt = ( (self.n_p*1e6) * (beta_i_p_sa*60) 
                       * (CvSa/CsSa)
                       ) 
        #organic condensation
        dlogCsOrg_dt = ( (self.n_p*1e6) * (beta_i_p_org*60) 
                        * ((CvOrg[:]/CsOrg[:]) 
                           - self.Co[:]*gamma_i[:]*(1./Cpart)*10**(self.dk/Dp)
                           )
                        ) 
        
        # use biscetion search to find proper derivative
        idx = bisect.bisect_left(self.time,t)

        dCvSa_dt = ((self.gasConcSa[idx]-self.gasConcSa[idx-1])
                    /(self.time[idx]-self.time[idx-1])
                    )
        dCvOrg_dt = ((self.gasConcOrg[idx,:]-self.gasConcOrg[idx-1,:])
                     /(self.time[idx]-self.time[idx-1])
                     )
        
        if self.particle_phase!='None':
            if isinstance(self.particle_phase, float):
                L = self.particle_phase
                P = self.particle_phase
                for i in range(self.numBins):
                    if i==0:
                        dlogCsOrg_dt[i] = dlogCsOrg_dt[i] - L
                    if ((i > 0) and (i!=self.numBins-1)):
                        dlogCsOrg_dt[i] = dlogCsOrg_dt[i] - L + P*CsOrg[i-1]/CsOrg[i]
                    if i==self.numBins-1:
                        dlogCsOrg_dt[i] = dlogCsOrg_dt[i] + P*CsOrg[i-1]/CsOrg[i]
                        
            elif isinstance(self.particle_phase, tuple):
                L = self.particle_phase[0]
                P = self.particle_phase[0]
                k = self.particle_phase[1]
                for i in range(self.numBins):   
                    if i<len(L):
                            dlogCsOrg_dt[i] = dlogCsOrg_dt[i] - L[i]
                    elif i>=k and i<len(L)+k:
                            dlogCsOrg_dt[i] = dlogCsOrg_dt[i]  + P[i-k]*CsOrg[i-k]/CsOrg[i]

            else:
                ValueError(self.particle_phase)
        
        return np.concatenate((np.array([dCvSa_dt]),
                               dCvOrg_dt,
                               dlogCsSa_dt,
                               dlogCsOrg_dt))
    
    def _solve(self,dt=0.05):
        """
        solver for set of differential equations
        
        Parameters
        ----------
        dt : float, optional
            time step [min] for internal solver calculation, default 0.05
        
        Notes
        ----------
        uses scipy.integrate solve_ivp, does not provide non-negativity 
        constraint, which is different to MATLAB version
        """
        # initial condition for differentials, vapor Cv at t=0 and logCs at t=0
        # while CsOrg ideally is 0, logCs cannot be inf
        # therefore logCsOrg approximated to a tiny value
        Cinit = np.concatenate((np.array([self.gasConcSa[0]]),
                                self.gasConcOrg[0,:],
                                np.log(np.array([self.cs_sa_init])),
                                -25*np.ones(self.numBins))
                               )


        # time interval for solution
        t0 = self.time[0] 
        tend = self.time[-1]

        sol = solve_ivp(self._dynVBS, [t0, tend], Cinit, 
                        t_eval=np.arange(t0,tend,dt),
                        method='BDF',
                        rtol=1e-12,atol=1e-12)
        ts = sol.t
        ys = sol.y

        return ts,ys
    
    def calc_vbs_dynamics(self):
        """
        calls VBS growth solver, calculates the system state variables for each 
        diameter including growth rates
            

        Returns:
        ----------
        tuple
            system state variables for each time/diameter step:
            diameter [nm], 
            total growth rate [nm h-1], 
            growth rate per vbs bin [nm h-1], 
            total condensed organic mass [ug m-3],
            condensed phase mass concentration per bin [ug m-3],
            vapor phase concentration [cm-3],
            vapor phase mass concentration per bin [ug m-3],
            mass flux to condensed phase per bin [ug m-3],
            saturation ratio per bin [ug m-3]
                    
        Notes
        ----------
        call for full output including all concentrations, fluxes, sat. ratio
        """
        t_prod, Cprod = self._solve()
            
        # calculates concentration and diameter evolution from solver solution
        CvSa_prod = Cprod[:1,:]
        CvOrg_prod = Cprod[1:self.numBins+1,:]
        CsSa_prod = np.exp(Cprod[self.numBins+1:self.numBins+2,:])
        CsOrg_prod = np.exp(Cprod[self.numBins+2:2*self.numBins+2,:])
        
        COA_prod = np.sum(CsOrg_prod,axis=0)
        Cpart_prod = COA_prod + self.interact_sa*CsSa_prod
        
        n_p_prod = self.n_p*np.ones(t_prod.shape[0])        
        Vp_prod = ((CsSa_prod[:]/self.rho_sa + COA_prod[:]/self.rho_org)
                   /((n_p_prod[:]*1e6)*1e9)
                   )
        Dp_prod = ((6/np.pi)*Vp_prod)**(1/3.) * 1e9

        if self.solve_activity==True:
            #calculate mass-weighted solvent carbon fraction
            fC_s = np.ones((self.numBins, t_prod.shape[0]))
            OtoC_s = np.ones((self.numBins, t_prod.shape[0]))
            for i in range(self.numBins):
                nC_Cs_weighted = np.zeros(t_prod.shape[0])
                nO_Cs_weighted = np.zeros(t_prod.shape[0])
                for j in range(self.numBins):
                    if j!=i:
                        nC_Cs_weighted[:] = (nC_Cs_weighted[:]
                                             +self.nC[j]*CsOrg_prod[j,:]
                                             )
                        nO_Cs_weighted[:] = (nO_Cs_weighted[:]
                                             +self.nO[j]*CsOrg_prod[j,:]
                                             )
                OtoC_s[i,:] = nO_Cs_weighted[:]/nC_Cs_weighted[:]
                fC_s[i,:] = 1./(1+nO_Cs_weighted[:]/nC_Cs_weighted[:])
            #calculate solute carbon fractions
            fC_i = 1./(1+self.nO[:]/self.nC[:])
            #OtoC_i = self.nO[:]/self.nC[:]
            # calculate activity coefficients
            gamma_prod = np.exp(-2*(self.nC[:,np.newaxis]+self.nO[:,np.newaxis])
                                *((fC_i[:,np.newaxis])**2+(fC_s[:,:])**2
                                  -2*fC_i[:,np.newaxis]*fC_s[:,:]
                                  )
                                *(self.bCO)*(690/self.T)
                                )
        else:
            gamma_prod = np.ones((self.numBins,t_prod.shape[0]))

        
        beta_i_p_org = self.coll_kernel_vp(self.dpOrg[:,np.newaxis],
                                           Dp_prod,
                                           self.rho_org,
                                           self.rho_org)

        beta_i_p_sa = self.coll_kernel_vdw(self.dp_sa,
                                           Dp_prod,
                                           self.rho_sa,
                                           self.rho_sa,
                                           hamaker=5.2e-20,
                                           diff_coeff_v=self.DvSa,
                                           method='sceats') 
        
        #kinetic H2SO4+H2O condensation
        dCsSa_dt = ( (self.n_p*1e6) * (beta_i_p_sa[:])
                    * CvSa_prod[:]
                    )
        ddpdt_sa = ( 2/(np.pi*self.rho_sa*(Dp_prod[:]*1e-9)**2*self.n_p*1e6)
                    *dCsSa_dt
                    )*3600
        #organic condensation
        dCsOrg_dt = ( (self.n_p*1e6) * (beta_i_p_org[:,:])   
                     * (CvOrg_prod[:,:] - gamma_prod[:,:]*CsOrg_prod[:,:]/Cpart_prod[:]*10**(self.dk/Dp_prod[:]) * self.Co[:,np.newaxis])
                     )
        if self.particle_phase!='None':
            if isinstance(self.particle_phase, float):
                L = self.particle_phase
                P = self.particle_phase
                for i in range(self.numBins):
                    if i==0:
                        dCsOrg_dt[i] = dCsOrg_dt[i] - L*CsOrg_prod[i]
                    if ((i > 0) and (i!=self.numBins-1)):
                            dCsOrg_dt[i] = dCsOrg_dt[i] - L*CsOrg_prod[i] + P*CsOrg_prod[i-1]
                    if i==self.numBins-1:
                        dCsOrg_dt[i] = dCsOrg_dt[i] + P*CsOrg_prod[i-1]
                        
            elif isinstance(self.particle_phase, tuple):
                L = self.particle_phase[0]
                P = self.particle_phase[0]
                k = self.particle_phase[1]
                for i in range(self.numBins):
                    if i<len(L):
                            dCsOrg_dt[i] = dCsOrg_dt[i] - L[i]*CsOrg_prod[i]
                    if i>=k and i<len(L)+k:
                            dCsOrg_dt[i] = dCsOrg_dt[i]  + P[i-k]*CsOrg_prod[i-k]

            else:
                ValueError(self.particle_phase)
            
        ddpdt_org = (2/(np.pi*self.rho_org*(Dp_prod[:]*1e-9)**2*self.n_p*1e6)
                     *dCsOrg_dt
                     )*3600
        
        
        Dp_prod = np.squeeze(Dp_prod)
        gr_bins = np.concatenate((ddpdt_sa,ddpdt_org),axis=0)
        gr_tot = np.sum(gr_bins,axis=0)
        
        return (Dp_prod,gr_tot,gr_bins,
                COA_prod,CvSa_prod,CvOrg_prod,CsSa_prod,CsOrg_prod,gamma_prod,t_prod)