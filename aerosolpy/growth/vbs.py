# -*- coding: utf-8 -*-

import numpy as np
import bisect
from scipy.integrate import solve_ivp
#import necessary base package functions and classes 
from aerosolpy.kinetics import AerosolKinetics


class VbsModel(AerosolKinetics):
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
    particle_reactions : string or float or 3-tuple
        whether to calculate particle phase reactions or not.
        default is 'None'. 
    particle_diffusion : string or float
        whether to take into account particle-phase diffusion limitations.
        float is diffusion coefficient in [cm2 s-1]
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
    
    particle_reactions and particle_diffusion should not be used simultaneously
    as the solution for the particle-phase diffusion assume a negligible 
    reactivity of the solute. Moreover, activity coefficients are set to 1
    if particle_diffusion is calculated
    """
    def __init__(self, time, vbs_traces, vbs_mass, vbs_logC,
                 sa_trace='None',
                 activity='unity',
                 particle_reactions='None',
                 particle_diffusion='None',
                 temp_kelvin=300,
                 **kwargs):
        
        super().__init__(temp_kelvin=temp_kelvin, **kwargs)
        
        # load time axis, needs to be equal for SA and VBS
        self.time = np.array(time)
        
        # fixed parameters, should be changed by experienced user only
        self.n_p = 1e3 
        self.dp_seed = 0.9

        ## properties of organics
        ## simulated within the VBS scheme
        self.dk = 4.8*(300./temp_kelvin) # Kelvin-diameter in nm
        self.rho_org = 1400 # equal for all VBS bins
        self.log_c0 = np.array(vbs_logC)
        self.c0 = 10**self.log_c0
        self.m_org = np.array(vbs_mass)
        self.dp_org = (((6/np.pi)*self.m_org*1.6605e-27
                        /self.rho_org)**(1/3.)
                       *1e9) # in nm
        
        # load VBS input data
        vbs_traces = np.array(vbs_traces)
        # extend gas-phase data if input is 0d
        if ((vbs_traces.ndim==1) or (vbs_traces.shape[0]==1)):
            self.time = np.arange(0,1000,10)
            vbs_traces = np.repeat(vbs_traces.flatten()[:,np.newaxis],
                                   len(self.time),axis=1).transpose()    
        # convert organic concentration to saturation mass concentration
        self.cv_org = vbs_traces[:,:]*1e6/6.022e23*self.m_org[:]*1e6
        self.n_bins = self.cv_org.shape[1]
        
        #sulfuric acid
        if  isinstance(sa_trace, str):
            if sa_trace=='None':
                self.rho_sa = 1800      # H2SO4+H2O density [kg m-3]
                self.m_sa = 98.        # vapor cluster mass in amu 
                self.dp_sa = (((6./np.pi)*self.m_sa*1.6605e-27
                               /self.rho_sa)**(1/3.)
                              *1e9) # size of SA vapor in nm
                self.diff_coeff_sa = self.diff_coeff_v(mv=98, diff_vol_v=51.96)
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
                           /self.rho_sa)**(1/3.)
                          *1e9) # size of SA vapor in nm
            self.diff_coeff_sa = self.diff_coeff_v(mv=134, diff_vol_v=73.42)
            # seed is made of initial condensed phase SA
            self.cs_sa_init = (np.pi/6.*(self.dp_seed*1e-9)**3
                               *(self.n_p*1e6)*self.rho_sa*1e9)
            self.interact_sa = 1
        
            # load H2SO4 input data 
            sa_trace = np.array(sa_trace)
            # extend input data if input is 0d
            if ((sa_trace.ndim==0) or (sa_trace.shape[0]==1)):
                sa_trace = np.repeat(sa_trace, len(self.time), axis=0)
            # convert SA concentration to saturation mass concentration
            self.cv_sa = sa_trace[:]*1e6/6.022e23*self.m_sa*1e6

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
        
        # particle-phase reaction flag
        self.particle_reactions = particle_reactions
        
        # particle-phase diffusion flag
        self.particle_diffusion = particle_diffusion

        
    
    def _dynVBS(self, t, c_comb):
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
        cv_sa = c_comb[:1] # gaseous sulfuric acid
        cv_org = c_comb[1:self.n_bins+1]
        # solution of the problem is achieved in the logCs space 
        # this avoids negative values. However, init values cannot be 0 
        # but must take a very small number if necessary
        cs_sa = np.exp(c_comb[1+self.n_bins:self.n_bins+2]) 
        #print(cs_sa)
        cs_org = np.exp(c_comb[2+self.n_bins:2*self.n_bins+2])

        # raoult term
        c_oa = np.sum(cs_org)
        c_part = c_oa + self.interact_sa*cs_sa
        
        #calculate mass-weighted organic activity coefficients
        if self.solve_activity==True:
            fc_s = np.ones(self.n_bins)
            for i in range(self.n_bins):
                nc_cs_weighted = 0
                no_cs_weighted = 0
                for j in range(self.n_bins):
                    if j!=i:
                        nc_cs_weighted = nc_cs_weighted+self.n_c[j]*cs_org[j]
                        no_cs_weighted = no_cs_weighted+self.n_o[j]*cs_org[j]
                fc_s[i] = 1./(1+no_cs_weighted/nc_cs_weighted)
            #calculate solute carbon fractions
            fc_i = 1./(1+self.n_o[:]/self.n_c[:])
            # calculate activity coefficients
            gamma_i = np.exp(-2*(self.n_c[:]+self.n_o[:])
                             *((fc_i[:])**2+(fc_s[:])**2-2*fc_i[:]*fc_s[:])
                             *(self.b_co)*(690/self.temp_kelvin)
                             )
            gamma_i = np.nan_to_num(gamma_i,nan=1.0)
        else:
            gamma_i = np.ones(self.n_bins)
        
        # particle volume,diameter and mass from condensed mass concentration
        vp = ((cs_sa/self.rho_sa + c_oa/self.rho_org)
              /(self.n_p*1e6*1e9)
              )
        dp = ((6/np.pi)*vp)**(1/3.) * 1e9
        

        beta_i_p_org = self.coll_kernel_vp(self.dp_org[:], dp,
                                           self.rho_org, self.rho_org)

        beta_i_p_sa = self.coll_kernel_vp_vdw(self.dp_sa, dp,
                                              self.rho_sa, self.rho_sa,
                                              hamaker=5.2e-20,
                                              diff_coeff_v=self.diff_coeff_sa,
                                              method='sceats') 

        #kinetic H2SO4+H2O condensation
        dlogcssa_dt = ( (self.n_p*1e6) * (beta_i_p_sa*60) 
                       * (cv_sa/cs_sa)
                       ) 
        #organic condensation
        dlogcsorg_dt = ( (self.n_p*1e6) * (beta_i_p_org*60) 
                        * ((cv_org[:]/cs_org[:]) 
                           - self.c0[:]*gamma_i[:]*(1./c_part)*10**(self.dk/dp)
                           )
                        )
        
        #modified condensation in case of particle-phase diffusion limitation
        if (self.particle_diffusion!='None'):
            if isinstance(self.particle_diffusion, float):
                diff_limit = (1./
                              (1+((dp*1e-9)**2*self.n_p*1e6
                                  *self.c0[:]*beta_i_p_org*60
                                  /(60*self.particle_diffusion*6e-3*c_part)
                                  )
                               )
                              )
                # activity coefficients are set to be 1
                dlogcsorg_dt = ((self.n_p*1e6) * (beta_i_p_org*60) * diff_limit
                                * ((cv_org[:]/cs_org[:]) 
                                   - self.c0[:]*(1./c_part)*10**(self.dk/dp)
                                   )
                                )
        else:
            dlogcsorg_dt = dlogcsorg_dt
        
        # use biscetion search to find proper derivative
        idx = bisect.bisect_left(self.time, t)

        dcvsa_dt = ((self.cv_sa[idx]-self.cv_sa[idx-1])
                    /(self.time[idx]-self.time[idx-1])
                    )
        dcvorg_dt = ((self.cv_org[idx,:]-self.cv_org[idx-1,:])
                     /(self.time[idx]-self.time[idx-1])
                     )
        
        if self.particle_reactions!='None':
            # model case of gradual breakdown of each bin into the next
            # higher volatility bin
            if isinstance(self.particle_reactions, float):
                L = self.particle_reactions
                P = self.particle_reactions
                for i in range(self.n_bins):
                    if i==0:
                        dlogcsorg_dt[i] = (dlogcsorg_dt[i]
                                           - L)
                    if ((i > 0) and (i!=self.n_bins-1)):
                        dlogcsorg_dt[i] = (dlogcsorg_dt[i]
                                           - L 
                                           + P*cs_org[i-1]/cs_org[i])
                    if i==self.n_bins-1:
                        dlogcsorg_dt[i] = (dlogcsorg_dt[i]
                                           + P*cs_org[i-1]/cs_org[i])
            
            # model case of breakdown of the first (least volatile) bins
            # into the bins shifted by k.
            elif isinstance(self.particle_reactions, tuple):
                L = self.particle_reactions[0]
                P = self.particle_reactions[0]
                k = self.particle_reactions[1]
                for i in range(self.n_bins):   
                    if i<len(L):
                            dlogcsorg_dt[i] = (dlogcsorg_dt[i]
                                               - L[i])
                    elif i>=k and i<len(L)+k:
                            dlogcsorg_dt[i] = (dlogcsorg_dt[i]
                                               + P[i-k]*cs_org[i-k]/cs_org[i])

            else:
                ValueError(self.particle_reactions)
            

        return np.concatenate((np.array([dcvsa_dt]),
                               dcvorg_dt,
                               dlogcssa_dt,
                               dlogcsorg_dt))
    
    def _solve(self, dt=0.05):
        """
        solver for set of differential equations
        
        Parameters
        ----------
        dt : float, optional
            time step [min] for internal solver calculation, default 0.05
        
        Notes
        ----------
        uses scipy.integrate solve_ivp, does not provide non-negativity 
        constraint, which is different to MATLAB version. 
        Uses log-space for solution, and therefore implicitly achieves non-
        neagtivity. 
        """
        # initial condition for differentials, 
        # vapor cv at t=0 and log_cs at t=0
        # while cs_org ideally is 0, log_cs cannot be inf
        # therefore log_cs_org approximated to a tiny value
        c_init = np.concatenate((np.array([self.cv_sa[0]]),
                                 self.cv_org[0,:],
                                 np.log(np.array([self.cs_sa_init])),
                                 -25*np.ones(self.n_bins)
                                 )
                                )

        #print(c_init)
        # time interval for solution
        t0 = self.time[0] 
        tend = self.time[-1]

        sol = solve_ivp(self._dynVBS, [t0, tend], c_init, 
                        t_eval=np.arange(t0, tend, dt),
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
            growth rate per vbs bin (including SA) [nm h-1],
            time of the step [min]
                    
        """
        t_prod, c_prod = self._solve()
        # calculates concentration and diameter evolution from solver solution
        cv_sa_prod = c_prod[:1,:]
        cv_org_prod = c_prod[1:self.n_bins+1,:]

        cs_sa_prod = np.exp(c_prod[self.n_bins+1:self.n_bins+2,:])
        cs_org_prod = np.exp(c_prod[self.n_bins+2:2*self.n_bins+2,:])

        c_oa_prod = np.sum(cs_org_prod,axis=0)
        c_part_prod = c_oa_prod + self.interact_sa*cs_sa_prod
        
        n_p_prod = self.n_p*np.ones(t_prod.shape[0])        
        vp_prod = ((cs_sa_prod[:]/self.rho_sa + c_oa_prod[:]/self.rho_org)
                   /((n_p_prod[:]*1e6)*1e9)
                   )
        dp_prod = ((6/np.pi)*vp_prod)**(1/3.) * 1e9

        if self.solve_activity==True:
            #calculate mass-weighted solvent carbon fraction
            fc_s = np.ones((self.n_bins, t_prod.shape[0]))

            for i in range(self.n_bins):
                nc_cs_weighted = np.zeros(t_prod.shape[0])
                no_cs_weighted = np.zeros(t_prod.shape[0])
                for j in range(self.n_bins):
                    if j!=i:
                        nc_cs_weighted[:] = (nc_cs_weighted[:]
                                             +self.n_c[j]*cs_org_prod[j,:]
                                             )
                        no_cs_weighted[:] = (no_cs_weighted[:]
                                             +self.n_o[j]*cs_org_prod[j,:]
                                             )

                fc_s[i,:] = 1./(1+no_cs_weighted[:]/nc_cs_weighted[:])
            #calculate solute carbon fractions
            fc_i = 1./(1+self.n_o[:]/self.n_c[:])

            # calculate activity coefficients
            gamma_prod = np.exp(-2*(self.n_c[:,np.newaxis]+self.n_o[:,np.newaxis])
                                *((fc_i[:,np.newaxis])**2+(fc_s[:,:])**2
                                  -2*fc_i[:,np.newaxis]*fc_s[:,:]
                                  )
                                *(self.b_co)*(690/self.temp_kelvin)
                                )
        else:
            gamma_prod = np.ones((self.n_bins,t_prod.shape[0]))

        
        beta_i_p_org = self.coll_kernel_vp(self.dp_org[:,np.newaxis],
                                           dp_prod,
                                           self.rho_org,
                                           self.rho_org)

        beta_i_p_sa = self.coll_kernel_vp_vdw(self.dp_sa,
                                              dp_prod,
                                              self.rho_sa,
                                              self.rho_sa,
                                              hamaker=5.2e-20,
                                              diff_coeff_v=self.diff_coeff_sa,
                                              method='sceats') 
        
        #kinetic H2SO4+H2O condensation
        dcssa_dt = ( (self.n_p*1e6) * (beta_i_p_sa[:])
                    * cv_sa_prod[:]
                    )
        ddpdt_sa = ( 2/(np.pi*self.rho_sa*(dp_prod[:]*1e-9)**2*self.n_p*1e6)
                    *dcssa_dt
                    )*3600
        #organic condensation
        if self.particle_diffusion!='None':
            if isinstance(self.particle_diffusion, float):
                diff_limit = (1./
                              (1+((dp_prod[:]*1e-9)**2*self.n_p*1e6
                                  *self.c0[:, np.newaxis]*beta_i_p_org[:,:]*60
                                  /(60*self.particle_diffusion*6e-3*c_part_prod[:])
                                  )
                               )
                              )

        else:
            diff_limit = 1
        
        
        dcsorg_dt = ( (self.n_p*1e6) * (beta_i_p_org[:,:]) * diff_limit
                     * (cv_org_prod[:,:] 
                        -(self.c0[:,np.newaxis] 
                          * gamma_prod[:,:]*cs_org_prod[:,:]/c_part_prod[:]
                          *10**(self.dk/dp_prod[:])
                          )                    
                        )
                     ) 
        
        if self.particle_reactions!='None':
            if isinstance(self.particle_reactions, float):
                L = self.particle_reactions
                P = self.particle_reactions
                for i in range(self.n_bins):
                    if i==0:
                        dcsorg_dt[i] = (dcsorg_dt[i]
                                        - L*cs_org_prod[i])
                    if ((i > 0) and (i!=self.n_bins-1)):
                            dcsorg_dt[i] = (dcsorg_dt[i] 
                                            - L*cs_org_prod[i] 
                                            + P*cs_org_prod[i-1])
                    if i==self.n_bins-1:
                        dcsorg_dt[i] = (dcsorg_dt[i] 
                                        + P*cs_org_prod[i-1])
                        
            elif isinstance(self.particle_reactions, tuple):
                L = self.particle_reactions[0]
                P = self.particle_reactions[0]
                k = self.particle_reactions[1]
                for i in range(self.n_bins):
                    if i<len(L):
                            dcsorg_dt[i] = (dcsorg_dt[i]
                                            - L[i]*cs_org_prod[i])
                    if i>=k and i<len(L)+k:
                            dcsorg_dt[i] = (dcsorg_dt[i]
                                            + P[i-k]*cs_org_prod[i-k])

            else:
                ValueError(self.particle_reactions)
            
        ddpdt_org = (2/(np.pi*self.rho_org*(dp_prod[:]*1e-9)**2*self.n_p*1e6)
                     *dcsorg_dt
                     )*3600
        
        
        dp_prod = np.squeeze(dp_prod)

        gr_bins = np.concatenate((ddpdt_sa,ddpdt_org),axis=0)
        gr_tot = np.sum(gr_bins,axis=0)
        
        return (dp_prod,gr_tot,gr_bins,t_prod)
    
    def growth_rate(self, dp):
        """
        calculates VBS growth rates at diameter dp and set-up conditions
        of the VbsModel
    
        Parameters
        ----------
        dp : array_like of float
            particle mass diameter in [nm], subtract 0.3 from mob. diameter
        
        Returns
        -------
        array_like
            growth rates of nanoparticles from VbsModel in [nm h-1]
        
        Notes
        -----
        Uses numpy.gradient to determine the growth rate from the model solution
        which is given in condensed mass produced at certain times
        """
        t_prod, c_prod = self._solve()
        
        cs_sa_prod = np.exp(c_prod[self.n_bins+1:self.n_bins+2,:])
        cs_org_prod = np.exp(c_prod[self.n_bins+2:2*self.n_bins+2,:])
        c_oa_prod = np.sum(cs_org_prod,axis=0)

        n_p_prod = self.n_p*np.ones(t_prod.shape[0])        
        vp_prod = ((cs_sa_prod[:]/self.rho_sa + c_oa_prod[:]/self.rho_org)
                   /((n_p_prod[:]*1e6)*1e9)
                   )
        dp_prod = ((6/np.pi)*vp_prod)**(1/3.) * 1e9
        dp_prod = np.squeeze(dp_prod)
        dpdt = np.gradient(dp_prod, t_prod/60.)

        idx = (np.abs(dp_prod - dp)).argmin()
        
        return dpdt[idx]