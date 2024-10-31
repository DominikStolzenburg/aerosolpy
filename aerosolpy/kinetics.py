# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 10:41:35 2024

@author: domin
"""
import numpy as np
from aerosolpy.mechanics import AerosolMechanics

class AerosolKinetics(AerosolMechanics):
    """
    base class for kinetics calculations
    
    Parameters
    ----------
    temp_kelvin : float, optional
        temperature in [K], default 296.15 K (23 deg C)
    
    pres_pha : float, optional
        pressure in [hPa],  default 1013.3 hPa
    """
    
    def __init__(self, temp_kelvin=296.15, pres_hpa=1013.3, **kwargs):

        super(AerosolKinetics, self).__init__(temp_kelvin=temp_kelvin,
                                              pres_hpa=pres_hpa,
                                              **kwargs)
   
    def coll_kernel_pp(self, di, dj, rhoi=1000, rhoj=1000):
        """
        calculates the coagulation kernel for particles of different size
        
        Parameters
        ----------
        di : array_like
            vapour molecule equivalent diameter in [nm]
        dj : array_like
            particle diameter of collision partner in [nm]
        rhoi : float, optional
            density of particle in [kg m-3],default 1000 (unit density)
        rhoj : float, optional
            density of particle in [kg m-3], default 1000 (unit density)
        alpha : float
            accomodation coefficient, between 0 and 1, dimless
        
        Returns
        -------
        array_like
            collision kernel in m3 s-1 yielding collision frequency when 
            multiplied by the two concentrations of the two collision partners
        
        Notes
        ----------
        calculations follow chapter 13.3.1.2 pp. 547-549 of [14]_
            
        References
        ----------
        .. [14] J.H. Seinfeld, S.N. Pandis, "Atmospheric Chemistry and Pysics: 
           From Air Pollution to Climate Change", John Wiley & Sons, Inc., 
           New Jersey, 2006
        
        """
        pi = 3.14159
        mi = (pi*rhoi/6.)*(di*1e-9)**3
        mj = (pi*rhoj/6.)*(dj*1e-9)**3
        ci = np.sqrt(8*1.380658e-23*self.temp_kelvin/(pi*mi))   
        cj = np.sqrt(8*1.380658e-23*self.temp_kelvin/(pi*mj)) 
        mfp_i = 8*self.diff_coeff_p(di)/(pi*ci)
        mfp_j = 8*self.diff_coeff_p(dj)/(pi*cj)
        gi = (np.sqrt(2)/(3*di*1e-9*mfp_i) * ( (di*1e-9 + mfp_i)**3 
              - ( (di*1e-9)**2 + mfp_j**2)**(3./2.) ) 
              - (di*1e-9)
              )
        gj = (np.sqrt(2)/(3*dj*1e-9*mfp_i) * ( (dj*1e-9 + mfp_i)**3 
              - ( (dj*1e-9)**2 + mfp_j**2)**(3./2.) ) 
              - (dj*1e-9)
              )
        corr1 = ((di*1e-9 + dj*1e-9) 
                 / (di*1e-9 + dj*1e-9 + 2*(gi**2+gj**2)**0.5)
                 )
        corr2 = ((8*(self.diff_coeff_p(di)+self.diff_coeff_p(dj)))
                 /((ci**2+cj**2)**0.5 * (di*1e-9 + dj*1e-9))
                 )
        corr = 1./(corr1+corr2)
        collKernel = (2*pi
                      * (self.diff_coeff_p(di)+self.diff_coeff_p(dj))
                      * (dj*1e-9+di*1e-9)
                      * corr
                      )
        return collKernel

    def coll_kernel_vp(self, dv, dp, rhov=1000, rhop=1000,
                       alpha=1, 
                       diff_coeff_v=None,
                       dynamic_regime='transition',
                       mv=None): 
        """
        coagulation kernel for vapour molecules with particles
        
        Parameters
        ----------
        dv : array_like
            vapour molecule equivalent diameter in [nm]
        dp : array_like
            particle diameter of collision partner in [nm]
        rhov : float, optional
            density of vapour molecule in [kg m-3], defaul 1000 (unit density)
        rhop : float, optional
            density of particle in [kg m-3], default 1000 (unit density)
        alpha : float, optional
            accomodation coefficient, between 0 and 1, dimless, default 1
        diff_coeff_v : float, optional
            diffusion coefficient of vapour in [m2 s-1], might be calculated 
            differently than diffusion coefficient of particle
        mv : float, optional
            molecular mass [in amu] of vapor, default None (calculate from
            diameter) 
        dynamic_regime : str, optional
            Knudsen number regime for calculations, default transition
            
        
        Returns
        -------
        array_like
            collision kernel in [m3 s-1] yielding collision frequency when 
            multiplied by the two concentrations of the two collision partners
        
        Notes
        -----
        mean free path and Knudsen number are assumed to follow [15]_, which 
        takes into account that vapour and particle size can be of similar
        magnitude, which is different to the classical approach of [16]_
        and [17]_
        
        The collison kernel can also be defined from the free molecular regime 
        collision kernel as basis and a correction factor gamma is used
        This is found by equating the kernels for the free molecular regime 
        with the kernel for the contiuum regime. This is done e.g. in [18]_, 
        but both approaches are equivalent.
                
            
        References
        ----------
        .. [15] K.E.J. Lehtinen, M. Kulmala, "A model for particle formation 
           and growth in the atmosphere with molecular resolution in size", 
           Atmos. CHem. Phys., vol. 3, pp. 251-257
        .. [16] N.A. Fuchs, A.G. Sutugin, "High dispersed aerosols", in Topics
           in Current Aerosol Research, G.M Hindy and J.R.Brock (eds.),
           Pergamon, New York, pp. 1-60, 1971
        .. [17] J.H. Seinfeld, S.N. Pandis, "Atmospheric Chemistry and Pysics: 
           From Air Pollution to Climate Change", John Wiley & Sons, Inc., 
           New Jersey, 2006
        .. [18] T. Nieminen et al., "Sub-10nm particle growth by vapor 
           condensation – effects of vapor molecule size and particle 
           thermal speed", Atmos. Chem. Phys., vol. 10, pp. 9773-9779, 2010

        """
        if mv is None:
            mv = (np.pi*rhov/6.)*(dv*1e-9)**3
        else: 
            mv = mv*1.66e-27
        mp = (np.pi*rhop/6.)*(dp*1e-9)**3
        cv = np.sqrt(8*1.380658e-23*self.temp_kelvin/(np.pi*mv))   
        cp = np.sqrt(8*1.380658e-23*self.temp_kelvin/(np.pi*mp))
        if diff_coeff_v is None:
            diff_coeff_v = self.diff_coeff_p(dv)
        else:
            diff_coeff_v = diff_coeff_v
            
        mfp_v_p = (3*(diff_coeff_v+self.diff_coeff_p(dp))
                   / np.sqrt(cv**2+cp**2)
                   )
        Kn = 2*mfp_v_p/((dv+dp)*1e-9)
        beta_m = ((1+Kn)
                  /(1 + 4/(3*alpha)*Kn + 0.337*Kn + 4/(3*alpha)*Kn**2)
                  )
        if dynamic_regime=='transition':
            return (2*np.pi*(diff_coeff_v+self.diff_coeff_p(dp))
                    *(dv*1e-9+dp*1e-9)*beta_m
                    )
        elif dynamic_regime=='kinetic':
            return np.pi/4. * (dv*1e-9+dp*1e-9)**2 * (cv**2+cp**2)**(1./2.)
        elif dynamic_regime=='continuum':
            return (2*np.pi*(diff_coeff_v+self.diff_coeff_p(dp))
                    *(dv*1e-9+dp*1e-9)
                    )
        else:
            raise ValueError(dynamic_regime,
                             "needs to be transition, kinetic or continuum"
                             )
    
    def condensation_sink(self, dp, dn_dlogdp):
        """
        calculates condensation sink according to [19]_
        
        Parameters
        ----------
        pp : array_like
            diameter discretization of size distribution in [nm]
        dn_dlogdp : array_like
            size distribution in dN/dlogDp
            
        Returns
        -------
        float
            condensation sink (CS) in [s-1]
        np.ndarray
            1-d array with size resolved CS contribution
            
        References
        ----------
        .. [19] Kulmala et al., "Measurement of the nucleation of atmospheric 
           aerosol particles", Nat. Protoc., vol. 7, pp. 1651-1667, 2012
        """
        dp = dp*1e-9

        kn = 2*self.mfp_v()/dp
        alpha = 1 #accommodation coefficient for H2SO4 assumed to  be 1
        beta = (kn+1)/(1 + 0.377*kn + 4/(3*alpha)*kn + 4/(3*alpha)*kn**2)
        
        dintvals = [dp[k]+0.5*(dp[k+1]-dp[k]) for k in range(len(dp)-1)]
        dintval0 = [dp[0]-0.5*(dp[1]-dp[0])]
        dintvalmax = [dp[-1]+0.5*(dp[-1]-dp[-2])]
        dintvals = dintval0+dintvals+dintvalmax
        dintvals_diff_log = np.array([(np.log10(dintvals[k+1])
                                       -np.log10(dintvals[k]))
                                      for k in range(len(dintvals)-1)])
        dn = dn_dlogdp[:]*dintvals_diff_log[:]
        
        dcs = 2*np.pi*self.diff_coeff_v()*dn[:]*beta[:]*dp[:]*1e2
        cs = np.sum(dcs)
        return cs, dcs

    def _vdW_potential(self, r, ri, rj, hamaker):
        """
        calculates van-der-Waals potential between two spheres of sizes 
        ri and rj at distance r
        
        Parameters
        ----------
        r : float 
            distance from between two particles in [m]
        ri : float 
            radius of first entity [m]
        rj : float
            radius of second entity [m]
        hamaker : float
            Hamaker constant in [j]
        
        Returns
        -------
        float
            potential at distance r
        """
        term_1 = (2*ri*rj)/(r**2-(ri+rj)**2)
        term_2 = (2*ri*rj)/(r**2-(ri-rj)**2)
        term_3 = np.log((r**2-(ri+rj)**2)/(r**2-(ri-rj)**2))
        phi = -hamaker/6 * (term_1 + term_2 + term_3) 
        return phi  
    
    def coll_kernel_vp_vdw(self, dv, dp, 
                           rhov=1000, rhop=1000, 
                           hamaker = 5.2e-20,
                           diff_coeff_v=None,
                           method='sceats',
                           dynamic_regime='transition'):
        """
        calculates the collision kernel inculuding collison enhancement 
        according to [20]_ and [21]_ (method='sceats') or alternatively [22]_
        (method='fuchs')
    
        Parameters
        ----------
        dv : float
            diameter of vapor molecule in [nm]
        dp : array_like of float
            diameter of particle in [nm]
        rhov : float, optional
            density of vapor [kg m-3], default 1000
        rhop : float, optional
            density of particle [kg m-3], default 1000
        hamaker : float, optional
            Hamaker constant in [J], default 5.2e-20 (H2SO4)
        diff_coeff_v : float, optional
            diffusivity of vapor [m2 s-1], default None
        method : str, optional
            method for kernel calculation, default 'sceats'
        dynamic_regime : str, optional
            Knudsen number regime for calculations, default transition
    
        Returns
        ----------
        array_like of float
            collision frequency in [m3 s-1]
    
        References
        ----------
        .. [20] M.G. Sceats, "Brownian Coagulation in Aerosols-The Role of Long 
           Range Forces", J. Coll. Interf. Sci., vol. 129, pp. 105-112, 1989
        .. [21] T.W. Chan and M. Mozurkewich, "Measurement of the coagulation 
           rate constant for sulfuric acid particles as a function of particle 
           size using tandem differential mobility analysis", J. Aersol Sci., 
           vol. 32, pp. 321-339, 2001
        .. [22] N.A. Fuchs and A. G. Sutugin, "Coagulation rate of highly 
           dispersed aerosols", J. Coll. Sci., vol 20., pp. 492-500, 1965 
        """
        dv = dv*1e-9
        dp = dp*1e-9
        k = 1.380658e-23 
        mv = (3.14159*rhov/6.)*dv**3
        mp = (3.14159*rhop/6.)*dp**3
        cv = np.sqrt(8*k*self.temp_kelvin/(np.pi*mv)) 
        cp = np.sqrt(8*k*self.temp_kelvin/(np.pi*mp)) 
        if diff_coeff_v is not None:
            diff_coeff_v = diff_coeff_v
        else: # if nothing specified calculate as if vapor is particle
            diff_coeff_v = self.diff_coeff_p(dv)
        diff_coeff_p = self.diff_coeff_p(dp*1e9)        
        
        if method=='sceats':
            a_prime = hamaker/(k*self.temp_kelvin) * (4*dp*dv)/((dp+dv)**2)
            E_inf = (1 
                     + ( np.sqrt(a_prime/3.) / (1+0.0151*np.sqrt(a_prime)) ) 
                     - 0.186*np.log(1+a_prime) - 0.0163*(np.log(1+a_prime)**3)
                     )
            E_0 = 1+0.0757*(np.log(1+a_prime))+0.0015*((np.log(1+a_prime))**3)
            kk = np.pi/4. * ((dp+dv)**2) * np.sqrt(cp**2+cv**2) * E_inf
            kd = 2 * np.pi * (dp+dv) * (diff_coeff_p+diff_coeff_v) * E_0
            kt = kk* (np.sqrt(1+(kk/(2*kd))**2) - (kk/(2*kd)))        
        elif method=='fuchs':
            rv = dv/2.
            rp = dp/2.
            r = np.logspace(np.log10((rv+rp)),-7,10000)   
            phi = self._vdW_potential(r, rv, rp, hamaker)
            b = r*np.sqrt(1+2*np.abs(phi)/(3*1.38e-23*self.temp_kelvin))
            b_crit = np.nanmin(b)
            E_inf = (b_crit/(rv+rp))**2 
            kk = np.pi/4. * ((dp+dv)**2) * np.sqrt(cp**2+cv**2) * E_inf
            a_prime = hamaker/(k*self.temp_kelvin) * (4*dp*dv)/((dp+dv)**2)
            E_0 = 1+0.0757*(np.log(1+a_prime))+0.0015*((np.log(1+a_prime))**3)
            kd = 2 * np.pi * (dp+dv) * (diff_coeff_p+diff_coeff_v) * E_0
            kt = kk* (np.sqrt(1+(kk/(2*kd))**2) - (kk/(2*kd)))   
        else:
            raise ValueError(method, "sceats or fuchs")
        
        if dynamic_regime=='transition':
            return kt
        elif dynamic_regime=='kinetic':
            return kk
        elif dynamic_regime=='continuum':
            return kd
        else:
            raise ValueError(dynamic_regime,"transition, kinetic or continuum")