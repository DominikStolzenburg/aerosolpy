# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 20:47:21 2024

@author: domin
"""
import numpy as np
from scipy.misc import derivative
#import necessary base package functions and classes 
from aerosolpy.kinetics import AerosolKinetics

class SulfuricAcid(AerosolKinetics):
    """
    calculates kinetic growth from H2SO4, H2O, NH3 following [26]_
    
    Parameters
    ----------
    temp_kelvin : float, optional
        temperature in [k], default 293.15 K
    rh : float, optional 
        relative humidity in percent, default 60
    mv_dry : float
        molecular mass of dry vapor [g mol-1], default 98 (H2SO4 monomer)
    mv_wet : float
        molecular mass of hydrated vapor [g mol-1], default 134
    
    Raises
    ----------
    aeropy.InputError
        rh needs to be between 0 and 100
    
    References
    ----------
    .. [26] Stolzenburg, D., et al., "Enhanced growth rate of atmospheric 
       particles from sulfuric acid", Atmos. Chem. Phys., vol. 20, 
       pp. 7359–7372, 2020
    """
    def __init__(self, temp_kelvin=293.15, rh=60, mv_dry=98., mv_wet=134.,
                 **kwargs):
        if rh>=0 and rh<=100:
            self.rh = rh
        else:
            raise ValueError(rh,'rh needs to be between 0 and 100')
        self.mv_dry = mv_dry
        self.mv_wet = mv_wet
        super(SulfuricAcid,self).__init__(temp_kelvin=temp_kelvin, **kwargs)
    
    def rho_h2so4(self, wt): 
        """
        calculates the density of a sulfuric acid water mixture
        
        Parameters
        ----------
        wt : float
            mass fraction of sulfuric acid

        Returns
        ----------
        float
            density in [kg m-3]
    
        Notes
        ----------
        Coefficients rho[j][i] and calculation taken from [27]_
    
        References
        ----------
        .. [27] C.E.L. Myhre, "Density and Surface Tension of Aqueous H2SO4 at 
           Low Temperature", J. Chem. Eng. Data, vol. 43, pp. 617-622, 1998
        """
        rho = [[999.8426,547.2659,526.2950e1,-621.3958e2,409.0293e3,
                -159.6989e4,385.7411e4,-580.8064e4,530.1976e4,
                -268.2616e4,576.4288e3],
               [334.5402e-4,-530.0445e-2,372.0445e-1,-287.7670,
                127.0854e1,-306.2836e1,408.3714e1,-284.4401e1,
                809.1053,0,0],
               [-569.1304e-5,118.7671e-4,120.1909e-3,-406.4638e-3,
                326.9710e-3,136.6499e-3,-192.7785e-3,0,0,0,0],
               [0,599.0008e-6,-414.8594e-5,111.9488e-4,-137.7435e-4,
                637.3031e-5,0,0,0,0,0],
               [0,0,119.7973e-7,360.7768e-7,-263.3585e-7,0,0,0,0,0,0]
               ]
        summe = 0
        for i in [0,1,2,3,4,5,6,7,8,9,10]:
            for j in [0,1,2,3,4]:
                summe = summe+rho[j][i]*(wt**i)*((self.temp_kelvin-273.15)**j)
        rho = summe
        return rho

    def wt(self, dp, model='mabnag'):
        """
        sulfuric acid mass fraction derived from diameter, relative humidity 
        and temperature
    
        Parameters
        ----------
        dp : array_like of float
            diameter in [nm]
        model : str
            underlying model used to calculate the mass fraction, default
            mabnag
        
        Returns
        ----------
        array_like of float
            sulfuric acid mass fraction in water solution, dimless, 
            between 0 and 1
    
        Notes
        ----------
        empirical fits to the SAWNUC or MABNAG results, only vaild between 1 
        and 10 nm function only defined for either 38 or 60 percent RH 
        and 278.15 and 293.15 K and at 5 percent RH for 293 K 
        (dry sheath flow conditions)
        
        Raises
        ----------
        ValueError
            when (rh, temp_kelvin) tuple of ap.growth.SulfuricAcid instance is 
            not pre-defined through mabnag or sawnuc calculations
        
        """
        if model=='mabnag':
            if self.rh==60 and self.temp_kelvin==278.15:
                wt = 0.6352-1.662*np.exp(-0.7394*dp)-8.404e-4*dp+0.2701/dp
            elif self.rh==60 and self.temp_kelvin==293.15:
                wt = 0.6373-1.642*np.exp(-0.7413*dp)-8.805e-4*dp+0.2618/dp
            elif self.rh==38 and self.temp_kelvin==278.15:
                wt = 0.8286-1.422*np.exp(-0.7040*dp)-6.733e-4*dp+0.0691/dp
            else:
                raise ValueError((self.rh,self.temp_kelvin),
                                 "(rh,temp_kelvin) not specified for sawnuc")
        
        elif model=='sawnuc':
            if self.rh==60 and self.temp_kelvin==293.15:
                wt = (1.5431e-02/(dp**1.8096e+01)
                      +3.6681e-01/(dp**4.4980e-01)
                      +3.0126e-01)
            elif self.rh==38 and self.temp_kelvin==293.15:
                wt = (2.6322e-02/(dp**1.2910e+01)
                      +2.5915e-01/(dp**5.3424e-01)
                      +4.3189e-01)
            elif self.rh==38 and self.temp_kelvin==278.15:
                wt = (2.6913e-01/(dp**5.1937e-01)
                      +1.1019e-02/(dp**1.4906e+01)
                      +4.1835e-01)
            elif self.rh==60 and self.temp_kelvin==278.15:
                wt = (3.7474e-01/(dp**4.4721e-01)
                      +1.5589e-02/(dp**1.3276e+01)
                      +2.9149e-01)
            elif self.rh==50 and self.temp_kelvin==285.65:
                wt = (3.1553e-01/(dp**4.8445e-01)
                      +2.0990e-02/(dp**1.3017e+01)
                      +3.6306e-01)
            else:
                raise ValueError((self.rh,self.temp_kelvin),
                                 "(rh,temp_kelvin) not specified for sawnuc")
        else:
            raise ValueError(model,
                             "model needs to be mabnag or sawnuc")
        return wt
    
    def gf(self, dp, model='mabnag'):
        """
        calculates sulfuric acid hygroscopic growth factor from models
    
        Parameters
        ----------
        dp : array_like of float
            diameter in [nm]
        model : str, optional
            model to be used for hygroscopicity, defaul mabnag
            
        Returns
        ----------
        array_like of float
            sulfuric acid hygroscopic growth factor
        """
        if model=='sawnuc':
            # model runs at 5% RH and 293.15 K, instrument conditions (dry DMA)
            wt_instr = (2.7735e-02/(dp**1.4384e+00)
                        +7.9929e-02/(dp**5.0060e-01)
                        +6.7605e-01)
            rho_instr = self.rho_h2so4(wt_instr)
        
        elif model=='mabnag':
            # model runs at 5% RH and 293.15 K, instrument conditions (dry DMA)
            wt_instr = 0.9996-1.223*np.exp(-1.1115*dp)-3.490e-5*dp-0.1268/dp
            rho_instr = self.rho_h2so4(wt_instr)
        
        else:
            raise ValueError(model, "model needs to be mabnag or sawnuc")   
        
        gf = (((wt_instr*rho_instr)
               /(self.wt(dp, model=model)
                 *self.rho_h2so4(self.wt(dp, model=model))
                 )
               )**(1/3.)
              )
        return gf
    
    def dln_rho_wt_ddp(self, dp, model='mabnag'):
        """
        calculates the size derivative of the product rho wt
    
        Parameters
        ----------
        dp : array_like of float
            diameter in [nm]
        model : str, optional
            model to be used for hygroscopicity, defaul mabnag
            
        Returns
        ----------
        array_like of float
            size derivative of log(w rho)
        """
        res = derivative(
                lambda x: np.log(self.wt(x,model=model)
                                 *self.rho_h2so4(self.wt(x,model=model))
                                 ),
            dp,dx=1e-6
                    )
        return res
    
    def diff_coeff_h2so4(self):
        """
        calculates the diffusion ceofficient of hydrated gas phase 
        sulfuric acid
    
        Returns
        ----------
        array_like of float
            diffusion coefficient in [m2 s-1]
        
        Notes
        ----------
        calculation according to [28]_
    
        References
        ----------
        .. [28] D.R.Hanson, F.Eisele, "Diffusion of H2SO4 in Humidified 
           Nitrogen: Hydrated H2SO4", J. Phys. Chem. A, vol. 104, pp. 
           1715-1719, 2000
        """
        p =  0.969/1.013 #could be changed to self.P/101.3
        pD0 = 0.094
        pD1 = 0.85*pD0
        pD2 = 0.76*pD0
        K1 = 0.13
        K2 = 0.016 
        Dv = ((1/p)
              * (self.temp_kelvin/298.)**1.75
              * (pD0+pD1*K1*(self.rh/100)+pD2*K1*K2*(self.rh/100)**2)
              / (1+K1*(self.rh/100)+K1*K2*(self.rh/100)**2)
              )
        Dv = Dv * 1e-4
        return Dv
    
    def growth_rate(self, dp, cv,
                    hamaker=5.2e-20,kernel='hard sphere',hydration='naive'):
        """
        calculates sulfuric acid growth rates at diameter dp and vapor 
        concentration Cp
    
        Parameters
        ----------
        dp : array_like of float
            particle mass diameter in [nm], subtract 0.3 from mob. diameter
        Cv : array_like of float
            vapor concentration in [cm-3]
        Hamaker : float
            Hamaker Constant
        kernel : str
            collision kernel formulation
        hydration : str
            treatment of H2SO4 hydration

        Returns
        ----------
        array_like
            growth rates of nanoparticles at the kinetic limit in [nm h-1]
    
        Raises
        ----------
        ValueError
            kernel and hydration only allow specific key words
        
        """    
        amu = 1.6605e-27
        
        # hydration approach set which diameters to use for collision frequency
        # calculations: 
        # 'naive' same hydration of clusters and particles
        # 'dry measurement' hygroscopic gf of 1.25 for all colliding particles
        # 'wet measurement' all particles are measured at their colliding dp
        # 'mabnag' or 'sawnuc' model water content of measured and colliding
        if hydration=='naive':
            rhov_coll = self.rho_h2so4(self.mv_dry/self.mv_wet)
            rhop_coll = rhov_coll
            dv_coll = (6.*self.mv_wet*amu/(3.14159*rhov_coll))**(1/3.)
            dp_coll = dp
            Vv_grow = (self.mv_wet*amu/rhov_coll)
        elif hydration=='dry measurement':
            rhov_coll = self.rho_h2so4(self.mv_dry/self.mv_wet)
            rhop_coll = self.rho_h2so4(0.62)
            dv_coll = (6.*self.mv_wet*amu/(3.14159*rhov_coll))**(1/3.)
            dp_coll = dp*1.25
            Vv_grow = (self.mv_dry*amu/self.rho_h2so4(1.0))
        elif hydration=='wet measurement':
            rhov_coll = self.rho_h2so4(self.mv_dry/self.mv_wet)
            rhop_coll = self.rho_h2so4(0.62)
            dv_coll = (6.*self.mv_wet*amu/(3.14159*rhov_coll))**(1/3.)
            dp_coll = dp
            Vv_grow = (self.mv_wet*amu/rhov_coll)
        elif hydration=='mabnag' or hydration=='sawnuc':
            rhov_coll = self.rho_h2so4(self.mv_dry/self.mv_wet)
            rhop_coll = self.rho_h2so4(self.wt(dp, model=hydration)) 
            dv_coll = (6.*self.mv_wet*amu/(3.14159*rhov_coll))**(1/3.)  
            dp_coll = dp*self.gf(dp,model=hydration)
        else:
            raise ValueError(hydration,
                             ("needs to be naive, dry measuement,"
                              "wet measurement, mabnag or sawnuc")
                             )

        if kernel=='hard sphere':
            kcoll = self.coll_kernel_vp(dv_coll*1e9, dp_coll, 
                                        rhov_coll, rhop_coll,
                                        diff_coeff_v=self.diff_coeff_h2so4()
                                        )

        elif ((kernel=='sceats') or (kernel=='fuchs')):
            kcoll = self.coll_kernel_vp_vdw(dv_coll*1e9, dp_coll,
                                            rhov_coll, rhop_coll,
                                            hamaker,
                                            diff_coeff_v=self.diff_coeff_h2so4()
                                            )
    
        else:
            raise ValueError(kernel, "needs to be hard sphere, sceats or fuchs")
    
        if (hydration=='naive' 
            or hydration=='dry measurement'
            or hydration=='wet measurement'):
            gr = (kcoll * Vv_grow * (cv*1e6) 
                  /((3.14159/2)*(dp*1e-9)**2)
                  )
            gr = gr * 1e9 * 3600
        if hydration=='mabnag' or hydration=='sawnuc':
            gr = (kcoll * (self.mv_dry*amu) * (cv*1e6)
                  /(self.rho_h2so4(self.wt(dp,model=hydration))
                    *1e-27
                    *(self.wt(dp,model=hydration))
                    *(3.14159/2)*(dp_coll)**2
                    *(1+((dp_coll)/3
                         *self.dln_rho_wt_ddp(dp,model=hydration))
                      )
                    )
                  )
            gr = gr*3600
        return gr