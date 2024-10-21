# -*- coding: utf-8 -*-

#import necessary base package functions and classes 
from aerosolpy.kinetics import AerosolKinetics

class KineticLimit(AerosolKinetics):
    """
    calculates kinetic growth for arbitrary vapor following [25]_
    
    Parameters
    ----------
    mv : float
        molecular mass of dry vapor [g mol-1]
    rhov : float, optional 
        density of vapor in [kg m-3]   
    temp_kelvin : float, optional
        temperature in [k], default 293.15 K
    
    References
    ----------
    .. [25] Nieminen, T. et al., "Sub-10 nm particle growth by vapor 
       condensation - effects of vapor molecule size and particle thermal 
       speed", Atmos. Chem. Phys., vol. 10, pp. 9773–9779, 2010
    """
    def __init__(self, mv, rhov, temp_kelvin=293.15, **kwargs):
        self.mv = mv
        self.rhov = rhov
        super(KineticLimit,self).__init__(temp_kelvin=temp_kelvin, **kwargs)
        
    def growth_rate(self, dp, cv,
                    diff_coeff_v = None, hamaker=5.2e-20,
                    kernel='hard sphere',
                    dynamic_regime='transition'):
        """
        calculates growth rates at diameter dp and vapor concentration Cv
    
        Parameters
        ----------
        dp : array_like of float
            particle mass diameter in [nm], subtract 0.3 from mob. diameter
        Cv : array_like of float
            vapor concentration in [cm-3]
        Hamaker : float, optional
            Hamaker Constant in [J]
        kernel : str, optional
            collision kernel formulation
        dynamic_regime : str, optional
            Knudsen number regime for calculations, default transition

        Returns
        ----------
        array_like
            growth rates of nanoparticles at the kinetic limit in [nm h-1]
    
        Raises
        ----------
        Value
            kernel and dynamic regime only allows specific key words
        
        """    
        amu = 1.6605e-27
        dv = (6.*self.mv*amu/(3.14159*self.rhov))**(1/3.)
        Vv = (self.mv*amu/self.rhov)
        
        if kernel=='hard sphere':
            kcoll = self.coll_kernel_vp(dv*1e9, dp, self.rhov, self.rhov,
                                        diff_coeff_v=diff_coeff_v, 
                                        dynamic_regime=dynamic_regime
                                        )

        elif ((kernel=='sceats') or (kernel=='fuchs')):
            kcoll = self.coll_kernel_vp_vdw(dv*1e9, dp, self.rhov, self.rhov,
                                            hamaker,
                                            diff_coeff_v=diff_coeff_v, 
                                            dynamic_regime=dynamic_regime
                                            )
    
        else:
            raise ValueError(kernel, "needs to be hard sphere, sceats or fuchs")
    
       
        gr = (kcoll * Vv * (cv*1e6) 
              /((3.14159/2)*(dp*1e-9)**2)
              )
        gr = gr * 1e9 * 3600
        return gr
    

