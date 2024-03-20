# -*- coding: utf-8 -*-

class Dma(AerosolMechanics):
    """
    a class for differential mobility analyzers
    
    Parameters
    ----------
    Q_a : float
        aerosol flow rate in [lpm]
    Q_sh : float
        sheath flow rate in [lpm]
    L : float
        classification length in [m]
    R_i : float
        inner electectrode radius in [m]
    R_a : float
        outer electrode radius in [m]
    f_sigma : float, optional
        non-ideal instrument broadening of transfer function, default 1
    penetration : float or array_like or callable, optional
        size-dependent penetration efficiency of the DMA, can be given either 
        as a diffusional length (single float), array (from measurement) or 
        callable function, default None, which results in penetration of 1 for
        all sizes
    
    
    Attributes
    ----------
    Q_a : float
        aerosol flow rate in [lpm]
    Q_sh : float
        sheath flow rate in [lpm]
    L : float
        classification length in [m]
    R_i : float
        inner electectrode radius in [m]
    R_a : float
        outer electrode radius in [m]
    f_sigma : float
        non-ideal instrument broadening of transfer function, default 1
    _penetration_func : callable, hidden
        function calculating the penetration efficiency for given diameter
    
    Raises
    ----------
    TypeError
        If the penetration parameter is neither float, nor numpy.ndarray,
        pandas.DataFrame, nor callable
    
    Notes
    ----------
    an instance of this class is a DMA with specified dimensions assuming 
    balanced flows 
    the class is a base class for MobilitySizer
    DMA theory is inferred according to [1]_
    
    References
    ----------
    .. [1] M.R. Stolzenburg, P.H. McMurry, "Equations Governing Single and 
       Tandem DMA Configurations and a New Lognormal Approximation to 
       the Transfer Function", Aerosol Sci. Tech., vol. 42, iss. 6, pp.
       421-432, 2008
    """
    def __init__(self,
                 Q_a,Q_sh,L,R_i,R_a,*args,
                 f_sigma=1,penetration=None,**kwargs):
        self.Q_a = Q_a/60000. #transform in m3 s-1
        self.Q_sh = Q_sh/60000. #transform in m3 s-1
        self.L = L
        self.R_i = R_i
        self.R_a = R_a
        self.f_sigma = f_sigma
        self.beta = self.Q_a/self.Q_sh
        self.g = (self.R_i/self.R_a)**2
        self.kappa = (self.L*self.R_a)/(self.R_a**2-self.R_i**2)
        self.I_gamma = ( ((1./4.*(1-self.g**2)*(1-self.g)**2) 
                          +(5./18.*(1-self.g**3)*(1-self.g)*np.log(self.g)) 
                          +(1./12.*(1-self.g**4)*(np.log(self.g))**2)
                          ) 
                        /((1-self.g)
                          *(-0.5*(1+self.g)*np.log(self.g)-(1-self.g))**2
                          )
                        )
        self.G_DMA = (4*((1+self.beta)**2/(1-self.g)) 
                      * (self.I_gamma+(1./(2*(1+self.beta)*self.kappa)**2))
                      )
        #input type dependent treatment of optional argument penetration
        if penetration is None:
            self._penetration_func = lambda x: 1.0
        elif np.isscalar(penetration): 
            self._penetration_func = (lambda x: 
                                      self.diff_loss(x,penetration/self.Q_a)
                                      )
        elif isinstance( penetration, np.ndarray):
            self._penetration_func  = scipy.interpolate.interp1d(
                                      penetration[:,0],
                                      penetration[:,1]
                                      )
        elif isinstance( penetration, pd.DataFrame):
            self._penetration_func  = scipy.interpolate.interp1d(
                                      penetration.iloc[:,0].tolist(),
                                      penetration.iloc[:,1].tolist()
                                      )
        elif callable( penetration):
            self._penetration_func = penetration
        else:
            raise TypeError(penetration,
                            "input parameter penetration must be scalar, "
                            "callable, numpy.ndarray or pandas.DataFrame"
                            )
        
        super(Dma,self).__init__(*args,**kwargs)

    def VtoZ(self,V):
        """
        calculates electrical mobility from set volatage at DMA
        
        Parameters
        ----------
        V : array_like
            volatage in [V] at DMA
        
        Returns
        ----------
        array_like
            electrical mobility in [m2 V-1 s-1]
        """
        Z = self.Q_sh*np.log(self.R_a/self.R_i)/(2*3.1415926*self.L*V)
        return Z
    
    def ZtoV(self,Z):
        """
        calculates voltage at DMA from electrical mobility
        
        Parameters
        ----------
        Z : array_like of float
            electrical mobility in [m2 V-1 s-1]
        
        Returns
        ----------
        array_like of float
            voltage in [V]
        """
        V = self.Q_sh*np.log(self.R_a/self.R_i)/(2*3.1415926*self.L*Z)
        return V

    def VtoD(self,V,i=1):
        """
        calculates particle diameter from set voltage at DMA
        
        Parameters
        ----------
        V : array_like
            volatage in [V] at DMA
        i : int, optional
            charging state of particle
        
        Returns
        ----------
        array_like
            particle diameter in [nm]
        
        """
        if not isinstance(i, int):
            raise TypeError(i, "charging state must be int")
        Z = self.VtoZ(V)
        D = self.ZtoD(Z,i=i)
        return D 

    def VtoD_approx(self,
                    V): 
        """
        calculates particle diameter from set voltage at DMA using the
        appoximation from [1]_ 
        
        Parameters
        ----------
        V : array_like
            volatage in [V] at DMA
        
        Returns
        ----------
        array_like
            particle diameter in [nm]
        
        Notes
        ----------
        Only valid for dp<5-10 nm and singly charged particles
        
        See also:
        ----------
        aeropy.AerosolMechanics.ZtoD_approx
        
        References:
        ----------
        .. [1] J.M. Maekelae et al., "Comparison of mobility equivalent 
           diameter with Kelvin-Thomson diameter using ion mobility data",
           J. Chem. Phys., vol. 105, pp.1562, 1996
        """
        Z = self.VtoZ(V)
        D = np.power(Z/(2.2458e-22),(1/-1.9956))
        return D*1e9
    
    def DtoV(self,dp,i=1):
        """
        calculates voltage from given diameter
        
        Parameters
        ----------
        dp : array_like of float
            diameter in [nm]
        i : int, optional
            charging state of particle
        
        Returns
        ----------
        array_like
            voltage in [V]
        """
        if not isinstance(i, int):
            raise TypeError(i, "charging state must be int")
        Z = self.DtoZ(dp)
        D = self.ZtoV(Z)
        return D 
   
    def _dimless_diff_coeff(self,
                            dp):
        """
        dimensionless diffusion coefficient
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        
        Returns
        ----------
        array_like
            dimensionless diffusion coefficent
        """
        pi = 3.14159
        D_tilde = 2*pi*self.L*self.diff_coeff(dp)/self.Q_sh
        return D_tilde
    
    def _sigma_theo(self,
                    dp):
        """
        computes theroetical width of DMA transferfunction
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        
        Returns
        ----------
        array_like
            dimensionless transferfunction width
        """
        sig = np.sqrt(self.G_DMA*self._dimless_diff_coeff(dp))
        return sig

    def sigma(self,dp):
        """
        computes actual width of DMA transferfunction 
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        
        Returns
        ----------
        array_like
            dimensionless transferfunction width
        
        Notes
        ----------
        takes into account instrument artefacts which lead to additional 
        broadening        
        """
        return self._sigma_theo(dp)*self.f_sigma 
    
    def _epsilon(self,x):
        """
        epsilon helper-function for transfer function integration
        
        Parameters
        ----------
        x : array_like
        
        Returns
        --------
        array_like
        """
        pi = 3.14159
        return (x*scipy.special.erf(x)) + (np.exp(-x**2)/(np.sqrt(pi)))

    def Zdimless_triangular_transferfunction(self,Zp_tilde):
        """
        calculates dimensionless triangular shaped transferfunction in 
        mobility space
        
        Parameters
        ----------
        Zp_tilde : array_like
            dimensionless electrical mobility Zp_tilde = Zp/Zp_prime
        Returns
        ----------
        array_like
            dimensionless transferfunction
        """
        Omega_ztilde = (1/(2*self.beta) 
                        * (np.absolute(Zp_tilde-(1+self.beta)) 
                           + np.absolute(Zp_tilde-(1-self.beta)) 
                           - np.absolute(Zp_tilde-1) 
                           - np.absolute(Zp_tilde-1) 
                           )
                        )
        return Omega_ztilde
    
    def Zdimless_diffusional_transferfunction(self,
                                              Zp_tilde,dp_prime):
        """
        calculates dimensionless diffusional transferfunction in 
        mobility space
        
        Parameters
        ----------
        Zp_tilde : array_like
            dimensionless electrical mobility Zp_tilde = Zp/Zp_prime
        dp_prime : float
            reference centroid diameter in [nm] corresponding to voltage at DMA
        Returns
        ----------
        array_like
            dimensionless transferfunction
        """
        Omega_ztilde = (self.sigma(dp_prime)/(np.sqrt(2)*self.beta) 
                        *(self._epsilon((Zp_tilde-(1+self.beta))
                                        /(np.sqrt(2)*self.sigma(dp_prime))) 
                          + self._epsilon((Zp_tilde-(1-self.beta))
                                          /(np.sqrt(2)*self.sigma(dp_prime))) 
                          - 2*self._epsilon((Zp_tilde-1)
                                            /(np.sqrt(2)*self.sigma(dp_prime))) 
                          )
                        )
        return Omega_ztilde
 
    def D_triangular_transferfunction(self,
                                      dp,dp_prime,i=1):
        """
        calculates triangular shaped transferfunction in diameter space
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        dp_prime : float
            reference centroid diameter in [nm] corresponding to voltage at DMA
        i : int, optional
            charging state of particle
            
        Returns
        ----------
        array_like
            transferfunction
        """
        if not isinstance(i, int):
            raise TypeError(i, "charging state must be int")
        Zp_tilde = self.DtoZ(dp,i=i)/self.DtoZ(dp_prime,i=i)
        return self.Zdimless_triangular_transferfunction(Zp_tilde)

    def D_diffusional_transferfunction(self,
                                       dp,dp_prime,i=1):
        """
        calculates diffusional transferfunction in diameter space
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        dp_prime : float
            reference centroid diameter in [nm] corresponding to voltage at DMA
        i : int, optional
            charging state of particle
            
        Returns
        ----------
        array_like
            transferfunction
        """
        if not isinstance(i, int):
            raise TypeError(i, "charging state must be int")
        Zp_tilde = self.DtoZ(dp,i=i)/self.DtoZ(dp_prime,i=i)
        Omega_d = self.Zdimless_diffusional_transferfunction(Zp_tilde,dp_prime)
        return Omega_d
    
    def D_diffusional_transferfunction_lognorm(self,
                                               dp,dp_prime):
        """
        calculates diffusional transferfunction in diameter space with log-
        normal approximation
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        dp_prime : float
            reference centroid diameter in [nm] corresponding to voltage at DMA
            
        Returns
        ----------
        array_like
            dimensionless transferfunction
        
        Notes
        ----------
        for singly charged particles only
        """
        dp_tilde = dp/dp_prime
        geomean = -(self.sigma(dp_prime)**2)/self.a_function(dp_prime)
        geodev = ((self.beta**2/6 + self.sigma(dp_prime)**2
                   *(1+2*self.sigma(dp_prime)**2)
                   )
                  /(self.a_function(dp_prime)**2)
                  )
        M0 = (self.beta)/self.a_function(dp_prime)
        Omega_d = (M0/(np.sqrt(2*3.14159*geodev))
                   * np.exp(-0.5*((np.log(dp_tilde)-geomean)**2)/geodev)
                   )
        return Omega_d
    
    def calc_transferfunction_limits(self,
                                     dp_prime,range_mult=3):
        """
        calculate left and right limits of transfer-function in diameter-space
        
        Parameters
        ----------
        dp_prime : float
            reference centroid diameter in [nm] corresponding to voltage at DMA
        range_mult : int
            range multiplier on how many sigmas should be included
        
        Returns
        ----------
        (float, float)
            tuple of low and high limits of transferfunction in [nm]
        
        Notes
        ----------
        Assumes triangular shape, which is good if range_mult is large
        """
        Zp_prime = self.DtoZ(dp_prime)
        D_lim = [self.ZtoD((1+range_mult*self.beta)*Zp_prime),
                 self.ZtoD((1-range_mult*self.beta)*Zp_prime)
                 ]
        return D_lim

    def pen_eff(self,dp):
        """
        calculates penetration efficiency of DMA
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        
        Returns
        ----------
        array_like
            penetration probability between 0 and 1 dimensionless
        """
        eta = self._penetration_func(dp)
        return eta
