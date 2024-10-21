# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 22:43:53 2024

@author: user_1
"""
import numpy as np

class SimpolVolatility():
    """
    class for calculating volatility through a group-contribution approach,
    the SIMPOL model [29]_
    
    Notes
    -----
    is currently not configured to handle nitrogen related functionality, 
    no full SIMPOL model. 
    
    References
    ----------
    .. [29] Pankow, J.F. and Asher, W.E., "SIMPOL.1: a simple group 
       contribution method for predicting vapor pressures and enthalpies of 
       vaporization of multifunctional organic compounds", Atmos. Chem. Phys.,
       vol. 8, pp. 2773–2796, 2008
    
    """
    def __init__(self):
        pass
        
    def _molecular_mass_organics(self, n_c, n_h, n_o):
        mol_mass = n_c*12.010700 + n_h*1.007940 + n_o*15.999400
        return mol_mass

    def _b_k(self, temp_kelvin, b1, b2, b3, b4):
        return b1/temp_kelvin + b2 + b3*temp_kelvin + b4*np.log(temp_kelvin)

    def _zeroeth(self, temp_kelvin):
        b1 = -4.26938e2
        b2 = 2.89223e-1
        b3 = 4.42057e-3
        b4 = 2.92846e-1
        return self._b_k(temp_kelvin, b1, b2, b3, b4)

    def _carbonnumber(self, n, temp_kelvin):
        b1 = -4.11248e2
        b2 = 8.96919e-1
        b3 = -2.48607e-3
        b4 = 1.40312e-1
        return n*self._b_k(temp_kelvin, b1, b2, b3, b4)

    def _hydroxyl(self, n, temp_kelvin):
        b1 = -7.25373e2
        b2 = 8.26326e-1
        b3 = 2.50957e-3
        b4 = -2.32304e-1
        return n*self._b_k(temp_kelvin, b1, b2, b3, b4)

    def _aldehyde(self, n, temp_kelvin):
        b1 = -7.29501e2
        b2 = 9.86017e-1
        b3 = -2.92664e-3
        b4 = 1.78077e-1
        return n*self._b_k(temp_kelvin, b1, b2, b3, b4)
    
    def _ketone(self, n, temp_kelvin):
        b1 = -1.37456e1
        b2 = 5.23486e-1
        b3 = 5.50298e-4
        b4 = -2.76950e-1
        return n*self._b_k(temp_kelvin, b1, b2, b3, b4)

    def _carboxylicacid(self, n, temp_kelvin):
        b1 = -7.98796e2
        b2 = -1.09436e0
        b3 = 5.24132e-3
        b4 = -2.28040e-1
        return n*self._b_k(temp_kelvin, b1, b2, b3, b4)

    def _hydroperoxide(self, n, temp_kelvin):
        b1 = -2.03387e1
        b2 = -5.48718e0
        b3 = 8.39075e-3
        b4 = 1.07884e-1
        return n*self._b_k(temp_kelvin, b1, b2, b3, b4)
    
    def _nonaromaticring(self, n, temp_kelvin):
        b1 = -8.72770e1
        b2 = 1.78059e0
        b3 = -3.07187e-3
        b4 = -1.04341e-1
        return n*self._b_k(temp_kelvin, b1, b2, b3, b4)

    def _carbonylperoxyacid(self, n, temp_kelvin):
        b1 = -8.38064e2
        b2 = -1.09600e0
        b3 = -4.24385e-4
        b4 = 2.81812e-1
        return n*self._b_k(temp_kelvin, b1, b2, b3, b4)
    
    def _peroxide(self, n, temp_kelvin):
        b1 = 1.50093e2
        b2 = 2.39875e-2
        b3 = -3.37969e-3
        b4 = 1.52789e-2
        return n*self._b_k(temp_kelvin, b1, b2, b3, b4)
    
    def log_c(self, temp_kelvin,
              n_c, n_h, n_o,
              nonaromaticring=0,
              hydroxyl=0,
              aldehyde=0,
              ketone=0,
              carboxylicacid=0,
              hydroperoxide=0,
              carbonylperoxyacid=0,
              peroxide=0):
        """
        Parameters
        ----------
        temp_kelvin : float
            temperature in [K].
        n_c : int
            number of carbon atoms in molecule.
        n_h : int
            number of hydrogen atoms in molecule.
        n_o : int
            number of oxygen atoms in molecule.
        nonaromaticring : int, optional
            number of nonaromaticrings in molecule. The default is 0.
        hydroxyl : int, optional
            number of hydroxyl groups in molecule. The default is 0.
        aldehyde : int, optional
            number of aldehyde groups in molecule. The default is 0.
        ketone : int, optional
            number of ketone groups in molecule. The default is 0.
        carboxylicacid : int, optional
            number of carboxylic acid groups in molecule. The default is 0.
        hydroperoxide : int, optional
            number of hydroperoxide groups in molecule. The default is 0.
        carbonylperoxyacid : int, optional
            number of carbonylperoxyacid groups in molecule. The default is 0.
        peroxide : int, optional
            number of peroxide groups in molecule. The default is 0.

        Returns
        -------
        log10_c : float
            log10 of saturation mass concentration in [ug m-3].

        """
        
        log10_p = (self._zeroeth(temp_kelvin)
                   + self._carbonnumber(n_c, temp_kelvin) 
                   + self._nonaromaticring(nonaromaticring, temp_kelvin) 
                   + self._hydroxyl(hydroxyl, temp_kelvin) 
                   + self._aldehyde(aldehyde, temp_kelvin) 
                   + self._ketone(ketone, temp_kelvin) 
                   + self._carboxylicacid(carboxylicacid, temp_kelvin) 
                   + self._hydroperoxide(hydroperoxide, temp_kelvin) 
                   + self._carbonylperoxyacid(carbonylperoxyacid, temp_kelvin) 
                   + self._peroxide(peroxide, temp_kelvin)
                   )
        
        p = np.power(10,log10_p) # in atm.
        p_pascal = p*101325
        c = (p_pascal*self._molecular_mass_organics(n_c,n_h,n_o)
             /(8.3144598*temp_kelvin)
             *1e6
             )
        log10_c = np.log10(c)
        return log10_c

class TwoDimVolatility():
    """
    volatility models for describing the 2D-volatility basis set, following
    the equations of [30]_
 
    Parameters
    ----------
    nC0 : float, optional
        carbon backbone reference, default 25
    bC : float, optional
        effect on volatility per carbon atom, default 0.475
    bO : float, optional
        effect on volatility per oxygen atom, default 2.3
    bCO : float, optional
        non-linearity parameter, default -0.3
    bN : float, optional
        effect on volatility per nitrogen atom, default 2.5
    model : str, optional
        using a pre-defined model to set nC0, bC, bO, bCO, bN. Unless None,
        makes all other definitions obsolete, default None
 
    References
    ----------
    .. [30] Donahue, N. M., "A two-dimensional volatility basis set: 1. 
       organic-aerosol mixing thermodynamics", Atmos. Chem. Phys., vol. 11,
       pp. 3303-3318, 2011
 
    """
    def __init__(self, n_c0=25, b_c=0.475, b_o=2.3, b_co=-0.3, b_n=2.5,
                 model=None):
        if model is None:
            self.n_c0 = n_c0
            self.b_c = b_c
            self.b_o = b_o
            self.b_co = b_co
            self.b_n = b_n
            self.model = None
        if model=='donahue':
            self.n_c0 = 25
            self.b_c = 0.475
            self.b_o = 2.3
            self.b_co = -0.3
            self.b_n = 2.5
            self.model = 'donahue'
        if model=='stolzenburg':
            self.n_c0 = 25
            self.b_c = 0.475
            self.b_o = 1.4
            self.b_co = -0.3
            self.b_n = 2.5
            self.model = 'stolzenburg'
        if model=='mohr':
            self.n_c0 = 25
            self.b_c = 0.475
            self.b_o = 0.2
            self.b_co = 0.9
            self.b_n = 2.5
            self.model = 'mohr'
        if model=='qiao':
            self.n_c0 = 25
            self.b_c = 0.475
            self.b_n = 2.5
            self.model = 'qiao'
    
    def _molecular_mass_organics(self, n_c, n_h, n_o, n_n):
        mol_mass = n_c*12.010700 + n_h*1.007940 + n_o*15.999400 + n_n*14.00670
        return mol_mass
 
    def log_c_300(self, n_c, n_o, n_h=0, n_n=0, dimers=False):
        """
        calculates the volatility at 300 K according to the set values of nC0,
        bC, bO, bCO, bN and the carbon, oxygen and nitrogen number of a 
        specific molecule
     
        Parameters
        ----------
        nC : array-like of int
            number of carbon atoms in molecule
        nO : array-like of int
            number of oxygen atoms in molecule
        nN : array-like of int, optional
            number of nitrogen atoms in molecule, default 0
        nH : array-like of int, optional
            number of hydrogen atoms in molecule, default 0
        dimers : bool, optional
            if dimer/monomer differntiation should be done. Default is False.
        
        Returns
        ----------
        array_like of float
            log10 c0 [mug m-3] of molecule at 300K
        
        Notes
        -----
        provides different pre-defined parameters to be used for different
        systems
        """     
        if (self.model=='stolzenburg') and (dimers==True): 
            # monomer/dimer differentiation for nC>10
            n_c = np.atleast_1d(n_c)
            n_o = np.atleast_1d(n_o)
            if np.isscalar(n_n):
                n_n = np.ones(len(n_c))*n_n
            else:
                n_n = np.atleast_1d(n_n)
            log_c_300 = np.zeros(len(np.atleast_1d(n_c)))
            dim = (n_c>10)
            mon = (n_c<=10)
            log_c_300[dim] = ((self.n_c0-n_c[dim])*self.b_c
                              -(n_o[dim]-3*n_n[dim])*(self.b_o-0.23)
                              -2*(n_c[dim]*(n_o[dim]-3*n_n[dim])
                                  /(n_c[dim]+n_o[dim]-3*n_n[dim])
                                  )*self.b_co
                              -n_n[dim]*self.b_n
                           )
            log_c_300[mon] = ((self.n_c0-n_c[mon])*self.b_c
                              -(n_o[mon]-3*n_n[mon])*self.b_o
                              -2*(n_c[mon]*(n_o[mon]-3*n_n[mon])
                                  /(n_c[mon]+n_o[mon]-3*n_n[mon])
                                  )*self.b_co
                              -n_n[mon]*self.b_n
                              )
        elif self.model=='qiao':
            h_to_c = n_h/n_c
            o_to_c = n_o/n_c           
            aox = (h_to_c>=-0.2*o_to_c+1.5)
            mox = (h_to_c<-0.2*o_to_c+1.5)

            log_c_300 = np.zeros(len(np.atleast_1d(n_c)))
            log_c_300[aox] = ((self.n_c0-n_c[aox])*self.b_c
                              -(n_o[aox]-3*n_n[aox])*1.4
                              -2*(n_c[aox]*(n_o[aox]-3*n_n[aox])
                                  /(n_c[aox]+n_o[aox]-3*n_n[aox])
                                  )*(-0.3)
                              -n_n[aox]*self.b_n
                              )
            log_c_300[mox] = ((self.n_c0-n_c[mox])*self.b_c
                              -(n_o[mox]-2*n_n[mox])*2.3
                              )
        else: #all other volatility models: no monomer/dimer differentiation        
            log_c_300 =  ((self.n_c0-n_c)*self.b_c
                          -(n_o-3*n_n)*self.b_o
                          -2*(n_c*(n_o-3*n_n)/(n_c+n_o-3*n_n))*self.b_co
                          -n_n*self.b_n
                          )
        return log_c_300
 
    def _deltaHvap(self,log_c_300,a=-5.7,b=129):
        """
        calculates the enthalpy of vaporization of a organic molecule assuming
        its volatility at 300 K is known
     
        Parameters
        ----------
        logC300 : array_like of float
            reference volatility in log10 c [mug m-3] at 300 K
        a : float, optional
            slope of deltaHvap-logC300 relation, default -5.7
        b : float, optional
            offset at 1 mug m-3 of the deltaHvap-volatility relation,
            default 129
     
        Returns
        ----------
        array_like of float
            deltaHvap in [kJ mol-1 K-1]
        """
        return a*log_c_300+b
        
 
    def log_c_300_to_log_c_temp(self, temp_kelvin, log_c_300,
                                a=-5.7, b=129):
        """
        calculates the volatility at a given temperature T assuming log c is
        known at 300 K
     
        Parameters
        ----------
        temp_kelvin : float
            temperature in [K]
        logC300 : array_like of float
            reference volatility in log10 c [mug m-3] at 300 K
        a : float, optional
            slope of deltaHvap-logC300 relation, default -5.7
        b : float, optional
            offset at 1 mug m-3 of the deltaHvap-volatility relation,
            default 129
     
        Returns
        ----------
        array_like of float
            log10 c0 [mug m-3] of molecule at specified T
     
        """
        dHvap = self._deltaHvap(log_c_300, a=a, b=b)
        return (log_c_300 
                + 1000*dHvap/(8.31446*np.log(10)) * (1./300 - 1./temp_kelvin)
                )
    
    def log_c(self, temp_kelvin,
              n_c, n_o, n_h=0, n_n=0, 
              a=-5.7, b=129, dimers=False
              ):
        """
        calculates the volatility at a given temperature for known molecular 
        composition
        
        Parameters
        ----------
        temp_kelvin : array_like of float
            temperature in [K]
        n_c : array like of int
            number of carbon atoms in molecule.
        n_o : array_like of int
            number of oxygen atoms in molecule.
        n_h : array_like of int, optional
            number of hydrogen atoms in molecule. Default is 0. 
        n_n : Array_like of int, optional
            number of nitrogen atoms in molecule. Default is 0. 
        a : float, optional
            slope of deltaHvap-logC300 relation, default -5.7
        b : float, optional
            offset at 1 mug m-3 of the deltaHvap-volatility relation,
            default 129
        dimers : bool, optional
            if dimer/monomer differntiation should be done. Default is False.
        
        Returns
        -------
        array_like of float
            log10 c0 [mug m-3] of molecule at temp_kelvin
            
        """
        log_c_300 = self.log_c_300(n_c, n_o, 
                                   n_n=n_n, n_h=n_h, dimers=dimers)
        
        log_c_t = self.log_c_300_to_log_c_temp(temp_kelvin, log_c_300,
                                               a=a, b=b)
        return log_c_t
    
    def log_n(self, temp_kelvin,
              n_c, n_o, n_h, n_n, 
              a=-5.7, b=129, dimers=False
              ):
        """
        Parameters
        ----------
        temp_kelvin : array_like of float
            temperature in [K]
        n_c : array like of int
            number of carbon atoms in molecule.
        n_o : array_like of int
            number of oxygen atoms in molecule.
        n_h : array_like of int, optional
            number of hydrogen atoms in molecule. Default is 0. 
        n_n : Array_like of int, optional
            number of nitrogen atoms in molecule. Default is 0. 
        a : float, optional
            slope of deltaHvap-logC300 relation, default -5.7
        b : float, optional
            offset at 1 mug m-3 of the deltaHvap-volatility relation,
            default 129
        dimers : bool, optional
            if dimer/monomer differntiation should be done. Default is False.
        
        Returns
        -------
        array_like of float
            log10 n0 [cm-3] of molecule at temp_kelvin
        """
        log_c_t = self.log_c(temp_kelvin, n_c, n_o, n_h=n_h, n_n=n_n, 
                             a=a, b=b, dimers=dimers)
        n_t = (10**(log_c_t)*1e-12
               /(self._molecular_mass_organics(n_c, n_h, n_o, n_n)*1.6605e-24))
        return np.log10(n_t)
        
    
    
    