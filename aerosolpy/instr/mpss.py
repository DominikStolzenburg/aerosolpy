# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import quad

from aerosolpy.mechanics import AerosolMechanics
from aerosolpy.instr.dma import Dma
from aerosolpy.instr.dma import DmaCylindrical
from aerosolpy.instr.cpc import Cpc

class Mpss(AerosolMechanics):
    """
    a class for mobility particle size spectrometers (DMPS/SMPS)
    providing data inversion basics such as inversion kernels
    
    Parameters
    ----------
    channels : array_like
        array of sampled sizes assuming singly charged particles
        if only voltages are known, use aerosolpy.instr.Dma.v_to_dp for
        coversion
    dma : aerosolpy.instr.Dma 
        DMA to be used within the Mpss
    cpc : aerosolpy.instr.Cpc,  optional
        CPC to be used within the Mpss, default None
    inlet_loss : float or array_like or callable, optional
        losses in the inlet lines of the Mpss, default None
        if float it is line length in [m] divided by flow in [m-3 s]
        if 2d (np.array or pd.DateFrame) first index gives size, second loss
    polarity : str, optional
        polarity of the particles sampled with the DMA
    volts_channels : bool, optional
        if channels are given in volts, induces automatic conversion using 
        aerosolpy.instr.dma.v_to_dp
    
    Raises:
    -------
    TypeError
        channels need to be of type int or array_like
        dma needs to be of type aerosolpy.instr.Dma
    
    See also:
    ---------
    aerosolpy.instr.Dma 
    """
    
    def __init__(self, channels, dma,
                 cpc=None, inlet_loss=None, polarity='neg', volts_channels=False,
                 **kwargs):

        if (isinstance(dma, Dma) or isinstance(dma, DmaCylindrical)):          
            self.mpss_dma = dma
        else:
            raise TypeError(dma, "dma needs to be of type aerosolpy.instr.Dma")
            
        if cpc is None:
            self.mpss_cpc = Cpc(activation=None)
        elif isinstance(cpc, Cpc):
            self.mpss_cpc = cpc
        else:
            raise TypeError(cpc, "cpc needs to be of type aerosolpy.instr.Cpc")
        
        #input type dependent treatment of optional argument penetration
        if inlet_loss is None:
            # use x*0+1.0 to appropriately shape the output. 
            self._inletloss_func = lambda x: x*0+1.0
        elif np.isscalar(inlet_loss): 
            self._inletloss_func = lambda x: self.diff_loss(x,inlet_loss)
                                    
        elif isinstance(inlet_loss, np.ndarray):
            self._inletloss_func  = interp1d(inlet_loss[:,0],
                                              inlet_loss[:,1])
        elif isinstance(inlet_loss, pd.DataFrame):
            self._inletloss_func  = interp1d(inlet_loss.iloc[:,0].tolist(),
                                             inlet_loss.iloc[:,1].tolist()
                                             )
        elif callable(inlet_loss):
            self._inletloss_func = inlet_loss
        else:
            raise TypeError(inlet_loss,
                            "inlet_loss must be scalar, callable, numpy.ndarray"
                            "or pandas.DataFrame"
                            )
        
        self.polarity = polarity       
        
        if volts_channels==True:
            self.channels = np.array(self.mpss_dma.v_to_dp(channels, i=1))
        else:
            self.channels = np.array(channels)
        
        #super takes arguments in order as given in the class definition
        super(Mpss,self).__init__(**kwargs)
       
    def tot_eff(self, dp):
        """
        calculates total instrument efficiency
            
        Parameters
        ----------
        dp : array_like of float
            particle diameter in [nm]            
            
        Returns
        -------
        array_like of float
            total detection efficiency, dimless, between 0 and 1
            
        See also
        --------
        aerosolpy.instr.Cpc.count_eff
        aerosolpy.instr.Dma.pen_eff
        """
        eta = (self.mpss_cpc.count_eff(dp)
               *self.mpss_dma.pen_eff(dp)
               *self._inletloss_func(dp)
               )
        return eta
        
    def tot_transfunc(self, ch, imax=5):
        """
        defines total transfer function of instrument channel k
            
        Extended Summary
        ----------------
        combines losses, DMA transferfunction and cpc activation 
            
        Parameters
        ----------
        ch : int or float
            channel number or centroid diameter of channel
        imax : int, optional
            maximum charging states to be considered, default 5
                
        Returns
        -------
        callable
            total transfer function of Mobility Sizer including measurement
            time and counter flow
                
        Raises
        ------
        ValueError
            channel specifier ch needs to be either int or float 
            imax needs to be positive int
            polarity needs to be either "pos" or "neg"
            
        Notes
        -----
        The returned function is the kernel function for every instrument
        inversion
        """
        #channel can be integer (channel number) or float (channel diameter)
        if not (isinstance(ch, int)  | isinstance(ch, float)):
            raise ValueError(ch, ("channel number must be positive int or "
                                  "channel diameter must be float")
                             )
        if isinstance(ch, int):
            dp_single = self.channels[ch]
        if isinstance(ch, float):
            dp_single = ch
            ch = np.where(self.channels==dp_single)[0][0]
       
        if not isinstance(imax, int):
            raise ValueError(imax, "max charging state must be int")
        dp_multiple = [dp_single]
        # calculation thresholds for multiple charge corrections
        dp_thresholds = [20.0,70.0,100.0,100.0]
        for n in range(imax-1):
            if dp_single>=dp_thresholds[n]:
                dp_multiple.append(self.zp_to_dp(self.dp_to_zp(dp_single), 
                                                 i=n+2
                                                 )
                                   )
        
        if self.polarity=='neg':
            char_states = list(range(1,len(dp_multiple)+1))
            char_states = [-c for c in char_states]
        elif self.polarity=='pos':
            char_states = list(range(1,len(dp_multiple)+1))
        else:
            raise ValueError(self.polarity,
                             "polarity must be either 'pos' or 'neg'")
        
        def transferfunction_of_channel(x):
            dma_trans = 0
            for n in range(len(dp_multiple)):
                dma_trans = (dma_trans 
                             + (self.mpss_dma.dp_diffusional_transfunc_lognorm(
                                x,dp_multiple[n])
                                *self.charge_prob(x, char_states[n])
                                )
                             )
            tf = dma_trans*self.tot_eff(x) #*self.Q_count*self.t_res[ch]*
            return tf
        
        return transferfunction_of_channel
            
    def c_expected(self, f_dNdlogDp, imax=1, xmin=1, xmax=1000):
        """
        calculates expected concentration in size channel from a known
        size-distribution 
            
        Parameters
        ----------
        f_dNdlogDp : callable
            size_distribution dNdlog/Dp as function of diameter dp with dp in 
            [nm]
        imax : int, optional
            maximum number of charges to be considered, default 1, 
            i.e. consider only singly charged particles
        xmin : float, optional
            minimum diameter [nm] considered for transfunc integration,
            default 1
        xmax : float, optional
            maximum diameter [nm] considered for transfunc integration,
            default 1000
            
        Returns
        -------
        array_like of float
            expected concentration in every size-channel 
        """
        if isinstance(self.mpss_dma, Dma):
            tf_shape = 'unity'
        elif isinstance(self.mpss_dma, DmaCylindrical):
            tf_shape = 'lognorm'
        else:
            TypeError(self.mpss_dma, 
                      ("for caclcualtion of expected signal, "
                       "aerosolpy.instr.Mpss.mpss_dma needs to by of Type "
                       "aerosolpy.instr.Dma or aerosolpy.instr.DmaCylindrical")
                      )
        c = []
        if self.polarity=='neg': pol=-1
        if self.polarity=='pos': pol=1
        if imax==1:
            for ch in range(len(self.channels)):
                d = self.channels[ch]
                integral = quad(lambda x: (f_dNdlogDp(x)
                                    *self.mpss_dma.dp_transfunc(x,d,shape=tf_shape)
                                    *self.charge_prob(x,pol)
                                    *self.tot_eff(x)/x/np.log(10)
                                           ),
                                xmin,
                                xmax
                                )[0] 
                c.append(integral)
        elif imax>=2:
            for ch in range(len(self.channels)):
                # tot_eff is incorporated into tot_transfunc
                integral = quad(lambda x: (f_dNdlogDp(x)
                                    *self.tot_transfunc(ch)(x)/x/np.log(10)
                                           ),
                                xmin,
                                xmax
                                )[0] 
                c.append(integral)
        else:
            raise ValueError(imax, "imax needs to be positive integer")
        return np.array(c)
    
    def n_expected(self, f_dNdlogDp, q_sample, t_res, 
                   imax=1, xmin=1, xmax=1000):
        """
        calculates number of counts in size channel from a known
        size-distribution 
            
        Parameters
        ----------
        f_dNdlogDp : callable
            size_distribution dNdlog/Dp as function of diameter dp with dp in 
            [nm]
        q_sample : float 
            sample flow in [lpm]
        t_res : array_like 
            time resolution of mpss channels
        imax : int, optional
            maximum number of charges to be considered, default 1, 
            i.e. consider only singly charged particles
        xmin : float, optional
            minimum diameter [nm] considered for transfunc integration,
            default 1
        xmax : float, optional
            maximum diameter [nm] considered for transfunc integration,
            default 1000
            
        Returns
        -------
        array_like of float
            expected number of counts in every size-channel for the 
            corresponding measurement time t_res and counter flow q_sample
        """
        q_sample = q_sample*1000/60. # for [cm3 s-1]
        if np.isscalar(t_res):
            t_res = np.array([t_res]*len(self.channels))

        n = (self.c_expected(f_dNdlogDp, imax=imax, xmin=xmin, xmax=xmax)[:]
             *q_sample*t_res[:])
        return n
        
    def std_inv(self, Craw, imax=1):  
        """
        standard point-by-point inversion of a MPSS response 
        using multiple charge correction
            
        Parameters
        ----------
        Craw : np.array
            instrument response as recorded in concentration [cm-3] 
        imax : int
            maximum number of charges to be considered for inversion      
                
        Returns
        -------
        tuple: (array_like of float, array_like of float)
            inverted size-distribution, dp in first column, dNdlog/Dp in second
            column
            
        Raises
        ------
        TypeError
            Mpss needs to have a Dma of Type aerosolpy.dma.DmaCylindrical
            
        Notes
        -----
        assumes constant size-distribution and efficiencies across the width
        of the transfer function
        """
        nsd = []
        if self.polarity=='neg': pol=-1
        if self.polarity=='pos': pol=1
        
        if not isinstance(self.mpss_dma, DmaCylindrical):
            raise TypeError(self.mpss_dma, 
                            ("Mpss needs to have a Dma of Type "
                             "aerosolpy.dma.DmaCylindrical")
                            )
        
        if imax==1:
            for k in range(len(Craw)):
                if (Craw[k]==0) or (self.tot_eff(self.channels[k])==0):
                    nsd.append(0)
                else:
                    beta = self.mpss_dma.q_a/self.mpss_dma.q_sh
                    # using Nraw would need conversion Nraw/(t_res*q_sample)
                    nsd.append(Craw[k]*self.a_function(self.channels[k])
                               /(beta
                                 *self.charge_prob(self.channels[k],pol)
                                 *self.tot_eff(self.channels[k])
                                 )
                                *np.log(10)
                               )
            nsd = np.array(nsd)
            
        elif imax>1:
            ch = len(self.channels)-1
            # sort channels, if they are not ascending iun diameter
            sort = np.argsort(self.channels)
            self.channels = self.channels[sort]
            Craw = Craw[sort]
            
            volts = self.mpss_dma.dp_to_v(self.channels)
            volts_max = np.max(volts)
            nsd = np.zeros(self.channels.shape)
            beta = self.mpss_dma.q_a/self.mpss_dma.q_sh
            
            # only correct channels with volts above volts_max/charge
            for i in range(1,imax+1):
                if len(self.channels[ch:])>1:
                    interpnsd = interp1d(self.channels[ch:],
                                         nsd[ch:],
                                         kind='cubic')
                while volts[ch] > volts_max/(i+1):
                    corr = 0
                    for j in range(2,i+1):
                        dp_inter = self.dma_mpss.v_to_dp(volts[ch], i=j)
                        a_inter = self.a_function(dp_inter)
                        fc_inter = self.charge_prob(dp_inter, pol*j)
                        eff_inter = self.tot_eff(dp_inter)
                        nsd_inter = interpnsd(dp_inter)
                        corr = corr + nsd_inter*fc_inter*eff_inter*beta/a_inter
                    
                    #calc nsd
                    dp1 = self.channels[ch]
                    a1 = self.a_function(dp1)
                    fc1 = self.charge_prob(dp1, pol)
                    eff1 = self.tot_eff(dp1)
                    # do not correct to negative nsd
                    if (Craw[ch]-corr)>=0:
                        nsd[ch] = (Craw[ch]-corr)*a1/(fc1*eff1*beta)
                    else:
                        nsd[ch] = 0
                    
                    ch = ch-1
            # rest of channels with volts[ch] < volts_max/(imax+1) 
            # where all charge corrections are considered
            while ch>=0:
                if len(self.channels[ch:])>1:
                    interpnsd = interp1d(self.channels[ch:],
                                         nsd[ch:],
                                         kind='cubic')
                corr = 0
                for j in range(2,imax+1):
                    dp_inter = self.dma_mpss.v_to_dp(volts[ch], i=j)
                    a_inter = self.a_function(dp_inter)
                    fc_inter = self.charge_prob(dp_inter, pol*j)
                    eff_inter = self.tot_eff(dp_inter)
                    nsd_inter = interpnsd(dp_inter)
                    corr = corr + nsd_inter*fc_inter*eff_inter*beta/a_inter
                
                #calc nsd
                dp1 = self.channels[ch]
                a1 = self.a_function(dp1)
                fc1 = self.charge_prob(dp1, pol)
                eff1 = self.tot_eff(dp1)
                if (Craw[ch]-corr)>=0:
                    nsd[ch] = (Craw[ch]-corr)*a1/(fc1*eff1*beta)
                else:
                    nsd[ch] = 0
                
                ch = ch-1
            nsd = nsd*np.log(10)
        
        return (self.channels, nsd)
        
        def kernel(self,discretizations=128,dj_low=1,dj_high=1000):
            """
            creates inversion kernel matrix for this instrument
            
            Parameters
            ----------
            discretizations : int, or np.array, optional
                number of discretizations for the inversion
                or np.ndarray of doimension 1 with pre-defined 
                diameters, default 128
            dj_low : float, optional
                lower diameter in [nm] for discretization, default 1
            dj_high : float, optional
                higher diameter in [nm] for discretization, default 1000
            
            Returns
            ----------
            array_like of float
                discretization diameters dj
            np.ndarray of float
                2d array fo kernel function evaluated at all discretizations
            
            Notes
            ----------
            kernel is discretized on logarithmic intervals if not
            specified
            """
            transfuncs = []
            I = len(self.channels)
            if np.isscalar(discretizations): 
                transfuncs = []
                J = discretizations
                I = len(self.channels)
                Delta = (np.log10(dj_high)-np.log10(dj_low))/(J-1) 
                dj = np.array([dj_low*10**(j*Delta) for j in range(J)])
                for k in range(len(self.channels)): 
                            transfuncs.append(self.tot_transfunc(self.channels[k]))
                Ainit = np.zeros((I,J))
                for i in range(I):
                    for j in range(J):
                        Ainit[i][j] = transfuncs[i](dj[j])*Delta
                return dj, Ainit
            elif (isinstance(discretizations,np.ndarray)):
                J = len(discretizations)
                for k in range(len(self.channels)): 
                            transfuncs.append(self.tot_transfunc(self.channels[k])) 
                Ainit = np.zeros((I,J))
                ldp = np.log10(discretizations)
                dp_lim = np.zeros(J+1)
                dp_lim[0] = ldp[0]-(ldp[1]-ldp[0])/2.
                dp_lim[1:-1] = (ldp[1:]+ldp[:-1])/2.
                dp_lim[-1] = ldp[-1]+(ldp[-1]-ldp[-2])/2.
                for i in range(I):
                    for j in range(J):
                        Ainit[i][j] = quad(
                                          lambda x: transfuncs[i](10**x),
                                          dp_lim[j],
                                          dp_lim[j+1])[0]
                return dp_lim, Ainit
            else:
                TypeError("discretizations needs to be scalar or np.ndarray")