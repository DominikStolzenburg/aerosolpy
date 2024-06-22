Introduction to aerosolpy
=========================

Purpose
-------

``aerosolpy`` is a collection of functions and classes which developed over 
the years since my PhD in aerosol physics with all kinds of utilities which
facilitate calculations related to aerosol physics and chemistry. 

``aerosolpy`` is not a package which serves a single puropse and by no means 
a compelete software tool which aims to provide everything to the user. 
It is really a collection of different things. There might be evry important 
functions missing and some functions might seem useless to another user. 

However, as it is structured currently it provides a python package which 
is useful especially for calculating basic mechanics and kinetics operations,
provides classes which describe typical aerosol sizing instruments, and will 
be extended to also include data inversion software for mobility particle size
spectrometers and a simulation package which facilitates aerosol growth 
calculations. 

Disclaimer
----------

I will not be liable for potential errors in these calculations. 
It is actually highly likely that there might be some mistakes. Please be
careful when using ``aerosolpy``.

Basics
------

Once installed import aerosolpy as follows, using the recommended abbreviation 
ap::

        import aerosolpy as ap
    
There are serval base functions provided within aerosolpy which include
typical unit conversions, time format conversions and math operations. 
They can be accessed within the ap namespace::

        ap.lpm_to_m3pers(1.5)
    
This for example converst the typical air flow unit litre per minute to SI
units, i.e. cubic meter per second. 

Another operation which was often needed when analyzing data with different 
software was the conversion of time stamps, especially going from the widely
used ``MATLAB`` to ``python``::

        ap.matlab2datetime(matlab_timestamp)
        
Whenever the need for a new basic function will appear in my future work I will
update them in ``aerosolpy``.

AerosolMechanics and AerosolKinetics
------------------------------------

Two major classes included in the basic ``aerosolpy`` package are
``ap.AerosolMechanics`` and ``ap.AerosolKinetics``. 

``AerosolMechanics`` was designed to provide calculations related to the 
motion of aerosol particles and vapor molecules within air. It is designed
as a class, that means for accessing its methods (the necessary calculations)
we need to create an instance of that class, an ``AerosolMechanics`` object::

        amech = ap.AerosolMechanics(temp_kelvin=273.15, press_hpa=1013)

The object defines all physics parameters of air at the given temperature and
pressure. It for example calculates the mean free path of air molecules and the
dynamic viscosity of air directly. Temperature and pressure of a given 
``ap.AerosolMechanics`` object can be changed during program execution, e.g.::

        amech.set_temp(293.15)
        
``AerosolKinetics`` was designed to provide kinetic calculations in the sense
of chemical kinetics, i.e. collision frequencies between aerosols and aerosols
and vapor molecules. As with ``AerosolMechanics`` it is implemented as class,
which inherits from ``ap.AerosolMechanics``, i.e. an object of ``AerosolKinetics``
also includes the methods from ``AerosolMechanics``::

        akin = ap.AerosolKinetics(temp_kelvin=273.15)
        akin.set_temp(293.15)

Its main purpose however is to provide additional functions such as the 
coagulation kernel of two aerosol particles::

        akin.coll_kernel_pp(10.0,100.0)
        
Subpackage: aerosolpy.instr
---------------------------

The ``aerosolpy.instr`` package provides additional classes which define 
objects related to typical aerosol measurement devices, currently implemented
is a differential mobility analyzer (DMA; ``ap.instr.Dma`` and its subclass 
``ap.instr.DmaCylindrical``), a condensation particle counter (CPC; 
``ap.instr.Cpc``) and a mobility particle size spectrometer 
(MPSS; ``ap.instr.Mpss``). Here the object-orientated programming
style of ``aerosolpy`` might be slightly more intuitive. 

An instance of the ``ap.instr.DmaCylindrical`` class just represents a 
cylindrical DMA with its geometry (radii of the cylinders, length of the 
aerosol extraction pathways), flow conditions (aerosol sample and sheath flows)
and other parameters. It the provides methods to e.g. convert voltages applied
to the central electrode to diameters classified assuming a certain charge on 
the particles::

        cy_dma = ap.instr.DmaCylindrical(q_a=1.5, q_sh=7.5, length=0.109, 
                                         r_i=0.025, r_o=0.033)
        cy_dma.v_to_dp(1000)
    
In a similar manner we can define CPCs using ``ap.instr.Cpc`` class with e.g.
predefined activation curves::

        tsi3772_cpc = ap.instr.Cpc(activation='tsi_3772_ag')

Mobility Particle Size Spectrometers build upon these two classes as every
MPSS typically consists of a DMA for particle size classification and a CPC
for particle detection (concentration measurement)::

        channel_sizes = np.array([10.,20.,30.,50.,100.])
        custom_mpss = ap.instr.Mpss(channel_sizes, cy_dma, cpc=tsi3772_cpc)

It then contains methods which can invert MPSS data::

        raw_conc_per_channel = np.array([2, 7.2, 8.9, 3., 0.])
        nsd = custom_mpss.std_inv(raw_conc_per_channel, imax=5)

where ``nsd`` is then the inverted size distribution taking into account 
multiple charges on aerosol particles up to ``imax=5``. 