# aerosolpy - version histroy

## v1.0.2

fixes scipy.misc.derivative being deprecated. uses numdifftools instead.n

## v1.0.0

`aerosolpy.growth.VbsModel` includes particle phase diffusion. 
no backward-compatibility to previous versions due to changes in method
definitions

## v0.4.0

includes subpackage `aerosolpy.growth` with classes for kinetically-limited
growth, H2SO4 growth, volatility calculation and VBS growth module

## v0.3.0

includes dayofyear_to_datetime function in time.py module

fixed bug in DMA class penetration loss calculation

## v0.2.0

includes subpackage `aerosolpy.instr` with classes for DMAs, CPC, and MPSS.

fixed an essential calculation error in 
`aerosolpy.AerosolMechanics._slipcorr_deriv` which was there since the 
beginning. 
    

## v0.1.0 

first stable version. only included aerosolpy

### v0.1.1 

some fixes, and updates on documentation