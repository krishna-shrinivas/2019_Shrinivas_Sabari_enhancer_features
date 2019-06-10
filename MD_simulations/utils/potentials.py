#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 19:45:35 2018

@author: krishna
"""
import numpy
# Define the tabulate bond potential

def harmonic(r, rmin, rmax, k,d1,d2):
    """
    
    HOOMD default table-bond potential parameters include
    r (spherical distance between the monomers), rmin and rmax
    (which dictate the range of the potential), and user defined parameters
    as listed below.
    
    Bonded potentials between two monomers of diameters
    d1 and d2 with a bond stiffness of k are estimated using the 
    harmonic/quadratic potential.
    
    A future version must include the possibility to also determine a bond
    size parameter which would premultiply the hard-core distance. For now,
    that value is hardcoded at 0.65
    
    """
    r0 = 0.65*(d1+d2);
    V = 0.5 * k * (r-r0)**2;
    F = -k*(r-r0);
    return (V, F)


# Define the tabulated shifted LJ potential
# along with the screened charge potential

def lj_shifted(r, rmin, rmax, epsilon,d1,d2,z1,z2,lb,kappa):

    """
    HOOMD default table-pair potential parameters include
    r (spherical distance between the monomers), rmin and rmax
    (which dictate the range of the potential), and user defined parameters
    as listed below.


    Non-bonded shifted Lennard Jones interactions between any two monomers
    The input parameters are the radii of the two types, the strength of
    interaction epsilon. The hard core repulsive potential has a sigma
    determined by tangential overlap length. The potentials are shifted by
    the cut-off length r_max = 2.5 sigma, thought in the future, the plan is to allow
    the zero-shift to happen at a different multiplier.
    
    A now deprecated feature included calculating charged interactions
    by employing the Debye-Huckel potential which requires the charges z1,z2 of
    the monomers, as well as the Bjeurrm length (lb) and inverse Debye screening
    length(k, of units inverse length).
    

    A future version must include the possibility to allow for different r_max
    and different sigma multiplier
    
    """



    sigma = 0.5*(d1+d2);
    r_max = 2.5*sigma;
    
    
    V_shift =  4 * epsilon * ( (sigma / r_max)**12 - (sigma / r_max)**6);
    V =   4 * epsilon * ( (sigma / r)**12 - (sigma / r)**6) - V_shift;
    F =  4 * epsilon / r * ( 12 * (sigma / r)**12 - 6 * (sigma / r)**6);
    
    
#    V_shift = lb * z1 * z2 * numpy.exp(-kappa * r_max) / r_max   + 4 * epsilon * ( (sigma / r_max)**12 - (sigma / r_max)**6);
#    V =  lb * z1 * z2 * numpy.exp(-kappa * r) / r + 4 * epsilon * ( (sigma / r)**12 - (sigma / r)**6) - V_shift;
#    F = lb * z1 * z2 * numpy.exp(-kappa * r) / r**2 + kappa * lb * z1 * z2 * numpy.exp(-kappa * r) / r + 4 * epsilon / r * ( 12 * (sigma / r)**12 - 6 * (sigma / r)**6);
    return (V, F)