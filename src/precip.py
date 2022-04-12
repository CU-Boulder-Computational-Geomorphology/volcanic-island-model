# -*- coding: utf-8 -*-

#From Roe 2005

'''
Author: Kevin Rozmiarek
email: kevin.rozmiarek@colorado.edu
'''
import scipy.stats as stats
import os

def file_gen(filename, years, mean, method, noise, option):
    """
    This code generates a temperature file for use in precip runs.

    Parameters:
        -filename
            The name of the generated file
        -years
            The number of years to generate values for
        -mean
            The average value of the parameter being generator
        -method
            The particular method to use to generate data
                -norm: Add noise to each value sampled from a normal distribustion
                -slope: Adds a linear change to the parameter over years. Then runs norm
        -option
            If needed, an option added to each method
                -norm: Adds a linear scalar on the noise
                -slope: The slope in change per year

    Returns:
        A text file with temperature as a single column with one time step per row
    """
    try:
        os.remove(filename)
    except OSError:
        pass
    
    temp = open(filename,"x")
    for i in range(years): #Main build loop
        
        if method == "norm":
            temp_num = norm_noise(mean, noise)
            
        if method == "slope":
            temp_num = increase(mean, option, i)
            temp_num = norm_noise(temp_num, noise) #Adding norm noise anyways
            
        file_write_value = temp_num
            
        temp.write(str(file_write_value) + "\n")
    temp.close()

#Example call: file_gen("./data/T_air.txt", 1000, -6, "slope", 1, 0.01)

###############################################################################

#file_gen utils

def norm_noise(number, scale): #Function that adds normal noise
    norm_value = stats.norm()
    value_report = float(number + norm_value.rvs(1)*scale)
    return(value_report)
    
def increase(number, slope, year): #Function the addes to a value an amount determined by the number of years that have gone by
    value_report = number + slope * year
    return (value_report)

###############################################################################


def getPrecip():
    '''
    Description here

        Parameters:
                STUFF

        Returns:
                STUFF
    '''
    
    precip = 0
    return precip

###############################################################################

#getPrecip utils

def esat(T, a, b): #Saturated vapor pressure from Clausius-Clapeyron relationship
    esat = 0
    return esat

def pz(z): #Pressure determined from lapse rate
    pz = 0
    return pz

def Tz(T0, gamma, z): #Temperature approximation by lapse rate
    Tz = 0
    return Tz

def qsat(esat, p): #Saturation-specific humidity aka mass of saturated water vapor per unit mass of the parcel (Wallace & Hobbs 1977)
    qsat = 0
    return qsat

def pqsat(p0, q0sat, z, a, b, gamma): #Mass of water vapor per unit volume in saturated air
    pqsat = 0
    return pqsat

def condes(w, pqsatdz): #Rate of condensation of water wapor
    condes = 0
    return condes

###############################################################################