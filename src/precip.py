# -*- coding: utf-8 -*-

#From Roe 2005

'''
Author: Kevin Rozmiarek
email: kevin.rozmiarek@colorado.edu
'''
import scipy.stats as stats
import os
from math import e

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
        -A text file with temperature as a single column with one time step per row
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


def getPrecip(wind_speed, temp, altitude, moisture_source_altitude):
    '''
    Description here

        Parameters:
                -wind_speed
                    the vertical wind speed, simplified in Roe to be equal to the horizontal wind speed
                -temp
                -altitude
                -moisture_source_altitude

        Returns:
                -A value for the mass of condensating water
    '''
    
    # Find height at each node
    # Find temperature from generated file WRITE IMPORT SCRIPT
    # Find wind for wind file NEED TO MAKE WIND FILE FOR SPACE
    
    z = moisture_source_altitude - altitude
    t = Tz(temp, z)
    
    precip = condes(wind_speed, t, z)
    return precip

###############################################################################

#getPrecip utils

def esat(T, a, b): #Saturated vapor pressure from Clausius-Clapeyron relationship
    esat = 6.112*e((a*T)/(b+T))
    return esat

def pz(z): #Pressure determined from lapse rate
    p0 = 101325 #Pa
    rho0 = 1.225 #kg*m-3
    g = 9.81 #m*s-2
    pz = p0*e(-(rho0*g*z)/p0)
    return pz

def Tz(T0, z): #Temperature approximation by lapse rate. Using dry adiabatic lapse rate from wiki
    Tz = T0 - 0.0098*z
    return Tz

def qsat(T, z): #Saturation-specific humidity aka mass of saturated water vapor per unit mass of the parcel (Wallace & Hobbs 1977)
    qsat = 0.622*(esat(T, 17.67, 243.5)/pz(z))
    return qsat

def pqsatdz(T, z): #Mass of water vapor per unit volume in saturated air
    rho0 = 1.225 #kg*m-3
    gamma = 4000 #m IN THE TRPOICS
    pqsatdz = rho0*qsat(T,z)*e(z/gamma)*(1/gamma)
    return pqsatdz

def condes(w, T, z): #Rate of condensation of water wapor
    condes = -w*pqsatdz(T,z)
    return condes

###############################################################################