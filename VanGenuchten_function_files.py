# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 17:14:47 2021

@author: Richa_Pc01
"""

"""
----------------------------defining required Van genuchten functions
"""

import numpy as np

def Van_Genuchten_and_Nelson_moisure( x , thetar, thetas , hs , ho , N, Ss):
    """
    calculation of soil moisture based on modified Van Genuchten and  Nelson (1985)
    thetar = residual moisture content, thetas = saturated soil moisture , h = presure head, 
    hs = air entry presure head , N = Van Genuchten paramter , Ss =  specific storage
    """
    y = np.empty([len(x)])
    tempI = x > ho
    m = 1 - 1/N
    beta = pow( abs(x[~tempI]/hs), N )
    betao = pow(ho/hs,N)
    y[~tempI] =  thetar + ( thetas - thetar )* pow( 1+ beta,-m) 
    y[tempI] =  thetar + ( thetas - thetar )*pow( 1 + betao, -m) + Ss*( x[tempI] - ho )   
    return y


def Van_Genuchten_and_Nelson_hydraulic_conductivity(x, Ksat, hs, N):
    """
    calculation of soil moisture based on modified Van Genuchten and  Nelson (1985)
    h = presure head, 
    hs = air entry presure head , N = Van Genuchten paramter , 
    Ksat = saturated hyraulic conductivity ,
    """
    y = np.empty([len(x)])
    tempI = x < 0
    m = 1 - 1/N
    beta = pow( abs(x[tempI]/hs), N )
    y[tempI] =Ksat * (  pow(1+beta,-5*m/2))* pow( pow(1 + beta, m) -  pow(beta,m) , 2 ) 
    y[~tempI] =  Ksat  
    return y

def Van_Genuchten_and_Nelson_specific_moisure_capacity( x , thetar, thetas, hs, ho, N, Ss ):
    """
    calculation of soil specific moisture capacity based on 
    modified Van Genuchten and Nelson (1985)
    thetar = residual moisture content, thetas = saturated soil moisture , h = presure head, 
    hs = air entry presure head , N = Van Genuchten paramter , % Ss =  specific storage
    """
    y = np.empty([len(x)])
    tempI = x > ho 
    m = 1 - 1/N 
    beta = pow( abs(x[~tempI]/hs), N )
    y[tempI] = Ss
    y[~tempI] = (N-1)*(thetas - thetar)* pow(abs(x[~tempI]),N-1)/( pow(abs(hs),N)*pow(1 + beta, m+1) )   
    return y

def Van_Genuchten_and_Nelson_C( x , thetar, thetas , hs , ho , N, Ss ):
    """
    calculation of soil specific moisture capacity based on 
    modified Van Genuchten and Nelson (1985)
    thetar = residual moisture content, thetas = saturated soil moisture , h = presure head, 
    hs = air entry presure head , N = Van Genuchten paramter , % Ss =  specific storage
    """
    if x > ho:
        y = Ss
    else:
        m = 1 - 1/N 
        beta = pow( abs(x/hs), N )
        y = (N-1)*(thetas - thetar)* pow(abs(x),N-1)/( pow(abs(hs),N)*pow(1+ beta, m+1) )   
    return y

