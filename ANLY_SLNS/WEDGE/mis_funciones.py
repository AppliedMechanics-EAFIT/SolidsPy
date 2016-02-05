# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 13:08:39 2015

@author: Grupo 063
"""
import numpy as np
   
def radianes(ang_grad):
    ang_rad=ang_grad*np.pi/180    
    return ang_rad      
    
def grados(ang_rad):
    ang_grad=ang_rad*180/np.pi    
    return ang_grad      
    
    
def tensor_polar(r,teta,f,beta,alfa):
    
    alfa=radianes(alfa) # para semi-espacio
    teta=radianes(teta) # paso teta a radianes
    beta=radianes(beta) # paso beta a radianes
    p=f*np.sin(beta)
    q=f*np.cos(beta)
    srp=-2.*p*np.cos(teta)/((2.*alfa+np.sin(2.*alfa))*r)
    srq=-2.*q*np.sin(teta)/((2.*alfa-np.sin(2.*alfa))*r)
    sigma=np.zeros((2,2))
    sigma[0,0]=srp+srq
    return sigma  
    
def tensor_cart(r,teta,f,beta):
    #################################### variables de entrada ####################
    # r=...........................radio (coordenada polar)
    # teta=........................Angulo en grados (coordenada polar)
    # f=...........................Fuerza ()
    # beta=....................... Angulo de la Fuerza en grados
    ##############################################################################
    teta=radianes(teta) # paso teta a radianes
    beta=radianes(beta) # paso beta a radianes
    #
    p=f*np.sin(beta)
    q=-f*np.cos(beta)
   
    sxp=-2.*p*np.cos(teta)**3./(np.pi*r)
    syp=-2.*p*np.cos(teta)*np.sin(teta)**2./(np.pi*r)
    txyp=-2.*p*np.cos(teta)**2.*np.sin(teta)/(np.pi*r)
    
    sxq=2.*q*np.cos(teta)**2.*np.sin(teta)/(np.pi*r)
    syq=2.*q*np.sin(teta)**3./(np.pi*r)
    txyq=2.*q*np.sin(teta)**2.*np.cos(teta)/(np.pi*r)
    #
    sigmaf=np.zeros((2,2))
    sigmaf[0,0]=sxp+sxq
    sigmaf[1,1]=syp+syq
    sigmaf[0,1]=txyp+txyq
    sigmaf[1,0]=txyp+txyq    
    return sigmaf       
    
def tensor_cart_m(r,teta,f,m,beta):    
    teta=radianes(teta) # paso teta a radianes
    beta=radianes(beta) # paso beta a radianes
    #
    p=f*np.sin(beta)
    q=-f*np.cos(beta)
# Carga Vertical   
    sxp=-2.*p*np.cos(teta)**3./(np.pi*r)
    syp=-2.*p*np.cos(teta)*np.sin(teta)**2./(np.pi*r)
    txyp=-2.*p*np.cos(teta)**2.*np.sin(teta)/(np.pi*r)
# CArga Horizontal   
    sxq=2.*q*np.cos(teta)**2.*np.sin(teta)/(np.pi*r)
    syq=2.*q*np.sin(teta)**3./(np.pi*r)
    txyq=2.*q*np.sin(teta)**2.*np.cos(teta)/(np.pi*r)
# Momento
    sxm=np.cos(teta)**3.*np.sin(teta)*(8.*m/(np.pi*r**2.))
    sym=np.cos(2.*teta)*np.sin(2.*teta)*(-2.*m/(np.pi*r**2.))
    txym=(3.*np.sin(teta)**2.*np.cos(teta)**2.-np.cos(teta)**4.)*(2.*m/(np.pi*r**2.))
# Tensor solución       
    sigmaf=np.zeros((2,2))
    sigmaf[0,0]=sxp+sxq+sxm
    sigmaf[1,1]=syp+syq+sym
    sigmaf[0,1]=txyp+txyq+txym
    sigmaf[1,0]=txyp+txyq+txym    
    return sigmaf
#
def cunia(x,y,phi,l,nu,E,S):
#
    K1=(np.cos(phi)/np.sin(phi))+nu*(np.sin(phi)/np.cos(phi))
    K2=(np.sin(phi)/np.cos(phi))+nu*(np.cos(phi)/np.sin(phi))
    ux=(S/E)*K1*(x-l*np.cos(phi))
    uy=-(S/E)*K2*y
    return ux , uy    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    