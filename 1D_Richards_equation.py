# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 15:06:44 2021

@author: Saurabh Kumar 
Code based on Dogan and Mortz 2002a,b 
"""
import numpy as np
import matplotlib.pyplot as plt
from VanGenuchten_function_files import Van_Genuchten_and_Nelson_moisure,Van_Genuchten_and_Nelson_hydraulic_conductivity,Van_Genuchten_and_Nelson_specific_moisure_capacity

# Dogan, A., & Motz, L. H. (2005). Saturated-unsaturated 3D groundwater model. I: Development. 
# Dogan, A., & Motz, L. H. (2005). Saturated-unsaturated 3D groundwater model. II: Verification and application

# test problem : 1D Infiltration and Redistribution test problem from Paniconi et al. (1991) 
# Paniconi, C., Aldama, A. A., & Wood, E. F. (1991). Numerical evaluation of iterative and noniterative methods for the solution of the nonlinear Richards equation


zmin = 0
zmax = 10  # in m
dz = 0.1 # in m
z = np.arange(zmin + dz/2,zmax+dz/2,dz)
tmin = 0; # starting time in hr
tmax = 32; # max time in hr
dt = 0.1; #in hr
t = np.arange(tmin,tmax+dt,dt) 

Ksat = 5   # m/h
thetas = 0.45 
thetar = 0.08
hs = -3    # air entry presure in m
ho = -0.19105548 # in cm/h
N = 3 # van Genucten Paramter 
Ss = 0.001 #  specific storage 

q = t/64 # flux in m/hr at the top boundary 
tnodes = len(t)
znodes = len(z)
# preallocate h and theta 
h = np.empty([znodes,tnodes])
theta=np.empty([znodes,tnodes])



h[:,0] = -z # initial presure condition 
# initial moisture condition
theta[:,0] = Van_Genuchten_and_Nelson_moisure( h[:,0] , thetar, thetas, hs, ho, N, Ss)


threshold = .001
MaxIterations = 50
Error = np.zeros([tnodes,MaxIterations])


RHS = np.empty([znodes])
A = np.empty([znodes,znodes])

#K = Van_Genuchten_and_Nelson_hydraulic_conductivity(h[:,0], Ksat, hs, N)

for n in range(tnodes-1):
    h[:,n+1] = h[:,n]
    theta [:,n+1] = theta[:,n]
    
    for m in range(MaxIterations):
        RHS = np.empty([znodes])
        A = np.zeros([znodes,znodes])
        K = Van_Genuchten_and_Nelson_hydraulic_conductivity(h[:,n+1], Ksat, hs, N)
        C = Van_Genuchten_and_Nelson_specific_moisure_capacity(h[:,n+1], thetar, thetas, hs, ho,N, Ss)
        # lower boundary condition 
        i = 0 ;
        #forming the A matrix and RHS vector
        A[i,i]=1;          
        RHS[i] = 0
        
        for i in range(1,znodes):
            h2  = h[i,n+1]
            K2  = K[i]
            K1Z = K[i-1] 
            theta_current = theta[i,n+1]
            theta_old = theta[i,n]
            # calculation of s
            s =(theta_current-theta_old)/dt
            # calculation of c           
            Kszminus = ( K2 + K1Z ) / 2
            CNZminusmean = -Kszminus/dz
            c = -CNZminusmean/dz;             
            # calculation of p1 
#            C = Van_Genuchten_and_Nelson_C( h2 , thetar, thetas , hs , ho , N, Ss )
            p1 = C[i]/dt
            # calculation of p2
            p2 = 0; # Sw*Ss/dt; since Sw=0
            if i < znodes-1:
                # calculation of g
                K3Z = K[i+1]
                Kszplus = ( K2 + K3Z )/2
                CNZplusmean = -Kszplus/dz
                g = -CNZplusmean/dz  
                # calculation of d          
                d = -(  c + g + p1 + p2)
                # forming the A matrix and RHS vector
                A[i,i-1] = c   # coeffiecient of H(k-1)
                A[i,i] = d      # coeffiecient of H(k)
                A[i,i+1] = g  # coeffiecient of H(k+1)
                RHS[i] = s - p1 *( h[i,n+1] + z[i] ) - p2*( h[i,n] + z[i])  # RHS at i
            elif i == znodes-1:
                d = -(  c +  p1 + p2)
                A[i,i-1] = c   # coeffiecient of H(k-1)
                A[i,i] = d      # coeffiecient of H(k)
                RHS[i] = s - p1 *( h[i,n+1] + z[i] ) - p2*( h[i,n] + z[i]) - q[n+1]/dz
           
            
        X_input = h[:,n+1]+z
        #X = np.linalg.inv(A).dot(RHS)
        X = np.linalg.solve(A,RHS)
        Error[n,m] =  np.sqrt(np.sum(pow(X-X_input,2)))
        # updating the presure head and theta
        h[:,n+1] = X-z
        theta[:,n+1] = Van_Genuchten_and_Nelson_moisure( h[:,n+1], thetar, thetas, hs, ho, N, Ss)
        if Error[n,m] <= threshold:
            break
            

# observed data     

h_10hr = [ -3.955696203, -3.987341772, -3.971518987,
-4.003164557,-4.034810127, -4.034810127,-4.034810127,-4.018987342,
-4.034810127,-4.003164557,-3.908227848,-3.75,
-3.528481013,-3.180379747,-2.800632911,-1.946202532,-0.949367089,-0.443037975]

Z_10hr = [ 9.95229981, 9.439537256, 8.926687764, 8.41392521, 7.919477674, 7.42497218, 6.930466685, 6.45424723, 5.959770715,
5.446892243,4.933897853, 4.439102564, 3.962506375, 3.467363333, 2.935532295, 1.981586452, 0.95411972, 0.422056846 ]


h_1hr = [-7.215189873, -8.876582278, -8.987341772,-8.481012658,
-8.022151899,-7.5,-6.993670886,-6.487341772,
-6.012658228,-5.506329114,-5,-4.493670886,-4.003164557,
-3.496835443,-3.006329114,-1.993670886,-0.996835443, -0.474683544]

Z_1hr = [ 9.976584597,9.485121946,8.954189271,8.477071452,
7.963410535,7.467948718,6.990830899,6.477083044,
5.981708165,5.449645291,4.954212454,4.42214958,
3.963375759,3.467942922,2.954224046,1.945043353,0.954206658,0.403799787]
# -------figures-----------------------------
plt.figure()
plt.plot(h[:,0],z,label='t = 0 hr')
plt.plot(h[:,10],z,label='t = 1 hr')
plt.plot(h_1hr,Z_1hr,'o',  label='t = 1 hr (o)')
plt.plot(h[:,20],z,label='t = 2 hr')
plt.plot(h[:,40],z,label='t = 4 hr') 
plt.plot(h[:,100],z,label='t = 10 hr')
plt.plot(h_10hr,Z_10hr,'o',  label='t = 10 hr (o)')
plt.plot(h[:,320],z,label='t = 32 hr')        
plt.ylabel('Elevation (m)')
plt.xlabel("Presure head(m)")
plt.legend()
plt.tight_layout()

plt.savefig('Richards_1D_plot_h_vs_elevation.png', dpi=200)


plt.figure()
plt.plot(theta[:,0],z,label='t = 0 hr')
plt.plot(theta[:,20],z,label='t = 2 hr')
plt.plot(theta[:,40],z,label='t = 4 hr') 
plt.plot(theta[:,100],z,label='t = 10 hr')
plt.plot(theta[:,320],z,label='t = 32 hr')  
plt.ylim([10,0])      
plt.xlim([0,0.5])      
plt.ylabel('Elevation (m)')
plt.xlabel(r"$\theta$")
plt.legend()
plt.tight_layout()
plt.savefig('Richards_1D_plot_theta_vs_elevation.png', dpi=200)










