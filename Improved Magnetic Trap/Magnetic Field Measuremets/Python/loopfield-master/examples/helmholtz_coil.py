#!/usr/bin/env python3

#
# Helmholtz coil plot example
#
import math
import loopfield as lf
import loopfield.plot as lfp
import numpy as np
import matplotlib.pylab as plt

def Ioffe(Zmin3):
    
    fieldIoffe = lf.Field(length_units = lf.m, current_units = lf.A, field_units = lf.T)

    coil_dx = 0.0e-3
    coil_dy = 0.0e-3
    
    #main part
    Sr0=1.1e-3 #radial spacing of windings
    Sz0=2.6e-3 #axial spacing of windings
    Rin0=2.5e-3 #inner radius
    
    Nz0=7 #windings in axial direction - 4
    Nr0=9 #windings in radial direction - 6
    Ii=-Iq #ioffe current - Ampere
    
    #secondary parts
    #3
    #Zmin3=14.5e-3 #coil axial distance
    # Nz=1 #windings in axial direction 
    Nr3=5 #windings in radial direction
    
    for R in  np.linspace(Rin0,Rin0+Sr0*(Nr3-1),Nr3):    
        fieldIoffe.addLoop(lf.Loop([ coil_dx, coil_dy+ Zmin3, 0.], [0, 1, 0], R, Ii))
    
        #2
        Zmin2=Zmin3+2e-3 #coil axial distance
        #Nz=1 #windings in axial direction 
        Nr2=4 #windings in radial direction
    
    for R in np.linspace(Rin0,Rin0+Sr0*(Nr2-1),Nr2):     
        fieldIoffe.addLoop(lf.Loop([ coil_dx,coil_dy+ Zmin2, 0.], [0, 1, 0], R, Ii))
    
    #1
    Zmin1=Zmin2+2e-3 #coil axial distance
    #Nz=1 #windings in axial direction 
    Nr1=8 #windings in radial direction
    
    for R in np.linspace(Rin0,Rin0+Sr0*(Nr1-1),Nr1):     
        fieldIoffe.addLoop(lf.Loop([coil_dx,coil_dy+ Zmin1, 0.], [0, 1, 0], R, Ii))
    
    Zmin0=Zmin1+2e-3 #coil axial distance
    
    for R in np.linspace(Rin0,Rin0+Sr0*(Nr0-1),Nr0): 
        for dz in np.linspace(0,Sz0*(Nz0-1),Nz0):
            fieldIoffe.addLoop(lf.Loop([ coil_dx,coil_dy+ Zmin0+dz, 0.], [0, 1, 0], R, Ii))
    
    return fieldIoffe

# field object
field = lf.Field(length_units = lf.m,
                 current_units = lf.A,
                 field_units = lf.T)

Sr = 1e-3       # radial spacing of windings
Sz = 2.6e-3     # axial spacing of windings
Rin = 15e-3     # inner radius
Zmin = 37e-3/2.  # lower coil axial distance 14
Nz = 4          # windings in axial direction - 4
Nr = 26         # windings in radial direction - 6
Iq = 1          # Quad current - Ampere

for R in np.linspace(Rin,Rin+Sr*(Nr-1),Nr):
    for dz in np.linspace(0,Sz*(Nz-1),Nz):
        field.addLoop(lf.Loop([ Zmin+dz, 0., 0.], [1, 0, 0], R, -Iq))
        field.addLoop(lf.Loop([-(Zmin+dz), 0., 0.], [1, 0, 0], R, Iq))

##Ioffe coil
num_points = 41
#X_Axis = np.linspace(-Zmin*2, Zmin*2, num_points)
Y_Axis = np.linspace(-20e-3,20e-3, num_points)
X_Field = np.array([field.evaluate([x, 0., 0.]) for x in X_Axis])


for z in np.linspace(0.5e-3,20e-3,40):
    fieldIoffe = Ioffe(z)
    
    # evaluate field at center of coil
        
    XIoffe_Field = np.array([fieldIoffe.evaluate([x, 0., 0.]) for x in X_Axis])
    
    #    plt.figure('X Axis')
    #    plt.plot(X_Axis, 1e4*np.linalg.norm(X_Field,axis = 1))
    #    plt.plot(X_Axis, 1e4*np.linalg.norm(XIoffe_Field,axis = 1))
    #    plt.plot(X_Axis, 1e4*np.linalg.norm(XIoffe_Field+X_Field, axis =1))
    #    plt.title('Abs Value of magnetic field VS distance along X Axis')
    #    plt.xlabel('X Dist(mm)')
    #    plt.ylabel('Abs Magnetic Field(Gauss)')
    #    plt.legend(['Helmholtz', 'Ioffe', 'Sum'])
       
    
    Y_Field = np.array([field.evaluate([0., y, 0.]) for y in Y_Axis])
    YIoffe_Field = np.array([fieldIoffe.evaluate([0., y, 0.]) for y in Y_Axis])
    
    plt.figure('Y Axis')
    plt.plot(Y_Axis*1e3, 1e4*np.linalg.norm(Y_Field,axis = 1))
    plt.plot(Y_Axis*1e3, 1e4*np.linalg.norm(YIoffe_Field,axis = 1))
    plt.plot(Y_Axis*1e3, 1e4*np.linalg.norm(YIoffe_Field+Y_Field, axis =1))
    plt.title('Abs Value of magnetic field VS distance along Y Axis')
    plt.xlabel('Y Dist(mm)')
    plt.ylabel('Abs Magnetic Field(Gauss)')
    plt.legend(['Helmholtz', 'Ioffe', 'Sum'])
    plt.savefig('dy{}'.format(z))

'''
print('Calculating plot...')
# function returns ratio of x-component to that at coil center
def x_ratio(B):
  return B[0] / Bc[0]

# create XY plot
min_x = -Zmin
max_x = Zmin
min_y = -Zmin
max_y = +Zmin
n_x = 101
n_y = 101 
plot = lfp.plotXY(field,
                  min_x, max_x, n_x,
                  min_y, max_y, n_y)

# add field lines
plot.fieldLines()

# add loop symbols
plot.loopSymbols(scale = 1.)

# add 1# error bound region
tol = 0.01
plot.region(x_ratio, [1.-tol, 1.+tol], color='red', alpha=0.5,
            label = ('Field error < %2.1f%%' % (100*tol)))

# add circled area hand-adjusted to fit in 1# error volume "octopus"
center_r = 3.2
plot.circle([0., 0.], radius = center_r, color='blue', alpha=0.5,
            label = ('r = %2.1f cm' % center_r))


# add text
plot.labels(title = '10cm Helmholtz Coil',
            xlabel = 'x (cm)', ylabel = 'y (cm)')

# save plot
plot.save('helmholtz_coil.png')
print('Plot written to "helmholtz_coil.png"')
'''