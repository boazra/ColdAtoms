#!/usr/bin/env python3

#
#
import loopfield as lf
import loopfield.plot as lfp
import numpy as np
import matplotlib.pylab as plt
import PyAstronomy.pyasl as pyasl


# Constants

MagneticDipoleK = 34.523

def quadFit(y,p,shift):
    return p[0]*(y-shift)**2+p[1]*(y-shift)+p[2]

def Bias(Nz, Nr, I,y0):
    field = lf.Field(length_units = lf.m, current_units = lf.A, field_units = lf.T)    
    Sr = 1.25e-3       # radial spacing of windings
    Sz = 2.5e-3     # axial spacing of windings
    Rin = 40e-3     # inner radius
    Zmin = 40e-3    # coil axial distance
    
    for R in np.linspace(Rin,Rin+Sr*(Nr-1),Nr):
        for dz in np.linspace(0,Sz*(Nz-1),Nz):
            field.addLoop(lf.Loop([ 0., Zmin+dz+y0, 0.], [0, -1, 0], R, I))
            field.addLoop(lf.Loop([0, -(Zmin+dz)+y0, 0.], [0, -1, 0], R, I))
    return field
        

def Ioffe(Zmin4,z0,Iq):
    
    fieldIoffe = lf.Field(length_units = lf.m, current_units = lf.A, field_units = lf.T)

    coil_dx = 0.0e-3
    coil_dy = 0.0e-3
    
    #main part
    Sr0=1.25e-3 #radial spacing of windings
    Sz0=2.5e-3 #axial spacing of windings
    Rin0=2.5e-3 #inner radius
    
    Nz0=6 #windings in axial direction - 4
    Nr0=7 #windings in radial direction - 6
    Ii=-Iq #ioffe current - Ampere
    
    #secondary parts
    #4
    Nr4 = 2
    
    for R in  np.linspace(Rin0,Rin0+Sr0*(Nr4-1),Nr4):    
        fieldIoffe.addLoop(lf.Loop([ coil_dx, coil_dy+ Zmin4, z0], [0, 1, 0], R, Ii))    
    
    #3
    Zmin3 = Zmin4 + Sz0
    #Zmin3=14.5e-3 #coil axial distance
    # Nz=1 #windings in axial direction 
    Nr3=4 #windings in radial direction
    
    for R in  np.linspace(Rin0,Rin0+Sr0*(Nr3-1),Nr3):    
        fieldIoffe.addLoop(lf.Loop([ coil_dx, coil_dy+ Zmin3, z0], [0, 1, 0], R, Ii))
    
    #2
    Zmin2=Zmin3+Sz0 #coil axial distance
    #Nz=1 #windings in axial direction 
    Nr2=5 #windings in radial direction

    for R in np.linspace(Rin0,Rin0+Sr0*(Nr2-1),Nr2):     
        fieldIoffe.addLoop(lf.Loop([ coil_dx,coil_dy+ Zmin2, z0], [0, 1, 0], R, Ii))

    #1
    Zmin1=Zmin2+Sz0 #coil axial distance
    #Nz=1 #windings in axial direction 
    Nr1=6 #windings in radial direction
    
    for R in np.linspace(Rin0,Rin0+Sr0*(Nr1-1),Nr1):     
        fieldIoffe.addLoop(lf.Loop([coil_dx,coil_dy+ Zmin1,z0], [0, 1, 0], R, Ii))
    
    Zmin0=Zmin1+Sz0 #coil axial distance
    
    for R in np.linspace(Rin0,Rin0+Sr0*(Nr0-1),Nr0): 
        for dz in np.linspace(0,Sz0*(Nz0-1),Nz0):
            fieldIoffe.addLoop(lf.Loop([ coil_dx,coil_dy+ Zmin0+dz, z0], [0, 1, 0], R, Ii))
    
    return fieldIoffe

def plot_trap(move_Ioffe_y,move_Ioffe_x,Current,save=False):
    plt.figure(0,(20,10))
    for i in Current:        
        for y in move_Ioffe_y:
            index_move_y = int(np.round(y/dy))
            if index_move_y > 0:
                Y_Axis_plot = Y_Axis[:-index_move_y]
            elif index_move_y < 0:
                Y_Axis_plot = Y_Axis[-index_move_y:]
            else:
                Y_Axis_plot = Y_Axis                
            for x in move_Ioffe_x:
                index_move_x = int(np.round(x/dx))    
                if index_move_y > 0:
                    if index_move_x == 0:
                        temp_Y_Field = Y_Field[:,index_move_y:]                    
                    else:
                        temp_Ioffe_Field = Ioffe_Field[index_move_x:,:-index_move_y]*i/30.0
                        temp_Y_Field = Y_Field[:-index_move_x,index_move_y:]                                  
                elif index_move_y < 0:
                    if index_move_x == 0:
                        temp_Y_Field = Y_Field[:,:index_move_y]                             
                    else:
                        temp_Y_Field = Y_Field[:-index_move_x,:index_move_y]                       
                else:
                    if index_move_x == 0:
                        temp_Y_Field = Y_Field                       
                    else:                         
                        temp_Y_Field = Y_Field[:-index_move_x,:]                  
                temp_Ioffe_Field = Ioffe_Field[int(abs(index_move_x)):,int(abs(index_move_y)):]*i/30.0
                sumField = 1e4*np.linalg.norm(temp_Ioffe_Field+temp_Y_Field,axis=-1)                              
                #plt.title('I = {0:03.2f} A, Ioffe at ({1:03.2f},{2:03.2f})mm'.format(i,x*1e3,(y0+y)*1e3)) 
                Extent=[-15-y*1e3,15,-10+x*1e3,10]
                plt.subplot(1,2,1)
                plt.cla()   
                plt.title('I = {0:03.2f} A, Ioffe at {1:03.2f}mm'.format(i,(y0+y)*1e3)) 
                plt.imshow(sumField,vmin = 0, vmax = 300, extent=Extent)
                plt.contour(sumField,list(range(0, 600, 15)), linewidths=1, colors = 'r', extent=Extent)   
                plt.subplot(1,2,2)                                         
                plt.cla()   
                plt.plot(Y_Axis_plot*1e3, 1e4*np.linalg.norm(temp_Y_Field[200,:,:],axis = 1))
                plt.plot(Y_Axis_plot*1e3, 1e4*np.linalg.norm(temp_Ioffe_Field[200,:,:],axis = 1))                              
                plt.plot(Y_Axis_plot*1e3,sumField[200,:])  
                plt.ylim([-10,150])
                plt.annotate('({0:02.2f},{1:02.2f})'.format(Y_Axis_plot[sumField[200,:].argmin()]*1e3,sumField[200,:].min()),(Y_Axis_plot[sumField[200,:].argmin()-38]*1e3,sumField[200,:].min().min()-7), fontsize = '14')
                plt.legend(['Quad', 'Ioffe', 'Sum'], fontsize = '18', loc = 2)
                if save:
                    filename = r'C:\Users\Boaz\Documents\GitHub\ColdAtoms\Improved Magnetic Trap\Magnetic Field Measuremets\Python\Ioffe movie\dx={1:05.2f}mm_dy={2:05.2f}mm_i{0:05}mA.jpg'.format(i*10,x*1.0e3,(y0+y)*1.0e3)                
                    plt.savefig(filename, dpi = 300, bbox_inches='tight', pad_inches=0.2)

'''
# field object
field = lf.Field(length_units = lf.m,
                 current_units = lf.A,
                 field_units = lf.T)

Sr = 1.25e-3       # radial spacing of windings
Sz = 2.5e-3     # axial spacing of windings
Rin = 12.5e-3     # inner radius
Zmin = 37.0e-3/2.0 # lower coil axial distance 14
Nz = 4          # windings in axial direction - 4
Nr = 26         # windings in radial direction - 6
Iq = 30.0         # Quad current - Ampere

for R in np.linspace(Rin,Rin+Sr*(Nr-1),Nr):
    for dz in np.linspace(0,Sz*(Nz-1),Nz):
        field.addLoop(lf.Loop([ Zmin+dz, 0., 0.], [1, 0, 0], R, -Iq))
        field.addLoop(lf.Loop([-(Zmin+dz), 0., 0.], [1, 0, 0], R, Iq))
'''

num_points_x = 401
num_points_y = 601
X_Axis = np.linspace(-10e-3,10e-3,num_points_x)
Y_Axis = np.linspace(-15e-3,15e-3, num_points_y)
dy = Y_Axis[1]-Y_Axis[0]
y0 = 14.5e-3
dx = X_Axis[1]-X_Axis[0]

'''
points = np.reshape(np.array([[x, y,0] for x in X_Axis for y in Y_Axis]),(-1,3))
Y_Field = np.reshape(field.evaluate(points),(num_points_x,num_points_y,3))
fieldIoffe = Ioffe(y0,0,Iq)
Ioffe_Field = np.reshape(fieldIoffe.evaluate(points),(num_points_x,num_points_y,3))
'''

'''
Y_Field = np.load(r'C:\Users\Boaz\Documents\GitHub\ColdAtoms\Improved Magnetic Trap\Magnetic Field Measuremets\Python\Ioffe movie\Y_Field.npy')
Ioffe_Field = np.load(r'C:\Users\Boaz\Documents\GitHub\ColdAtoms\Improved Magnetic Trap\Magnetic Field Measuremets\Python\Ioffe movie\Ioffe_Field.npy')
'''
font = {'family' : 'DejaVu Sans',        
        'size'   : 24}
plt.rc('font', **font)
plt.ioff()

move_Ioffe_y =[-3.5e-3] #(np.arange(41)-20)*0.25e-3
move_Ioffe_x = [2e-3] #(np.arange(41))*0.2e-3
Current = [30] #np.linspace(0,30,61)   
plot_trap(move_Ioffe_y,move_Ioffe_x,Current,True)
plot_trap( [-10e-3],[0],[30],True)


    
