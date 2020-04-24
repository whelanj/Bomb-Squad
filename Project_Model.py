# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 16:06:27 2020

@author: jwhel
"""

# Documentation

# The University of Texas at Austin - Spring 2020
# ME 337G - Nuclear Safety and Security - Dr. HAAS, Derek 
# Team 7 - Bomb Squad - INANC, Ece Shelby WHELAN, Jack  

# This code uses the Bateman equation to solve for the concentrations of
# the daughter products of Xe-140, a fission fragment of U-235.

## Constants 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

#Place holders for velocity terms
LeakVelocity = 1 #[m/s]
SeepVelocity = 10 #[m/s]
VentVelocity = 100 #[m/s]

#Reactivity Constants including placeholders [1/m]
XeReactivity = 0
CsReactivity = np.log(50)
BaReactivity = np.log(40)
LaReactivity = np.log(30)
CeReactivity = np.log(20)

#placeholder time value
t = 1
#Reactivity Decay Terms : exp([1/m]*[m/s]*[s])
RD_Xe_Leak = np.exp(XeReactivity*LeakVelocity*t)
RD_Xe_Seep = np.exp(XeReactivity*SeepVelocity*t)
RD_Xe_Vent = np.exp(XeReactivity*VentVelocity*t)

RD_Cs_Leak = np.exp(CsReactivity*LeakVelocity*t)
RD_Cs_Seep = np.exp(CsReactivity*SeepVelocity*t)
RD_Cs_Vent = np.exp(CsReactivity*VentVelocity*t)

RD_Ba_Leak = np.exp(BaReactivity*LeakVelocity*t)
RD_Ba_Seep = np.exp(BaReactivity*SeepVelocity*t)
RD_Ba_Vent = np.exp(BaReactivity*VentVelocity*t)

RD_La_Leak = np.exp(LaReactivity*LeakVelocity*t)
RD_La_Seep = np.exp(LaReactivity*SeepVelocity*t)
RD_La_Vent = np.exp(LaReactivity*VentVelocity*t)

RD_Ce_Leak = np.exp(CeReactivity*LeakVelocity*t)
RD_Ce_Seep = np.exp(CeReactivity*SeepVelocity*t)
RD_Ce_Vent = np.exp(CeReactivity*VentVelocity*t)

#timespan adjusts the time for which the model will examine, also affects frames in gif
timespan = 1000 # [s]


HL_XE = 14 # Half-life of Xe-140, [s]
HL_CS = 64 # Half-life of Cs-140, [s]
HL_BA = 13*24*60*60 # Half-life of Ba-140, [s]
HL_LA = 40*60*60 # Half-life of La-140, [s]

LL_XE = np.log(2)/HL_XE # Decay constant of Xe-140, [1/s]
LL_CS = np.log(2)/HL_CS # Decay constant of Cs-140, [1/s]
LL_BA = np.log(2)/HL_BA # Decay constant of Ba-140, [1/s]
LL_LA = np.log(2)/HL_LA # Decay constant of La-140, [1/s]
LL_CE = 0

LL = [LL_XE, LL_CS, LL_BA, LL_LA, LL_CE]

## Compuate Coefficients

coefficients = []
for NUM_COEFF in range(0,5):
    row = []
    for j in range(0,NUM_COEFF+1):
        NUMERATOR = 1
        DENOMINATOR = 1
        for k in range(0,NUM_COEFF+1):
            NUMERATOR *= LL[k]
            if k is not j:
                DENOMINATOR *= LL[k]-LL[j]
        row.append(NUMERATOR/DENOMINATOR)
    while len(row) < 5:
        row.append(0)
    coefficients.append(row)
#    print(row)
COEFF = np.array(coefficients)


## Compute Radionuclide Density
## NN[0] is the starting amounts of each radionuclide
Initial_amount_fission_product = 1
NN = [[Initial_amount_fission_product,0,0,0,0]]
T = 0
for i in range(0,int(timespan/2)):
    row = []
    T += 2 # Increment time, [s]
    for j in range(0,5):
        SUMM = 0
        for k in range(0,j+1):
            SUMM += COEFF[j][k]*np.exp((-1)*LL[k]*T) #this is where we would multiply by the relevant reactivity decay term
        if j < 4:
            row.append(NN[0][0]*SUMM/LL[j])
    row.append((NN[i-1][4]+LL[3]*(row[3]-NN[i-1][3])))    
        
    NN.append(row)
    
Xe = []
Cs = []
Ba = []
La = []
Ce = []
for x in NN:
    Xe.append(x[0])
    Cs.append(x[1])
    Ba.append(x[2])
    La.append(x[3])
    Ce.append(x[4])
  
NN = np.array(NN)

time = [*range(0,timespan+1,2)]
plt.figure()

plt.plot(time,Xe,label = 'Xe-140')

plt.plot(time,Cs,label = 'Cs-140')

plt.plot(time,Ba,label = 'Ba-140')

plt.plot(time,La,label = 'La-140')

plt.plot(time,Ce,label = 'Ce-140')

plt.yscale('log')
plt.xlabel('Time, [s]')
plt.ylabel('Nuclear Density, [# of nuclei/cm^3]')
plt.title('Nuclear Density of Xe-140 Daughter Products')
plt.legend()
plt.grid()
plt.savefig('Model_Graph.png')

# create the scatter plot.
fig, ax = plt.subplots()
fig.set_tight_layout(True)
    
# Query the figure's on-screen size and DPI. 
# Note that when saving the figure to a file, 
# we need to provide a DPI for that separately.
print('fig size: {0} DPI, size in inches {1}'.format(fig.get_dpi(), fig.get_size_inches()))
    
# Plot the set of scatter points

#plt.axis([0, 200, 0, 1])
plt.yscale('log')
plt.ylabel('Nuclear Density, [# of nuclei/cm^3]')
plt.title('Nuclear Density of Xe-140 Daughter Products')
plt.grid()
legend = plt.legend()
plt.xlim(0,timespan)

lineXe, = ax.plot(time[0:2],Xe[0:2],label = 'Xe-140', c = 'g')
lineCs, = ax.plot(time[0:2],Cs[0:2],label = 'Cs-140', c = 'c')
lineBa, = ax.plot(time[0:2],Ba[0:2],label = 'Ba-140', c = 'b')
lineLa, = ax.plot(time[0:2],La[0:2],label = 'La-140', c = 'r')
lineCe, = ax.plot(time[0:2],Ce[0:2],label = 'Ce-140', c = 'm')
# This function updates the plot when it is called to generate
# each frame requested by the FuncAnimation method
def update(i):
    label = 'Timestep {0} [s]'.format(i*2-2)
    print(label)
    # Update the axes (with a new xlabel). Return a tuple of
    # "artists" that have to be redrawn for this frame.
    ax.set_xlabel(label)

    plt.legend()
    lineXe.set_xdata(time[0:i+1])
    lineXe.set_ydata(Xe[0:i+1])
    lineCs.set_xdata(time[0:i+1])
    lineCs.set_ydata(Cs[0:i+1])
    lineBa.set_xdata(time[0:i+1])
    lineBa.set_ydata(Ba[0:i+1])
    lineLa.set_xdata(time[0:i+1])
    lineLa.set_ydata(La[0:i+1])
    lineCe.set_xdata(time[0:i+1])
    lineCe.set_ydata(Ce[0:i+1])

    return lineXe,lineCs,lineBa,lineLa,lineCe,ax

# Create an animation object from the created figure that includes

anim = FuncAnimation(fig, update, frames=len(Xe)+1, interval=100, repeat_delay = 10000)

# save a gif of the animation using the writing package from magick
anim.save('Model_Gif.gif', dpi=72, writer='imagemagick')
