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
NN = [[1,0,0,0,0]]
T = 0
for i in range(0,100):
    row = []
    T += 2 # Increment time, [s]
    for j in range(0,5):
        SUMM = 0
        for k in range(0,j+1):
            SUMM += COEFF[j][k]*np.exp((-1)*LL[k]*T)
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

time = [*range(0,201,2)]
# =============================================================================
# plt.figure()
# 
# plt.plot(time,Xe,label = 'Xe-140')
# 
# plt.plot(time,Cs,label = 'Cs-140')
# 
# plt.plot(time,Ba,label = 'Ba-140')
# 
# plt.plot(time,La,label = 'La-140')
# 
# plt.plot(time,Ce,label = 'Ce-140')
# 
# plt.yscale('log')
# plt.xlabel('Time, [s]')
# plt.ylabel('Nuclear Density, [# of nuclei/cm^3]')
# plt.title('Nuclear Density of Xe-140 Daughter Products')
# plt.legend()
# plt.grid()
# plt.savefig('Model_Graph.png')
# =============================================================================

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
plt.xlim(0,200)

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
    
#    ax.plot(time[0:i+1],Xe[0:i+1],label = 'Xe-140', c = 'g')
#    
#    ax.plot(time[0:i+1],Cs[0:i+1],label = 'Cs-140', c = 'c')
#    
#    ax.plot(time[0:i+1],Ba[0:i+1],label = 'Ba-140', c = 'b')
#    
#    ax.plot(time[0:i+1],La[0:i+1],label = 'La-140', c = 'r')
#    
#    ax.plot(time[0:i+1],Ce[0:i+1],label = 'Ce-140', c = 'm')
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
# =============================================================================
# line.set_xdata(x[0:i])    
# line.set_ydata(y_values(x[0:i],a,b))
# =============================================================================
    return lineXe,lineCs,lineBa,lineLa,lineCe,ax

# Create an animation object from the created figure that includes

anim = FuncAnimation(fig, update, frames=len(Xe)+1, interval=100, repeat_delay = 5000)

# save a gif of the animation using the writing package from magick
anim.save('Model_Gif.gif', dpi=72, writer='imagemagick')
