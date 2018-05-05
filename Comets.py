"""
Steven Clark
Ashley Piccone

This code animates the motion of bodies orbitting a center.

Set usePlanets on line 69 to true to show our solar system. set to false and will solve for N random bodies.

NOTE: The code will not actually finish running due to the save animation function, so you will have to ctrl-c it after the animation stops
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random 
from mpl_toolkits.mplot3d import Axes3D

# function to calculate evolution of bodies using discretized equations of motion
def orb_mech(pos,vel):
    for ii in range(0, steps-1):
        
        r = np.sqrt(pos[0][ii]**2 + pos[1][ii]**2 + pos[2][ii]**2)
        
        vel[0][ii+1] = dt * (-mu*pos[0][ii] / r**3) + vel[0][ii]
        vel[1][ii+1] = dt * (-mu*pos[1][ii] / r**3) + vel[1][ii]
        vel[2][ii+1] = dt * (-mu*pos[2][ii] / r**3) + vel[2][ii]
        
        pos[0][ii+1] = vel[0][ii]*dt + pos[0][ii]
        pos[1][ii+1] = vel[1][ii]*dt + pos[1][ii]
        pos[2][ii+1] = vel[2][ii]*dt + pos[2][ii]
        
        # plotting energy makes sure it remains constant- bodies are in orbit
        # used for debugging
        #energy[ii] = m * (vel[0][ii]**2 + vel[1][ii]**2 + vel[2][ii]**2) / 2 - mu* m / r
              
    vel[0] = vel[0][:-1]
    vel[1] = vel[1][:-1]
    vel[2] = vel[2][:-1]
    pos[0] = pos[0][:-1]
    pos[1] = pos[1][:-1]
    pos[2] = pos[2][:-1]
    #energy = energy[:-1]
    
    return pos, vel

# define a body from initial conditions
def Planet(xi,yi,zi,vxi,vyi,vzi,steps):
    pos = [np.zeros(steps),np.zeros(steps),np.zeros(steps)]
    vel = [np.zeros(steps),np.zeros(steps),np.zeros(steps)]
    
    pos[0][0], pos[1][0], pos[2][0] = xi,yi,zi
    vel[0][0], vel[1][0], vel[2][0] = vxi,vyi,vzi
       
    pos, vel = orb_mech(pos,vel)
    
    return pos, vel

# time steps and dt, these are the necessary values for energy to remain constant, found through plotting energy
dt = 10e-3
finaltime = 10
steps = int(finaltime / dt)
time = np.linspace(0, finaltime, steps)
time = time[:-1] # the last value of the position arrays cause issues, keep lengths constant

G = 6.674e-11 # gravitational constants
M = 1.989e30 # mass of sun kg 
Rsun = 695700000 # radius of sun m
mu = 1
m = 5.972e24 / M # mass of earth / sun

usePlanets = True # decides whether to use N bodies or the data for the planets
N = 3

# initial conditions for the planets in our solar system, from Mercury to Neptune
# data from https://sites.temple.edu/math5061/files/2016/12/final_project.pdf
planetposx = [-4.6e10, -1.0748e11, -1.47095e11, -2.0662e11, -7.4052e11,-1.35255e12,-2.7413e12, -4.44445e12]
planetposy = [0,0,0,0,0,0,0,0]
planetposz = [0,0,0,0,0,0,0,0]
planetvelx = [0,0,0,0,0,0,0,0]
planetvely = [-58980,-35260,-30300,-26500,-13720,-10180,-7110,-5500]
planetvelz = [0,0,0,0,0,0,0,0]

planetmass = [0.33011e24, 4.8675e24, 5.972e24, 6.4171e23, 1898.19e24, 568.34e24, 86.813e24, 102.413e24]

# normalize planets to the Earth's position and size
nplanetposx = np.zeros(len(planetposx))
nplanetposy = np.zeros(len(planetposy))
nplanetposz = np.zeros(len(planetposz))
nplanetvelx = np.zeros(len(planetvelx))
nplanetvely = np.zeros(len(planetvely))
nplanetvelz = np.zeros(len(planetvelz))
nplanetmass = np.zeros(len(planetmass))
sunsize=0
normindex = 2


for ii in range(len(planetposx)):
  nplanetposx[ii] = (planetposx[ii]-sunsize)/ planetposx[normindex] 
  nplanetvely[ii] = planetvely[ii] / planetvely[normindex]
  nplanetmass[ii] = planetmass[ii] / planetmass[normindex]



# create position and velocity arrays for bodies
pos_arr = []
vel_arr = []
tp = []
tv = []

# creates our solar system
if(usePlanets == True):
  visual_masses = [2,4,5,3,18,14,12,12]
  sunsize = 1500 
  for ii in range(0, len(nplanetposx)):
    x = nplanetposx[ii]
    y = nplanetposy[ii]
    z = nplanetposz[ii]
    vx = nplanetvelx[ii]
    vy = nplanetvely[ii]
    vz = nplanetvelz[ii]
    tp, tv = Planet(x,y,z,vx,vy,vz,steps)
    pos_arr.append(tp)
    vel_arr.append(tv)

# creates N planets with random initial conditions and random sizes
elif(usePlanets == False):
  visual_masses = []
  sunsize = 150
  for ii in range(0,N):
    x = random.randint(1,11)/10
    y = random.randint(1,11)/10
    z = random.randint(1,11)/10
    vx = random.randint(1,11)/10
    vy = random.randint(1,11)/10
    vz = random.randint(1,11)/10
    tp, tv = Planet(x,y,z,vx,vy,vz,steps)
    pos_arr.append(tp)
    vel_arr.append(tv)
    visual_masses.append(random.randint(1,11))

else: print("Something went horribly wrong") #If usePlanets is neither true nor false

#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=18000)

fig = plt.figure()
ax = Axes3D(fig)

# creates lines for animation of bodies
line_array = []
for jj in range(0,len(pos_arr)):
    line, = ax.plot(pos_arr[jj][0],pos_arr[jj][1],pos_arr[jj][2],'o-', markersize = visual_masses[jj],zorder = 1)
    line_array.append(line,)

# animation process
def init():
    for kk in range(0,len(pos_arr)):
        line_array[kk].set_data(np.ma.array(pos_arr[kk][0],mask = True), np.ma.array(pos_arr[kk][1],mask = True))
        line_array[kk].set_3d_properties(np.ma.array(pos_arr[kk][2]))
    return line_array

def animate(i):
    for kk in range(0,len(pos_arr)):
        #These three lines plot with line to the center, but also a big dot at the center
        #thisx = [0,pos_arr[kk][0][i]] 
        #thisy = [0,pos_arr[kk][1][i]]
        #thisz = [0,pos_arr[kk][2][i]]

        #These ones have not dot in center or line connecting to sun
        thisx = pos_arr[kk][0][i]
        thisy = pos_arr[kk][1][i]
        thisz = pos_arr[kk][2][i]
        line_array[kk].set_data(thisx,thisy)
        line_array[kk].set_3d_properties(thisz)
    return line_array

xmax, ymax, zmax = 0,0,0
# loop below finds absolute maximums for x,y,z
for ii in range(0, len(pos_arr)):
  if(max([max(pos_arr[ii][0]), np.abs(min(pos_arr[ii][0]))]) > xmax):
    xmax = max([max(pos_arr[ii][0]), np.abs(min(pos_arr[ii][0]))])

  if(max([max(pos_arr[ii][1]), np.abs(min(pos_arr[ii][1]))]) > ymax):
    ymax = max([max(pos_arr[ii][1]), np.abs(min(pos_arr[ii][1]))])

  if(max([max(pos_arr[ii][2]), np.abs(min(pos_arr[ii][2]))]) > zmax):
    zmax = max([max(pos_arr[ii][2]), np.abs(min(pos_arr[ii][2]))])

ani = animation.FuncAnimation(fig,animate,interval = finaltime,blit=True,init_func = init, frames=steps)

ax.set_zlim(-zmax,zmax)
plt.xlim(-xmax,xmax)
plt.ylim(-ymax,ymax)

# create the sun 
ax.scatter(0,0,0,'o', c = 'y', s = sunsize, zorder = 2)


plt.show()
#ani.save("im.mp4", writer=writer)

#Some debugging plots below
#ax.plot(pos_arr[0][0],pos_arr[0][1],pos_arr[0][2])
#plt.show()
#
#plt.clf()
#plt.plot(time, energy, label='Straight Line?')
#plt.legend(loc='best')
#plt.title("Energy vs. Time")
#plt.show()
