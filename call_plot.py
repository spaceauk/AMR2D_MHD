import numpy as np
import matplotlib.pyplot as plt
from os import listdir
import sys

data_path="./data/"
plot_path="./plots/"
allfiles=listdir(data_path)
allfiles_sorted=sorted(allfiles)

nx=int(sys.argv[1])
ny=int(sys.argv[2])
plotvar=int(sys.argv[3])
ICtype=sys.argv[4]
if plotvar<0 or plotvar>7:
    print("Error as wrong variable type selected!")
    print("Default variable of density taken instead.")
    plotvar=0
    # sys.exit()
print("Begin plotting on Python script: nx="+str(nx)+" ny="+str(ny))

def choose_plotvar(x):
    return {
            0:'Density',
            1:'U_vel',
            2:'V_vel',
            3:'W_vel',
            4:'Pressure',
            5:'Bx',
            6:'By',
            7:'Bz',
            }.get(x,'Unknown')

if plotvar==5 and ICtype=="VST":
    pvar='T'
else:
    pvar=choose_plotvar(plotvar)

frame_no=0
for fil_name in allfiles_sorted:
    fname=data_path+fil_name
    data=np.loadtxt(fname)

# extract x and y coordinates
sorted_indices = np.argsort(data[:,0])
sorted_x = data[:,0][sorted_indices]
sorted_y = data[:,1][sorted_indices]
x = data[:, 0]
y = data[:, 1]

# create plot for mesh
fig, ax = plt.subplots()
ax.set_aspect("equal")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_title("Mesh Grid")
ax.grid(True)

for i in range(ny):
    start = i * nx +1
    end = start + nx 
    ax.plot(x[start:end], y[start:end], "r-", lw=0.5)
    
for i in range(nx):
    ax.plot(x[i::nx], y[i::nx], "r-", lw=0.5)

# Linear plot
plt.figure()
plt.subplot(1,2,1)
sorted_indices = np.argsort(data[:,0])
sorted_X = data[:,0][sorted_indices]
sorted_Z = data[:,plotvar+4][sorted_indices]
plt.plot(sorted_X,sorted_Z,".")
plt.grid()
plt.subplot(1,2,2)
sorted_indices = np.argsort(data[:,1])
sorted_Y = data[:,1][sorted_indices]
sorted_Z = data[:,plotvar+4][sorted_indices]
plt.plot(sorted_Y,sorted_Z,".")
plt.grid()
#plt.figure()
#plt.scatter(data[:,0],data[:,1])
#plt.grid()
plt.show()

