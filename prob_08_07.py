# start from Example 8.5 program odesim.py
import numpy as np
import matplotlib.pyplot as plt

# constants
g     = 3.71	# m s^-2
m     = 1.0	# kg
rho   = 0.20	# kg m^-3
C     = 0.47	# unitless
R     = 0.08    # m
h     = 0.001   # seconds
theta = 45.0*(np.pi/180) # radians
v0    = 100.0	# m s^-1
const = (rho*C*np.pi*R**2)/(2.0*m)

# define the equations of motion
def f(r,const):
    x   = r[0]
    y   = r[1]
    vx  = r[2]
    vy  = r[3]
    fx  = vx
    fy  = vy
    fvx = -const*vx*np.sqrt(vx**2+vy**2)
    fvy = -g-const*vy*np.sqrt(vx**2+vy**2)
    return np.array([fx,fy,fvx,fvy],float)

# containers for output
r = np.array([0.0,0.0,v0*np.cos(theta),v0*np.sin(theta)],float)
xpoints = []
ypoints = []

# use fourth-order Runge-Kutta
while r[1]>=0:
    k1 = h*f(r,const)
    k2 = h*f(r+0.5*k1,const)
    k3 = h*f(r+0.5*k2,const)
    k4 = h*f(r+k3,const)
    r += (k1+2*k2+2*k3+k4)/6
    xpoints.append(r[0])
    ypoints.append(r[1])

# make plot for part (b)
p1 = plt.figure(1)
plt.plot(xpoints,ypoints)
plt.xlabel("x [m]")
plt.ylabel('y [m]')
p1.show()

# try different values of m
p2 = plt.figure(2)
for m in 0.4536*np.array([18,24,32,42]):
    const = (rho*C*np.pi*R**2)/(2.0*m)
    r = np.array([0.0,0.0,v0*np.cos(theta),v0*np.sin(theta)],float)
    xpoints = []
    ypoints = []

    # use fourth-order Runge-Kutta
    while r[1]>=0:
        k1 = h*f(r,const)
        k2 = h*f(r+0.5*k1,const)
        k3 = h*f(r+0.5*k2,const)
        k4 = h*f(r+k3,const)
        r += (k1+2*k2+2*k3+k4)/6
        xpoints.append(r[0])
        ypoints.append(r[1])

    plt.plot(np.array(xpoints)/1000,np.array(ypoints)/1000,label='m = '+str(m)+' kg')
    
    print('Max Horizontal distance (' + str(m) + ' kg): ' + str(xpoints[-1]/1000) + ' km')
    print('Max Vertical distance (' + str(m) + ' kg): ' + str(max(ypoints)/1000) + ' km')

plt.xlabel("x [km]")
plt.ylabel('y [km]')
plt.legend()
p2.show()
