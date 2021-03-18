from astropy.constants import G, M_sun, M_earth #(m3 / (kg s2) and kg)
from numpy import array, sqrt, arange
import matplotlib.pyplot as plt
   
#function that calculates the radius with x and y         
def r_func(x,y):
    return sqrt(x**2 + y**2)

#function calculating the RHS of all 4 1st order ODEs
def f(r, v, t):
    x = r[0]
    y = r[1]
    v_x = v[0]
    v_y= v[1]
    radius = r_func(x,y)
    fx = v_x
    fv_x = -GM*x/radius**3
    fy = v_y
    fv_y = -GM*y/radius**3
    
    return array([fx, fy], float), array([fv_x, fv_y], float)

#function appending each part of vector r and v to the correct list
def append_all(r, v):
    x_points.append(r[0])
    y_points.append(r[1])
    vx_points.append(v[0])
    vy_points.append(v[1])
    
#Using astropy to define constants
GM = G.value*M_sun.value
m = M_earth.value

x_points = []
y_points = []
vx_points = []
vy_points = []

#Definition of intial values for all variables
x_0 = 1.4710e11 #intial x posirion (m)
y_0 = 0 #inital y position (m)
vx_0 = 0 #initial x velocity (m/s)
vy_0 = 3.0287e4 #intial y velocity (m/s)

r = array([x_0, y_0], float)  #array containing x and y position
v = array([vx_0, vy_0], float) #array containing x and y velocity
a = 0 #intial time (s)
b = 31536000*2 #final time (s) (3 years)
#N = 100
h = 5
tpoints = arange(a,b,h) 

for t in tpoints:
    append_all(r,v)
    k1 = h*f(r,v,t)
    k2 = h*f(r+0.5*k1[0], v+0.5*k1[1], t+0.5*h)
    r += k2[0]
    v += k2[1]

#Plot of x position and y position in function of time
plt.plot(tpoints, x_points, color = "k", label = "Position en x")
plt.plot(tpoints, y_points, color = "r", label = "Position en y")
plt.xlabel("temps (s)")
plt.ylabel("Position (m)")
plt.legend()
plt.show()

#Plot of x position in function of y position
plt.plot(x_points, y_points, "k")
plt.plot(0,0,"ko")
plt.xlabel("Position en x (m)")
plt.ylabel("Position en y (m)")
plt.show()