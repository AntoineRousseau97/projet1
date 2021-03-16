from astropy.constants import G, M_sun #(m3 / (kg s2) and kg)
from numpy import array, sqrt, arange
import matplotlib.pyplot as plt
           
def r_func(x,y):
    return sqrt(x**2 + y**2)

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

def append_all(r, v):
    x_points.append(r[0])
    y_points.append(r[1])
    vx_points.append(v[0])
    vy_points.append(v[1])
    

GM = G.value*M_sun.value

x_points = []
y_points = []
vx_points = []
vy_points = []

x_0 = 1.4710e11 #m
y_0 = 0 #m
vx_0 = 0 #m/s
vy_0 = 3.0287e4 #m/s 

r = array([x_0, y_0], float)  
v = array([vx_0, vy_0], float)
a = 0
b = 31536000*3 # 3 years
h = 3600 # 1 hour in seconds
tpoints = arange(a,b,h)
v_int = v + 0.5*h*f(r,v,0)[1]

for t in tpoints:
    append_all(r,v)
    r += h*v_int
    k = h*f(r,v,t)[1]
    v = v_int + 0.5*k
    v_int += k

    
plt.plot(tpoints, x_points, "r")
plt.show()
