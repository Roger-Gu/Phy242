import numpy as np
import math
# from scipy.integrate import dblquad
from integeration import dblquad
import matplotlib.pyplot as plt

epsilon_0 = 8.8541878128*(10**-12)
k = 1 / (4* np.pi * epsilon_0)

# define a point class
class Point:
    def __init__(self, x: float, y: float, z: float):
        self.x = x
        self.y = y
        self.z = z

# calculate x-distance between two points
def disX(p1: Point, p2: Point):
    return p1.x-p2.x

# calculate y-distance between two points
def disY(p1: Point, p2: Point):
    return p1.y-p2.y

# calculate z-distance between two points
def disZ(p1: Point, p2: Point):
    return p1.z-p2.z

# calculate distance between two points
# reqiure: p1 != p2
def dis(p1: Point, p2: Point):
    return math.sqrt((p1.x-p2.x)**2+(p1.y-p2.y)**2+(p1.z-p2.z)**2)

# get the field at point p2, generated by charge q as p1, ignore constants pi and epsilon_0 
def pointField(p1, p2, q):
    Ex = q*disX(p2, p1)/(dis(p1, p2)**3)
    Ey = q*disY(p2, p1)/(dis(p1, p2)**3)
    Ez = q*disY(p2, p1)/(dis(p1, p2)**3)
    return Ex, Ey, Ez

# the integrand of the result for x-component
def inte_x(y, x, p):
    q = Point(x,y,0)
    return disX(p, q)/(dis(q, p)**3)

# the x-component of E at p from the disk with raduis R, ignoring the constants
def Ex(R, p):
    return dblquad(inte_x, -R, R, lambda x: -math.sqrt(R**2-x**2), lambda x: math.sqrt(R**2-x**2), args=(p,))

# the integrand of the result for y-component
def inte_y(y, x, p):
    q = Point(x,y,0)
    return disY(p, q)/(dis(q, p)**3)

# the y-component of E at p from the disk with raduis R, ignoring the constants
def Ey(R, p):
    return dblquad(inte_y, -R, R, lambda x: -math.sqrt(R**2-x**2), lambda x: math.sqrt(R**2-x**2), args=(p,))

# the integrand of the result for z-component
def inte_z(y, x, p):
    q = Point(x,y,0)
    return disZ(p, q)/(dis(q, p)**3)

# the z-component of E at p from the disk with raduis R, ignoring the constants
def Ez(R, p):
    return dblquad(inte_z, -R, R, lambda x: -math.sqrt(R**2-x**2), lambda x: math.sqrt(R**2-x**2), args=(p,))

# get the field at p generated by disk with radius R and density delta
# require: p is not on the disk
def diskField(p: Point, R, delta, slient = False):
    if (p.z == 0 and dis(p,Point(0,0,0)) <= R):
        print("Error: the point is on the disk")
        return 0
    C = delta * k
    E_x = Ex(R, p)
    E_y = Ey(R, p)
    E_z = Ez(R, p)
    E_size = math.sqrt(E_x**2 + E_y**2 + E_z**2) * C
    E_x *= C
    E_y *= C
    E_z *= C
    # print the result 
    if not slient:
        print("The Radius is {}, the field at point ({},{},{}) is {}x, {}y, {}z, with size {}".format(R, p.x, p.y, p.z, E_x, E_y, E_z, E_size))
    return E_x, E_y, E_z, E_size

# Q3
# try some examples
# p = (1,0,0), R = 1, delta = 1
delta = 1
p = Point(1,0,0)
R = 1
diskField(p, R, delta)
# p = (1,1,0), R = 1, delta = 1
delta = 1
p = Point(1,1,0)
R = 1
diskField(p, R, delta)
# p = (1,0,1), R = 1, delta = 1
delta = 1
p = Point(1,0,1)
R = 1
diskField(p, R, delta)
# p = (0,1,1), R = 1, delta = 1
delta = 1
p = Point(0,1,1)
R = 1
diskField(p, R, delta)
# p = (0,0,1), R = 2, delta = 1
delta = 1
p = Point(0,0,1)
R = 2
diskField(p, R, delta)
# p = (1,1,1), R = 1, delta = 80
delta = 80
p = Point(1,1,1)
R = 1
diskField(p, R, delta)
# p = (0,0,-1), R = 2, delta = 1
delta = 1
p = Point(0,0,-1)
R = 2
diskField(p, R, delta)


# Q4 Compare with the formula on central axis, with height h
# for the sake of simplicity, set the constant k and delta to 1
k = 1
delta = 1
def centralField(R, delta, h):
    E_z = 2*np.pi*k*delta*(1-h/math.sqrt(h**2+R**2))
    print("The Radius is {}, the field at point ({},{},{}) is {}, {}, {}".format(R, 0, 0, h, 0, 0, E_z))
    return E_z

# h = 10, R = 1
h = 10
p = Point(0,0,h)
R = 1
__, __, __, E1 = diskField(p, R, delta)
Ez2 = centralField(R, delta, h)
print("the relative difference is {}".format(abs((E1-Ez2)/Ez2)))

# h = 9999, R = 1
h = 9999
p = Point(0,0,h)
R = 1
__, __, __, E1 = diskField(p, R, delta)
Ez2 = centralField(R, delta, h)
print("the relative difference is {}".format(abs((E1-Ez2)/Ez2)))

# plot symmetry
ax = plt.figure().add_subplot(projection='3d')
# Make the grid
x, y, z = np.meshgrid(np.arange(-10, 10, 2),
                      np.arange(-10, 10, 2),
                      np.arange(-10, 10 ,4))

# wrapper for diskField, but take constant delta=1, R =5, and x,y,z instead of a point p
def disField1(x,y,z):
    p = Point(x,y,z)
    return diskField(p,5,1,True)

u,v,w = [], [], []
# plot symmetry
ax = plt.figure().add_subplot(projection='3d')
# Make the grid
x, y, z = np.meshgrid(np.arange(-10, 10, 2),
                      np.arange(-10, 10, 2),
                      np.arange(-10, 10 ,4))

# wrapper for diskField, but take constant delta=1, R =5, and x,y,z instead of a point p
def disField1(x,y,z):
    p = Point(x,y,z)
    return diskField(p,5,1,True)

u,v,w = [], [], []
for i in np.arange(len(x)):
    u.append([])
    v.append([])
    w.append([])
    for j in np.arange(len(x[0])):
        u[-1].append([])
        v[-1].append([])
        w[-1].append([])
        for k in np.arange(len(x[0][0])):
            # after investigation on the meshgrid function, it very wierdly is in the form y-x-z loop
            E_x, E_y, E_z, __ = disField1(x[i][j][k],y[i][j][k],z[i][j][k])
            u[-1][-1].append(E_x)
            v[-1][-1].append(E_y)
            # Ez is way too large in comparison to the other two, divide by constant just to show symmetry
            w[-1][-1].append(E_z/2)

ax.quiver(x, y, z, u, v, w, length=0.7, normalize=False)
plt.show()
