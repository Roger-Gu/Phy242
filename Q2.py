import numpy as np
import math
from scipy.integrate import dblquad
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

# the integrand of the result for x-component
def inte_x(y, x, p):
    q = Point(x,y,0)
    return disX(p, q)/(dis(q, p)**3)

# the integrand of the result for y-component
def inte_y(y, x, p):
    q = Point(x,y,0)
    return disY(p, q)/(dis(q, p)**3)

# the integrand of the result for z-component
def inte_z(y, x, p):
    q = Point(x,y,0)
    return disZ(p, q)/(dis(q, p)**3)

# the x-component of E at p from the rectangle of (a,b), ignoring the constants
def Ex(a, b, p):
    return dblquad(inte_x, -a/2, a/2, lambda x: -b/2, lambda x: b/2, args=(p,))

# the y-component of E at p from the rectangle of (a,b), ignoring the constants
def Ey(a, b, p):
    return dblquad(inte_y, -a/2, a/2, lambda x: -b/2, lambda x: b/2, args=(p,))

# the z-component of E at p from the rectangle of (a,b), ignoring the constants
def Ez(a, b, p):
    return dblquad(inte_z, -a/2, a/2, lambda x: -b/2, lambda x: b/2, args=(p,))

# get the field at p generated from the rectangle of (a,b) and density delta
# require: p is not on the rectangle, a,b >0
def rectangleField(p: Point, a, b, delta, slient = False):
    if (p.z == 0 and p.x <= a/2 and p.x >= -a/2 and p.y <= b/2 and p.y >= -b/2):
        print("Error: the point is on the rectangle")
        return 0
    C = delta * k
    E_x = Ex(a, b, p)[0]
    E_y = Ey(a, b, p)[0]
    E_z = Ez(a, b, p)[0]
    E_size = math.sqrt(E_x**2 + E_y**2 + E_z**2) * C
    E_x *= C
    E_y *= C
    E_z *= C
    # print the result 
    if not slient:
        print("The rectangle is {}*{}, the field at point ({},{},{}) is {}x, {}y, {}z, with size {}".format(a, b, p.x, p.y, p.z, E_x, E_y, E_z, E_size))
    return E_x, E_y, E_z, E_size

# Q3
# try some examples
# p = (4,3,0), a = 1, b = 1, delta = 1
delta = 1
p = Point(4,3,0)
a = 1
b = 1
rectangleField(p, a, b, delta)

# p = (4,3,3), a = 1, b = 1, delta = 1
delta = 1
p = Point(4,3,3)
a = 1
b = 1
rectangleField(p, a, b, delta)

# p = (4,0,3), a = 1, b = 1, delta = 1
delta = 1
p = Point(4,0,3)
a = 1
b = 1
rectangleField(p, a, b, delta)

# p = (0,3,7), a = 1, b = 1, delta = 1
delta = 1
p = Point(0,3,7)
a = 110
b = 123
rectangleField(p, a, b, delta)

# Q4 Compare with the formula on central axis, with height h
# for the sake of simplicity, set the constant delta to 1
delta = 1
def planeField(delta, h):
    E_z = delta/(2*epsilon_0)
    print("The field at point ({},{},{}) is {}, {}, {}".format(0, 0, h, 0, 0, E_z))
    return E_z

# large plane
h = 1
p = Point(0,0,h)
a = 10000
b = 10000
__, __, __, E1 = rectangleField(p, a, b, delta)
Ez2 = planeField(delta, h)
print("the relative difference is {}".format(abs((E1-Ez2)/Ez2)))

# smaller plane
h = 1
p = Point(0,0,h)
a = 100
b = 100
__, __, __, E1 = rectangleField(p, a, b, delta)
Ez2 = planeField(delta, h)
print("the relative difference is {}".format(abs((E1-Ez2)/Ez2)))

# plot symmetry
ax = plt.figure().add_subplot(projection='3d')
# Make the grid
x, y, z = np.meshgrid(np.arange(-10, 10, 2),
                      np.arange(-10, 10, 2),
                      np.arange(-10, 10 ,4))

# wrapper for diskField, but take constant delta=1, a = 5, b=7, and x,y,z instead of a point p
def recField1(x,y,z):
    p = Point(x,y,z)
    return rectangleField(p,5,7,1,True)

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
            E_x, E_y, E_z, __ = recField1(x[i][j][k],y[i][j][k],z[i][j][k])
            u[-1][-1].append(E_x)
            v[-1][-1].append(E_y)
            # Ez is way too large in comparison to the other two, divide by constant just to show symmetry
            w[-1][-1].append(E_z)


ax.quiver(x, y, z, u, v, w, length=0.7, normalize=False)

plt.show()