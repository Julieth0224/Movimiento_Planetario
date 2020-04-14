
import sys
sys.path.insert(0, '../../')

import numpy as np
import matplotlib.pyplot as plt

import particle.particle as pt
import solver.solver as sl
import forces.forces as fr

GM = 4 * (np.pi)**2

def Newton_again(state, params):
    x, y, vx, vy, t = state
    r = (x**2 + y**2)**0.5
    ax = (-GM*x)/(r**3)
    ay = (-GM*y)/(r**3)
    return vx, vy, ax, ay, 1

def area(deltat, x, y, vx, vy):
    a = (1/2) * deltat * (x * vy - y * vx)
    return a

def integrate(obj):
    cont = 0
    xpos, ypos, tpos, areas = [], [], [], []
    tc = 0
    for i in range(200000):
        xc, yc, vxc, vyc, tc = obj.objs.get_state()
        ar = round(area(deltat, xc, yc, vxc, vyc), 5)
        areas.append(ar)
        xpos.append(xc)
        ypos.append(yc)
        tpos.append(tc)
        obj.do_step()
        xc1, yc1, vxc1, vyc1, tc1 = obj.objs.get_state()
        if yc < 0 and yc1 >= 0:
            cont += 1
    print("Orbitas", cont)
    return xpos, ypos, tpos, areas


sim_params = pt.grav_
deltat = 0.01/256

m1 = "Euler"
m2 = "Euler-Cromer"
m3 = "Midpoint"

num_method = m2 # m1 # m3

# k = 2 * np.pi
m, x0, y0, v0, a0 = 1., 1.5, .0, 1.5, 90. # v0 = k

earth = pt.Particle(x0, y0, v0, a0)
earth_force = fr.Forces(Newton_again, sim_params)
earth.set_force(earth_force)
euler = sl.Solver(earth, num_method, deltat)
xvac, yvac, tvac, area = integrate(euler)


print("Area Ini:", area[0], ", Area Fin:", area[-1])

# Muestra la orbita circular o eliptica
# fig, ax = plt.subplots()
# ax.plot(xvac, yvac, '--', c = 'purple', label=num_method)
#
# ax.set(xlabel='x (a.u.)', ylabel='y (a.u.)')

# Muestra grafica area vs tiempo
fig, ax = plt.subplots()
ax.plot(tvac, area, '--', c = 'purple', label=num_method)

ax.set(xlabel='t (yr)', ylabel="A (a.u.**2)")


ax.grid()
fig = plt.gcf()
plt.legend()
fig.set_size_inches(10., 5.)
plt.show()
