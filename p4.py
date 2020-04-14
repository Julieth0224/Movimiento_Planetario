
import sys
sys.path.insert(0, '../../')

import numpy as np
import matplotlib.pyplot as plt

import particle.particle as pt
import solver.solver as sl
import forces.forces as fr

GM = 4 * (np.pi)**2

def dist(a1, a2):
    dis = abs(a2 - a1)
    return dis

def Newton_again(state, params):
    x, y, vx, vy, t = state
    r = (x**2 + y**2)**0.5
    ax = (-GM*x)/(r**3)
    ay = (-GM*y)/(r**3)
    return vx, vy, ax, ay, 1

def integrate(obj):
    cont = 0
    xpos, ypos, tpos = [], [], []
    tc = 0
    for i in range(200000):
        xc, yc, vxc, vyc, tc = obj.objs.get_state()
        xpos.append(xc)
        ypos.append(yc)
        tpos.append(tc)
        obj.do_step()
        xc1, yc1, vxc1, vyc1, tc1 = obj.objs.get_state()
        if yc < 0 and yc1 >= 0:
            cont += 1
    print("Orbitas", cont)
    return xpos, ypos, tpos

def ellipse_data(xpos, ypos):
    xmin = min(xpos)
    xmax = max(xpos)
    ymin = min(ypos)
    ymax = max(ypos)

    sema = dist(xmin, xmax) / 2
    seme = dist(ymin, ymax) / 2

    exc = np.sqrt(sema**2 - seme**2) / sema

    p = np.sqrt(sema ** 3)

    return sema, seme, exc, p


m = 1.
sim_params = pt.grav_
deltat = 0.01/1024

m1 = "Euler"
m2 = "Euler-Cromer"
m3 = "Midpoint"

num_method = m2 # m1 # m3

info = []

x0 = 1.
vy0 = 1.
for i in range(200):
    print("cond ini, x0 & vy0 =", round(x0, 2))
    v0 = vy0
    m, x0, y0, v0, a0 = 1., x0, .0, v0, 90.,

    earth = pt.Particle(x0, y0, v0, a0)
    earth_force = fr.Forces(Newton_again, sim_params)
    earth.set_force(earth_force)
    euler = sl.Solver(earth, num_method, deltat)
    xvac, yvac, tvac = integrate(euler)

    sema, seme, exc, p = ellipse_data(xvac, yvac)

    print("Sema =", sema, ", Seme =", seme)
    print("Exc =", exc, ", p =", p)
    print()

    info.append((round(x0, 2), round(vy0, 2), round(sema, 5), round(p, 5)))

    x0 += 0.008
    vy0 += 0.008

semapos = []
ppos = []
for k in info:
    semapos.append(k[2])
    ppos.append(k[3])


fig, ax = plt.subplots()
ax.plot(ppos, semapos, '--', label=num_method)

ax.set(xlabel='sema (a.u.)', ylabel='p (yr)')


ax.grid()
plt.legend()
plt.show()
