
import sys
sys.path.insert(0, '../../')

import numpy as np
import matplotlib.pyplot as plt

import particle.particle as pt
import solver.solver as sl
import forces.forces as fr

GM = 4 * (np.pi)**2


def init_c(x, y, vx):
    v = np.sqrt(GM / np.sqrt(x**2 + y**2))
    vy = np.sqrt(v**2 - vx**2)
    return float(v), vy

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
        xc, yc, _, _, tc = obj.objs.get_state()
        xpos.append(xc)
        ypos.append(yc)
        tpos.append(tc)
        obj.do_step()
        xc1, yc1, _, _, tc1 = obj.objs.get_state()
        if yc < 0 and yc1 >= 0:
            cont += 1
    print("Orbitas", cont)
    return xpos, ypos, tpos


v0, vy = init_c(1, 0, 0)
sim_params = pt.grav_
m, x0, y0, v0, a0 = 1., 1., .0, v0, 90.
deltat = 0.01/664 # 0.0001

m1 = "Euler"
m2 = "Euler-Cromer"
m3 = "Midpoint"

num_method = m2 # m1 # m3

earth = pt.Particle(x0, y0, v0, a0)
earth_force = fr.Forces(Newton_again, sim_params)
earth.set_force(earth_force)
euler = sl.Solver(earth, num_method, deltat)
xvac, yvac, tvac = integrate(euler)


fig, ax = plt.subplots()
ax.plot(xvac, yvac, '--', label=num_method)

ax.set(xlabel='x (a.u.)', ylabel='y (a.u.)',
       title= 'Metodo: ' + num_method)
ax.grid()

plt.legend()
plt.show()
