import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, '../')

import particle.particle as pt
import forces.forces as fr
import solver.solver as sol
import animator.animator as ani

GM = 4*(np.pi**2)

def Newton_again(state, params):
    x, y, vx, vy, t = state
    # G, M = params
    r = (x**2 + y**2)**0.5
    ax = (-GM*x)/(r**3)
    ay = (-GM*y)/(r**3)
    return vx, vy, ax, ay, 1

def integrate(obj):
    xpos, ypos, tpos = [], [], []
    tc = 0
    while tc < 2:
        xc, yc, _, _, tc = obj.objs.get_state()
        xpos.append(xc)
        ypos.append(yc)
        tpos.append(tc)
        obj.do_step()

    return xpos, ypos, tpos

sim_params = pt.grav_
deltat = 0.01/1024

m1 = "Euler"
m2 = "Euler-Cromer"
m3 = "Midpoint"

num_method = m2 # m1 # m3

m, x0, y0, v0, a0 = 1., 1.5, 0., 1.5, 90. # x0 = 1.0, v0 = 1.0

earth = pt.Particle(x0, y0, v0, a0)
f1 = fr.Forces(Newton_again, sim_params)
earth.set_force(f1)
i1 = sol.Solver(earth, num_method, deltat)
xvac, yvac, tvac = integrate(i1)

m2, x02, y02, v02, a02 = 1., 1.6, 0., 1.6, 90. # x0 = 1.1, v0 = 1.1

mars = pt.Particle(x02, y02, v02, a02)
f2 = fr.Forces(Newton_again, sim_params)
mars.set_force(f2)
i2 = sol.Solver(mars, num_method, deltat)
xvac2, yvac2, tvac2 = integrate(i2)

m3, x03, y03, v03, a03 = 1., 1.7, 0., 1.7, 90. # x0 = 1.2, v0 = 1.2

mercury = pt.Particle(x03, y03, v03, a03)
f3 = fr.Forces(Newton_again, sim_params)
mercury.set_force(f3)
i3 = sol.Solver(mercury, num_method, deltat)
xvac3, yvac3, tvac3 = integrate(i3)

m4, x04, y04, v04, a04 = 1., 1.8, 0., 1.8, 90. # x0 = 1.3, v0 = 1.3

venus = pt.Particle(x04, y04, v04, a04)
f4 = fr.Forces(Newton_again, sim_params)
venus.set_force(f4)
i4 = sol.Solver(venus, num_method, deltat)
xvac4, yvac4, tvac4 = integrate(i4)

m5, x05, y05, v05, a05 = 1., 1.9, 0., 1.9, 90. # x0 = 1.4, v0 = 1.4

Jupiter = pt.Particle(x05, y05, v05, a05)
f5 = fr.Forces(Newton_again, sim_params)
Jupiter.set_force(f5)
i5 = sol.Solver(Jupiter, num_method, deltat)
xvac5, yvac5, tvac5 = integrate(i5)


# Muestra las orbitas de cinco planetas

sx, sy = 0, 0
ex, ey = 1.5, 0. # 1.0, 0.
x2, y2 = 1.6, 0 # 1.1, 0.
x3, y3 = 1.7, 0 # 1.2, 0.
x4, y4 = 1.8, 0 # 1.3, 0.
x5, y5 = 1.9, 0 # 1.4, 0.

fig, ax = plt.subplots()
ax.plot(xvac, yvac, c = 'blue', label='Earth')
ax.plot(xvac2, yvac2, c = 'red', label='Earth')
ax.plot(xvac3, yvac3, c = 'green', label='Earth')
ax.plot(xvac4, yvac4, c = 'black', label='Earth')
ax.plot(xvac5, yvac5, c = 'pink', label='Earth')

ax.plot(sx, sy, 'o-', c = 'orange', label='Sun')
ax.plot(ex, ey, 'o-', c='blue')
ax.plot(x2, y2, 'o-', c='red')
ax.plot(x3, y3, 'o-', c='green')
ax.plot(x4, y4, 'o-', c='black')
ax.plot(x5, y5, 'o-', c='pink')

ax.set(xlabel='x (a.u.)', ylabel='y (a.u.)')
ax.grid()

plt.legend()
plt.show()
