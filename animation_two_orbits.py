import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, '../')

import particle.particle as pt
import forces.forces as fr
import solver.solver as sol
import animator.animator as ani

GM = 4*(np.pi**2)

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
    xpos, ypos, tpos = [], [], []
    tc = 0
    while tc < 2:
        xc, yc, _, _, tc = obj.objs.get_state()
        xpos.append(xc)
        ypos.append(yc)
        tpos.append(tc)
        obj.do_step()

    return xpos, ypos, tpos

v0, vy = init_c(1, 0, 0)
sim_params = pt.grav_
m, x0, y0, v0, a0 = 1., 1., 0., 5., 70
deltat = .0001

m1 = "Euler"
m2 = "Euler-Cromer"
m3 = "Midpoint"

num_method = m2

earth = pt.Particle(x0, y0, v0, a0)
earth_force = fr.Forces(Newton_again, sim_params)
earth.set_force(earth_force)
euler = sol.Solver(earth, num_method, deltat)
xvac, yvac, tvac = integrate(euler)

m2, x02, y02, v02, a02 = 1., 1.5, 0., 5., 70

mars = pt.Particle(x02, y02, v02, a02)
f2 = fr.Forces(Newton_again, sim_params)
mars.set_force(f2)
i2 = sol.Solver(mars, num_method, deltat)
xvac2, yvac2, tvac2 = integrate(i2)

sx, sy = 0, 0
ex, ey = 1, 0.
mx, my = 1.5, 0

# Muestra las orbitas de los dos planetas
fig, ax = plt.subplots()
ax.plot(xvac, yvac, c = 'blue', label='Earth')
ax.plot(xvac2, yvac2, c = 'red', label='Mars')
ax.plot(sx, sy, 'o-', c = 'orange', label='Sun')
ax.plot(ex, ey, 'o-', c='blue')
ax.plot(mx, my, 'o-', c='red')

ax.set(xlabel='x (a.u.)', ylabel='y (a.u.)',
       title='Projectile motion with different drags')
ax.grid()

plt.legend()
plt.show()

# Crea la animacion de los dos planetas alrededor del sol
earth = (xvac, yvac)
mars = (xvac2, yvac2)
sun = ([0], [0])
anime = ani.Animator((earth, sun, mars))
anime.setup_anime()
anime.run_anime()
