import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, '../')

import particle.particle as pt
import forces.forces as fr
import solver.solver as sol


C_DRAG = .5 * pt.grav_


def analytic_solution(tf, v0, a0):
    """This expression assumes x0 = 0 and y0 = 0 and m = 1"""
    v0x = v0 * np.cos(np.radians(a0))
    v0y = v0 * np.sin(np.radians(a0))
    xa = v0x * (1. - np.exp(-C_DRAG * tf)) / C_DRAG
    ya = (v0y + pt.grav_ / C_DRAG) * (1. - np.exp(-C_DRAG * tf)) \
       - pt.grav_ * tf
    return xa, ya / C_DRAG


def free_falling(state, params):
    xp, yp, vxp, vyp, _ = state
    grav = params
    axp, ayp = 0., -grav
    return vxp, vyp, axp, ayp, 1.


def linear_drag(state, params):
    xp, yp, vxp, vyp, _ = state
    mass, grav, drag = params
    axp = -drag * vxp / mass
    ayp = -grav - drag * vyp / mass
    return vxp, vyp, axp, ayp, 1.


def quadratic_drag(state, params):
    xp, yp, vxp, vyp, _ = state
    mass, grav, drag = params
    vtp = np.sqrt(vxp**2 + vyp**2)
    axp = -drag * vtp * vxp / mass
    ayp = -grav - drag * vtp * vyp / mass
    return vxp, vyp, axp, ayp, 1.


def integrate(obj):
    xpos, ypos, tpos = [], [], []

    while True:
        xc, yc, _, _, tc = obj.objs.get_state()
        if yc < 0: break
        xpos.append(xc)
        ypos.append(yc)
        tpos.append(tc)
        obj.do_step()

    return xpos, ypos, tpos


# BEGINNING-OF-EXECUTION
deltat, num_method = .005, "Midpoint"
m, x0, y0, v0, a0 = 1., 0., .0, 1., 75

# create free-falling Particle
sim_params = pt.grav_    # gravity
vacuum = pt.Particle(x0, y0, v0, a0)
the_force = fr.Forces(free_falling, sim_params)
vacuum.set_force(the_force)
intgr = sol.Solver(vacuum, num_method, deltat)
xvac, yvac, tvac = integrate(intgr)

# create falling Particle with linear drag
sim_params = m, pt.grav_, C_DRAG    # mass, gravity, drag
linear = pt.Particle(x0, y0, v0, a0)
the_force = fr.Forces(linear_drag, sim_params)
linear.set_force(the_force)
intgr = sol.Solver(linear, num_method, deltat)
xlin, ylin, tlin = integrate(intgr)

# create falling Particle with quadratic drag
sim_params = m, pt.grav_, C_DRAG    # mass, gravity, drag
qdatic = pt.Particle(x0, y0, v0, a0)
the_force = fr.Forces(quadratic_drag, sim_params)
qdatic.set_force(the_force)
intgr = sol.Solver(qdatic, num_method, deltat)
xqua, yqua, tqua = integrate(intgr)

# generate exact evolution for linear drag
rana = [analytic_solution(t, v0, a0) for t in tlin]
xana, yana = map(list, zip(*rana))

# generate plots
fig, ax = plt.subplots()
ax.plot(xvac, yvac, label='Free falling')
ax.plot(xlin, ylin, label='Linear drag')
ax.plot(xqua, yqua, label='Quadratic drag')
ax.plot(xana, yana, '--', label='Analytic linear')

ax.set(xlabel='x (a.u.)', ylabel='y (a.u.)',
       title='Projectile motion with different drags')
ax.grid()

plt.legend()
plt.show()
# END-OF-EXECUTION
