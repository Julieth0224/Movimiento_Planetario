import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, '../')

import particle.particle as pt
import forces.forces as fr
import solver.solver as sol

def analytic_solution(tf, v0, a0):
    """This expression assumes x0 = 0 and y0 = 0 and m = 1"""
    v0x = v0 * np.cos(np.radians(a0))
    v0y = v0 * np.sin(np.radians(a0))
    xa = v0x * (1. - np.exp(-pt.drag_ * tf)) / pt.drag_
    ya = (v0y + pt.grav_ / pt.drag_) * (1. - np.exp(-pt.drag_ * tf)) \
       - pt.grav_ * tf
    return xa, ya / pt.drag_

def resistive_falling(state, params):
    xp, yp, vxp, vyp, _ = state
    mass, grav, drag = params
    axp = -drag * vxp / mass
    ayp = -grav - drag * vyp / mass
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
deltat = .01
m, x0, y0, v0, a0 = 1., 0., .0, 1., 45
sim_params = m, pt.grav_, pt.drag_ # mass, gravity, drag

# create Particle and set force
ball = pt.Particle(x0, y0, v0, a0, m)
ball_force = fr.Forces(resistive_falling, sim_params)
ball.set_force(ball_force)
print("Projectile:", ball)


# ########## do time evolution ##########
# create an Euler solver to evolve Particle
euler = sol.Solver(ball, "Euler", deltat)
print("Euler", euler)
xeul, yeul, teul = integrate(euler)

# create an Euler-Cromer solver to evolve Particle
ball.set_state_angle(x0, y0, v0, a0, 0.)
euler_cromer = sol.Solver(ball, "Euler-Cromer", deltat)
print("Euler-Cromer:", euler_cromer)
xeuc, yeuc, teuc = integrate(euler_cromer)

# create an Midpoint solver to evolve Particle
ball.set_state_angle(x0, y0, v0, a0, 0.)
midpoint = sol.Solver(ball, "Midpoint", deltat)
print("Midpoint", midpoint)
xmid, ymid, tmid = integrate(midpoint)

# generate exact evolution
rana = [analytic_solution(t, v0, a0) for t in teul]
xana, yana = map(list, zip(*rana))
# ########## end time evolution ##########

# generate plots
fig, ax = plt.subplots()
ax.plot(xeul, yeul, '--', label='Euler')
ax.plot(xeuc, yeuc, '--', label='Euler-Cromer')
ax.plot(xmid, ymid, '--', label='Midpoint')
ax.plot(xana, yana, label='Analytic solution')

ax.set(xlabel='x (a.u.)', ylabel='y (a.u.)',
       title='Projectile motion with drag')
ax.grid()

plt.legend()
plt.show()
# END-OF-EXECUTION
