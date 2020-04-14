import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, '../')

import particle.particle as pt
import forces.forces as fr


def resistive_falling(state, params):
    xp, yp, vxp, vyp, _ = state
    mass, grav, drag = params
    axp = -drag * vxp / mass
    ayp = -grav - drag * vyp / mass
    return vxp, vyp, axp, ayp, 1.


# BEGINNING-OF-EXECUTION
deltat = .01
m, x0, y0, v0, a0 = 1., 0., .0, 1., 45
sim_params = m, pt.grav_, pt.drag_ # mass, gravity, drag

ball = pt.Particle(x0, y0, v0, a0, m)
ball_force = fr.Forces(resistive_falling, sim_params)
ball.set_force(ball_force)
print("Projectile:", ball)

# do time evolution
xpos, ypos, tpos = [], [], []
while True:
    xc, yc, _, _, tc = ball.get_state()
    if yc < 0: break
    xpos.append(xc)
    ypos.append(yc)
    tpos.append(tc)
    ball.euler_step(deltat)

# generate plots
fig, ax = plt.subplots()
ax.plot(xpos, ypos, '--', label='numerical')

ax.set(xlabel='x (a.u.)', ylabel='y (a.u.)',
       title='Projectile motion with drag')
ax.grid()

plt.legend()
plt.show()
# END-OF-EXECUTION
