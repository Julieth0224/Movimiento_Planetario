import sys
sys.path.insert(0, '../')

#from particle import particle as pt
import particle.particle as pt
import matplotlib.pyplot as plt
import numpy as np


def analytic_solution(tf, v0, a0):
    """This expression assumes x0 = 0 and y0 = 0 and m = 1"""
    v0x = v0 * np.cos(np.radians(a0))
    v0y = v0 * np.sin(np.radians(a0))
    xa = v0x * (1. - np.exp(-pt.drag_ * tf)) / pt.drag_
    ya = (v0y + pt.grav_ / pt.drag_) * (1. - np.exp(-pt.drag_ * tf)) \
       - pt.grav_ * tf
    return xa, ya / pt.drag_

# BEGINNING-OF-EXECUTION
deltat = 0.01
x0, y0, v0, a0 = 0., .0, 1., 45

ball = pt.Particle(x0, y0, v0, a0)
print("Earth:", ball)

# do time evolution
xpos, ypos, tpos = [], [], []
while True:
    xc, yc, _, _, tc = ball.get_state()
    if yc < 0: break
    xpos.append(xc)
    ypos.append(yc)
    tpos.append(tc)
    ball.midpoint_step(deltat)

# compare to algebraic solution
rana = [analytic_solution(t, v0, a0) for t in tpos]
xana, yana = map(list, zip(*rana))

xerr = [abs(a - b) for a, b in zip(xana, xpos)]
yerr = [abs(a - b) for a, b in zip(yana, ypos)]

# generate plots
fig, ax = plt.subplots()
ax.plot(tpos, xerr, '--', label='|x_ana - x_pos|')
ax.plot(tpos, yerr, '-.', label='|y_ana - y_pos|')
#ax.plot(xpos, ypos, '--', label='numerical')
#ax.plot(xana, yana, '-.', label='algebraic')

ax.set(xlabel='x (a.u.)', ylabel='y (a.u.)',
       title='Projectile motion. Method: Euler')
ax.grid()

plt.legend()
plt.show()
# END-OF-EXECUTION
