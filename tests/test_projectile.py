import sys
sys.path.insert(0, '../')

#from particle import particle as pt
import particle.particle as pt
import matplotlib.pyplot as plt
import numpy as np


def algebraic_solution(tf, x0, y0, v0, a0):
    v0x = v0 * np.cos(np.radians(a0))
    v0y = v0 * np.sin(np.radians(a0))
    xa = x0 + v0x * tf
    ya = y0 + v0y * tf - 0.5 * pt.grav_ * tf**2
    return xa, ya

# BEGINNING-OF-EXECUTION
deltat = 0.01
x0, y0, v0, a0 = 0., .5, 1., 45

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
    ball.step(deltat)

# compare to algebraic solution
rana = [algebraic_solution(t, x0, y0, v0, a0) for t in tpos]
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
