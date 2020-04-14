
import sys
sys.path.insert(0, '../../')

import numpy as np
import matplotlib.pyplot as plt

import particle.particle as pt
import solver.solver as sl
import forces.forces as fr

GMe = 20

A = 0.
A2 = 0.3538 # 0.16 # 0.18
print("Con A =", A)
print("Luego de una orbita, A =", A2)


def init_c(x, y, vx):
    v = np.sqrt(GMe / np.sqrt(x**2 + y**2))
    vy = np.sqrt(v**2 - vx**2)
    return float(v), vy

def Newton_again(state, params):
    x, y, vx, vy, t = state
    mass, gu, A = params
    r = (x**2 + y**2)**0.5
    ax = (-gu*x)/(r**3) + A * np.sqrt(vx**2 + vy**2) * x
    ay = (-gu*y)/(r**3) - A * np.sqrt(vx**2 + vy**2) * y
    return vx, vy, ax, ay, 1

def mechanic_energy(v, m, x, y):
    emt = (1/2) * m * (v**2) - GMe * m / np.sqrt(x**2 + y**2)
    return emt


v0, vy = init_c(1, 0, 0)
m, x0, y0, v0, a0 = 1., 1., .0, v0, 90.
sim_params = m, GMe, A
deltat = 0.01/128

m1 = "Euler"
m2 = "Euler-Cromer"
m3 = "Midpoint"

num_method = m2 # m1 # m3

satellite = pt.Particle(x0, y0, v0, a0)
satellite_force = fr.Forces(Newton_again, sim_params)
satellite.set_force(satellite_force)
euler = sl.Solver(satellite, num_method, deltat)
xpos, ypos, tpos, em, vel = [], [], [], [], []
cont = 0

for i in range(200000):
    xc, yc, vxc, vyc, tc = satellite.get_state()
    v = np.sqrt(vxc**2 + vyc**2)
    vel.append(v)
    emt = mechanic_energy(v, m, xc, yc)
    em.append(round(emt, 2))
    xpos.append(xc)
    ypos.append(yc)
    tpos.append(tc)
    euler.do_step()
    xc1, yc1, vxc1, vyc1, tc1 = satellite.get_state()
    if yc < 0 and yc1 >= 0:
        cont += 1
        if cont == 1:
            A = A2
            sim_params = m, GMe, A
            satellite_force = fr.Forces(Newton_again, sim_params)
            satellite.set_force(satellite_force)

print("Method:", num_method, ", dt:", deltat)
print("Energia Mec Tot In:", em[0], ", Energia Mec Tot Fin:", em[-1])
print("Velocidad Tot In:", vel[0], ", Velocidad Tot Fin:", vel[-1])
print("Orbitas:", cont)

print(A)

# Grafica Orbitas
fig, ax = plt.subplots()
ax.plot(xpos, ypos, '--', c = 'purple', label=num_method)

ax.set(xlabel='x (a.u.)', ylabel='y (a.u.)')
ax.grid()

# Grafica energia vs tiempo
# fig, ax = plt.subplots()
# ax.plot(tpos, em, '--', c = 'purple',label=num_method)
#
# ax.set(xlabel='t (h)', ylabel='emt')
# ax.grid()

# Grafica velocidad vs tiempo
# fig, ax = plt.subplots()
# ax.plot(tpos, vel, '--', c = 'purple',label=num_method)
#
# ax.set(xlabel='t (h)', ylabel='v(E.U / h)')
# ax.grid()

plt.legend()
plt.show()
