import sys
sys.path.insert(0, '../../')

import numpy as np
import matplotlib.pyplot as plt

import particle.particle as pt
import solver.solver as sl
import forces.forces as fr

GM = 4 * (np.pi)**2

def Newton_again(state, params):
    x, y, vx, vy, t = state
    r = (x**2 + y**2)**0.5
    ax = (-GM*x)/(r**3)
    ay = (-GM*y)/(r**3)
    return vx, vy, ax, ay, 1

def integrate(obj):
    cont = 0
    xpos, ypos, tpos, vpos = [], [], [], []
    tc = 0
    for i in range(200000):
        xc, yc, vxc, vyc, tc = obj.objs.get_state()
        v = np.sqrt(vxc**2 + vyc**2)
        vpos.append(v)
        xpos.append(xc)
        ypos.append(yc)
        tpos.append(tc)
        obj.do_step()
        xc1, yc1, vxc1, vyc1, tc1 = obj.objs.get_state()
        if yc < 0 and yc1 >= 0:
            cont += 1
    print("Orbitas", cont)
    return xpos, ypos, tpos, vpos


sim_params = pt.grav_
deltat = 0.01/1024

m1 = "Euler"
m2 = "Euler-Cromer"
m3 = "Midpoint"

num_method = m2 # m1 # m3

m, x0, y0, v0, a0 = 1., 1.5, 0., 1.5, 90.

earth = pt.Particle(x0, y0, v0, a0)
earth_force = fr.Forces(Newton_again, sim_params)
earth.set_force(earth_force)
euler = sl.Solver(earth, num_method, deltat)
xvac, yvac, tvac, v = integrate(euler)


#Calculamos las velocidades minimas y maximas con x0 = 1.5, v0 = 1.5

vmin = min(v)
tmin = tvac[v.index(vmin)]
xmin = xvac[v.index(vmin)]
ymin = yvac[v.index(vmin)]

vmax = max(v)
tmax = tvac[v.index(vmax)]
xmax = xvac[v.index(vmax)]
ymax = yvac[v.index(vmax)]


print("Vmin:", round(vmin, 4), ", with t = ", round(tmin, 4), ", in the point x =", round(xmin, 4), "and y =", round(ymin, 4))
print("Vmax:", round(vmax, 4), ", with t = ", round(tmax, 4), ", in the point x =", round(xmax, 4), "and y =", round(ymax, 4))



# Muestra grafica de velocidad con tiempo

fig, ax = plt.subplots()
ax.plot(tvac, v, '--', label=num_method)

ax.set(xlabel='t (yr)', ylabel="v (a.u./yr)")

ax.grid()


plt.legend()
plt.show()
