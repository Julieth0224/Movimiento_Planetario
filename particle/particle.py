# particle.py

import numpy as np


grav_ = 1.
drag_ = 4.


class Particle:

    def __init__(self, x0=0., y0=0., v0=0., alpha0=0., m0=1., t0=0.):
        self.m, self.t = m0, t0
        self.x, self.y = x0, y0
        self.vx = v0 * np.cos(np.radians(alpha0))
        self.vy = v0 * np.sin(np.radians(alpha0))
        self.force = None
        
    def __str__(self):
        strng = "particle state\n"
        strng += "m = {}, t = {}\n".format(self.m, self.t)
        strng += "r = ({:.4f}, {:.4f})\n".format(self.x, self.y)
        strng += "v = ({:.4f}, {:.4f})\n".format(self.vx, self.vy)
#        strng += f"v = ({self.vx:.4f}, {self.vy:.4f})\n"
        return strng

    def get_state(self):
        return self.x, self.y, self.vx, self.vy, self.t
    
    def set_state(self, *state):
        self.x, self.y, self.vx, self.vy, self.t = state

    def set_state_angle(self, *state):
        self.x, self.y, v0, a0, self.t = state
        self.vx = v0 * np.cos(np.radians(a0))
        self.vy = v0 * np.sin(np.radians(a0))
    
    def set_force(self, netforce):
        self.force = netforce


if __name__ == "__main__":
    print(__name__)
    earth = Particle()
    print("Earth:", earth)
    mars = Particle(1., 2., 3., 4.)
    print("Mars:", mars)
