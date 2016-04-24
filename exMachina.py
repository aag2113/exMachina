"""
Recreation of the title screen animation from the DVD menu of ExMachina (http://www.imdb.com/title/tt0470752/).
Works well when paired with music that is unsettlingly serene. 
For the full effect please have a fully sentient robot locked in a nearby room that is primarily interested in your demise.

Author: Adam Gall
email: adam@allgall.net
website: allgall.net
license: BSD

-------
This project was heavily influenced by a tutorial written by Jake Vanderplas.

The original tutorial is here: https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/

-------
Animation of Elastic collisions with Gravity

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""
import numpy as np
from scipy.spatial.distance import pdist, squareform
from matplotlib import collections as mc
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

class ParticleBox:
    def __init__(self,
                 init_state = [[1, 0, 0, -1],
                               [-0.5, 0.5, 0.5, 0.5],
                               [-0.5, -0.5, -0.5, 0.5]],
                 bounds = [-3.2, 3.2, -1.8, 1.8],
                 size = 0.04,
                 M = 0.05,
                 G = 9.8):
        self.init_state = np.asarray(init_state, dtype=float)
        self.M = M * np.ones(self.init_state.shape[0])
        self.size = size
        self.state = self.init_state.copy()
        self.time_elapsed = 0
        self.bounds = bounds
        self.G = G
        self.lines = []


    def step(self, dt):
        self.time_elapsed += dt
        
        self.state[:, :2] += dt * self.state[:, 2:]
        self.find_lines()

        

    def find_lines(self):
        # Clear out current lines
        self.lines = []

        # find pairs of nearby particles
        D = squareform(pdist(self.state[:, :2]))
        ind1, ind2 = np.where(D < .5)
        unique = (ind1 < ind2)
        ind1 = ind1[unique]
        ind2 = ind2[unique]

        # Define a line between "near" pairs
        for i1, i2 in zip(ind1, ind2):
            p1 = (self.state[i1, 0], self.state[i1, 1])
            p2 = (self.state[i2, 0], self.state[i2, 1])
            self.lines.append([p1, p2])

        # Determine if a particle is hitting a boundary
        crossed_x1 = (self.state[:, 0] < self.bounds[0] + self.size)
        crossed_x2 = (self.state[:, 0] > self.bounds[1] - self.size)
        crossed_y1 = (self.state[:, 1] < self.bounds[2] + self.size)
        crossed_y2 = (self.state[:, 1] > self.bounds[3] - self.size)

        self.state[crossed_x1, 0] = self.bounds[1] - self.size
        self.state[crossed_x2, 0] = self.bounds[0] + self.size

        self.state[crossed_y1, 1] = self.bounds[3] - self.size
        self.state[crossed_y2, 1] = self.bounds[2] + self.size


# Initialize animation
def init():
    global box, lns
    particles.set_data([], [])
    lns.set_segments([])
    return particles, lns

# perform animation
def animate(i):
    global box, dt, ax, fig, lns
    box.step(dt)

    ms = int(fig.dpi * 2 * box.size * fig.get_figwidth()
             / np.diff(ax.get_xbound())[0])
    
    # update pieces of the animation
    particles.set_data(box.state[:, 0], box.state[:, 1])
    particles.set_markersize(ms)
    lns.set_segments(box.lines[:])
    return particles, lns

if __name__ == '__main__':

    # set up initial state
    init_state = -0.5 + np.random.random((100, 4))
    init_state[:, 0] *= 6.4
    init_state[:, 1] *= 3.6
    init_state[:, 2:] *= .5

    box = ParticleBox(init_state, size=0.02)
    box.find_lines()
    dt = 1. / 30 # 30fps


    # set up figure and animation
    fig = plt.figure()
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                         xlim=(-3.2, 3.2), ylim=(-1.8, 1.8), axisbg='black')

    # particles holds the locations of the particles
    particles, = ax.plot([], [], 'wo', ms=6)

    # lns is the line collection
    lns = mc.LineCollection([], linewidth=1, color='w')


    ax.add_collection(lns)

    ani = animation.FuncAnimation(fig, animate, frames=600,
                                  interval=10, blit=False, init_func=init)

    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    #ani.save('exMachina.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
    ani.save('exMachina.mp4', fps=30)

    plt.show()