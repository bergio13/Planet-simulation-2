from matplotlib import style
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from astropy.time import Time
from astroquery.jplhorizons import Horizons

sim_start_date = "2022-01-01"     # simulating a solar system starting from this date
sim_duration = 10 * 365                # (int) simulation duration in days
m_earth = 5.9722e24 / 1.98847e30  # Mass of Earth relative to mass of the sun
m_moon = 7.3477e22 / 1.98847e30

# Define the objects: the planets and the sun
class Object:              
    def __init__(self, name, rad, color, distance, velocity):
        self.name = name
        self.distance = np.array(distance, dtype=float)
        self.velocity = np.array(velocity, dtype=float)
        self.xs = []
        self.ys = []
        self.plot = ax.scatter(distance[0], distance[1], color=color, s=rad**2, edgecolors=None, zorder=10)
        self.line, = ax.plot([], [], color=color, linewidth=0.5, linestyle='--')

class SolarSystem:
    def __init__(self, sun):
        self.sun = sun
        self.planets = []
        self.time = None
        self.timestamp = ax.text(.03, .94, 'Date: ', color='darkviolet', transform=ax.transAxes, family='fantasy')
        
    def add_planet(self, planet):
        self.planets.append(planet)
     
    # Evolve the trajectories   
    def evolve(self):           
        dt = 1.0
        self.time += dt
        plots = []
        lines = []
        for planet in self.planets:
            planet.distance += planet.velocity * dt
            acc = -2.959e-4 * planet.distance / np.sum(planet.distance**2)**(3./2)  # in units of AU/day^2
            planet.velocity += acc * dt
            planet.xs.append(planet.distance[0])
            planet.ys.append(planet.distance[1])
            planet.plot.set_offsets(planet.distance[:2])
            planet.line.set_xdata(planet.xs)
            planet.line.set_ydata(planet.ys)
            plots.append(planet.plot)
            lines.append(planet.line)
        self.timestamp.set_text('Date: ' + self.time.iso[:10])
        return plots + lines + [self.timestamp]

plt.style.use('dark_background')
fig = plt.figure(figsize=[7, 7])
ax = plt.axes([0., 0., 1., 1.], xlim=(-1.8, 1.8), ylim=(-1.8, 1.8))
ax.set_aspect('equal')
ax.axis('off')
ss = SolarSystem(Object("Sun", 28, 'gold', [0, 0, 0], [0, 0, 0]))
ss.time = Time(sim_start_date)
ss.time.format = 'jd'
colors = ['dimgray', 'lightgray', 'blue', 'orangered']
sizes = [0.38, 0.95, 1., 0.53]
names = ['Mercury', 'Venus', 'Earth', 'Mars']
texty = [.47, .73, 1, 1.5]
for i, nasaid in enumerate([1, 2, 3, 4]):  # The 1st, 2nd, 3rd, 4th planet in solar system
    obj = Horizons(id=nasaid, location="@sun", epochs=ss.time, id_type='id').vectors()
    ss.add_planet(Object(nasaid, 20 * sizes[i], colors[i], 
                         [np.double(obj[xi]) for xi in ['x', 'y', 'z']], 
                         [np.double(obj[vxi]) for vxi in ['vx', 'vy', 'vz']]))
    ax.text(0, - (texty[i] + 0.1), names[i], color=colors[i], zorder=1000, ha='center', fontsize='small', family='monospace')
    
def animate(i):
    return ss.evolve()

anim = animation.FuncAnimation(fig, animate, repeat=False, frames=sim_duration, blit=True, interval=10,)
plt.show()