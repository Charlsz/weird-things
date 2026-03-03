# Real-time gravity simulation using Newton's law
import numpy as np
import matplotlib.pyplot as plt
import time

# Gravitational constant (m^3 kg^-1 s^-2)
G = 6.67430e-11

# Body class
class Body:
	def __init__(self, mass, position, velocity):
		self.mass = mass
		self.position = np.array(position, dtype=float)
		self.velocity = np.array(velocity, dtype=float)

def compute_gravitational_force(body1, body2):
	r_vec = body2.position - body1.position
	r_mag = np.linalg.norm(r_vec)
	if r_mag == 0:
		return np.zeros(2)
	force_mag = G * body1.mass * body2.mass / r_mag**2
	force_vec = force_mag * r_vec / r_mag
	return force_vec

# Simulation parameters
dt = 0.01  # time step (seconds)
total_time = 10  # seconds

# Example: Earth and Moon (scaled for visualization)
earth = Body(mass=5.972e24, position=[0, 0], velocity=[0, 0])
moon = Body(mass=7.348e22, position=[384400000, 0], velocity=[0, 1022])
bodies = [earth, moon]

# Set up plot
fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_xlim(-5e8, 5e8)
ax.set_ylim(-5e8, 5e8)
earth_plot, = ax.plot([], [], 'bo', markersize=10, label='Earth')
moon_plot, = ax.plot([], [], 'go', markersize=5, label='Moon')
ax.legend()

def update_plot():
	earth_plot.set_data(earth.position[0], earth.position[1])
	moon_plot.set_data(moon.position[0], moon.position[1])
	plt.draw()
	plt.pause(0.001)

# Simulation loop
steps = int(total_time / dt)
for _ in range(steps):
	# Compute forces
	force_on_earth = compute_gravitational_force(earth, moon)
	force_on_moon = compute_gravitational_force(moon, earth)

	# Update velocities
	earth.velocity += force_on_earth / earth.mass * dt
	moon.velocity += force_on_moon / moon.mass * dt

	# Update positions
	earth.position += earth.velocity * dt
	moon.position += moon.velocity * dt

	update_plot()
	time.sleep(dt)

plt.show()
