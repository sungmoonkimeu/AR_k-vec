import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Create a figure and a 3D axis
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Define the parameters for the sphere
u = np.linspace(0, np.pi/2, 100)
v = np.linspace(0, np.pi/2, 100)

# Create a meshgrid
u, v = np.meshgrid(u, v)

# Define the radius of the sphere
r = 1.0

# Calculate the coordinates of the points on the sphere
x = r * np.sin(u) * np.cos(v)
y = r * np.sin(u) * np.sin(v)
z = r * np.cos(u)

# Plot the 1/8 sphere
ax.plot_surface(x, y, z, color='b', alpha=0.6)

# Define the vertices of the square
square_vertices = np.array([
    [0.5, 0.5, 0.5],  # Vertex 1
    [0.5, -0.5, 0.5], # Vertex 2
    [-0.5, -0.5, 0.5],# Vertex 3
    [-0.5, 0.5, 0.5]  # Vertex 4
])

# Define the edges of the square
square_edges = [
    [square_vertices[0], square_vertices[1]],
    [square_vertices[1], square_vertices[2]],
    [square_vertices[2], square_vertices[3]],
    [square_vertices[3], square_vertices[0]]
]

# Plot the square on the sphere
for edge in square_edges:
    xs, ys, zs = zip(*edge)
    ax.plot(xs, ys, zs, color='r')

# Set axis labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Remove the grid
ax.grid(False)

# Show the plot
plt.show()
