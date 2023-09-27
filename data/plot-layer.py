import matplotlib.pyplot as plt
import numpy as np
import sys

plt.axis('equal')
plt.axis('off')

if len(sys.argv) > 1:
	n = int(sys.argv[1])
else:
	n = 1

# Open the file for reading
with open('hat.path', 'r') as file:
	# Read the first line
	line = file.readline()
	
	# Continue reading lines until the end of the file is reached
	while line:
		x = []
		y = []
		num_vertices = int(line)
		for i in range(num_vertices):
			line = file.readline()
			line_data = line.strip().split()
			x.append(float(line_data[0]))
			y.append(float(line_data[1]))
			z = float(line_data[2])
		if z == n * 0.08:
			plt.plot(x, y, color='#6699cc', linewidth=0.5)
		line = file.readline()
plt.savefig('Layer ' + str(n) + '.svg')

