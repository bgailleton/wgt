import numpy as np
import pywtt as pw
import cppwtt as cw
import math
import matplotlib.pyplot as plt





class FastFlood(object):
	"""
		Python/c++ version of Fastflood
		This is the python side, handling the set up, persistence, monitoring and visualisation of the data
		It needs a graph object (at the moment a SMgraph) and a RasterData object holding the dem data and neighbourer.

		B.G. 

	"""
	def __init__(self, dem, graph):
		'''
		Construction the fastflood object and init a bunch of internal variables.
		Arguments:
			- dem: a RasterData object (DEM + neighbourer + other stuff)
			- graph: a SMG graph objectß
		B.G - 2022
		
		'''
		
		# mandatory python stuff
		super(FastFlood, self).__init__()
		
		# dem is the RasterData object holding the base topography
		self.dem = dem
		
		# The SMG graph
		self.graph = graph
		
		# Water height, initialised to barely nothing
		self.hw = np.full_like(dem.data, 1e-6)

		# Precipitation rates in m/s
		self.precipitations = np.full_like(dem.data, 1e-4)

		# Manning coefficient
		self.manning = 0.033

		# Holds the divergence of Q (calculated in the c++ side)
		self.diff = None

		# Variables helping with animated visualisations 
		self.fig = None
		self.ax = None
		self.cbar = None
		self.hwim = None


	def run(self, dt = 1e-3):
		'''
		Run one iteration of the fastflood model with a given time step dt
		B.G.
		'''
		# Calculating water surface (topo + water height)
		topohw = self.dem.data + self.hw
		# The graph computes receivers and fill the topography
		filled = self.graph.compute_graph_multi_filled(topohw,self.dem.neighbourer)
		# A first correction of water height by filling the local minima 
		self.hw += (filled - topohw)
		# Calculating the difference of water height
		self.diff = cw.run_multi_fastflood_static(self.graph, self.dem.neighbourer, self.hw, self.dem.data, self.manning, self.precipitations)
		# Applying the diff to the water height
		self.hw += self.diff * dt
		# Just making sure we do not have negative water height (not possible)
		self.hw[self.hw<0] = 0
		# Done


	def plot_hw(self, figsize = None, vmin = None, vmax = None):
		'''
		Set up interactive plot for water height
		Arguments:
			- figsize: size of the figure [X,Y] in (erk) inches
			- vmin,vmax: the min max of the water height to plot (<vmin will be transparent)

		B.G.
		'''

		# Formatting the vmin and vmax and defaulting them to 60 to 98th percentile if None
		if(vmin is None):
			vmin = np.percentile(self.hw,60)
		if(vmax is None):
			vmax = np.percentile(self.hw,98)

		# Getting the base figure (HS with extent and everything)
		self.fig, self.ax = self.dem.get_basemap(figsize)
		# Getting rid of what's bellow vmin
		toplot = np.copy(self.hw.reshape(self.dem.reshape))
		toplot[toplot<vmin] = np.nan
		# Plotting water height
		self.hwim = self.ax.imshow(toplot, extent = self.dem.extent, vmin = vmin, vmax = vmax, cmap = 'Blues', alpha = 0.8)
		# Plotting the colorbar
		self.cbar = plt.colorbar(self.hwim, label = r'$H_w$ (m)')
		# Updating the display
		self.fig.canvas.draw()

	def update_hw(self, vmin = None, vmax = None):
		'''
		Updating the hw interatctive minute.
		See above for details
		'''
		# Formatting vmin/vmax
		if(vmin is None):
			vmin = np.percentile(self.hw,60)
		if(vmax is None):
			vmax = np.percentile(self.hw,98)
		# Removing < vmin
		toplot = np.copy(self.hw.reshape(self.dem.reshape))
		toplot[toplot<vmin] = np.nan
		# updating the imshow data
		self.hwim.set_data(toplot)
		# updating the range (colorbar gets updated automatically)
		self.hwim.set_clim(vmin = vmin,vmax = vmax)
		# Update the visual
		self.fig.canvas.draw()

	def plot_A(self, figsize = None, vmin = 0, vmax = 1.5):
		filled = self.graph.compute_graph_multi_filled(self.dem.data + self.hw,self.dem.neighbourer)
		A = self.graph.get_DA_proposlope(self.dem.neighbourer,filled)
		# plt.ioff()
		self.fig, self.ax = self.dem.get_basemap(figsize)
		self.hwim = self.ax.imshow(np.log10(A).reshape(self.dem.reshape), extent = self.dem.extent, vmin = vmin, vmax = vmax, cmap = 'Blues', alpha = 0.5)
		self.cbar = plt.colorbar(self.hwim)
		# plt.ion()
		self.fig.canvas.draw()
		# plt.pause(0.1)

	def update_A(self):
		filled = self.graph.compute_graph_multi_filled(self.dem.data + self.hw,self.dem.neighbourer)
		A = self.graph.get_DA_proposlope(self.dem.neighbourer,filled)
		self.hwim.set_data(np.log10(A).reshape(self.dem.reshape))
		self.fig.canvas.draw()








