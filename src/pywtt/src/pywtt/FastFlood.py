import numpy as np
import numexpr as ne
import pywtt as pw
import cppwtt as cw
import math
import matplotlib.pyplot as plt
# import cv2

# def fill_python(marker: np.ndarray, mask: np.ndarray, radius: int = 1):
#     """Iteratively expand the markers white keeping them limited by the mask during each iteration.
#     :param marker: Grayscale image where initial seed is white on black background.
#     :param mask: Grayscale mask where the valid area is white on black background.
#     :param radius Can be increased to improve expansion speed while causing decreased isolation from nearby areas.
#     :returns A copy of the last expansion.
#     Written By Semnodime.
#     """
    
#     kernel = np.ones(shape=(radius * 2 + 1,) * 2, dtype=np.uint8)

#     while True:
#       expanded = cv2.dilate(src=marker, kernel=kernel)
#       cv2.bitwise_and(src1=expanded, src2=mask, dst=expanded)

#       # Termination criterion: Expansion didn't change the image at all
#       if (marker == expanded).all():
#           return expanded
#       marker = expanded





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

		# Holds the divergence of Q (calculated inplace in the c++ side, but init on the python side)
		self.Qwin = np.zeros_like(self.dem.data)
		self.Qwout = np.zeros_like(self.dem.data)
		self._Sw = None
		self._susmt = np.zeros_like(self.dem.data)

		# Variables helping with animated visualisations 
		self.fig = None
		self.ax = None
		self.cbar = None
		self.hwim = None

		self.mask = np.asarray(self.dem.neighbourer.get_mask_array(), dtype = bool)


	def run(self, dt = 1e-3, graph_mode = 'multi_opti', force_flood = False, depressions = "cordonnier_fill", slope_mode = "manning"):
		'''
		Run one iteration of the fastflood model with a given time step dt
		B.G.
		'''
		# Calculating water surface (topo + water height)
		topohw = self.dem.data + self.hw
		if(depressions == "cordonnier_fill"):
			filled = self.graph.compute_graph_v6("fill",topohw,self.dem.neighbourer)
		elif(depressions == "cordonnier_carve"):
			filled = self.graph.compute_graph_v6("carve",topohw,self.dem.neighbourer)
		elif(depressions == "priority_flood"):
			filled = self.graph.compute_graph_v6_PQ(topohw,self.dem.neighbourer)
		elif(depressions == "none"):
			filled = self.graph.compute_graph_v6_nodep(topohw,self.dem.neighbourer)


		# initilising internal array ifnot done yet (it needs the graph to have been calculated at least once so it is unsafe to do it at build time)
		if(self._Sw is None):
			self._Sw = np.zeros(self.graph.get_rec_array_size())

		# A first correction of water height by filling the local minima 


		self.hw += (filled - topohw)
		self.hw[self.hw<0] = 1e-6

		self._Sw.fill(1e-6)
		self._susmt.fill(0)

		if(force_flood):
			pQwin = self.Qwin.copy()

		cw.compute_Sw_sumslopes(self.graph, self.dem.neighbourer, self.hw, self.dem.data, self._Sw, self._susmt)
		
		if(slope_mode == "manning"):
			# Calculating the difference of water height
			cw.run_multi_fastflood_static(self.graph, self.dem.neighbourer, self.hw, self.dem.data, self.manning, self.precipitations,self.Qwin,self.Qwout, self._Sw, self._susmt)
		elif(slope_mode == "prop"):
			self.Qwin = self.graph.get_DA_proposlope(self.dem.neighbourer,filled)
			cw.run_multi_fastflood_static_ext_Qwin(self.graph, self.dem.neighbourer, self.hw, self.dem.data, self.manning,self.Qwin,self.Qwout, self._Sw)
		elif(slope_mode == "precipiton_like"):
			cw.run_multi_fastflood_static_precipitonlike(self.graph, self.dem.neighbourer, self.hw, self.dem.data, self.manning, self.precipitations,self.Qwin,self.Qwout, self._Sw, self._susmt)


		if(force_flood):
			np.maximum(self.Qwin,pQwin,out=self.Qwin)

		# Applying the diff to the water height
		self.hw += (self.Qwin-self.Qwout)/(self.dem.dx * self.dem.dy) * dt
		# print(np.min(self.Qwin), np.min(self.Qwout))


		# Just making sure we do not have negative water height (not possible)
		self.hw[self.hw<0] = 1e-6
		self.hw[~self.mask] = 0
		# Done

	def run_single_flow(self, dt = 1e-3, no_numexpr = False, depressions = "cordonnier_fill"):
		'''
		'''

		topohw = self.dem.data + self.hw
		if(depressions == "cordonnier_fill"):
			filled = self.graph.compute_graph_v6_SS("fill",topohw,self.dem.neighbourer)
		elif(depressions == "cordonnier_carve"):
			filled = self.graph.compute_graph_v6_SS("carve",topohw,self.dem.neighbourer)
		elif(depressions == "priority_flood"):
			filled = self.graph.compute_graph_v6_PQ(topohw,self.dem.neighbourer)
		elif(depressions == "none"):
			filled = self.graph.compute_graph_v6_nodep(topohw,self.dem.neighbourer)

		self.hw += (filled - topohw)
		self.hw[self.hw<0] = 1e-6
		Sw = np.full_like(self.hw,1e-6)
		self.Qwin.fill(0)
		self.Qwout.fill(0)
		self.Qwin = self.graph.get_DA_SS(self.dem.neighbourer, filled) * self.precipitations
		Sr = self.graph.get_Sreceivers()
		dx = self.graph.get_dx_array()
		if(no_numexpr):
			Sw = (filled - filled[Sr])/dx
			tmask = (Sw>0)
			# Sw[Sw <= 0] = 1e-6
			# self.Qwout[tmask] = dx[tmask] * 1/self.manning * np.power(self.hw[tmask],5/3) * np.sqrt(Sw[tmask])
			self.Qwout[tmask] = rd.dx * 1/self.manning * np.power(self.hw[tmask],5/3) * np.sqrt(Sw[tmask])
			# self.Qwout[tmask] = 1 * 1/self.manning * np.power(self.hw[tmask],5/3) * np.sqrt(Sw[tmask])


			self.hw += (self.Qwin-self.Qwout)/(self.dem.dx * self.dem.dy) * dt
			self.hw[self.hw<=0 | ~self.mask] = 1e-6
		else:
			recZ = filled[Sr]
			Sw = ne.evaluate("(filled - recZ)/dx")
			Sw = ne.evaluate("where(Sw<=0,1e-6,Sw)")
			self.Qwout = ne.evaluate("dx * 1/manning * hw**(5/3) * sqrt(Sw)", local_dict = {'manning': self.manning,'hw':self.hw, 'Sw':Sw, 'dx':dx})
			self.hw = ne.evaluate("hw + ((Qwin-Qwout)/(dx * dy) * dt)" , local_dict = {'Qwin': self.Qwin,'Qwout': self.Qwout,'hw':self.hw,'dx': self.dem.dx,'dy': self.dem.dy, 'dt':dt})
			self.hw = ne.evaluate("where(((hw<0) | (mask == False)),1e-6, hw)", local_dict = {'hw':self.hw, 'mask': self.mask})

		





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

	def plot_Qw(self, figsize = None, vmin = 0, vmax = 1.5, out = True):
		# filled = self.graph.compute_graph_multi_filled(self.dem.data + self.hw,self.dem.neighbourer)
		# A = self.graph.get_DA_proposlope(self.dem.neighbourer,filled)
		# plt.ioff()
		self.fig, self.ax = self.dem.get_basemap(figsize)
		toplot = self.Qwout if (out) else self.Qwin
		self.hwim = self.ax.imshow(np.log10(toplot).reshape(self.dem.reshape), extent = self.dem.extent, vmin = vmin, vmax = vmax, cmap = 'Blues', alpha = 0.5)
		self.cbar = plt.colorbar(self.hwim)
		# plt.ion()
		self.fig.canvas.draw()
		# plt.pause(0.1)

	def update_Qw(self, out = True):
		# filled = self.graph.compute_graph_multi_filled(self.dem.data + self.hw,self.dem.neighbourer)
		# A = self.graph.get_DA_proposlope(self.dem.neighbourer,filled)
		toplot = self.Qwout if (out) else self.Qwin
		self.hwim.set_data(np.log10(toplot).reshape(self.dem.reshape))
		self.fig.canvas.draw()










