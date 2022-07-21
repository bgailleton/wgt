import numpy as np
import pywtt as pw
import cppwtt as cw
import matplotlib.pyplot as plt
import pywtt.plot_helper as ph
import math


class RasterData(object):
	"""
Raster data is a wrapper class containing spatial, raster-like data and its geographical context.
B.G.
	"""

	def __init__(self):
		
		super(RasterData, self).__init__()

		self.data = None
		self.cpp = None
		self.neighbourer = None
		self.extent = None
		self.reshape = None
		self.dx = None
		self.dy = None
		self.dxy = None
		self.nx = None
		self.ny = None
		self.nxy = None



	def D2(self):
		return self.data.reshape(self.reshape)


	def get_basemap(self, figsize = None):

		HS = self.neighbourer.get_HS(self.data).reshape(self.reshape)
		fig,ax = plt.subplots(figsize = figsize)
		ax.imshow(HS, extent = self.extent, vmin =0, vmax = 1, cmap = 'gray')
		ph.fix_map_labels(fig, ax, self)
		ax.set_xlabel('X (km)')
		ax.set_ylabel('Y (km)')
		return fig,ax







def RasterData_from_file(fname):
	"""
	Load a raster file like tif, asc, or anything supported by rasterio into a RasterData file
	"""
	out = pw.raster_loader.load_raster(fname)
	topo = out['array'].ravel().astype(np.float64)
	neighbourer = cw.D8N(out["nx"],out["ny"],out["dx"],out["dy"],out["x_min"],out["y_min"])
	
	rd = RasterData()
	print(topo.dtype)
	rd.data = topo
	rd.neighbourer = neighbourer
	rd.extent = [out['x_min'],out['x_max'],out['y_max'], out['y_min']]
	rd.reshape =[ out['ny'], out['nx']]
	rd.dx = out['dx']
	rd.dy = out['dy']
	rd.dxy = math.sqrt(out['dx']**2 + out['dy']**2)
	rd.nx = out['nx']
	rd.ny = out['ny']
	rd.nxy = rd.nx * rd.ny
	rd.cpp = cw.numvecf64(rd.data)

	return rd






































# end of file