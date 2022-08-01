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
		self.crs = None



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
	rd.crs = out['crs']


	return rd

def CreateRasterData(nx = 200, ny = 200, dx = 30, dy = 30, x_min = 0, y_min = 0 ):
	'''
	Wrapper function creating a RasterData object from scratch, for example to hold other data later or to generate white noise
	'''

	rd = RasterData()
	neighbourer = cw.D8N(nx,ny,dx,dy,x_min,y_min)
	rd.neighbourer = neighbourer
	x_max = (nx + 1) * dx + x_min
	y_max = (ny + 1) * dy + y_min
	rd.extent = [x_min,x_max,y_max, y_min]
	rd.reshape =[ny,nx]
	rd.dx = dx
	rd.dy = dy
	rd.dxy = math.sqrt(dx**2 + dy**2)
	rd.nx = nx
	rd.ny = ny
	rd.nxy = rd.nx * rd.ny
	rd.data = np.zeros(rd.nxy)
	return rd



def fill_raster(rd,method = "priority_flood"):
	'''
	return a filled version of the RasterData, assuming it is topography of course
	'''
	smg = cw.smgraph(rd.nxy,8)
	topo = np.copy(rd.data).ravel()
	ret = None
	if("priority_flood" == method):
		ret = np.array(smg.just_fill_it(topo, rd.neighbourer).reshape(rd.reshape))
	elif("cordonnier" == method):
		smg.init_graph_v6(rd.neighbourer)
		ret = np.array(smg.compute_graph_v6('fill',topo, rd.neighbourer).reshape(rd.reshape))

	return ret





































# end of file