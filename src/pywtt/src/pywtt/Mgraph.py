import numpy as np
import pywtt as pw
import cppwtt as cw
import math






class Graph(object):
	"""docstring for Mgraph"""
	def __init__(self, dem, method = 'carve', topology = 'multi', OMP = False, n_proc = 4):
		super(Graph, self).__init__()
		if(topology == 'multi'):
			self.graph = cw.mgraph(dem.nxy)
			self.multi = True
		
		if(OMP):
			self.graph.compute_graph_OMP(method,dem.neighbourer,dem.data.ravel(), n_proc)
		else:
			self.graph.compute_graph(method,dem.neighbourer,dem.data.ravel())



	def update(self,dem,method = 'carve', replace_topo = None):
		if(replace_topo is None):
			return self.graph.update_graph(method,dem.neighbourer,dem.data.ravel())
		else:
			return self.graph.update_graph(method,dem.neighbourer,replace_topo)


	def drainage_area(self, dem):
		if(self.multi):
			return np.array(self.graph.get_DA_proposlope(dem.neighbourer,dem.data))













































#End of File
