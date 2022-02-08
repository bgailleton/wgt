// First function registering basic topoplot
const registerTopoPlot = function(pyPlotter){
	pyPlotter.topoPlot = pyodide.runPython(`
def topoPlot(dataCauldron):
	dataCauldron = dataCauldron.to_py()
	import matplotlib.pyplot as plt
	import numpy as np
	from js import document

	# First getting all the data from the Cauldron
	topo = np.asarray(dataCauldron["topo"]).reshape(dataCauldron["ny"],dataCauldron["nx"])
	HS = np.asarray(dataCauldron["HS"]).reshape(dataCauldron["ny"],dataCauldron["nx"])
	#pXriv = np.asarray(dataCauldron["Xriv"]);
	#pYriv = np.asarray(dataCauldron["Yriv"]);
	#size = np.asarray(dataCauldron["Ariv"])

	fig, ax = plt.subplots()

	cb = ax.imshow(topo, cmap = 'gist_earth', extent = dataCauldron["extents"], aspect = "auto")
	ax.imshow(HS, cmap = 'gray', extent = dataCauldron["extents"], aspect = "auto", alpha = 0.65)

	#size = size/size.max()
	#size = (size * (5 - 2)) + 2

	#ax.scatter(pXriv,pYriv, c = 'b', s = size, lw = 0, alpha = 0.1)

	plt.colorbar(cb, label="Elevation")

	ax.set_xlabel("Easting")
	ax.set_ylabel("Northing")


	def create_root_element(self):
		return document.getElementById('pyplotdiv')
	# proxy = create_proxy(create_root_element2)

	##override create_root_element method of canvas by one of the functions above
	fig.canvas.create_root_element = create_root_element.__get__(create_root_element, fig.canvas.__class__)
	#fig.canvas.create_root_element = create_root_element

	fig.canvas.show()



topoPlot
`)

}


const registerPlotters = function(pyplotter){
	registerTopoPlot(pyPlotter);
}











// old function, keepign for legacy
function plot_topo()
{	

	const tttopo = mg.get_topo_array()
	const tttHS = mg.get_HS()
	window.ttopo = new Float32Array(mg.nx * mg.ny)
	window.ttHS = new Float32Array(mg.nx * mg.ny)
	for (var i = 0; i < ttopo.length; i++) {
      ttopo[i] = tttopo.get(i);
      ttHS[i] = tttHS.get(i);
  }

  window.nriv = mg.getNriv()
  const txriv = mg.getXriv()
  const tyriv = mg.getYriv()
  const tAriv = mg.getAriv()
	window.ttxriv = new Float32Array(window.nriv)
	window.ttyriv = new Float32Array(window.nriv)
	window.textents = new Float32Array(4)
	window.trA = new Float32Array(window.nriv)

	for (var i = 0; i < window.nriv; i++) {
		window.ttxriv[i] = txriv.get(i);
		window.ttyriv[i] = tyriv.get(i);
		window.trA[i] = tAriv.get(i);
	}

	for (var i = 0; i < 4; i++) {
		window.textents[i] = mg.extents.get(i)	
	}



  window.Xriv = window.ttxriv;
  window.Yriv = window.ttyriv;

	window.nx = mg.nx
	window.ny = mg.ny
	window.dx = mg.dx
	window.dy = mg.dy

	console.log("n rivers::" +window.nriv )
	console.log("extents::" + window.textents )

	pyodide.runPython(`
import matplotlib.pyplot as plt
import numpy as np
from pyodide import create_proxy
from js import document
from js import ttopo, ttHS, nx, ny, dx, dy, Xriv, Yriv, textents, trA
#import io, base64

#Loading the topography

topo = np.asarray(ttopo.to_py()).reshape(ny,nx)
HS = np.asarray(ttHS.to_py()).reshape(ny,nx)

pXriv = np.asarray(Xriv.to_py());
pYriv = np.asarray(Yriv.to_py());

fig, ax = plt.subplots()

#ax.imshow([[0,1],[2,6]], cmap = 'gist_earth')
cb = ax.imshow(topo, cmap = 'gist_earth', extent = textents.to_py(), aspect = "auto")
ax.imshow(HS, cmap = 'gray', extent = textents.to_py(), aspect = "auto", alpha = 0.65)

size = np.asarray(trA.to_py())
size = size/size.max()
size = (size * (5 - 2)) + 2

ax.scatter(pXriv,pYriv, c = 'b', s = size, lw = 0, alpha = 0.1)

plt.colorbar(cb, label="Elevation")

ax.set_xlabel("Easting")
ax.set_ylabel("Northing")


def create_root_element2(self):
	return document.getElementById('pyplotdiv')
proxy = create_proxy(create_root_element2)

##override create_root_element method of canvas by one of the functions above
fig.canvas.create_root_element = create_root_element2.__get__(
create_root_element2, fig.canvas.__class__)

fig.canvas.show()


	`);

	// document.getElementById("pyplotfigure").src=pyodide.globals.img_str
	delete tttopo
	delete tttHS

	window.ttopo = null
	window.ttHS = null

};








































// End of file