// This file host the functions managing the python plotters.
// apart from the initialiser, these are not called directly, but from other analyser
// this really only manage the writing/initialisation/registering of the python function inside JS

// Initialise the figure and create a div holding it
const initpyPlot = function(){
	addToLog("Creating the figure holder")
	pyodide.runPython(`
import matplotlib.pyplot as plt
from js import document
import numpy as np

global_holder = {}

fig, ax = plt.subplots()
ax.set_xlabel("Easting")
ax.set_ylabel("Northing")
def create_root_element(self):
	return document.getElementById('pyplotdiv')
# proxy = create_proxy(create_root_element2)

##override create_root_element method of canvas by one of the functions above
fig.canvas.create_root_element = create_root_element.__get__(create_root_element, fig.canvas.__class__)
#fig.canvas.create_root_element = create_root_element
fig.canvas.show()

		`)
	addToLog("Figure ready")
}

// First function registering basic topoplot
const registerTopoPlot = function(pyPlotter){
	pyPlotter.topoPlot = pyodide.runPython(`
def topoPlot(dataCauldron):
	dataCauldron = dataCauldron.to_py()

	# First getting all the data from the Cauldron
	topo = np.asarray(dataCauldron["topo"]).reshape(dataCauldron["ny"],dataCauldron["nx"])
	HS = np.asarray(dataCauldron["HS"]).reshape(dataCauldron["ny"],dataCauldron["nx"])
	cb = ax.imshow(topo, cmap = 'gist_earth', extent = dataCauldron["extents"], aspect = "auto", vmin = dataCauldron["seaLvl"])
	ax.imshow(HS, cmap = 'gray', extent = dataCauldron["extents"], aspect = "auto", alpha = dataCauldron["alphaHS"])

	#size = size/size.max()
	#size = (size * (5 - 2)) + 2

	#ax.scatter(pXriv,pYriv, c = 'b', s = size, lw = 0, alpha = 0.1)

	plt.colorbar(cb, label="Elevation")

	ax.set_xlabel("Easting")
	ax.set_ylabel("Northing")

	fig.canvas.show()


topoPlot
`)

}

// Basic figure replotter
const registerToporePlot = function(pyPlotter){
	pyPlotter.toporePlot = pyodide.runPython(`
def toporePlot(dataCauldron):
	dataCauldron = dataCauldron.to_py()
	fig.clear()
	ax = fig.add_subplot()

	# First getting all the data from the Cauldron
	topo = np.asarray(dataCauldron["topo"]).reshape(dataCauldron["ny"],dataCauldron["nx"])
	HS = np.asarray(dataCauldron["HS"]).reshape(dataCauldron["ny"],dataCauldron["nx"])
	cb = ax.imshow(topo, cmap = 'gist_earth', extent = dataCauldron["extents"], aspect = "auto", vmin = dataCauldron["seaLvl"])
	ax.imshow(HS, cmap = 'gray', extent = dataCauldron["extents"], aspect = "auto", alpha = dataCauldron["alphaHS"])

	#size = size/size.max()
	#size = (size * (5 - 2)) + 2

	#ax.scatter(pXriv,pYriv, c = 'b', s = size, lw = 0, alpha = 0.1)

	plt.colorbar(cb, label="Elevation")

	ax.set_xlabel("Easting")
	ax.set_ylabel("Northing")

	fig.canvas.show()


toporePlot
`)

}

// FUnction managing the plotting of rivers
const registerExtractRivers = function(pyPlotter){
	pyPlotter.riverPlot = pyodide.runPython(`
def riverPlot(X,Y,A, minsize, maxsize):

	if("scatter_riv" in global_holder.keys()):
		global_holder["scatter_riv"].remove()
		fig.canvas.show()

	
	size = np.asarray(A.to_py())
	#print("4")
	#size = np.log(size)
	size = size/size.max()
	size = (size) * (maxsize - minsize) + minsize
	global_holder["scatter_riv"] = ax.scatter(np.asarray(X.to_py()), np.asarray(Y.to_py()), s = size, lw = 0, c = "blue")
	
	fig.canvas.show()

riverPlot`)
}


// Load and initialise the plotting engine
const initPlotters = function(pyplotter){
	// Iniplot register the python plotting function
	initpyPlot()
	// Registered the different python functions
	registerTopoPlot(pyPlotter);
	registerExtractRivers(pyPlotter);
	registerToporePlot(pyPlotter);
	addToLog("Registered python plotters")

	// Once this is ready, I can add the events
	document.getElementById('chooser_analysis').addEventListener('change',displayTheRightAnalysis);
	document.getElementById('extract_river').addEventListener('click',extract_and_plot_river);

	// finding the new div holding matplotlibs' figure
	allids = document.querySelectorAll('*[id]')
	for(el of allids){
		if(el.id.includes('styles')){continue;}
		if(el.id.includes('matplotlib') && idplotdiv === 'pyplotdiv'){
			idplotdiv = el.id;
		}
		else if(el.id.includes('matplotlib')){
			var comp =document.getElementById(idplotdiv);
			if(el.compareDocumentPosition(comp) & Node.DOCUMENT_POSITION_CONTAINED_BY){
				idplotdiv = el.id;
			}

		}
	}
	console.log("selected is " + idplotdiv)

	// Applying style to the figue
	document.getElementById(idplotdiv).style.position = "absolute";
	document.getElementById(idplotdiv).style.left = "0px";
	addToLog("Figure ready to plot data.")

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