
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


// import loam from "./node_modules/loam/lib/loam";
document.querySelector("#statuspan").innerHTML = "Initialising ..."
document.querySelector("#statuspan").style.color = "red"	
// document.querySelector("#loader").style.color = "red"

let mg = null
let pyodide = null

async function load_piodide_and_initialise_wasm() {
  pyodide = await loadPyodide({
  indexURL: "https://cdn.jsdelivr.net/pyodide/v0.18.1/full/",
})
  await pyodide.loadPackage(["matplotlib", "numpy"]);

  // await loam.initialize("./","/node_modules/gdal-js/");
  await loam.initialize('', 'https://unpkg.com/gdal-js@2.0.0/');
  
  document.querySelector("#statuspan").innerHTML = "Ready!"
	document.querySelector("#statuspan").style.color = "green"
	document.querySelector("#loader").style.display = "inline-block"
	// pyodide.runPython("print('YOLOBENGYOLOBENG')")

  // console.log(mg)

}

// First I need to load pyodide engine
load_piodide_and_initialise_wasm()

const load_DEM = async function(){
	const file = document.querySelector('#geotiff-file').files[0];
	document.querySelector("#loader").style.display = "none"

	document.querySelector("#statuspan").innerHTML = "Loading file ..."
	document.querySelector("#statuspan").style.color = "orange"

	let ds = await loam.open(file)
	document.querySelector("#statuspan").innerHTML = "Converting to right format ..."

	ds = await ds.convert(['-ot','Float32', '-r', 'cubic', '-of', 'ENVI'])

	document.querySelector("#statuspan").innerHTML = "Reading info ..."

	const nx = await ds.width()
	const ny = await ds.height()
	const dxy = await ds.transform()
	const dx = Math.abs(dxy[1])
	const dy = Math.abs(dxy[5])
	const xmin = Math.abs(dxy[0])
	const ymin = Math.abs(dxy[3]) - dy * (ny + 1)

	document.querySelector("#statuspan").innerHTML = "Loading bytes in memory ..."
	console.log("Final should have " + nx*ny + " elements, dx,dy are " + dx + " " + dy)
	let Module = await createModule()

	let array = await ds.bytes()
	// console.log(array)
	array = array.buffer
	// console.log(array)
	array = new Float32Array(array)

	var vec = new Module.VectorFloat();

	vec.resize(array.length,0)
  for (var i = 0; i < array.length; i++) {
		vec.set(i,array[i]);
  }

  mg = new Module.MinimalGraph(vec)
	document.querySelector("#statuspan").innerHTML = "Feeding webassembly ..."
	mg.set_dimensions(Number(nx), Number(ny), Number(nx * ny), Number(dx), Number(dx), Number(xmin), Number(ymin))
	mg.set_default_boundaries("4edges")
	let datvec = new Module.VectorFloat()
	document.querySelector("#statuspan").innerHTML = "Computing graph ..."
	mg.compute("carve")

	document.querySelector("#statuspan").innerHTML = "Ready!"
	document.querySelector("#statuspan").style.color = "green"
	document.querySelector("#plotter_topo").style.display = "inline-block"


	document.getElementById("plotfigtopo").onclick = plot_topo()
	
}

document.getElementById('geotiff-file').onchange =  async function () {
    load_DEM();
};





