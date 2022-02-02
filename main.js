
function plot_topo(thispy, this_mg)
{	

	const tttopo = mg.get_topo_array()
	window.ttopo = new Float32Array(mg.nx * mg.ny)
	for (var i = 0; i < ttopo.length; i++) {
      ttopo[i] = tttopo.get(i);
  } 

	window.nx = mg.nx
	window.ny = mg.ny


	thispy.runPython(`
import matplotlib.pyplot as plt
import numpy as np
from js import document
from js import ttopo, nx, ny

#Loading the topography

topo = np.asarray(ttopo.to_py()).reshape(ny,nx)

fig, ax = plt.subplots()

ax.imshow(topo, cmap = 'gist_earth')

# ordinary function to create a div
def create_root_element1(self):
	div = document.createElement('div')
	document.body.appendChild(div)
	return div
def create_root_element2(self):
	return document.getElementById('figure1')

#override create_root_element method of canvas by one of the functions above
fig.canvas.create_root_element = create_root_element1.__get__(
create_root_element1, fig.canvas.__class__)

fig.canvas.show()


	`);
};


// import loam from "./node_modules/loam/lib/loam";
document.querySelector("#statuspan").innerHTML = "Initialising ..."
document.querySelector("#statuspan").style.color = "red"	
// document.querySelector("#loader").style.color = "red"

let mg = null
let pyodide = null

async function load_piodide_and_initialise_wasm() {
  pyodide = await loadPyodide({
  indexURL: "https://cdn.jsdelivr.net/pyodide/dev/full/",
})
  await pyodide.loadPackage(["matplotlib", "numpy"]);

  // await loam.initialize("./","/node_modules/gdal-js/");
  await loam.initialize('', 'https://unpkg.com/gdal-js@2.0.0/');
  
  document.querySelector("#statuspan").innerHTML = "Ready!"
	document.querySelector("#statuspan").style.color = "green"
	document.querySelector("#loader").style.display = "inline-block"
	pyodide.runPython("print('YOLOBENGYOLOBENG')")

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
	document.querySelector("#statuspan").innerHTML = "Loading bytes in memory ..."
	console.log("Final should have " + nx*ny + " elements, dx,dy are " + dx + " " + dy)
	let Module = await createModule()

	let array = await ds.bytes()
	// console.log(array)
	array = array.buffer
	// console.log(array)
	array = new Float32Array(array)

	var vec = new Module.VectorFloat();
	var test = new Module.test_class(vec)
	test.mock_function(vec)

	vec.resize(array.length,0)
	// vec.resize(100,0)
  for (var i = 0; i < array.length; i++) {
  // for (var i = 0; i < 100; i++) {
      vec.set(i,array[i]);
  }
	console.log("0:" + array.length)
  mg = new Module.MinimalGraph(vec)
	document.querySelector("#statuspan").innerHTML = "Feeding webassembly ..."
	mg.set_dimensions(Number(nx), Number(ny), Number(nx * ny), Number(dx), Number(dx))
	mg.set_default_boundaries("4edges")
	let datvec = new Module.VectorFloat()
	// mg.mock_function(vec)
	mg.compute_graph("none")


	// array = Float32Array.from(array)
	// console.log(array)
	// console.log(`dimensions: nx ${nx} ny ${ny}`)
	// console.log(array)
	// console.log(ds)
	console.log("here")


	document.querySelector("#statuspan").innerHTML = "Ready!"
	document.querySelector("#statuspan").style.color = "green"
	document.querySelector("#plotter_topo").style.display = "inline-block"
	document.getElementById("plotfigtopo").onclick = plot_topo(pyodide)
	
}

document.getElementById('geotiff-file').onchange =  async function () {
    load_DEM();
};





