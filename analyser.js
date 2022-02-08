function cArrayFloat64(size) {
    var offset = mgModule._malloc(size * 8);
    mgModule.HEAPF64.set(new Float64Array(size), offset / 8);
    return {
        "data": mgModule.HEAPF64.subarray(offset / 8, offset / 8 + size),
        "offset": offset
    }
}

function cArrayFloat32(size) {
    var offset = mgModule._malloc(size * 4);
    mgModule.HEAPF32.set(new Float32Array(size), offset / 4);
    return {
        "data": mgModule.HEAPF32.subarray(offset / 4, offset / 4 +  size),
        "offset": offset
    }
}

function cArrayFloat32_(array) {
    var offset = mgModule._malloc(array.length * 4);
    mgModule.HEAPF32.set(array, offset / 4);
    return offset;
}

function cArrayFloat32FromOffset(offset,size) {
		array = new Float32Array(size);
    mgModule.HEAPF32.set(array, offset / 4);
    return  mgModule.HEAPF32.subarray(offset / 4, offset / 4 +  size);
}

async function load_piodide_and_initialise_wasm() {
  pyodide = await loadPyodide({
  indexURL: "https://cdn.jsdelivr.net/pyodide/v0.18.1/full/",
})
  await pyodide.loadPackage(["matplotlib", "numpy"]);

  // pyodide.runPython(await (await fetch("./pyscripts.py")).text());

  // await loam.initialize("./","/node_modules/gdal-js/");
  await loam.initialize('', 'https://unpkg.com/gdal-js@2.0.0/');
  
  document.querySelector("#statuspan").innerHTML = "Ready!"
	document.querySelector("#statuspan").style.color = "green"
	// document.querySelector("#loader").style.display = "inline-block"
	// document.querySelector("#chooser").style.display = "inline-block"

	// pyodide.runPython("print('YOLOBENGYOLOBENG')")

	// test_jsobj(meobj)
	registerPlotters(pyPlotter)

	displayChooser()

	checker.ready = true

}




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

	dataCauldron.nx = nx;
	dataCauldron.ny = ny;
	dataCauldron.dxy = dxy;
	dataCauldron.dx = dx;
	dataCauldron.dy = dy;
	dataCauldron.xmin = xmin;
	dataCauldron.ymin = ymin;
	dataCauldron.extents = [xmin, xmin +  (nx + 1) * dx,ymin, ymin + (ny + 1) * dy];

	document.querySelector("#statuspan").innerHTML = "Loading bytes in memory ..."
	console.log("Final should have " + nx*ny + " elements, dx,dy are " + dx + " " + dy)
	mgModule = await createModule()

	let array = await ds.bytes()
	// console.log(array)
	array = array.buffer
	// console.log(array)
	// array = new Float32Array(array)
	// let ptrCtopo = cArrayFloat32_(array)
	dataCauldron.topo = new Float32Array(array)
	// console.log(array)

	var vec = new mgModule.VectorFloat();

	// vec.resize(array.length,0)
 //  for (var i = 0; i < array.length; i++) {
	// 	vec.set(i,array[i]);
 //  }

  mg = new mgModule.MinimalGraph(vec)
	document.querySelector("#statuspan").innerHTML = "Feeding webassembly ..."
	// console.log(array.data)
	// console.log(Number(array.data.length))
	// mg.ingest_topo_from_C(ptrCtopo, nx * ny)
	mg.set_dimensions(Number(nx), Number(ny), Number(nx * ny), Number(dx), Number(dx), Number(xmin), Number(ymin))
	mg.ingest_topo(dataCauldron.topo)


	mg.set_default_boundaries("4edges")
	let datvec = new mgModule.VectorFloat()
	
	document.querySelector("#statuspan").innerHTML = "Computing graph ..."
	
	mg.compute_graph(paramCauldron.local_minima)
	document.querySelector("#statuspan").innerHTML = "Ready!"
	document.querySelector("#statuspan").style.color = "green"

	checker.loaded = true

	displayChooser();
	
}



const hillshade = async function(){
	console.log("hillshade")
	dataCauldron.HS = mg.get_HS()
	// dataCauldron.HS = cArrayFloat32FromOffset(offset, dataCauldron.nx * dataCauldron.ny)
	console.log(dataCauldron.HS)


}


// const 







