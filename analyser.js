// Bunch of utility functions to translate from JS to C++, ignore
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


// Interesting functions start here

// Load pyodide in the process
async function load_piodide_and_initialise_wasm() {

	// Loading pyodide
	addToLog("Initialising python webassembly engine")
  pyodide = await loadPyodide({
  indexURL: "https://cdn.jsdelivr.net/pyodide/v0.18.1/full/",
})
  await pyodide.loadPackage(["matplotlib", "numpy"]);
  addToLog("Python engine ready!")

  addToLog("Loading GDAL")
  await loam.initialize('', 'https://unpkg.com/gdal-js@2.0.0/');
  addToLog("GDAL Loaded")

  // Now initialising the plotters
	initPlotters(pyPlotter)

	// and displaying the chooser
	displayChooser()

	// checker fed
	checker.ready = true

}



// Function loading a DEM
const load_DEM = async function(){

  addToLog("Loading DEM:")

	const file = document.querySelector('#geotiff-file').files[0];
	document.querySelector("#loader").style.display = "none"
  addToLog("file is: " + file)

  addToLog("Converting format to binary representation of float32")
	let ds = await loam.open(file)
	ds = await ds.convert(['-ot','Float32', '-r', 'cubic', '-of', 'ENVI'])

	addToLog("Converted!")


	// No loading all the parameter in the cauldron
	const minZ = Number(document.querySelector("#seaLvlRemove").value)
	dataCauldron.seaLvl = minZ
	addToLog("Everything bellow " + minZ + " will be nodata.")

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
	dataCauldron.alphaHS = Number(document.querySelector("#alphaHSAtLoad").value)

	addToLog("Cauldron has all the parameters loaded")

	addToLog('Computing DEM ...')
	mgModule = await createModule()
	let array = await ds.bytes()
	array = array.buffer
	dataCauldron.topo = new Float32Array(array)
	var vec = new mgModule.VectorFloat();
  mg = new mgModule.MinimalGraph(vec)
	mg.set_dimensions(Number(nx), Number(ny), Number(nx * ny), Number(dx), Number(dx), Number(xmin), Number(ymin))
	mg.ingest_topo(dataCauldron.topo)
	mg.set_default_boundaries("4edges")
	mg.remove_seas(dataCauldron.seaLvl);
	checker.demLoaded = true
	checker.Acalc = false
		
	if(document.querySelector("#defaultCarvingCheck").value){
		addToLog("Computing graph info and resoving depressions ...")
		mg.compute_graph(paramCauldron.local_minima)
		checker.preproc = true
	}

	addToLog("DEM loaded and ready to go")

	displayChooser();
	
}


const extract_and_plot_river = async function(){
	// hideAn();
	if(checker.demLoaded === false){
		addToLog("Cannot extract rivers if no topography ingested. Please load or create DEM first.")
	}
	else if( checker.preproc === false){
		addToLog("Cannot extract river network if the dem is not preprocessed for flow routines")
	}
	else
	{
		const Ath = Number(document.querySelector("#AthRiverExtraction").value);
		console.log("Ath is ", Ath)
		dataCauldron.minsize = 1 
		dataCauldron.maxsize = 5
		dataCauldron.area_threshold = Number(Ath);
		// Checking if the number is too small
		if(dataCauldron.area_threshold < dataCauldron.dx * dataCauldron.dy * 20)

		if(checker.Acalc === false){
			mg.calculate_area();
			checker.Acalc = true;
		}
		mg.d_sources(dataCauldron.area_threshold);
		await mg.compute_river_nodes();

		const tempN = mg.getNriv();
		dataCauldron.nriv = null;
		dataCauldron.nriv = tempN;
		const tempX = mg.getXriv();
		dataCauldron.xriv = null;
		dataCauldron.xriv = tempX;
		console.log(tempX)
		const tempY = mg.getYriv();
		dataCauldron.yriv = null;
		dataCauldron.yriv = tempY;
		const tempA = mg.getAriv();
		dataCauldron.Ariv = null;
		dataCauldron.Ariv = tempA;
		// ok extracting rivers here
		pyPlotter.riverPlot(dataCauldron.xriv,dataCauldron.yriv,dataCauldron.Ariv, 0.1, 2);
	}
	// displayChooser();
}


const hillshade = async function(){
	console.log("hillshade")
	dataCauldron.HS = mg.get_HS()
	// dataCauldron.HS = cArrayFloat32FromOffset(offset, dataCauldron.nx * dataCauldron.ny)
	// console.log(dataCauldron.HS)
}

const replotTopoHs = async function(){
	// hideAn();
	dataCauldron.alphaHS = Number(document.querySelector("#alphaHSAtLoad_replot").value);
	await pyPlotter.toporePlot(dataCauldron);
	// displayChooser();
}



// const 







