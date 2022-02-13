addToLog("Starting the initialisation")

let mg = null;
let pyodide = null;
let mgModule = null;
let dataCauldron = {sirius : "yolo yolo beng beng"}
let paramCauldron = {
	area_threshold: Number(2e6),
	local_minima: "carve",

}
let pyPlotter = {}
let checker = {
	ready: false,
	demLoaded: false,
}
let idplotdiv = 'pyplotdiv'


document.getElementById('OKbuttonloaderPost').addEventListener('click', dealWithLoaderPost);

// First I need to load pyodide engine
load_piodide_and_initialise_wasm()

document.getElementById('OKbuttonloader').onclick =  async function () {

	if(document.getElementById('geotiff-file').files.length === 0){
		addToLog("No File Selected, Cannot load!")
		return;
	}

	await load_DEM();
	await hillshade();
	console.log(dataCauldron)
	await pyPlotter.topoPlot(dataCauldron)
};




















































// #end of file


