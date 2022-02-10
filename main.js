// Starting up the engine:
document.querySelector("#statuspan").innerHTML = "Initialising ..."
document.querySelector("#statuspan").style.color = "red"	

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



// First I need to load pyodide engine
load_piodide_and_initialise_wasm()

document.getElementById('geotiff-file').onchange =  async function () {
    await load_DEM();
    await hillshade();
    console.log(dataCauldron)
    await pyPlotter.topoPlot(dataCauldron)
};



















































// #end of file


