// This file contains all the function managing page display. For example when shall we show or not a page or another

const displayChooser = function(){
	// Display the main choosing platform for 
	document.querySelector("#chooser").style.display = "inline-block";
	document.getElementById('chooser_analysis').value = "NaA"
}

const hideChooser = function(){
	document.querySelector("#chooser").style.display = "none";
}

const hideAn = function(){
	document.querySelectorAll('.analyser').forEach((el)=>{el.style.display='none'});
}


const displayTheRightAnalysis = function(event){

	// First I am hiding all the other analyser
	document.querySelectorAll('.analyser').forEach((el)=>{el.style.display='none'});
	// Then checking which one is chosen
	const analysis = event.target.value

	if(analysis === "NaA")
	{
		displayChooser();
		return
	}

	if(analysis === "loadDEM")
	{
		if(checker.demLoaded === true){
			document.querySelector("#loaderPost").style.display = "inline-block";
			return;
		}
		document.querySelector('#seaLvlRemove').value = 0;
		document.querySelector("#loader").style.display = "inline-block";
  	addToLog("Displaying raster loader")
		return;
	}

	if(analysis === 'replotDEM'){
		if(checker.demLoaded === true){
			document.querySelector("#rePlotBase").style.display = "inline-block";
			return;
		}
		else{
			addToLog("Need to load a DEM before replotting it yo.")
			displayChooser();
		}
	}

	if(analysis === "riverExt")
	{
  	addToLog("Displaying River Extractor")
		document.querySelector("#river_extraction").style.display = "inline-block";
		return;
	}
	
}


document.querySelector("#valueOfalphaHSAtLoad").innerHTML = document.querySelector("#alphaHSAtLoad").value

document.querySelector("#alphaHSAtLoad").addEventListener('input',(event) => {
	let val = new Intl.NumberFormat('en-uk', { minimumFractionDigits: 2, maximumFractionDigits: 2 }).format(document.querySelector("#alphaHSAtLoad").value);
	document.querySelector("#valueOfalphaHSAtLoad").innerHTML = val;
} ); 

document.querySelector("#valueOfalphaHSAtLoad_replot").innerHTML = document.querySelector("#alphaHSAtLoad_replot").value

document.querySelector("#alphaHSAtLoad_replot").addEventListener('input',(event) => {
	let val = new Intl.NumberFormat('en-uk', { minimumFractionDigits: 2, maximumFractionDigits: 2 }).format(document.querySelector("#alphaHSAtLoad_replot").value);
	document.querySelector("#valueOfalphaHSAtLoad_replot").innerHTML = val;
} );


async function addToLog(txt){
	p = document.createElement("p");
	p.innerHTML = txt
	document.querySelector("#logTop").appendChild(p)

	var div = document.querySelector("#logTop");
	div.scrollTop = div.scrollHeight - div.clientHeight;

}


function dealWithLoaderPost(){
	console.log("dfljasldfgjlkjsd")
	hideAn();
	displayChooser();
}





function OKBtn(){
	addToLog("Done with action, back to menu")
	hideAn();
	displayChooser();
}





































// end of file