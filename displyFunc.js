// This file contains all the function managing page display. For example when shall we show or not a page or another

const displayChooser = function(){
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
	console.log("STUFF")
	// First I am hiding all the other analyser
	document.querySelectorAll('.analyser').forEach((el)=>{el.style.display='none'});
	// Then checking which one is chosen
	const analysis = event.target.value

	if(analysis == "NaA")
	{
		displayChooser();
		return
	}

	if(analysis == "loadDEM")
	{
		document.querySelector('#seaLvlRemove').value = 0;
		document.querySelector("#loader").style.display = "inline-block";
		return;
	}

	if(analysis == "riverExt")
	{
		document.querySelector("#river_extraction").style.display = "inline-block";
		return;
	}
	
}


document.querySelector("#valueOfalphaHSAtLoad").innerHTML = document.querySelector("#alphaHSAtLoad").value

document.querySelector("#alphaHSAtLoad").addEventListener('change',(event) => {
	document.querySelector("#valueOfalphaHSAtLoad").innerHTML = document.querySelector("#alphaHSAtLoad").value
	console.log("DFKPDSJFLSKHDJGHJKFGJKH")
} ); 


function updateStatus(txt, col){
	document.querySelector("#statuspan").innerHTML = txt
	document.querySelector("#statuspan").style.color = col
}












































// end of file