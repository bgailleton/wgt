// This file contains all the function managing page display. For example when shall we show or not a page or another

const displayChooser = function(){
	document.querySelector("#chooser").style.display = "inline-block";
}

const hideChooser = function(){
	document.querySelector("#chooser").style.display = "none";
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
		document.querySelector("#loader").style.display = "inline-block";
		return;
	}
	
}














































// end of file