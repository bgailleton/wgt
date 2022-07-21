import numpy as np
from matplotlib import pyplot as plt


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "helvetica"




def fix_map_labels(fig, ax, dem, ndec = 2):
	'''
	Take a fig ax reprensenting a map and apply various filters to embellish the axis labels
	So far only set them to kilometers
	'''

	xlabels = list(ax.get_xticks())
	xlabelstr = []
	xlabelsnew = []
	for lab in xlabels:
		if(lab > dem.extent[0] and lab < dem.extent[1]):
			xlabelsnew.append(lab)
			xlabelstr.append(round(lab/1000, ndec))
	ax.set_xticks(xlabelsnew)
	ax.set_xticklabels(xlabelstr)

	ylabels = list(ax.get_yticks())
	ylabelstr = []
	ylabelsnew = []
	for lab in ylabels:
		if(lab > dem.extent[3] and lab < dem.extent[2]):
			ylabelsnew.append(lab)
			ylabelstr.append(round(lab/1000, ndec))
	ax.set_yticks(ylabelsnew)
	ax.set_yticklabels(ylabelstr)

































# Edn of file