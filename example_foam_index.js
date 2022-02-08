/* global loam */
// import loam
// import {yoloyolobengbeng} from './testmpl.js';
// console.log("Loading python enjiugkuygkuygjhgine ...");

// async function main() {
//   let pyodide = await loadPyodide({
//     indexURL: "https://cdn.jsdelivr.net/pyodide/v0.18.1/full/",
//   });
//   return pyodide;
// }
// const pyodide = main();
// pyodide.runPython(`print('baft'`);

function yoloyolobengbeng(thispy)
{
	thispy.runPython(`
import matplotlib.pyplot as plt
import numpy as np
from js import document

fig, ax = plt.subplots()

X = np.arange(0,10,0.01)
Y = np.sin(X) + np.random.rand(*X.shape)

ax.plot(X,Y, color = "red")

ax.set_xlabel("Bite")
ax.set_ylabel("LOL")

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

// Use the locally built version of loam, with a CDN copy of GDAL from unpkg.
loam.initialize('', 'https://unpkg.com/gdal-js@2.0.0/');

const EPSG4326 =
    'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]';

function displayInfo() {
    const file = document.querySelector('#geotiff-file').files[0];
    const displayElem = document.getElementById('gdalinfo');

    // Clear display text
    displayElem.innerText = '';
    // Use Loam to get GeoTIFF metadata
    loam.open(file).then((ds) => {
        return Promise.all([ds.width(), ds.height(), ds.count(), ds.wkt(), ds.transform()]).then(
            ([width, height, count, wkt, geoTransform]) => {
                displayElem.innerText +=
                    'Size: ' + width.toString() + ', ' + height.toString() + '\n';
                displayElem.innerText += 'Band count: ' + count.toString() + '\n';
                displayElem.innerText += 'Coordinate system:\n' + wkt + '\n';

                const cornersPx = [
                    [0, 0],
                    [width, 0],
                    [width, height],
                    [0, height],
                ];
                const cornersGeo = cornersPx.map(([x, y]) => {
                    return [
                        // http://www.gdal.org/gdal_datamodel.html
                        geoTransform[0] + geoTransform[1] * x + geoTransform[2] * y,
                        geoTransform[3] + geoTransform[4] * x + geoTransform[5] * y,
                    ];
                });

                loam.reproject(wkt, EPSG4326, cornersGeo).then((cornersLngLat) => {
                    displayElem.innerText += 'Corner Coordinates:\n';
                    cornersLngLat.forEach(([lng, lat], i) => {
                        displayElem.innerText +=
                            '(' +
                            cornersGeo[i][0].toString() +
                            ', ' +
                            cornersGeo[i][1].toString() +
                            ') (' +
                            lng.toString() +
                            ', ' +
                            lat.toString() +
                            ')\n';
                    });
                });
            }
        );
    });
}

document.getElementById('geotiff-file').onchange = function () {
    displayInfo();
};

document.getElementById('yoloer').onclick = async function() {
	console.log("yolo")
	// let pyodide = main();
	let pyodide = await initPyodide;
	yoloyolobengbeng(pyodide);
};



