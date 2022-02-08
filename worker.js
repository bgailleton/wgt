var initialized = false;

var GDALOpen,
    GDALClose,
    GDALGetRasterCount,
    GDALTranslate,
    GDALTranslateOptionsNew,
    GDALTranslateOptionsFree,
    CSLCount;
    
// Set up Module object for gdal.js to populate. Emscripten sets up its compiled
// code to look for a Module object in the global scope. If found, it reads runtime
// configuration from the existing object, and then further populates that object
// with other helpful functionality (e.g. ccall() and cwrap(), which are used in
// the onRuntimeInitialized callback, below).
var Module = {
    'print': function(text) { console.log('stdout: ' + text); },
    'printErr': function(text) { console.log('stderr: ' + text); },
    // Optimized builds contain a .js.mem file which is loaded asynchronously;
    // this waits until that has finished before performing further setup.
    'onRuntimeInitialized': function() {
        // Initialize GDAL
        Module.ccall('GDALAllRegister', null, [], []);

        // Set up JS proxy functions
        // Note that JS Number types are used to represent pointers, which means that
        // any time we want to pass a pointer to an object, such as in GDALOpen, which in
        // C returns a pointer to a GDALDataset, we need to use 'number'.
        GDALOpen = Module.cwrap('GDALOpen', 'number', ['string']);
        GDALClose = Module.cwrap('GDALClose', 'number', ['number']);

        // GDALGetRasterCount = Module.cwrap('GDALGetRasterCount', 'number', ['number']);
        // // Params:
        // //  1. Output path
        // //  2. Pointer to a GDALDataset
        // //  3. Pointer to a GDALTranslateOptions
        // //  4. Int to use for error reporting
        // // Returns a pointer to a new GDAL Dataset
        // GDALTranslate = Module.cwrap('GDALTranslate', 'number', ['string', 'number', 'number', 'number']);
        // // Params: array of option strings as to gdal_translate; pointer to a struct that should be null.
        // GDALTranslateOptionsNew = Module.cwrap('GDALTranslateOptionsNew', 'number', ['number', 'number']);
        // GDALTranslateOptionsFree = Module.cwrap('GDALTranslateOptionsFree', 'number', ['number']);

        // CSLCount = Module.cwrap('CSLCount', 'number', ['number']);
        initialized = true;
    }
};







