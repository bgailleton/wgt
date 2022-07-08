//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef RASTERDATA_HPP
#define RASTERDATA_HPP

// STL imports
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <ctime>
#include <fstream>
#include <queue>
#include <stack>
#include <iostream>
#include <numeric>
#include <cmath>
#include <initializer_list>
#include <thread>
#include<stdlib.h>
#include<ctime>

#include "chonkutils.hpp"
#include "graph.hpp"
#include "wapart.hpp"
#include "PerlinNoise.hpp"

#include "npy.hpp"
#ifdef __EMSCRIPTEN__
	#include <emscripten.h>
#else
	#include <pybind11/pybind11.h>
	#include <pybind11/stl.h>
	#include <pybind11/numpy.h>
#endif

#ifdef __EMSCRIPTEN__
#else
	namespace py = pybind11;
#endif





// Small class emulating some aspect of a vector vector from its pointer
template<class T>
class pvector
{
public:
	std::vector<T>* data;

	pvector(){};
	pvector(std::vector<T>& dat){this->data = &dat;}
	pvector(std::vector<T>* dat){this->data = dat;}
	size_t size(){return this->data->size();}
	T operator [](int i) const    {return (*this->data)[i];}
  T & operator [](int i) {return (*this->data)[i];}
  void emplace_back(T& tin){this->data->emplace_back(tin);}
  void push_back(T& tin){this->data->push_back(tin);}
  void shrink_to_fit(){this->data->shrink_to_fit();}
  void clear(){this->data->clear();}

}


template<class Neighbourer, class Holder>
class RasterData
{
public:

	Holder data;

	RasterData(){}
	RasterData(Holder& ho){this->data = ho;}

	T operator [](int i) const    {return this->data[i];}
  T & operator [](int i) {return this->data[i];}

};







































#endif