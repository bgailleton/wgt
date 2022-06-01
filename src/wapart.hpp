//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef wapart_HPP
#define wapart_HPP

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

#include <glm/vec2.hpp> // glm::vec3
#include <glm/vec3.hpp> // glm::vec3

#include <glm/glm.hpp>
#include "chonkutils.hpp"
#include "graph.hpp"

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


template<class f_t, class i_t>
class wapart2
{
public:
	i_t i;
	i_t rec;
	f_t S = 0;
	f_t volume = 1.;
	f_t sediments = 0.;

	Graph* graph;
	wapart2(){;};
	wapart2(Graph& g, i_t i)
	{
		this->graph = &g;
		this->i = i;
	}

	bool compute_sloperec(std::vector<float>& topo)
	{
		auto neighbours = this->graph->get_neighbours(this->i);
		int i_n = -1;
		std::vector<float> slopes;slopes.reserve(8);
		std::vector<int> recs;recs.reserve(8);
		for(auto n:neighbours)
		{
			if(topo[n.node] < topo[this->i])
			{
				recs.emplace_back(n.node);
				slopes.emplace_back((topo[this->i] - topo[n.node])/n.distance);
			}
		}

		if(slopes.size() == 0)
			return false;

		int randpos = floor(static_cast <float> (rand()) / (static_cast <float> (RAND_MAX))  * slopes.size());
		this->rec = recs[randpos];
		this->S = slopes[randpos];

		return true;
	}

	void move(){this->S = 0;this->i = this->rec;}

	void mass_transfer(f_t k_err, f_t exp_err, f_t k_dep, f_t exp_dep, f_t factor)
	{
		// f_t equ = this->S * k_eq;
		// f_t delta = this->sediments - equ;
		// this->sediments = this->sediments - delta;
		// // std::cout <<"|" << delta << "|" << this->S << "||";
		// if(this->sediments < 0)
		// {
		// 	this->sediments = 0;
		// }

		// this->graph->topography[this->i] += delta * factor; 

		f_t err = k_err * std::pow(this->S,exp_err);
		f_t dep = k_dep * this->sediments * std::exp(-this->sediments/exp_dep);

		if(dep > sediments)
			dep = this->sediments;

		this->sediments -= dep;
		this->sediments += err;
		this->graph->topography[this->i] += (dep - err) * factor;

	}

	void mass_transfer_SPL(float K, float m, float n, float L, float factor)
	{
		f_t err = K * std::pow(this->S,n) * std::pow(this->graph->area[this->i],m);
		f_t dep = this->sediments/(L * this->graph->area[this->i]/this->graph->cellarea);
		dep = (this->sediments >= dep) ? dep : this->sediments;
		this->sediments -= dep;
		this->sediments += err;
		if(this->sediments < 0)
			throw std::runtime_error("Galgabu");
		this->graph->topography[i] += (dep - err) * factor;
	}

	void fill()
	{
		auto neighbours = this->graph->get_neighbours(this->i);
		int i_n = -1;
		float higher = 0;
		for(auto n:neighbours)
		{
			if(this->graph->topography[n.node] > higher)
			{
				higher = this->graph->topography[n.node];
			}
		}
		this->graph->topography[this->i] = higher;
	}





};



template<class f_t, class i_t>
class wapart
{
public:

	i_t i = 0;
	glm::vec2 pos;
	glm::vec2 speed = glm::vec2(0.0);
	f_t volume = 1.;
	f_t sediment = 0.;
	f_t density = 1;
	f_t c_eq = 0;

	Graph* graph;

	wapart(){;};
	wapart(Graph& g, i_t i)
	{
		this->graph = &g;
		this->i = i;
		int row,col;
		this->graph->rowcol_from_node_id(i, row, col);
		pos.x = row;
		pos.y = col;
	};

	void set_pos(i_t ni)
	{
		int row,col;
		this->graph->rowcol_from_node_id(i, row, col);
		pos.x = row;
		pos.y = col;
	}

	void compute_speed(f_t dt, f_t friction, f_t scale)
	{
		glm::vec3 n = this->surfaceNormal(scale);
		this->speed += dt * glm::vec2(n.x, n.z)/(this->volume * this->density);
		this->pos += dt * this->speed;
		this->speed *= (1-dt*friction);
		this->i = this->graph->nodeid_from_row_col(this->pos.x, this->pos.y);
	}

	void compute_c_eq(float dz)
	{
		this->c_eq = this->volume * glm::length(this->speed) * dz;
		if(this->c_eq < 0) this->c_eq = 0;
	}

	void compute_ED(f_t dt, f_t depositionRate, f_t k_err)
	{
		f_t cdiff = k_err * (this->c_eq  - this->sediment);
		f_t gougne = dt * depositionRate * cdiff;
		this->sediment += gougne;
		this->graph->topography[this->i] -= gougne * this->volume;
	}

	void evaporate(f_t dt, f_t evaporation_rate){this->volume *= 1.0 - dt * evaporation_rate;}
  
  bool check_if_in_dem()
  {
  	if(this->pos.x < 0 ||this->pos.y < 0 || this->pos.x > this->graph->nx - 1 || this->pos.y > this->graph->ny - 1)
  		return false;
  	return true;
  }
	
	glm::vec3 surfaceNormal(float scale)
	{
	  /*
	    Note: Surface normal is computed in this way, because the square-grid surface is meshed using triangles.
	    To avoid spatial artifacts, you need to weight properly with all neighbors.
	  */

	  glm::vec3 n = glm::vec3(0.15) * glm::normalize(glm::vec3(scale*(this->graph->topography[this->i]-this->graph->topography[this->graph->get_bottom_index(this->i)]), 1.0, 0.0));  //Positive X
	  n += glm::vec3(0.15) * glm::normalize(glm::vec3(scale*(this->graph->topography[this->graph->get_top_index(this->i)]-this->graph->topography[this->i]), 1.0, 0.0));  //Negative X
	  n += glm::vec3(0.15) * glm::normalize(glm::vec3(0.0, 1.0, scale*(this->graph->topography[this->i]-this->graph->topography[this->graph->get_right_index(this->i)])));    //Positive Y
	  n += glm::vec3(0.15) * glm::normalize(glm::vec3(0.0, 1.0, scale*(this->graph->topography[this->graph->get_left_index(this->i)]-this->graph->topography[this->i])));  //Negative Y

	  //Diagonals! (This removes the last spatial artifacts)
	  n += glm::vec3(0.1) * glm::normalize(glm::vec3(scale*(this->graph->topography[this->i]-this->graph->topography[this->graph->get_topright_index(this->i)])/sqrt(2), sqrt(2), scale*(this->graph->topography[this->i]-this->graph->topography[this->graph->get_topright_index(this->i)])/sqrt(2)));    //Positive Y
	  n += glm::vec3(0.1) * glm::normalize(glm::vec3(scale*(this->graph->topography[this->i]-this->graph->topography[this->graph->get_topleft_index(this->i)])/sqrt(2), sqrt(2), scale*(this->graph->topography[this->i]-this->graph->topography[this->graph->get_topleft_index(this->i)])/sqrt(2)));    //Positive Y
	  n += glm::vec3(0.1) * glm::normalize(glm::vec3(scale*(this->graph->topography[this->i]-this->graph->topography[this->graph->get_bottomright_index(this->i)])/sqrt(2), sqrt(2), scale*(this->graph->topography[this->i]-this->graph->topography[this->graph->get_bottomright_index(this->i)])/sqrt(2)));    //Positive Y
	  n += glm::vec3(0.1) * glm::normalize(glm::vec3(scale*(this->graph->topography[this->i]-this->graph->topography[this->graph->get_bottomleft_index(this->i)])/sqrt(2), sqrt(2), scale*(this->graph->topography[this->i]-this->graph->topography[this->graph->get_bottomleft_index(this->i)])/sqrt(2)));    //Positive Y

	  return n;
	}




};

























































#endif