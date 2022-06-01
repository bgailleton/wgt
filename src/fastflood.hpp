//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef fastflood_HPP
#define fastflood_HPP

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

class FastFlood
{
	public:

		// A graph object
		Graph* graph;
		std::vector<float> precipitations;
		std::vector<float> hw;
		std::vector<float> base_topo;
		float rho = 1000;
		float g = 9.81;
		float manning = 0.033;
		float dt = 1e-1;
		float tolerance = 1e-6;


		// Constructor
		FastFlood(){}
		FastFlood(Graph& g){this->graph = &g;}

		void update_graph(Graph& g){this->graph = &g;}

		void set_dt(float dt){this->dt = dt;}

		// Set the precipitation raster
		void set_uniform_precipitation(float val)
		{
			this->precipitations = std::vector<float>(graph->nnodes,val);
		}

		void init()
		{
			this->hw = std::vector<float>(this->graph->nnodes_t,0.);
			this->base_topo = std::vector<float>(this->graph->nnodes_t,0.);
			this->base_topo = std::vector<float>(this->graph->topography);
		}

		void run()
		{
			// this->graph->topography = this->get_full_ht();
			// this->graph->compute_graph("fill");
			// this->hw = On_gaussian_blur(1,this->hw,this->graph->nx, this->graph->ny);
			this->fill_pits();
			// first calculating Qin
			std::vector<float> Qin = this->graph->accumulate_downstream_var(this->graph->cellarea, this->precipitations);
			float maxs = 0;
			for(int i = this->graph->nnodes - 1 ; i >= 0; --i)
			{
				int tnode = this->graph->stack[i];
				int rec = this->graph->receivers[tnode];

				if(this->graph->can_flow_even_go_there(tnode) == false || this->graph->can_flow_out_there(tnode) || this->get_ht(tnode) <= this->get_ht(rec) )
					continue;

				float tS = this->get_Sw(tnode);
				float Qout = this->graph->distance2receivers[tnode] * 1/this->manning * std::pow(this->hw[tnode],5/3) * std::pow(std::abs(tS),0.5);
				float dQ = Qin[tnode] - Qout;
				float dhw = dQ * this->dt / this->graph->cellarea;
				// if(dhw > 10)
				// 	std::cout << Qin[tnode] << "|" << Qout << std::endl;
				if(Qin[tnode] > maxs)
					maxs = Qin[tnode];

				this->hw[tnode] = fmax(dhw + this->hw[tnode],0);
			}
			// std::cout << "Max Qin was " << maxs << std::endl;
		}


		float get_ht(int tnode){return this->base_topo[tnode] + this->hw[tnode];}
		float get_Sw(int tnode){return (this->get_ht(tnode) - this->get_ht(this->graph->receivers[tnode]) ) / this->graph->distance2receivers[tnode];}
		std::vector<float> get_full_ht()
		{
			std::vector<float> fullht(this->graph->nnodes_t,this->graph->NDV);
			for(int i =0; i<this->graph->nnodes; ++i)
			{
				fullht[i] = this->get_ht(i);
			}
			return fullht;
		}

		void ingest_hw_other_dim(std::vector<float>& ohw, Graph& ograph)
		{
			this->init();
			for(int i = 0; i<this->graph->nnodes; ++i)
			{
				float tX, tY;
				this->graph->XY_from_nodeid(i,tX,tY);
				int j = ograph.nodeid_from_XY(tX,tY);
				this->hw[i] = ohw[j];
			}
		}

		void reinit_boundaries()
		{
			for(int i = 0; i<this->graph->nnodes; ++i)
			{
				if(this->graph->can_flow_out_there(i) || this->graph->can_flow_even_go_there(i) == false)
					this->hw[i] = 0;
			}
		}

		void fill_pits()
		{
			std::vector<float> ntopo(this->base_topo);
			for(int i=0; i<this->graph->nnodes; ++i)
				ntopo[i] = this->get_ht(i);
			this->graph->topography = ntopo;

			this->graph->compute_graph("fill");
			for(int i=0; i<this->graph->nnodes; ++i)
			{
				
				if(this->graph->can_flow_out_there(i))
					continue;

				int rec = this->graph->receivers[i];
				if(this->get_ht(i) <= this->get_ht(rec))
				{

					this->hw[i] += 1e-4 + this->get_ht(rec) - this->get_ht(i);
				}
			}
		}

#ifdef __EMSCRIPTEN__
#else
		py::array get_hw_np(){return py::array(this->hw.size(), this->hw.data());}
#endif

};























#endif