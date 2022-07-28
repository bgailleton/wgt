//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
/*
fastflood.hpp -> c++ side of the fastflood algorithm
Older attemps at single flow directions are through a fully fledged c++ object
Multiple flow are functions called and monitored from python
B.G.
*/
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
#include "mgraph.hpp"
#include "SMgraph.hpp"
#include "npy.hpp"
#include "numvec.hpp"


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


// Multiple flow version of Fastflood
// Serial version
// 1) calculates hydrolic slope, sqrt and sum for the weighting
// 2) Calculate simultaneously from upstream to downstream the Qin, Qout and dhw
// returns the water diff
template<class Neighbourer_t,class topo_t, class T, class out_t>
void run_multi_fastflood_static(SMgraph& graph, Neighbourer_t& neighbourer, topo_t& thw, topo_t& ttopography, 
	T manning, topo_t& tprecipitations, topo_t& tQin, topo_t& tQout, topo_t& tSw, topo_t& tsumslopes)
{
	// preformat the fluxes
	// (does no cost a lot of performance, just makes sure everything is vector-like through pointer magic if needed)
	auto hw = format_input(thw);
	auto topography = format_input(ttopography);
	auto precipitations = format_input(tprecipitations);
	auto Qin = format_input(tQin);
	auto Qout = format_input(tQout);
	auto Sw = format_input(tSw);
	auto sumslopes = format_input(tsumslopes);


	// First iteration to calculates the hydraulic slope in a SMG manner (should move that part to SMG btw)
	for(int i = 0; i< neighbourer.nnodes; ++i)
	{
		// Taking advantage of this loop to reinit discharge to 0
		Qin[i] = 0;
		Qout[i] = 0;
	}
	

	// Now i have my hydraulic slope precalculated
	// I can iterate through the reverse stack (upstream to downstream)
	// this main loop calculates Qi Qout/Qout and propagate Qi downstream
	// std::vector<int> gabun(neighbourer.nnodes,0);
	// double QW_OUT = 0, WEIGHTS = 0;
	int totnrecs = 0,nnorec = 0;
	for(int i = graph.nnodes-1; i>=0; --i)
	{
		// Getting the next upstreamest node
		int node = graph.stack[i];
		// gabun[node] += 1;
		// if nodata -> ignore
		if(neighbourer.is_active(node) == false)
		{
			Qin[node] = 0;
			continue;
		}
		
		// Accumulating local disacharge (precipitation * area)
		Qin[node] += precipitations[node] * neighbourer.cellarea;

		// getting list of receivers
		auto recs = graph.get_receiver_link_indices(node,neighbourer);
		if(recs.size() == 0)
		{
			nnorec += 1;
			continue;
		}
		// 	throw std::runtime_error("norecs where should");
		// std::cout << recs.size() << '|';
		totnrecs += recs.size();
		// double sumsl = 0;
		// for(auto rec:recs)
		// 	sumsl += std::pow(Sw[rec.second],2);

		// I'll need to store the maximum slope and the sum of sts
		double maxslope = 1e-6;
		double sumst = 0;
		// double sumw = 0;
		for(auto rec:recs)
		{
		// std::cout << "FF::DEBUG::2.51::" << rec.first <<"|" << rec.second << std::endl;
			if(rec.first < 0 || rec.first >= neighbourer.nnodes)
				continue;
			double weight = Sw[rec.second]/sumslopes[node];
			// double weight = std::pow(Sw[rec.second],2)/sumsl;
			// sumw += weight;
			if(weight > 1)
				std::cout << "W>1::" << Sw[rec.second] << "|" << sumslopes[node] << std::endl;
			// WEIGHTS+= weight;

			if(neighbourer.is_active(rec.first))
			{
				Qin[rec.first] += Qin[node] * weight;
			}
			// else
			// {
			// 	QW_OUT += Qin[node] * weight;
			// }
			
			sumst += std::pow(Sw[rec.second],2); // POW 2 because Sw it is already sqrted
			if(Sw[rec.second] > maxslope)
				maxslope = Sw[rec.second];
		}

		// if(std::abs(sumw - 1) >1e-2 )
		// 	throw std::runtime_error("SUMW ERROR::" + std::to_string(sumw));

		Qout[node] = neighbourer.dx * 1/manning * 1/recs.size() * std::pow(hw[node],5/3) * sumst/maxslope;// / std::sqrt(2); // Note that smst is the sum of slopes and maxslope is squarerooted
		
		if(std::isfinite(Qout[node]) == false)
		{
			std::cout << recs.size() << "|" << std::pow(hw[node],5/3) << "|" << maxslope << "|" << sumst << "|node " << node <<"|Srec" << graph.Sreceivers[node] << "|" << (topography[node] + hw[node] - topography[graph.Sreceivers[node]] - hw[graph.Sreceivers[node]])  << std::endl;
		}
	}

	std::cout << "N_no recs = " << nnorec << std::endl;
	// std::cout << "FLUXES OUT = " << QW_OUT << " and weights " << WEIGHTS << " TOTNRECS " << totnrecs << " vs " << graph.isrec.size() << std::endl;
}

// Multiple flow version of Fastflood
// Serial version
// 1) calculates hydrolic slope, sqrt and sum for the weighting
// 2) Calculate simultaneously from upstream to downstream the Qin, Qout and dhw
// returns the water diff
template<class Neighbourer_t,class topo_t, class T, class out_t>
void run_multi_fastflood_static_ext_Qwin(SMgraph& graph, Neighbourer_t& neighbourer, topo_t& thw, topo_t& ttopography, 
	T manning, topo_t& tprecipitations, topo_t& tQin, topo_t& tQout, topo_t& tSw, topo_t& tsumslopes)
{
	// preformat the fluxes
	// (does no cost a lot of performance, just makes sure everything is vector-like through pointer magic if needed)
	auto hw = format_input(thw);
	auto topography = format_input(ttopography);
	auto precipitations = format_input(tprecipitations);
	auto Qin = format_input(tQin);
	auto Qout = format_input(tQout);
	auto Sw = format_input(tSw);
	auto sumslopes = format_input(tsumslopes);


	// First iteration to calculates the hydraulic slope in a SMG manner (should move that part to SMG btw)
	for(int i = 0; i< neighbourer.nnodes; ++i)
	{
		// Taking advantage of this loop to reinit discharge to 0
		Qout[i] = 0;
	}
	

	// Now i have my hydraulic slope precalculated
	// I can iterate through the reverse stack (upstream to downstream)
	// this main loop calculates Qi Qout/Qout and propagate Qi downstream
	// std::vector<int> gabun(neighbourer.nnodes,0);
	// double QW_OUT = 0, WEIGHTS = 0;
	int totnrecs = 0,nnorec = 0;
	for(int i = graph.nnodes-1; i>=0; --i)
	{
		// Getting the next upstreamest node
		int node = graph.stack[i];
		// gabun[node] += 1;
		// if nodata -> ignore
		if(neighbourer.is_active(node) == false)
		{
			Qin[node] = 0;
			continue;
		}
		

		// getting list of receivers
		auto recs = graph.get_receiver_link_indices(node,neighbourer);
		if(recs.size() == 0)
		{
			nnorec += 1;
			continue;
		}
		// 	throw std::runtime_error("norecs where should");
		// std::cout << recs.size() << '|';
		totnrecs += recs.size();
		// double sumsl = 0;
		// for(auto rec:recs)
		// 	sumsl += std::pow(Sw[rec.second],2);

		// I'll need to store the maximum slope and the sum of sts
		double maxslope = 1e-6;
		double sumst = 0;
		// double sumw = 0;
		for(auto rec:recs)
		{
		// std::cout << "FF::DEBUG::2.51::" << rec.first <<"|" << rec.second << std::endl;
			if(rec.first < 0 || rec.first >= neighbourer.nnodes)
				continue;
			
			sumst += std::pow(Sw[rec.second],2); // POW 2 because Sw it is already sqrted
			if(Sw[rec.second] > maxslope)
				maxslope = Sw[rec.second];
		}

		Qout[node] = neighbourer.dx * 1/manning * 0.5/recs.size() * std::pow(hw[node],5/3) * sumst/maxslope;// / std::sqrt(2); // Note that smst is the sum of slopes and maxslope is squarerooted
		// Qout[node] = neighbourer.dx * 1/manning * 0.5 * std::pow(hw[node],5/3) * sumst/maxslope;// / std::sqrt(2); // Note that smst is the sum of slopes and maxslope is squarerooted
		
		if(std::isfinite(Qout[node]) == false)
		{
			std::cout << recs.size() << "|" << std::pow(hw[node],5/3) << "|" << maxslope << "|" << sumst << "|node " << node <<"|Srec" << graph.Sreceivers[node] << "|" << (topography[node] + hw[node] - topography[graph.Sreceivers[node]] - hw[graph.Sreceivers[node]])  << std::endl;
		}
	}

	std::cout << "N_no recs = " << nnorec << std::endl;
	// std::cout << "FLUXES OUT = " << QW_OUT << " and weights " << WEIGHTS << " TOTNRECS " << totnrecs << " vs " << graph.isrec.size() << std::endl;
}


template<class Neighbourer_t,class topo_t, class T, class out_t>
void compute_Sw_sumslopes(SMgraph& graph, Neighbourer_t& neighbourer, topo_t& thw, topo_t& ttopography, topo_t& tSw, topo_t& tsumslopes)
{

	auto hw = format_input(thw);
	auto topography = format_input(ttopography);
	auto Sw = format_input(tSw);
	auto sumslopes = format_input(tsumslopes);

	// First iteration to calculates the hydraulic slope in a SMG manner (should move that part to SMG btw)
	for(size_t i = 0; i< graph.isrec.size(); ++i)
	{
		
		int j = graph.links[i*2];
		int k = graph.links[i*2 + 1];
		
		if(j < 0)
			continue;

		if(graph.isrec[i] == false)
			std::swap(j,k);

		Sw[i] = std::sqrt(std::max(std::abs(topography[j] + hw[j] - topography[k] - hw[k])/neighbourer.get_dx_from_isrec_idx(i), 1e-6));
		sumslopes[j] += Sw[i];

	}

}


template<class Neighbourer_t,class topo_t, class T, class out_t>
out_t run_multi_fastflood_static_OMP(SMgraph& graph, Neighbourer_t& neighbourer, topo_t& thw, topo_t& ttopography, T manning, topo_t& tprecipitations, T dt)
{
	// init the fluxes
	auto hw = format_input(thw);
	auto topography = format_input(ttopography);
	auto precipitations = format_input(tprecipitations);
	std::vector<double> diff(hw.size(),0.), Qin(hw.size(),0.), sumslopes(hw.size(),0.), Sw(graph.isrec.size(),1e-5);
	std::vector<int> nrecs(hw.size(),0);

	py::gil_scoped_release release;

	// std::cout << "FF::DEBUG::1" << std::endl;
	#pragma omp parallel for num_threads(4)
	for(int i = 0; i< neighbourer.nnodes; ++i)
	{

		// cannot be a neighbour anyway, abort
		if(neighbourer.can_flow_even_go_there(i) == false)
			continue;

		// double this_topo = topography[i];

		int n = neighbourer.get_id_right_SMG(i);
		int tn = neighbourer.get_right_index(i);
		if(n >=0 && n<neighbourer.nnodes * 4 && tn >= 0 && tn < neighbourer.nnodes)
		{
			Sw[n] = std::max(std::sqrt(std::abs(topography[i] + hw[i] - topography[tn] - hw[tn])/neighbourer.dx), 1e-6);
			sumslopes[tn] += Sw[n];

		}
		
		n = neighbourer.get_id_bottomright_SMG(i);
		tn = neighbourer.get_bottomright_index(i);
		if(n >=0 && n<neighbourer.nnodes * 4 && tn >= 0 && tn < neighbourer.nnodes)
		{
			Sw[n] = std::max(std::sqrt(std::abs(topography[i] + hw[i] - topography[tn] - hw[tn])/neighbourer.dxy), 1e-6);
			sumslopes[tn] += Sw[n];

		}

		n = neighbourer.get_id_bottom_SMG(i);
		tn = neighbourer.get_bottom_index(i);
		if(n >=0 && n<neighbourer.nnodes * 4 && tn >= 0 && tn < neighbourer.nnodes)
		{
			Sw[n] = std::max(std::sqrt(std::abs(topography[i] + hw[i] - topography[tn] - hw[tn])/neighbourer.dy), 1e-6);
			sumslopes[tn] += Sw[n];

		}

		n = neighbourer.get_id_bottomleft_SMG(i);
		tn = neighbourer.get_bottomleft_index(i);
		if(n >=0 && n<neighbourer.nnodes * 4 && tn >= 0 && tn < neighbourer.nnodes)
		{
			Sw[n] = std::max(std::sqrt(std::abs(topography[i] + hw[i] - topography[tn] - hw[tn])/neighbourer.dxy), 1e-6);
			sumslopes[tn] += Sw[n];

		}
	}
	std::cout << "FF::DEBUG::2" << std::endl;
	for(auto v:Sw)
	{
		if(std::isfinite(v) == false)
			throw std::runtime_error("YOLO");
	}

	for(int i = graph.nnodes-1; i>=0; --i)
	{
		// std::cout << "FF::DEBUG::2.1::" << i << std::endl;
		int node = graph.stack[i];

		// std::cout << "FF::DEBUG::2.2" << std::endl;
		if(neighbourer.can_flow_even_go_there(node) == false)
			continue;
		
		// std::cout << "FF::DEBUG::2.3" << std::endl;
		Qin[node] += precipitations[node] * neighbourer.cellarea;
		// std::cout << "FF::DEBUG::2.4::" << node << std::endl;

		

		auto recs = graph.get_receiver_link_indices(node,neighbourer);
		// std::cout << "FF::DEBUG::2.5" << std::endl;
		
		for(auto rec:recs)
		{
		// std::cout << "FF::DEBUG::2.51::" << rec.first <<"|" << rec.second << std::endl;
			if(rec.first < 0 || rec.first >= neighbourer.nnodes)
				continue;

			Qin[rec.first] += Qin[node] * Sw[rec.second]/sumslopes[rec.first];

			++nrecs[rec.first];
		}
		// std::cout << "FF::DEBUG::2.6" << std::endl;

	}



	// std::cout << "FF::DEBUG::3" << std::endl;


	#pragma omp parallel for num_threads(4)
	for(int i=0; i<neighbourer.nnodes; ++i)
	{
		diff[i] =  Qin[i] - 1/manning * std::pow(hw[i],5/3) * std::sqrt(sumslopes[i])/nrecs[i];
	}

	// std::cout << "FF::DEBUG::4" << std::endl;
	py::gil_scoped_acquire acquire;

	return format_output(diff);

}




// ==================================================
// ==================================================
// ==================================================
// OLDER TESTS BELLOW
// ==================================================
// ==================================================
// ==================================================




template<class Neighbourer_t,class topo_t, class T, class out_t>
out_t run_multi_fastflood_static_old(MGraph<T>& graph, Neighbourer_t& neighbourer, topo_t& thw, topo_t& ttopography, T manning, topo_t& tprecipitations, T pcoeff, T dt)
{
	// init the fluxes
	auto hw = format_input(thw);
	auto topography = format_input(ttopography);
	auto precipitations = format_input(tprecipitations);
	std::vector<double> diff(hw.size(),0.), Qin(hw.size(),0.);

	// From upstream to downstream
	for(int i = graph.nnodes-1; i>=0; --i)
	{
		int node = graph.stack[i];
		if(neighbourer.is_active(node) == false)
			continue;

		if(graph.receivers.size() > 0)
		{

			Qin[node] += precipitations[node] * neighbourer.cellarea;

			std::vector<T> weights(graph.receivers[node].size(),0.), slopes(graph.receivers[node].size());
			T sumslopes = 0, fact = 0;
			if(pcoeff == -1)
			{
				for(int j = 0;j < graph.receivers[node].size(); ++j)
				{
					int rec = graph.receivers[node][j];
					slopes[j] = std::sqrt((topography[node] + hw[node] - topography[rec] - hw[rec])/graph.distance2receivers[node][j]);
					if(slopes[j] <= 0)
						slopes[j] = 1e-5;
					sumslopes += slopes[j];
				}

				for(int j = 0;j < graph.receivers[node].size(); ++j)
				{
					int rec = graph.receivers[node][j];
					weights[j] = slopes[j]/sumslopes;
					fact += weights[j] * std::sqrt(slopes[j]);
					Qin[rec] += Qin[node] * weights[j];
				}
			}

			diff[node] =  Qin[node] - 1/manning * std::pow(hw[node],5/3) * fact;
		}

		else
		{
			auto neighbours = neighbourer.get_neighbours_only_id(node);
			T topomax = topography[node];
			for(auto k:neighbours)
			{
				if(topography[k] > topomax)
					topomax = topography[k];
			}
			diff[node] = (topomax + 1e-3)/dt;
		}

	}

	return format_output(diff);

}


// class FastFlood
// {
// 	public:

// 		// A graph object
// 		Graph* graph;
// 		std::vector<float> precipitations;
// 		std::vector<float> hw;
// 		std::vector<float> base_topo;
// 		float rho = 1000;
// 		float g = 9.81;
// 		float manning = 0.033;
// 		float dt = 1e-1;
// 		float tolerance = 1e-6;
// 		bool allow_neg = true;
// 		float hydrcof = 1e-2;
// 		float hydrex = 0.5;
// 		float Qwinth = 0;


// 		// Constructor
// 		FastFlood(){}
// 		FastFlood(Graph& g){this->graph = &g;}
// 		void set_allow_neg(bool tn){this->allow_neg = tn;}

// 		void update_graph(Graph& g){this->graph = &g;}

// 		void set_dt(float dt){this->dt = dt;}

// 		// Set the precipitation raster
// 		void set_uniform_precipitation(float val)
// 		{
// 			this->precipitations = std::vector<float>(graph->nnodes,val);
// 		}

// 		void init()
// 		{
// 			this->hw = std::vector<float>(this->graph->nnodes_t,0.);
// 			this->base_topo = std::vector<float>(this->graph->nnodes_t,0.);
// 			this->base_topo = std::vector<float>(this->graph->topography);
// 		}

// 		void run()
// 		{
			
// 			std::vector<float> ntopo(this->base_topo);
// 			for(int i=0; i<this->graph->nnodes; ++i)
// 				ntopo[i] = this->get_ht(i);
// 			this->graph->topography = ntopo;
// 			this->graph->compute_graph("carve");
// 			this->graph->topography = this->base_topo;

// 			// first calculating Qin
// 			std::vector<float> Qin = this->graph->accumulate_downstream_var(this->graph->cellarea, this->precipitations);
// 			float maxs = 0;
// 			for(int i = this->graph->nnodes - 1 ; i >= 0; --i)
// 			{
// 				int tnode = this->graph->stack[i];
// 				int rec = this->graph->receivers[tnode];

// 				if(this->graph->can_flow_even_go_there(tnode) == false || this->graph->can_flow_out_there(tnode) )
// 					continue;

// 				if(this->Qwinth > 0 && this->Qwinth > Qin[tnode] )
// 					continue;

// 				float tS = std::fmax(this->get_Sw(tnode), 0);
// 				float Qout = this->graph->distance2receivers[tnode] * 1/this->manning * std::pow(this->hw[tnode],5/3) * std::pow(std::abs(tS),0.5);
// 				float dQ = Qin[tnode] - Qout;
// 				float dhw = dQ * this->dt / this->graph->cellarea;
// 				// if(dhw > 10)
// 				// 	std::cout << Qin[tnode] << "|" << Qout << std::endl;
// 				if(Qin[tnode] > maxs)
// 					maxs = Qin[tnode];

// 				if(allow_neg)
// 					this->hw[tnode] = fmax(dhw + this->hw[tnode],0);
// 				else
// 				{
// 					if(dhw > 0)
// 						this->hw[tnode] = fmax(dhw + this->hw[tnode],0);

// 				}
// 			}

// 		}

// 		void run_analitically(float increment = 1e-3)
// 		{
// 			std::vector<float> ntopo(this->base_topo);
// 			for(int i=0; i<this->graph->nnodes; ++i)
// 				ntopo[i] = this->get_ht(i);
// 			this->graph->topography = ntopo;
// 			this->graph->compute_graph("carve");
// 			this->graph->topography = this->base_topo;
// 			std::vector<float> Qin = this->graph->accumulate_downstream_var(this->graph->cellarea, this->precipitations);

// 			for(int i = 0 ; i < this->graph->nnodes; ++i)
// 			{
// 				int tnode = this->graph->stack[i];
// 				int rec = this->graph->receivers[tnode];

// 				if(this->graph->can_flow_even_go_there(tnode) == false || this->graph->can_flow_out_there(tnode) )
// 					continue;

// 				float tS = std::fmax(this->get_Sw(tnode), 0);
// 				float Qout = this->graph->distance2receivers[tnode] * 1/this->manning * std::pow(this->hw[tnode],5/3) * std::pow(std::abs(tS),0.5);
// 				while(Qout < Qin[tnode])
// 				{
// 					this->hw[tnode] += increment;
// 					float tS = std::fmax(this->get_Sw(tnode), 0);
// 					Qout = this->graph->distance2receivers[tnode] * 1/this->manning * std::pow(this->hw[tnode],5/3) * std::pow(std::abs(tS),0.5);
// 				}
// 			}
// 		}

// 		std::vector<float> get_full_Sw()
// 		{
// 			std::vector<float> fSw(this->graph->nnodes,0.);

// 			for(int i=0; i<this->graph->nnodes; ++i)
// 				fSw[i] = this->get_Sw(i);

// 			return fSw;
// 		}


// 		void flood_in()
// 		{
// 			std::vector<float> ntopo(this->base_topo);
// 			for(int i=0; i<this->graph->nnodes; ++i)
// 				ntopo[i] = this->get_ht(i);
// 			this->graph->topography = ntopo;
// 			this->graph->compute_graph("carve");
// 			this->graph->topography = this->base_topo;

// 			// first calculating Qin
// 			std::vector<float> Qin = this->graph->accumulate_downstream_var(this->graph->cellarea, this->precipitations);
// 			float maxs = 0;
// 			for(int i = this->graph->nnodes - 1 ; i >= 0; --i)
// 			{
// 				int tnode = this->graph->stack[i];
// 				int rec = this->graph->receivers[tnode];

// 				if(this->graph->can_flow_even_go_there(tnode) == false || this->graph->can_flow_out_there(tnode) )
// 					continue;

// 				this->hw[tnode] = std::fmax((this->hydrcof * std::pow(Qin[tnode],this->hydrex))/this->graph->cellarea , this->hw[tnode] );
// 			}
// 		}

// 		void set_Qwinth(float Qwinth){this->Qwinth = Qwinth;};

// 		void flood_in_exp()
// 		{
// 			std::vector<float> ntopo(this->base_topo);
// 			for(int i=0; i<this->graph->nnodes; ++i)
// 				ntopo[i] = this->get_ht(i);
// 			this->graph->topography = ntopo;
// 			this->graph->compute_graph("carve");
// 			this->graph->topography = this->base_topo;

// 			// first calculating Qin
// 			std::vector<float> Qin = this->graph->accumulate_downstream_var(this->graph->cellarea, this->precipitations);
// 			float maxs = 0;
// 			for(int i = this->graph->nnodes - 1 ; i >= 0; --i)
// 			{
// 				int tnode = this->graph->stack[i];
// 				int rec = this->graph->receivers[tnode];

// 				if(this->graph->can_flow_even_go_there(tnode) == false || this->graph->can_flow_out_there(tnode) )
// 					continue;

// 				this->hw[tnode] = std::fmax((this->hydrcof * std::pow(Qin[tnode],this->hydrex))/this->graph->cellarea , this->hw[tnode] );
// 			}
// 		}

// 		void transient()
// 		{
// 			std::vector<float> ntopo(this->base_topo);
// 			for(int i=0; i<this->graph->nnodes; ++i)
// 				ntopo[i] = this->get_ht(i);
// 			this->graph->topography = ntopo;
// 			this->graph->compute_graph("carve");
// 			this->graph->topography = this->base_topo;

// 			// first calculating Qin
// 			std::vector<float> Qin(this->precipitations);
// 			float maxs = 0;
// 			for(int i = this->graph->nnodes - 1 ; i >= 0; --i)
// 			{
// 				int tnode = this->graph->stack[i];
// 				int rec = this->graph->receivers[tnode];

// 				if(this->graph->can_flow_even_go_there(tnode) == false || this->graph->can_flow_out_there(tnode) )
// 					continue;

// 				float tS = std::fmax(this->get_Sw(tnode), 0);
// 				float Qout = this->graph->distance2receivers[tnode] * 1/this->manning * std::pow(this->hw[tnode],5/3) * std::pow(std::abs(tS),0.5);
// 				float dQ = Qin[tnode] - Qout;
// 				float dhw = dQ * this->dt / this->graph->cellarea;
// 				// if(dhw > 10)
// 				// 	std::cout << Qin[tnode] << "|" << Qout << std::endl;
// 				if(Qin[tnode] > maxs)
// 					maxs = Qin[tnode];

// 				Qin[rec] += Qout;

// 				if(allow_neg)
// 					this->hw[tnode] = fmax(dhw + this->hw[tnode],0);
// 				else
// 				{
// 					if(dhw > 0)
// 						this->hw[tnode] = fmax(dhw + this->hw[tnode],0);

// 				}
// 			}
// 		}

// 		float get_ht(int tnode){return this->base_topo[tnode] + this->hw[tnode];}
// 		float get_Sw(int tnode){return (this->get_ht(tnode) - this->get_ht(this->graph->receivers[tnode]) ) / this->graph->distance2receivers[tnode];}
// 		float get_S(int tnode){return (this->base_topo[tnode] - this->base_topo[this->graph->receivers[tnode]] ) / this->graph->distance2receivers[tnode];}
// 		std::vector<float> get_full_ht()
// 		{
// 			std::vector<float> fullht(this->graph->nnodes_t,this->graph->NDV);
// 			for(int i =0; i<this->graph->nnodes; ++i)
// 			{
// 				fullht[i] = this->get_ht(i);
// 			}
// 			return fullht;
// 		}

// 		void set_hydremp(float hydrcof, float hydrex)
// 		{
// 			this->hydrcof = hydrcof;
// 			this->hydrex = hydrex;
// 		}

// 		void ingest_hw_other_dim(std::vector<float>& ohw, Graph& ograph)
// 		{
// 			this->init();
// 			this->graph->calculate_area();
// 			for(int i = 0; i<this->graph->nnodes; ++i)
// 			{
// 				if(this->graph->can_flow_even_go_there(i) == false || this->graph->area[i] < this->graph->dx * this->graph->dy * 500 )
// 					continue;

// 				int row,col;
// 				this->graph->rowcol_from_node_id(i,row,col);
// 				float percrow = float(row)/this->graph->ny;
// 				float perccol = float(col)/this->graph->nx;
// 				int j = ograph.nodeid_from_row_col(floor(percrow * ograph.ny),floor(perccol * ograph.nx));
// 				this->hw[i] = ohw[j];
// 			}
// 			On_gaussian_blur(1,this->hw, this->graph->nx, this->graph->ny);
// 		}

// 		void diffuse_hw(float r){On_gaussian_blur(r,this->hw, this->graph->nx, this->graph->ny);}
// 		void multiply_hw_by(float factor)
// 		{
// 			for (auto& v: this->hw)
// 			{
// 				v*=factor;
// 			}
// 		}

// 		void reinit_boundaries()
// 		{
// 			for(int i = 0; i<this->graph->nnodes; ++i)
// 			{
// 				if(this->graph->can_flow_out_there(i) || this->graph->can_flow_even_go_there(i) == false)
// 					this->hw[i] = 0;
// 			}
// 		}

// 		void fill_pits()
// 		{
// 			std::vector<float> ntopo(this->base_topo);
// 			for(int i=0; i<this->graph->nnodes; ++i)
// 				ntopo[i] = this->get_ht(i);
// 			this->graph->topography = ntopo;

// 			this->graph->compute_graph("fill");
// 			for(int i=0; i<this->graph->nnodes; ++i)
// 			{
				
// 				if(this->graph->can_flow_out_there(i) || this->graph->can_flow_even_go_there(i) == false)
// 					continue;

// 				int rec = this->graph->receivers[i];
// 				if(this->get_ht(i) <= this->get_ht(rec))
// 				{

// 					this->hw[i] += 1e-4 + this->get_ht(rec) - this->get_ht(i);
// 				}
// 			}
// 			this->graph->topography = this->base_topo;
// 		}

// #ifdef __EMSCRIPTEN__
// #else
// 		py::array get_hw_np(){return py::array(this->hw.size(), this->hw.data());}
// #endif

// };























#endif