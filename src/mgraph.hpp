//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef mgraph_HPP
#define mgraph_HPP

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
#include <omp.h>

// local includes 
// -> General routines and data structures
#include "chonkutils.hpp"
// -> Depression solvers
#include "cordonnier_versatile_2019.hpp"
// -> lightweigth numpy export
#include "npy.hpp"

#include "neighbourer.hpp"
#include "numvec.hpp"

// If compiling for web (using emscripten)
#ifdef __EMSCRIPTEN__
	#include <emscripten.h>
// -> Else I assume you are compiling for python (for c++ only use I'll make a separate file, just need to remvove these lines)
#else
	#include <pybind11/pybind11.h>
	#include <pybind11/stl.h>
	#include <pybind11/numpy.h>
#endif

#ifdef __EMSCRIPTEN__
#else
	namespace py = pybind11;
#endif


// the main calss managing the graph
template<class T>
class MGraph
{
public:

// ------------------------------------------------

//	                             	              __
//                                             / _)
//                                    _.----._/ /
//   CLASS MEMBERS                   /         /
//                                __/ (  | (  |
//                               /__.-'|_|--|_|

// ------------------------------------------------


	// General informations about the graph
	// #-> number of nodes (integer and unsized int for loops or list innit).
	int nnodes = 0;
	size_t nnodes_t = 0;

	// Steepest descent receiver
	std::vector<int> Sreceivers;
	// Steepest descent donors
	std::vector<std::vector<int> > Sdonors;
	// Steepest descent dx
	std::vector<T> Sdistance2receivers;


	// Multiple flow reveivers
	std::vector<std::vector<int> > receivers;
	// Multiple flow donors
	std::vector<std::vector<int> > donors;
	// Multiple flow dx
	std::vector<std::vector<T> > distance2receivers;


		// #->stack: topological order from downstream to upstream direction for general (sorted)
	std::vector<size_t> stack,Sstack;




	// ------------------------------------------------

	//	                             	              __
	//                                             / _)
	//                                    _.----._/ /
	//   CONSTRUCTOR METHODS             /         /
	//                                __/ (  | (  |
	//                               /__.-'|_|--|_|

	// ------------------------------------------------

	// Constructors
	// #->Empty constructor
	MGraph(){;};
	MGraph(int nnodes){this->nnodes = nnodes; this->nnodes_t = nnodes;};

	void yolomp(int n_proc)
	{
		int N = 10000000;

		std::vector<double> yabul(N);

		auto t1 = high_resolution_clock::now();
		



		#pragma omp parallel for default(shared) num_threads(2)
		for(int i = 0; i<N; ++i)
		{
			for(int j =0; j<100; ++j)
				yabul[i]+= j * i;
		}

		auto t2 = high_resolution_clock::now();
		duration<double, std::milli> ms_double = t2 - t1;
		double time_para = ms_double.count();


		t1 = high_resolution_clock::now();
		for(int i = 0; i<N; ++i)
		{
			for(int j =0; j<100; ++j)
				yabul[i]+= j * i;
		}

		t2 = high_resolution_clock::now();
		ms_double = t2 - t1;
		double time_no_para = ms_double.count();
		std:: cout << "OpenMP was " << time_para << " and vanilla was " << time_no_para << std::endl;
	}

	void yonolomp(int n_proc)
	{
		std::vector<double> yabul(this->nnodes * 10);
		#pragma omp parallel for num_threads(n_proc)
		for(int i = 0; i<this->nnodes; ++i)
			yabul[i]+= -1e-4 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(2e-4)));
	}


	// ------------------------------------------------

	//	                             	              __
	//                                             / _)
	//                                    _.----._/ /
	//   GRAPH CREATION METHODS          /         /
	//                                __/ (  | (  |
	//                               /__.-'|_|--|_|

	// All the methods related to accessing and calculating neighbours
	// ------------------------------------------------

	template<class Neighbourer_t, class topo_t, class out_t>
	out_t compute_graph(std::string depression_solver, Neighbourer_t& neighbourer, topo_t& ttopography)
	{
		// std::cout << "ASDFdsfa->1" << std::endl;
		auto topography = format_input(ttopography);

		// auto t1 = high_resolution_clock::now();
		this->compute_graph_both_v2(neighbourer,topography);
		// auto t2 = high_resolution_clock::now();
		// duration<double, std::milli> ms_double = t2 - t1;
		// double time_compute_graph_both_v2 = ms_double.count();
		// std::cout << "ASDFdsfa->2" << std::endl;
		
		// t1 = high_resolution_clock::now();
		this->compute_TO_SF_stack_version();
		// t2 = high_resolution_clock::now();
		// ms_double = t2 - t1;
		// double time_compute_TO_SF_stack_version1 = ms_double.count();

		// std::cout << "ASDFdsfa->3" << std::endl;
		std::vector<double> faketopo = to_vec(topography);

		if(depression_solver != "none")
		{
			// t1 = high_resolution_clock::now();

			this->solve_depressions(depression_solver, neighbourer, topography);
			// t2 = high_resolution_clock::now();
			// ms_double = t2 - t1;
			// double time_solve_depressions = ms_double.count();

			// t1 = high_resolution_clock::now();
			this->compute_TO_SF_stack_version();
			// t2 = high_resolution_clock::now();
			// ms_double = t2 - t1;
			// double time_compute_TO_SF_stack_version2 = ms_double.count();

			// t1 = high_resolution_clock::now();
			if(depression_solver == "fill")
				this->fill_topo(1e-3,neighbourer,faketopo);
			else
				this->carve_topo(1e-3,neighbourer,faketopo);
			// this->carve_topo(1e-3,neighbourer,faketopo);
			// t2 = high_resolution_clock::now();
			// ms_double = t2 - t1;
			// double time_carve_topo = ms_double.count();



			// t1 = high_resolution_clock::now();
			// this->receivers.clear();
			// this->donors.clear();
			// this->distance2receivers.clear();
			// this->receivers = std::vector<std::vector<int> >(this->nnodes);
			// this->donors = std::vector<std::vector<int> >(this->nnodes);
			// this->distance2receivers = std::vector<std::vector<T> >(this->nnodes);
			// neighbourer.rebuild_mgraph_after_solving_depression(Sdonors, Sreceivers, receivers,donors, Sdistance2receivers, distance2receivers, faketopo);
			this->update_receivers_v2(faketopo);
			// t2 = high_resolution_clock::now();
			// ms_double = t2 - t1;
			// double time_rebuild_mgraph = ms_double.count();


			// t1 = high_resolution_clock::now();
			// this->compute_MF_topological_order();
			this->compute_MF_topological_order_insort(faketopo);
			// t2 = high_resolution_clock::now();
			// ms_double = t2 - t1;
			// double time_compute_MF_topological_order = ms_double.count();

		// std::cout << "ASDFdsfa->4" << std::endl;
		// std::cout << "ASDFdsfa->5" << std::endl;

			// std::cout << "time_compute_graph_both_v2::" << time_compute_graph_both_v2 << std::endl;
			// std::cout << "time_compute_TO_SF_stack_version1::" << time_compute_TO_SF_stack_version1 << std::endl;
			// std::cout << "time_solve_depressions::" << time_solve_depressions << std::endl;
			// std::cout << "time_compute_TO_SF_stack_version2::" << time_compute_TO_SF_stack_version2 << std::endl;
			// std::cout << "time_carve_topo::" << time_carve_topo << std::endl;
			// std::cout << "time_rebuild_mgraph::" << time_rebuild_mgraph << std::endl;
			// std::cout << "time_compute_MF_topological_order::" << time_compute_MF_topological_order << std::endl;
		}


		return format_output(faketopo);
	}



	template<class Neighbourer_t, class topo_t, class out_t>
	out_t compute_graph_OMP(std::string depression_solver, Neighbourer_t& neighbourer, topo_t& ttopography, int n_proc)
	{
		// std::cout << "ASDFdsfa->1" << std::endl;
		auto topography = format_input(ttopography);

		// auto t1 = high_resolution_clock::now();
		this->compute_graph_both_v2_OMP(neighbourer,topography, n_proc);
		// auto t2 = high_resolution_clock::now();
		// duration<double, std::milli> ms_double = t2 - t1;
		// double time_compute_graph_both_v2 = ms_double.count();
		// std::cout << "ASDFdsfa->2" << std::endl;
		
		// t1 = high_resolution_clock::now();
		this->compute_TO_SF_stack_version();
		// t2 = high_resolution_clock::now();
		// ms_double = t2 - t1;
		// double time_compute_TO_SF_stack_version1 = ms_double.count();

		// std::cout << "ASDFdsfa->3" << std::endl;
		std::vector<double> faketopo = to_vec(topography);
		if(depression_solver != "none")
		{
			// t1 = high_resolution_clock::now();

			this->solve_depressions(depression_solver, neighbourer, topography);
			// t2 = high_resolution_clock::now();
			// ms_double = t2 - t1;
			// double time_solve_depressions = ms_double.count();

			// t1 = high_resolution_clock::now();
			this->compute_TO_SF_stack_version();
			// t2 = high_resolution_clock::now();
			// ms_double = t2 - t1;
			// double time_compute_TO_SF_stack_version2 = ms_double.count();



			// t1 = high_resolution_clock::now();
			if(depression_solver == "fill")
				this->fill_topo(1e-3,neighbourer,faketopo);
			else
				this->carve_topo(1e-3,neighbourer,faketopo);
			// this->carve_topo(1e-3,neighbourer,faketopo);
			// t2 = high_resolution_clock::now();
			// ms_double = t2 - t1;
			// double time_carve_topo = ms_double.count();



			// t1 = high_resolution_clock::now();
			// this->receivers.clear();
			// this->donors.clear();
			// this->distance2receivers.clear();
			// this->receivers = std::vector<std::vector<int> >(this->nnodes);
			// this->donors = std::vector<std::vector<int> >(this->nnodes);
			// this->distance2receivers = std::vector<std::vector<T> >(this->nnodes);
			// neighbourer.rebuild_mgraph_after_solving_depression(Sdonors, Sreceivers, receivers,donors, Sdistance2receivers, distance2receivers, faketopo);
			this->update_receivers_v2(faketopo);
			// t2 = high_resolution_clock::now();
			// ms_double = t2 - t1;
			// double time_rebuild_mgraph = ms_double.count();


			// t1 = high_resolution_clock::now();
			// this->compute_MF_topological_order();
			this->compute_MF_topological_order_insort(faketopo);
			// t2 = high_resolution_clock::now();
			// ms_double = t2 - t1;
			// double time_compute_MF_topological_order = ms_double.count();

		// std::cout << "ASDFdsfa->4" << std::endl;
		// std::cout << "ASDFdsfa->5" << std::endl;

			// std::cout << "time_compute_graph_both_v2::" << time_compute_graph_both_v2 << std::endl;
			// std::cout << "time_compute_TO_SF_stack_version1::" << time_compute_TO_SF_stack_version1 << std::endl;
			// std::cout << "time_solve_depressions::" << time_solve_depressions << std::endl;
			// std::cout << "time_compute_TO_SF_stack_version2::" << time_compute_TO_SF_stack_version2 << std::endl;
			// std::cout << "time_carve_topo::" << time_carve_topo << std::endl;
			// std::cout << "time_rebuild_mgraph::" << time_rebuild_mgraph << std::endl;
			// std::cout << "time_compute_MF_topological_order::" << time_compute_MF_topological_order << std::endl;
		}
		return format_output(faketopo);
	}


	template<class Neighbourer_t, class topo_t, class out_t>
	out_t update_graph(std::string depression_solver, Neighbourer_t& neighbourer, topo_t& ttopography)
	{
		// std::cout << "ASDFdsfa->1" << std::endl;
		auto topography = format_input(ttopography);
		// auto t1 = high_resolution_clock::now();
		// this->compute_graph_both_v2(neighbourer,topography);
		this->update_receivers_v2(topography);
		this->update_Srecs_from_recs(topography, neighbourer);
		this->recompute_SF_donors_from_receivers();

		// auto t2 = high_resolution_clock::now();
		// duration<double, std::milli> ms_double = t2 - t1;
		// double time_compute_graph_both_v2 = ms_double.count();
		// std::cout << "ASDFdsfa->2" << std::endl;
		
		// t1 = high_resolution_clock::now();
		this->compute_TO_SF_stack_version();
		// t2 = high_resolution_clock::now();
		// ms_double = t2 - t1;
		// double time_compute_TO_SF_stack_version1 = ms_double.count();

		// std::cout << "ASDFdsfa->3" << std::endl;
		std::vector<double> faketopo = to_vec(topography);
		if(depression_solver != "none")
		{
			// t1 = high_resolution_clock::now();

			this->solve_depressions(depression_solver, neighbourer, topography);
			// t2 = high_resolution_clock::now();
			// ms_double = t2 - t1;
			// double time_solve_depressions = ms_double.count();

			// t1 = high_resolution_clock::now();
			this->compute_TO_SF_stack_version();
			// t2 = high_resolution_clock::now();
			// ms_double = t2 - t1;
			// double time_compute_TO_SF_stack_version2 = ms_double.count();



			// t1 = high_resolution_clock::now();
			if(depression_solver == "fill")
				this->fill_topo(1e-3,neighbourer,faketopo);
			else
				this->carve_topo(1e-3,neighbourer,faketopo);
			// t2 = high_resolution_clock::now();
			// ms_double = t2 - t1;
			// double time_carve_topo = ms_double.count();



			// t1 = high_resolution_clock::now();
			// this->receivers.clear();
			// this->donors.clear();
			// this->distance2receivers.clear();
			// this->receivers = std::vector<std::vector<int> >(this->nnodes);
			// this->donors = std::vector<std::vector<int> >(this->nnodes);
			// this->distance2receivers = std::vector<std::vector<T> >(this->nnodes);
			// neighbourer.rebuild_mgraph_after_solving_depression(Sdonors, Sreceivers, receivers,donors, Sdistance2receivers, distance2receivers, faketopo);
			this->update_receivers_v2(faketopo);
			// t2 = high_resolution_clock::now();
			// ms_double = t2 - t1;
			// double time_rebuild_mgraph = ms_double.count();

			int Ndef = 0;
			for(int i =0 ;i< this->nnodes; ++i)
			{
				int row, col;
				neighbourer.rowcol_from_node_id(i,row,col);
				if(row <2 || col < 2 || row >= neighbourer.ny - 2 || col >= neighbourer.nx - 2)
					continue;

				if(this->donors[i].size() + this->receivers[i].size() != 8)
					Ndef++;

			}
			if(Ndef > 0)
				std::cout << "SOMETHING WRONG HERE" << std::endl;


			// t1 = high_resolution_clock::now();
			// this->compute_MF_topological_order();
			this->recompute_MF_topological_order_insort(faketopo);
			// t2 = high_resolution_clock::now();
			// ms_double = t2 - t1;
			// double time_compute_MF_topological_order = ms_double.count();

		// std::cout << "ASDFdsfa->4" << std::endl;
		// std::cout << "ASDFdsfa->5" << std::endl;

			// std::cout << "time_compute_graph_both_v2::" << time_compute_graph_both_v2 << std::endl;
			// std::cout << "time_compute_TO_SF_stack_version1::" << time_compute_TO_SF_stack_version1 << std::endl;
			// std::cout << "time_solve_depressions::" << time_solve_depressions << std::endl;
			// std::cout << "time_compute_TO_SF_stack_version2::" << time_compute_TO_SF_stack_version2 << std::endl;
			// std::cout << "time_carve_topo::" << time_carve_topo << std::endl;
			// std::cout << "time_rebuild_mgraph::" << time_rebuild_mgraph << std::endl;
			// std::cout << "time_compute_MF_topological_order::" << time_compute_MF_topological_order << std::endl;
		}
		return format_output(faketopo);
	}

	template<class Neighbourer_t, class topo_t>
	void compute_graph_both_v2(Neighbourer_t& neighbourer, topo_t& ttopography)
	{
		auto topography = format_input(ttopography);
		// Initialising the graph dimesions for the receivers
		// All of thenm have the graph dimension
		this->Sreceivers.clear();
		this->Sdonors.clear();
		this->Sdistance2receivers.clear();
		this->Sdonors = std::vector<std::vector<int> >(this->nnodes);
		this->Sreceivers = std::vector<int>(this->nnodes_t, -1);
		this->Sdistance2receivers = std::vector<T>(this->nnodes_t, -1);
		std::vector<T> SS(this->nnodes_t, 0.);

		this->receivers.clear();
		this->donors.clear();
		this->distance2receivers.clear();
		this->receivers = std::vector<std::vector<int> >(this->nnodes);
		this->donors = std::vector<std::vector<int> >(this->nnodes);
		this->distance2receivers = std::vector<std::vector<T> >(this->nnodes);


		for (int i=0; i < this->nnodes; ++i)
			this->Sreceivers[i] = i;

		neighbourer.build_mgraph(SS, Sdonors, Sreceivers, receivers,donors, Sdistance2receivers, distance2receivers, topography);
		
		this->recompute_SF_donors_from_receivers();
	}


	template<class Neighbourer_t, class topo_t>
	void compute_graph_both_v2_OMP(Neighbourer_t& neighbourer, topo_t& ttopography, int n_proc)
	{
		// Initialising the graph dimesions for the receivers
		// All of thenm have the graph dimension
		auto topography = format_input(ttopography);
		this->Sreceivers.clear();
		this->Sdonors.clear();
		this->Sdistance2receivers.clear();
		this->Sdonors = std::vector<std::vector<int> >(this->nnodes);
		this->Sreceivers = std::vector<int>(this->nnodes_t, -1);
		this->Sdistance2receivers = std::vector<T>(this->nnodes_t, -1);
		std::vector<T> SS(this->nnodes_t, 0.);

		this->receivers.clear();
		this->donors.clear();
		this->distance2receivers.clear();
		this->receivers = std::vector<std::vector<int> >(this->nnodes);
		this->donors = std::vector<std::vector<int> >(this->nnodes);
		this->distance2receivers = std::vector<std::vector<T> >(this->nnodes);


		for (int i=0; i < this->nnodes; ++i)
			this->Sreceivers[i] = i;

		neighbourer.build_mgraph_OMP(SS, Sdonors, Sreceivers, receivers,donors, Sdistance2receivers, distance2receivers, topography, n_proc);
		
		this->recompute_SF_donors_from_receivers();
	}
	



	void recompute_SF_donors_from_receivers()
	{
		// Initialising the graph dimesions for the donors
		// All of thenm have the graph dimension
		this->Sdonors = std::vector<std::vector<int> >(this->nnodes);

		for(int i=0; i < this->nnodes; ++i)
		{
			// SF so rid == i cause there is only 1 rec
			int trec = this->Sreceivers[i];
			if(trec == i)
				continue;
			this->Sdonors[trec].emplace_back(i);
		}
	}

	template<class topo_t>
	void update_receivers(topo_t& ttopography)
	{
		auto topography = format_input(ttopography);

		std::vector<std::vector<int> > nreceivers(this->receivers.size());
		std::vector<std::vector<int> > ndonors(this->receivers.size());
		std::vector<std::vector<T> > ndistance2receivers(this->distance2receivers);
		for(auto&v:nreceivers) {v.reserve(8);}
		for(auto&v:ndonors) {v.reserve(8);}
		for(auto&v:ndistance2receivers) {v.reserve(8);}

		for(int i = 0; i<this->nnodes; ++i)
		{
			for(int j=0; j<int(this->receivers[i].size()); ++j)
			{

				if(topography[i] >= topography[this->receivers[i][j]])
				{
					nreceivers[i].emplace_back(this->receivers[i][j]);
					ndonors[this->receivers[i][j]].emplace_back(i);
					ndistance2receivers[i].emplace_back(this->distance2receivers[i][j]);
				}
				else
				{
					nreceivers[this->receivers[i][j]].emplace_back(i);
					ndonors[i].emplace_back(this->receivers[i][j]);
					ndistance2receivers[this->receivers[i][j]].emplace_back(this->distance2receivers[i][j]);
				}
			}
		}

		this->receivers = std::move(nreceivers);
		this->donors = std::move(ndonors);
		this->distance2receivers = std::move(ndistance2receivers);

	}


	template<class topo_t>
	void update_receivers_v2(topo_t& ttopography)
	{
		auto topography = format_input(ttopography);

		std::vector<bool> need_redodon(this->nnodes,false);
		for(int i = 0; i<this->nnodes; ++i)
		{
			std::vector<int> tokeep;
			for(int j=0; j<int(this->receivers[i].size()); ++j)
			{

				if(topography[i] >= topography[this->receivers[i][j]])
				{

					tokeep.emplace_back(j);
				}
				else
				{
					need_redodon[this->receivers[i][j]] = true;
					need_redodon[i] = true;
					this->receivers[this->receivers[i][j]].emplace_back(i);
					this->distance2receivers[this->receivers[i][j]].emplace_back(this->distance2receivers[i][j]);
				}

			}

			if(tokeep.size() < this->receivers[i].size())
			{
				std::vector<int> trecs;
				std::vector<T> tdist;

			
				for(auto j:tokeep)
				{
					trecs.emplace_back(this->receivers[i][j]);
					tdist.emplace_back(this->distance2receivers[i][j]);
				}

				this->receivers[i] = trecs;
				this->distance2receivers[i] = tdist;
			}
		}

		for(int i=0 ; i<this->nnodes; ++i)
		{
			if(need_redodon[i])
			{
				this->donors[i].clear();
			}
		}

		for(int i=0 ; i<this->nnodes; ++i)
		{
			for(auto v:this->receivers[i])
			{
				if(need_redodon[v])
				{
					this->donors[v].emplace_back(i);
				}
			}
		}


	}

	template<class topo_t, class Neighbourer_t>
	void update_Srecs_from_recs(topo_t& ttopography, Neighbourer_t& neighbourer )
	{
		auto topography = format_input(ttopography);


		for(int i=0;i<this->nnodes; ++i)
		{
			this->Sreceivers[i] = i;

			if(neighbourer.is_active(i) == false)
				continue;

			int srec = -1;
			T dist = -1;
			T smax = -1e3;
			for(int j=0; j<this->receivers[i].size(); ++j)
			{
				T tsle = (topography[i] - topography[this->receivers[i][j]])/this->distance2receivers[i][j];
				if(tsle > smax)
				{
					smax = tsle;
					srec = this->receivers[i][j];
					dist =this->distance2receivers[i][j];
				}
			}
			if(srec != -1)
			{
				this->Sreceivers[i] = srec;
				this->Sdistance2receivers[i] = dist ;
			}
		
		}

	}

	template<class Neighbourer_t, class topo_t>
	void solve_depressions(std::string& depression_solver , Neighbourer_t& neighbourer, topo_t& ttopography)
	{
		auto topography = format_input(ttopography);
		LMRerouter depsolver(neighbourer, topography, this->Sreceivers,  this->Sdistance2receivers, this->Sstack);


		if(depsolver.npits > 0)
		{

			depsolver.update_receivers(depression_solver,neighbourer, topography, this->Sreceivers, this->Sdistance2receivers, this->Sstack);
			
			this->recompute_SF_donors_from_receivers();
		}
	}


	
	/// this function enforces minimal slope 
	template<class Neighbourer_t, class topo_t>
	void carve_topo(T slope, Neighbourer_t& neighbourer, topo_t& ttopography)
	{
		auto topography = format_input(ttopography);

		for(int i=this->nnodes-1; i >= 0; --i)
		{
			int node  = this->Sstack[i];
			if(neighbourer.can_flow_out_there(node) || neighbourer.can_flow_even_go_there(node) == false)
				continue;
			int rec = this->Sreceivers[node];
			T dz = topography[node] - topography[rec];

			if(dz <= 0)
			{
				T d2rec = this->Sdistance2receivers[node];
				topography[rec] = topography[node] - slope * d2rec;
			}

		}
	}

	/// this function enforces minimal slope 
	template<class Neighbourer_t, class topo_t>
	void fill_topo(T slope, Neighbourer_t& neighbourer, topo_t& ttopography)
	{
		auto topography = format_input(ttopography);
		for(int i=0; i < this->nnodes; ++i)
		{
			int node  = this->Sstack[i];
			if(neighbourer.can_flow_out_there(node) || neighbourer.can_flow_even_go_there(node) == false)
				continue;

			int rec = this->Sreceivers[node];
			T dz = topography[node] - topography[rec];

			if(dz <= 0)
			{
				T d2rec = this->Sdistance2receivers[node];
				topography[node] = topography[rec] + slope * d2rec;
			}
		}
	}

	// Return an array of receivers of a given node
	// Return array in a DirectedNeighbour structure: node ID receiver ID and distance
	Neighbour<int,T> get_steepest_receivers(int i)
	{
		return Neighbour<int,T>(this->Sreceivers[i],this->Sdistance2receivers[i]);
	}


	// Return an array of donors of a given node
	// Return array in a DirectedNeighbour structure: node ID donor ID and distance
	std::vector<Neighbour<int,T> > get_steepest_donors(int i)
	{
		// Preprocessing the output to the right size (number of donors for the given node)
		std::vector<Neighbour<int,T> > output;
		output.reserve(this->Sdonors[i].size());

		//Going through all the donor ID 
		for(int j = 0; j < this->Sdonors[i].size(); ++j)
		{
			// And gathering the node ID,donor ID and distance to i of the donor node
			output.emplace_back(this->Sdonors[i][j], this->Sdistance2receivers[this->Sdonors[i][j]]);
		}

		// AANd we are done
		return output;
	}


	// ------------------------------------------------

	//	                             	            __
	//                                             / _)
	//                                    _.----._/ /
	//   TOPOLOGICAL ORDER METHODS       /         /
	//                                __/ (  | (  |
	//                               /__.-'|_|--|_|

	// All the methods related to accessing and calculating neighbours
	// ------------------------------------------------
	void compute_TO_SF_stack_version()
	{
		// Initialising the stack
		this->Sstack.clear();
		// reserving the amount of stuff
		this->Sstack.reserve(this->nnodes_t);

		// The stack container helper
		std::stack<size_t, std::vector<size_t> > stackhelper;
		// std::vector<bool> isdone(this->nnodes_t,false);
		// going through all the nodes
		for(int i=0; i<this->nnodes; ++i)
		{
			// if they are base level I include them in the stack
			if(this->Sreceivers[i] == i)
			{
				stackhelper.emplace(i);
			}

			// While I still have stuff in the stack helper
			while(stackhelper.empty() == false)
			{
				// I get the next node and pop it from the stack helper
				int nextnode = stackhelper.top();stackhelper.pop();
				// std::cout << stackhelper.size() << "->" << nextnode << std::endl;

				// // I emplace it in the stack
				// if(isdone[nextnode] == true)
				//	throw std::runtime_error("node-duplicate");


				// isdone[nextnode] = true;
				this->Sstack.emplace_back(nextnode);

				// as well as all its donors which will be processed next
				for( int j = 0; j < this->Sdonors[nextnode].size(); ++j)
				{
					stackhelper.emplace(this->Sdonors[nextnode][j]);
				}

			}

		}

		if(this->Sstack.size() != this->nnodes_t)
		{
			std::cout << "SStack error, "<< this->Sstack.size() << "|" << this->nnodes << " checking for nans..." << std::endl;
			throw std::runtime_error("SStack error: should be " + std::to_string(this->nnodes_t) + " but is " + std::to_string(this->Sstack.size()));
		}

	}

	


	void compute_MF_topological_order()
	{

		this->stack.clear();
		this->stack.reserve(this->nnodes);

	  std::vector<T> vis(this->nnodes,0), parse(this->nnodes,-1), ndons(this->nnodes,0);
	  
	  int nparse = -1;
	  int nstack = -1;
	  for(int i=0; i<this->nnodes;++i)
	  {
	  	ndons[i] = this->donors[i].size();
	  }

	  // we go through the nodes
	  for(int i=0; i<this->nnodes;++i)
	  {
	    // when we find a "summit" (ie a node that has no donors)
	    // we parse it (put it in a stack called parse)
	    // auto dons = this->get_donors_ID(i);
	    if (ndons[i] == 0)
	    {
	      nparse =  nparse+1;
	      parse[nparse] = i;
	    }
	    // we go through the parsing stack
	    while (nparse > -1)
	    {
	      int ijn = parse[nparse];
	      nparse = nparse-1;
	      nstack = nstack+1;

	      this->stack.emplace_back(ijn);
	      // std::cout << ijn << "|";

	      auto recs = this->receivers[ijn];

	      for(int ijk=0; ijk < int(recs.size()); ++ijk)
	      {
	        int ijr = recs[ijk];
	        vis[ijr] = vis[ijr]+1;

	        // if the counter is equal to the number of donors for that node we add it to the parsing stack
	        if (vis[ijr] == ndons[ijr])
	        {
	          nparse=nparse+1;
	          parse[nparse]=ijr;
	        }
	        
	      }
	    } 
	  }

	  std::reverse(this->stack.begin(), this->stack.end());
	  
	  // std::cout << "TO full finished with " << this->topological_order.size() << " versus " << this->isize << std::endl;;
	}

	template<class topo_t>
	void compute_MF_topological_order_insort(topo_t& ttopography)
	{
		auto topography = format_input(ttopography);

		auto yolo = sort_indexes(topography);
		this->stack = std::move(yolo);

	}

	template<class topo_t>
	void recompute_MF_topological_order_insort(topo_t& ttopography)
	{
		auto topography = format_input(ttopography);

		for(std::size_t j = 1; j < this->stack.size(); ++j)
    {
    	int oval = this->stack[j];
      T otopo = topography[oval];
      int i = j-1;

      while(i >= 0 && topography[this->stack[i]] > otopo)
      {
         this->stack[i+1] = this->stack[i];
         --i;
      }
      this->stack[i+1] = oval;


    }
	
	}





	template<class Neighbourer_t, class topo_t, class out_t>
	out_t get_DA_proposlope(Neighbourer_t& neighbourer, topo_t& ttopography)
	{
		auto topography = format_input(ttopography);

		std::vector<double> DA(neighbourer.nnodes,0.);

		for(int i = neighbourer.nnodes - 1; i>=0; --i)
		{
			int node = this->stack[i];
			DA[node] += neighbourer.cellarea;

			if(neighbourer.is_active(node))
			{
				std::vector<double> slopes(this->receivers[node].size());
				double sumslopes = 0;
				for(int j = 0;j < this->receivers[node].size(); ++j)
				{
					int rec = this->receivers[node][j];
					slopes[j] = (topography[node] - topography[rec])/this->distance2receivers[node][j];
					if(slopes[j] <= 0)
						slopes[j] = 1e-5;
					sumslopes += slopes[j];
				}

				for(int j = 0;j < this->receivers[node].size(); ++j)
				{
					int rec = this->receivers[node][j];
					DA[rec] += DA[node] * slopes[j]/sumslopes;
				}

			}
		}

		return format_output(DA);
	}

	template<class Neighbourer_t, class topo_t, class out_t>
	out_t get_DA_SS(Neighbourer_t& neighbourer, topo_t& ttopography)
	{
		auto topography = format_input(ttopography);
		std::vector<double> DA(neighbourer.nnodes,0.);
		for(int i = neighbourer.nnodes - 1; i>=0; --i)
		{
			int node = this->stack[i];
			DA[node] += neighbourer.cellarea;

			if(neighbourer.is_active(node))
			{
				int Srec = this->Sreceivers[node];
				DA[Srec] += DA[node];
			}
		}

		return format_output(DA);
	}


	template<class Neighbourer_t>
	std::vector<std::vector<int> > get_rowcol_receivers(int row, int col,  Neighbourer_t& neighbourer)
	{
		int node = neighbourer.nodeid_from_row_col(row,col);
		std::vector<std::vector<int> > out_receivers;
		auto recs = this->receivers[node];
		for(auto r: recs)
		{
			int trow,tcol;
			neighbourer.rowcol_from_node_id(r,trow,tcol);
			out_receivers.emplace_back(std::vector<int>{trow,tcol});
		}
		return out_receivers;
	}

	template<class Neighbourer_t>
	std::vector<std::vector<int> > get_rowcol_donors(int row, int col,  Neighbourer_t& neighbourer)
	{
		int node = neighbourer.nodeid_from_row_col(row,col);
		std::vector<std::vector<int> > out_donors;
		auto recs = this->donors[node];
		for(auto r: recs)
		{
			int trow,tcol;
			neighbourer.rowcol_from_node_id(r,trow,tcol);
			out_donors.emplace_back(std::vector<int>{trow,tcol});
		}
		return out_donors;
	}

	template<class Neighbourer_t>
	std::vector<int> get_rowcol_Sreceivers(int row, int col,  Neighbourer_t& neighbourer)
	{
		int node = neighbourer.nodeid_from_row_col(row,col);
		std::vector<int> out_receivers;
		int trow,tcol;
		neighbourer.rowcol_from_node_id(this->Sreceivers[node],trow,tcol);
		out_receivers = std::vector<int>{trow,tcol};
		
		std::cout << "Srec is " << this->Sreceivers[this->Sreceivers[node]] << std::endl;
		return out_receivers;
	}

};






























#endif