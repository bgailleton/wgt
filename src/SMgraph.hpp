//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef smgraph_HPP
#define smgraph_HPP


/*
SMGraph stands for Static Multiple flow Graph.
It is a very memory efficient way to represent a graph and manipulate it.
It stores less informations than a fully fledged MGraph and therefore needs more processing to 
retrieve the same level of information (e.g. donors, distances, ...)
But is way more efficient in term of building/updating speed and extremely low in memory usage
It is useful for different cases, for example when updating and building the graph needs to be 
done multiple times and represents a significant part of the computation time (e.g. LEMs); or cases where
memory use becomes limitting.
It has a limitation though, the nodes need a fixed and even number of neighbours.

B.G. 2022
*/



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
// -> The neighbourer classes
#include "neighbourer.hpp"
// -> Manages operability with numpy types
#include "numvec.hpp"



class SMgraph
{

// Everything goes public, more straighforward
public:

	// Number of nodes in the graph
	int nnodes;
	// Number of neighbours by nodes
	int n_neighbours;
	// bool vector for each link: true is receiver direction and false is donor
	// The meaning of the index depends on the neighbourer
	std::vector<bool> isrec;
	std::vector<int> links;
	// Single graph receivers (optional)
	std::vector<int> Sreceivers,nSdonors,Sdonors;
	// Topological order and Single graph topological order (optional)
	std::vector<size_t> stack, Sstack;

	// Single graph distance to receivers
	std::vector<double> Sdistance2receivers;

	std::vector<double> SS;


	// default constructor
	SMgraph(){};
	SMgraph(int nnodes, int n_neighbours){this->nnodes = nnodes; this->n_neighbours = n_neighbours ;}


	void _allocate_vectors()
	{
		this->isrec = std::vector<bool>(int(this->nnodes * this->n_neighbours/2), false);
		this->links = std::vector<int>(int(this->nnodes * this->n_neighbours), 0);
		this->Sreceivers = std::vector<int>(this->nnodes,-1);
		this->Sstack = std::vector<size_t>(this->nnodes,0);
		for(int i=0;i<this->nnodes; ++i)
			this->Sreceivers[i] = i;
		this->Sdistance2receivers = std::vector<double >(this->nnodes,-1);
		this->SS = std::vector<double>(this->nnodes,0.);
	}
	void _reallocate_vectors()
	{

		for(int i=0;i<this->nnodes; ++i)
		{
			this->Sreceivers[i] = i;
			this->Sdistance2receivers[i] = 0;
			this->SS[i] = 0;
		}
	}

	template<class Neighbourer_t,class topo_t>
	void update_Mrecs(topo_t& topography, Neighbourer_t& neighbourer)
	{
		for(size_t i = 0; i<this->isrec.size(); ++i)
		{
			int from = this->links[i*2];
			int to = this->links[i*2 + 1];
			
			if(neighbourer.is_in_bound(from) == false || neighbourer.is_in_bound(to) == false)
				continue;

			if(topography[from] > topography[to])
				this->isrec[i] = true;
			else
				this->isrec[i] = false;
		}
	}
	template<class Neighbourer_t,class topo_t>
	void update_some_Mrecs(topo_t& topography, Neighbourer_t& neighbourer, std::vector<int>& some)
	{
		for(size_t j = 0; j<some.size(); ++j)
		{
			int node = some[j];
			auto ilinks = neighbourer.get_ilinks_from_node(node);
			for(auto i: ilinks)
			{
				int from = this->links[i*2];
				int to = this->links[i*2 + 1];
				
				if(neighbourer.is_in_bound(from) == false || neighbourer.is_in_bound(to) == false)
					continue;

				if(topography[from] > topography[to])
					this->isrec[i] = true;
				else
					this->isrec[i] = false;
			}
		}
	}

	template<class Neighbourer_t,class topo_t>
	void update_recs(topo_t& topography, Neighbourer_t& neighbourer)
	{
		for(size_t i = 0; i<this->isrec.size(); ++i)
		{
			int from = this->links[i*2];
			int to = this->links[i*2 + 1];
			
			if(neighbourer.is_in_bound(from) == false || neighbourer.is_in_bound(to) == false)
			{
				continue;
			}

			// if(neighbourer.is_active(from) == false || neighbourer.is_active(to) == false )
			// 	continue;

			double dx = neighbourer.get_dx_from_isrec_idx(i);
			double slope = (topography[from] - topography[to])/dx;
			// std::cout << from << "|" << to << "|";

			if(slope>0)
			{
				this->isrec[i] = true;
				if(this->SS[from]<slope)
				{
					this->Sreceivers[from] = to;
					this->Sdistance2receivers[from] = dx;
					this->SS[from] = slope;
				}
			}
			else
			{
				this->isrec[i] = false;
				slope = std::abs(slope);
				if(this->SS[to]<slope)
				{
					this->Sreceivers[to] = from;
					this->Sdistance2receivers[to] = dx;
					this->SS[to] = slope;
				}
			}

		}


	}

	template<class out_t>
	out_t test_Srecs()
	{
		std::vector<int> OUT(this->nnodes,0);
		for(int i=0; i< this->nnodes;++i)
		{
			if(i !=  this->Sreceivers[i])
				++OUT[i];
		}

		return format_output(OUT);
	}




	template<class Neighbourer_t>
	void init_graph_v6(Neighbourer_t& neighbourer)
	{
		// Allocate vectors
		this->_allocate_vectors();
		neighbourer.fill_links_SMG(this->links);
	}

	template<class Neighbourer_t>
	void reinit_graph_v6(Neighbourer_t& neighbourer)
	{
		// Allocate vectors
		this->_reallocate_vectors();
	}

	py::array_t<int,1> get_Sreceivers(){return py::array_t<int,1>(this->Sreceivers.size(),this->Sreceivers.data() ) ;} 
	py::array_t<double,1> get_dx_array(){return py::array_t<double,1>(this->Sdistance2receivers.size(),this->Sdistance2receivers.data() ) ;} 

	template<class Neighbourer_t,class topo_t, class out_t>
	out_t compute_graph_v6(std::string depression_solver, topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		// std::cout << "DEBUGGRAPH6::1" << std::endl;
		auto topography = format_input(ttopography);
		this->reinit_graph_v6(neighbourer);
		// std::cout << "DEBUGGRAPH6::2" << std::endl;
		this->update_recs(topography,neighbourer);
		// std::cout << "DEBUGGRAPH6::3" << std::endl;
		
		this->compute_SF_donors_from_receivers();
		// std::cout << "DEBUGGRAPH6::4" << std::endl;
		
		this->compute_TO_SF_stack_version();
		// std::cout << "DEBUGGRAPH6::5" << std::endl;

		std::vector<double> faketopo(to_vec(topography));
		
		LMRerouter_II depsolver;
		// std::cout << "DEBUGGRAPH6::prerun" << std::endl;
		bool need_recompute = depsolver.run(depression_solver, faketopo, neighbourer, this->Sreceivers, this->Sdistance2receivers, this->Sstack, this->links);
		// std::cout << "DEBUGGRAPH6::postrun" << std::endl;

		if(need_recompute)
		{
		
			this->recompute_SF_donors_from_receivers();
		
			this->compute_TO_SF_stack_version();

			if(depression_solver == "carve")
				this->carve_topo_v2(1e-5, neighbourer, faketopo);
					
			this->compute_MF_topological_order_insort(faketopo);
			this->update_Mrecs(faketopo,neighbourer);


			return format_output(faketopo);
		}
		else
		{
			// std::cout << "nodep" << std::endl;
			this->compute_MF_topological_order_insort(faketopo);
			return format_output(faketopo);	
		}

	}

	template<class Neighbourer_t,class topo_t, class out_t>
	out_t compute_graph_v6_SS(std::string depression_solver, topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		// std::cout << "DEBUGGRAPH6::1" << std::endl;
		auto topography = format_input(ttopography);
		this->reinit_graph_v6(neighbourer);
		// std::cout << "DEBUGGRAPH6::2" << std::endl;
		this->update_recs(topography,neighbourer);
		// std::cout << "DEBUGGRAPH6::3" << std::endl;
		
		this->compute_SF_donors_from_receivers();
		// std::cout << "DEBUGGRAPH6::4" << std::endl;
		
		this->compute_TO_SF_stack_version();
		// std::cout << "DEBUGGRAPH6::5" << std::endl;

		std::vector<double> faketopo(to_vec(topography));
		
		LMRerouter_II depsolver;
		// std::cout << "DEBUGGRAPH6::prerun" << std::endl;
		bool need_recompute = depsolver.run(depression_solver, faketopo, neighbourer, this->Sreceivers, this->Sdistance2receivers, this->Sstack, this->links);
		// std::cout << "DEBUGGRAPH6::postrun" << std::endl;

		if(need_recompute)
		{
		
			this->recompute_SF_donors_from_receivers();
		
			this->compute_TO_SF_stack_version();

			if(depression_solver == "carve")
				this->carve_topo_v2(1e-5, neighbourer, faketopo);
		
			// std::vector<int> to_recompute = this->carve_topo_v2(1e-4,neighbourer,faketopo);
			// this->fill_topo_v2(1e-5,neighbourer,faketopo);
			// std::cout << "sizetorecom::" << to_recompute.size() <<  std::endl;
			
			// this->compute_MF_topological_order_insort(faketopo);
		
			// this->update_some_Mrecs(faketopo,neighbourer,to_recompute);
			// this->update_recs(faketopo,neighbourer);


			return format_output(faketopo);
		}
		else
		{
			// std::cout << "nodep" << std::endl;
			// this->compute_MF_topological_order_insort(faketopo);
			return format_output(faketopo);	
		}

	}

	template<class Neighbourer_t,class topo_t, class out_t>
	out_t compute_graph_v6_nodep(topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		auto topography = format_input(ttopography);
		this->reinit_graph_v6(neighbourer);
		this->update_recs(topography,neighbourer);
		
		this->compute_SF_donors_from_receivers();
		
		this->compute_TO_SF_stack_version();

		std::vector<double> faketopo(to_vec(topography));
		this->compute_MF_topological_order_insort(faketopo);

	
		return format_output(faketopo);	

	}

	template<class Neighbourer_t,class topo_t, class out_t>
	out_t compute_graph_v6_PQ( topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		// std::cout << "DEBUGGRAPH6::1" << std::endl;
		auto topography = format_input(ttopography);

		std::vector<double> faketopo = neighbourer.PriorityFlood_Wei2018(topography);

		this->reinit_graph_v6(neighbourer);
		// std::cout << "DEBUGGRAPH6::2" << std::endl;
		this->update_recs(faketopo,neighbourer);
		// std::cout << "DEBUGGRAPH6::3" << std::endl;
		
		this->compute_SF_donors_from_receivers();
		// std::cout << "DEBUGGRAPH6::4" << std::endl;
		
		this->compute_TO_SF_stack_version();

		this->compute_MF_topological_order_insort(faketopo);
		// std::cout << "DEBUGGRAPH6::5" << std::endl;

		return format_output(faketopo);	


	}



	template<class Neighbourer_t,class topo_t, class out_t>
	out_t compute_graph(std::string depression_solver, topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		
		ocarina timer;


		timer.tik();

		auto topography = format_input(ttopography);
		add_noise_to_vector(topography,-1e-8,1e-8);

		this->isrec = std::vector<bool>(int(this->nnodes * this->n_neighbours/2), false);
		this->Sreceivers = std::vector<int>(this->nnodes,-1);
		this->Sstack = std::vector<size_t>(this->nnodes,0);
		for(int i=0;i<this->nnodes; ++i)
			this->Sreceivers[i] = i;
		this->Sdistance2receivers = std::vector<double >(this->nnodes,-1);
		
		this->SS = std::vector<double>(this->nnodes,0.);
		// timer.tok("Initiation");

		timer.tik();
		neighbourer.build_smgraph_SS_only_SS(topography, this->Sreceivers, this->Sdistance2receivers, this->SS);
		// timer.tok("Build_SS");

		timer.tik();
		this->compute_SF_donors_from_receivers();
		this->compute_TO_SF_stack_version();

		// timer.tok("SS_stuff");

		timer.tik();
		this->solve_depressions( depression_solver, neighbourer, topography);
		// timer.tok("Depression solver");

		timer.tik();

		std::vector<double> faketopo = to_vec(topography);

		this->compute_TO_SF_stack_version();

		for(int i=0; i < this->nnodes; ++i)
		{
			if(neighbourer.is_active(i) && i == this->Sreceivers[i])
				throw std::runtime_error("critical error");
		}

		if(depression_solver == "carve")
		{
			this->carve_topo(1e-4,neighbourer,faketopo);

		}
		else if (depression_solver == "fill")
			this->fill_topo(1e-4,neighbourer,faketopo);


		// timer.tok("Process topo post depression");
		for(int i=0; i < this->nnodes; ++i)
		{
			if(neighbourer.is_active(i) && faketopo[i] <= faketopo[this->Sreceivers[i]])
				throw std::runtime_error("critical error in carvefill");
		}


		
		timer.tik();
		neighbourer.build_smgraph_only_MF(faketopo, this->isrec);



		this->compute_MF_topological_order_insort(faketopo);
		
		// timer.tok("MF stuffs");


		return format_output(faketopo);
	}

	template<class Neighbourer_t, class topo_t, class out_t>
	out_t compute_graph_multi_filled( topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		
		auto topography = format_input(ttopography);
		// ocarina timer;
		// timer.tik();

		// // auto topography = format_input(ttopography);
		// add_noise_to_vector(topography,-1e-6,1e-6);

		this->isrec = std::vector<bool>(int(this->nnodes * this->n_neighbours/2), false);

		// timer.tik();
		// std::vector<double> faketopo = neighbourer.fill_barne_2014(topography);
		std::vector<double> faketopo = neighbourer.PriorityFlood_Wei2018(topography);

		// add_noise_to_vector(faketopo,-1e-6,1e-6);
		// timer.tok("filled depression");

		// timer.tik();
		neighbourer.build_smgraph_only_MF(faketopo, this->isrec);

		this->compute_MF_topological_order_insort(faketopo);
		
		// timer.tok("MF stuffs");


		return format_output(faketopo);
	}
	template<class Neighbourer_t, class topo_t, class out_t>
	out_t just_fill_it( topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		auto topography = format_input(ttopography);
		std::vector<double> faketopo = neighbourer.PriorityFlood_Wei2018(topography);
		return format_output(faketopo);
	}


	template<class Neighbourer_t, class topo_t, class out_t>
	out_t update_graph_multi_filled( topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		
		auto topography = format_input(ttopography);
		ocarina timer;
		timer.tik();

		// auto topography = format_input(ttopography);
		// add_noise_to_vector(topography,-1e-6,1e-6);

		// this->isrec = std::vector<bool>(int(this->nnodes * this->n_neighbours/2), false);

		timer.tik();
		// std::vector<double> faketopo = neighbourer.fill_barne_2014(topography);
		std::vector<double> faketopo = neighbourer.PriorityFlood_Wei2018(topography);

		// add_noise_to_vector(faketopo,-1e-6,1e-6);
		// timer.tok("filled depression");

		timer.tik();
		neighbourer.build_smgraph_only_MF(faketopo, this->isrec);

		this->compute_MF_topological_order_insort(faketopo);
		
		// timer.tok("MF stuffs");


		return format_output(faketopo);
	}

	template<class Neighbourer_t, class topo_t, class out_t>
	out_t compute_graph_single_filled( topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		
		auto topography = format_input(ttopography);
		ocarina timer;
		timer.tik();

		// auto topography = format_input(ttopography);
		add_noise_to_vector(topography,-1e-6,1e-6);

		std::vector<double> faketopo = neighbourer.PriorityFlood_Wei2018(topography);
		
		this->Sreceivers = std::vector<int>(this->nnodes,-1);
		this->Sstack = std::vector<size_t>(this->nnodes,0);
		for(int i=0;i<this->nnodes; ++i)
			this->Sreceivers[i] = i;
		this->Sdistance2receivers = std::vector<double >(this->nnodes,-1);
		
		this->SS = std::vector<double>(this->nnodes,0.);

		neighbourer.build_smgraph_SS_only_SS(topography, this->Sreceivers, this->Sdistance2receivers, SS);

		this->compute_SF_donors_from_receivers();

		this->compute_TO_SF_stack_version();


		return format_output(faketopo);
	}


	template<class Neighbourer_t, class topo_t, class out_t>
	out_t compute_graph_multi_filled_old( topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		
		auto topography = format_input(ttopography);
		ocarina timer;
		timer.tik();

		// auto topography = format_input(ttopography);
		add_noise_to_vector(topography,-1e-6,1e-6);

		this->isrec = std::vector<bool>(int(this->nnodes * this->n_neighbours/2), false);

		timer.tik();
		std::vector<double> faketopo = neighbourer.fill_barne_2014(topography);
		// std::vector<double> faketopo = neighbourer.PriorityFlood_Wei2018(topography);
		
		// add_noise_to_vector(faketopo,-1e-6,1e-6);
		// timer.tok("filled depression");

		timer.tik();
		neighbourer.build_smgraph_only_MF(faketopo, this->isrec);

		this->compute_MF_topological_order_insort(faketopo);
		
		// timer.tok("MF stuffs");


		return format_output(faketopo);
	}



	template<class Neighbourer_t, class topo_t, class out_t>
	out_t compute_graph_multi_filled_OMP( topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		py::gil_scoped_release release;
		
		ocarina timer;
		timer.tik();

		auto topography = format_input(ttopography);
		add_noise_to_vector(topography,-1e-6,1e-6);

		this->isrec = std::vector<bool>(int(this->nnodes * this->n_neighbours/2), false);

		// timer.tik();
		std::vector<double> faketopo = neighbourer.fill_barne_2014(topography);
		// timer.tok("filled depression");

		// timer.tik();
		neighbourer.build_smgraph_only_MF_OMP(faketopo, this->isrec);

		this->compute_MF_topological_order_insort(faketopo);
		
		// timer.tok("MF stuffs");

		py::gil_scoped_acquire acquire;

		return format_output(faketopo);
	}

	template<class Neighbourer_t,class topo_t, class out_t>
	out_t compute_graph_OMP(std::string depression_solver, topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		// py::gil_scoped_release release;
		// // std::cout << "made it here -2" << std::endl;
		// auto topography = format_input(ttopography);
		// // std::cout << "made it here -1,9" << std::endl;
		// add_noise_to_vector(topography,-1e-6,1e-6);
		// // std::cout << "made it here -1" << std::endl;

		// this->isrec = std::vector<bool>(int(this->nnodes * this->n_neighbours/2), false);
		// this->Sreceivers = std::vector<int>(this->nnodes,-1);
		// for(int i=0;i<this->nnodes; ++i)
		// 	this->Sreceivers[i] = i;
		// this->Sdistance2receivers = std::vector<double >(this->nnodes,-1);
		
		// std::vector<double> SS(this->nnodes,0.);
		// neighbourer.build_smgraph_SS_only_SS_OMP(topography, this->Sreceivers, this->Sdistance2receivers, SS);
		// this->compute_SF_donors_from_receivers();
		// this->compute_TO_SF_stack_version();

		// this->solve_depressions( depression_solver, neighbourer, topography);
		// std::vector<double> faketopo = to_vec(topography);

		// this->compute_TO_SF_stack_version();

		// if(depression_solver == "carve")
		// {
		// 	this->carve_topo(1e-3,neighbourer,faketopo);

		// }
		// else if (depression_solver == "fill")
		// 	this->fill_topo(1e-3,neighbourer,faketopo);


		// neighbourer.build_smgraph_only_MF_OMP(faketopo, this->isrec);
		
		// this->compute_MF_topological_order_insort(faketopo);
		// py::gil_scoped_acquire acquire;
		
		// return format_output(faketopo);
		return format_output(ttopography);
	}

	template<class Neighbourer_t,class topo_t, class out_t>
	out_t update_graph(std::string depression_solver, topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		auto topography = format_input(ttopography);

		for(int i=0;i<this->nnodes; ++i)
			this->Sreceivers[i] = i;

		this->SS = std::vector<double>(this->nnodes,0.);

		neighbourer.build_smgraph_SS_only_SS(topography, this->Sreceivers, this->Sdistance2receivers, SS);
		this->compute_SF_donors_from_receivers();
		
		this->compute_TO_SF_stack_version();
		
		this->solve_depressions( depression_solver, neighbourer, topography);
		
		this->compute_TO_SF_stack_version();

		std::vector<double> faketopo = to_vec(topography);

		if(depression_solver == "carve")
			this->carve_topo(1e-3,neighbourer,faketopo);
		else if (depression_solver == "fill")
			this->fill_topo(1e-3,neighbourer,faketopo);
		
		neighbourer.build_smgraph_only_MF(faketopo, this->isrec);
		
		this->compute_MF_topological_order_insort(faketopo);
	
		return format_output(faketopo);
	}




	template<class Neighbourer_t, class topo_t>
	void compute_graph_single_solo( topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		
		auto topography = format_input(ttopography);

		this->Sreceivers = std::vector<int>(this->nnodes,-1);
		this->Sstack = std::vector<size_t>(this->nnodes,0);
		for(int i=0;i<this->nnodes; ++i)
			this->Sreceivers[i] = i;

		this->Sdistance2receivers = std::vector<double >(this->nnodes,-1);
		
		this->SS = std::vector<double>(this->nnodes,0.);

		neighbourer.build_smgraph_SS_only_SS(topography, this->Sreceivers, this->Sdistance2receivers, SS);

		this->compute_SF_donors_from_receivers();

		this->compute_TO_SF_stack_version();

	}

	template<class Neighbourer_t, class topo_t>
	void update_graph_single_solo( topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		
		auto topography = format_input(ttopography);

		for(int i=0;i<this->nnodes; ++i)
		{
			this->Sreceivers[i] = i;		
			this->SS[i] = 0;
			this->Sdistance2receivers[i] = -1;
		}
		
		neighbourer.build_smgraph_SS_only_SS(topography, this->Sreceivers, this->Sdistance2receivers, this->SS);

		this->recompute_SF_donors_from_receivers();

		this->compute_TO_SF_stack_version();
	}

	template<class Neighbourer_t, class topo_t>
	void compute_graph_multi_solo( topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		
		auto topography = format_input(ttopography);
		this->isrec = std::vector<bool>(int(this->nnodes * this->n_neighbours/2), false);
		neighbourer.build_smgraph_only_MF(topography, this->isrec);
		this->compute_MF_topological_order_insort(topography);
	}

	template<class Neighbourer_t, class topo_t>
	void update_graph_multi_solo( topo_t& ttopography, Neighbourer_t& neighbourer)
	{		
		auto topography = format_input(ttopography);
		neighbourer.build_smgraph_only_MF(topography, this->isrec);
		this->compute_MF_topological_order_insort(topography);
	}


	template<class Neighbourer_t, class topo_t>
	void compute_graph_both_solo( topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		
		auto topography = format_input(ttopography);
		this->isrec = std::vector<bool>(int(this->nnodes * this->n_neighbours/2), false);
		this->Sreceivers = std::vector<int>(this->nnodes,-1);
		this->Sstack = std::vector<size_t>(this->nnodes,0);
		this->Sdistance2receivers = std::vector<double >(this->nnodes,-1);
		this->SS = std::vector<double>(this->nnodes,0.);

		for(int i=0;i<this->nnodes; ++i)
			this->Sreceivers[i] = i;
		
		neighbourer.build_smgraph_only_MF(topography, this->isrec);
		neighbourer.buildSrecfromIsrec(topography, this->isrec, this->Sreceivers, this->SS);
		this->compute_SF_donors_from_receivers();
		this->compute_TO_SF_stack_version();
		this->compute_MF_topological_order_insort(topography);
	}

	template<class Neighbourer_t, class topo_t>
	void update_graph_both_solo( topo_t& ttopography, Neighbourer_t& neighbourer)
	{		
		auto topography = format_input(ttopography);

		for(int i=0;i<this->nnodes; ++i)
		{
			this->SS[i] = 0;
			this->Sreceivers[i] = i;
		}
		
		neighbourer.build_smgraph_only_MF(topography, this->isrec);
		neighbourer.buildSrecfromIsrec(topography, this->isrec, this->Sreceivers, this->SS);

		this->compute_TO_SF_stack_version();
		this->compute_MF_topological_order_insort(topography);
	}


	template<class Neighbourer_t,class topo_t>
	std::vector<double> integrated_simple_fill(topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		auto topography = format_input(ttopography);

		std::vector<bool> isdep(this->nnodes, false);
		int ndep = 0;
		for(int i=0; i<this->nnodes;++i)
		{
			int node = this->Sstack[i];

			if(neighbourer.is_active(node) && this->Sreceivers[node] == node)
			{
				// std::cout << node << std::endl;
				isdep[node] = true;
				++ndep;
			}

			isdep[node] = isdep[this->Sreceivers[node]];
		}
		std::cout << "There are " << ndep << " depressions vs " << this->nnodes << std::endl;

		std::vector<double> faketopo = to_vec(topography);

		// std::cout << "WABUNDEBUG1 " << std::endl;
		neighbourer.fill_connected_component(isdep, ndep, faketopo, this->Sreceivers);
		// std::cout << "WABUNDEBUG2" << std::endl;

		return faketopo;

	}


	template<class Neighbourer_t,class topo_t, class out_t>
	out_t compute_graph_v4(topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		auto topography = format_input(ttopography);
		std::vector<double> faketopo = neighbourer.PriorityFlood_Wei2018(topography);
		this->compute_graph_both_solo(faketopo,neighbourer);
		return format_output(faketopo);
	}

	template<class Neighbourer_t,class topo_t, class out_t>
	out_t compute_graph_v4_simple(topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		auto topography = format_input(ttopography);
		bool worked = true;
		this->compute_graph_both_solo(topography,neighbourer);

		// std::vector<double> faketopo(to_vec(topography));
		std::vector<double> faketopo(this->nnodes);
		for(int i=0; i< this->nnodes; ++i)
			faketopo[i] = topography[i];
		std::vector<bool> filled = neighbourer.simple_fill(faketopo,this->Sreceivers);

		this->compute_graph_both_solo(faketopo,neighbourer);
		// neighbourer.build_smgraph_only_MF_mask(faketopo, this->isrec,filled);
		// neighbourer.buildSrecfromIsrec_mask(faketopo, this->isrec, this->Sreceivers, this->SS,filled);
		// this->compute_SF_donors_from_receivers();
		// this->compute_TO_SF_stack_version();
		// this->compute_MF_topological_order_insort(faketopo);

		for(int i = 0; i < this->nnodes; ++i)
		{
			if(neighbourer.is_active(i) && i == this->Sreceivers[i])
				worked = false;
		}


		if(worked == false)
		{
			std::cout << "Backing up" << std::endl;
			auto out = this->compute_graph_v4<Neighbourer_t,std::vector<double>, py::array >(faketopo,neighbourer);
			return format_output(out);
		}


		return format_output(faketopo);
	}

	template<class Neighbourer_t,class topo_t, class out_t>
	out_t compute_graph_v4_SF(topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		auto topography = format_input(ttopography);
		this->compute_graph_single_solo(topography,neighbourer);
		std::vector<double> faketopo = this->integrated_simple_fill(topography,neighbourer);
		this->recompute_SF_donors_from_receivers();
		this->compute_TO_SF_stack_version();
		return format_output(faketopo);
	}
	template<class Neighbourer_t,class topo_t, class out_t>
	out_t update_graph_v4_SF(topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		auto topography = format_input(ttopography);
		this->update_graph_single_solo(topography,neighbourer);
		std::vector<double> faketopo = this->integrated_simple_fill(topography,neighbourer);
		this->recompute_SF_donors_from_receivers();
		this->compute_TO_SF_stack_version();
		return format_output(faketopo);
	}

	template<class Neighbourer_t,class topo_t, class out_t>
	out_t update_graph_v4(topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		auto topography = format_input(ttopography);
		this->update_graph_single_solo(topography,neighbourer);
		std::vector<double> faketopo = this->integrated_simple_fill(topography,neighbourer);
		this->update_graph_multi_solo(faketopo,neighbourer);
		return format_output(faketopo);
	}


	template<class Neighbourer_t,class topo_t, class out_t>
	out_t compute_graph_v5(std::string method,  topo_t& ttopography, Neighbourer_t& neighbourer)
	{

		ocarina timer;
		auto topography = format_input(ttopography);
		this->compute_graph_single_solo(topography,neighbourer);
		timer.tik();
		this->solve_depressions(method,neighbourer,topography);
		// this->recompute_SF_donors_from_receivers();
		this->compute_TO_SF_stack_version();
		std::vector<double> faketopo = std::vector<double>(to_vec(topography));
		this->fill_topo(1e-3,neighbourer,faketopo);
		add_noise_to_vector(topography, -1e-6, 1e-6);

		timer.tok("Fulldepsolve::");
		this->compute_graph_multi_solo(faketopo,neighbourer);

		return format_output(faketopo);

	}

	template<class Neighbourer_t,class topo_t, class out_t>
	out_t update_graph_v5(std::string method,  topo_t& ttopography, Neighbourer_t& neighbourer)
	{

		ocarina timer;
		auto topography = format_input(ttopography);
		timer.tik();
		add_noise_to_vector(topography, -1e-6, 1e-6);
		timer.tok("random");
		this->update_graph_single_solo(topography,neighbourer);
		timer.tik();
		this->solve_depressions(method,neighbourer,topography);
		this->recompute_SF_donors_from_receivers();
		this->compute_TO_SF_stack_version();
		std::vector<double> faketopo = to_vec(topography);
		this->fill_topo(1e-4,neighbourer,faketopo);
		timer.tok("Fulldepsolve::");
		this->update_graph_multi_solo(faketopo,neighbourer);

		return format_output(faketopo);

	}




	template<class topo_t>
	void compute_MF_topological_order_insort(topo_t& ttopography)
	{
		auto topography = format_input(ttopography);


		auto yolo = sort_indexes(topography);
		this->stack = std::move(yolo);

	}

	void compute_SF_donors_from_receivers()
	{
		// Initialising the graph dimesions for the donors
		// All of thenm have the graph dimension
		this->Sdonors = std::vector<int>(this->nnodes * 8,-1);
		this->nSdonors = std::vector<int>(this->nnodes,0);

		for(int i=0; i < this->nnodes; ++i)
		{
			// SF so rid == i cause there is only 1 rec
			int trec = this->Sreceivers[i];
			if(trec == i)
				continue;

			this->Sdonors[trec * 8  + this->nSdonors[trec]] = i;
			this->nSdonors[trec] += 1;
		}

	}

	void recompute_SF_donors_from_receivers()
	{

		for(int i=0; i < this->nnodes; ++i)
		{
			for(int j=0; j<8; ++j)
				this->Sdonors[i * 8 + j] = -1;
			this->nSdonors[i] = 0;
		}

		for(int i=0; i < this->nnodes; ++i)
		{
			// SF so rid == i cause there is only 1 rec
			int trec = this->Sreceivers[i];
			if(trec == i)
				continue;
			this->Sdonors[trec * 8  + this->nSdonors[trec]] = i;
			this->nSdonors[trec] += 1;
		}

	}

	template<class Neighbourer_t, class topo_t>
	void solve_depressions(std::string& depression_solver , Neighbourer_t& neighbourer, topo_t& ttopography)
	{
		auto topography = format_input(ttopography);

		
		LMRerouter<int, double,  Neighbourer_t, topo_t> depsolver(neighbourer, topography, this->Sreceivers,  this->Sdistance2receivers, this->Sstack);


		if(depsolver.npits > 0)
		{
			std::cout << "updating deps" << std::endl;
			depsolver.update_receivers(depression_solver,neighbourer, topography, this->Sreceivers,  this->Sdistance2receivers, this->Sstack);
			this->recompute_SF_donors_from_receivers();
		}
	}

	void compute_TO_SF_stack_version()
	{
		// The stack container helper
		std::stack<size_t, std::vector<size_t> > stackhelper;
		// std::vector<bool> isdone(this->nnodes,false);
		// going through all the nodes
		int istack = 0;
		for(int i=0; i<this->nnodes; ++i)
		{
			// if they are base level I include them in the stack
			if(this->Sreceivers[i] == i)
			{
				stackhelper.emplace(i);
				// ++istack;
			}

			// While I still have stuff in the stack helper
			while(stackhelper.empty() == false)
			{
				// I get the next node and pop it from the stack helper
				int nextnode = stackhelper.top();stackhelper.pop();
				this->Sstack[istack] = nextnode;
				++istack;

				// as well as all its donors which will be processed next
				for( int j = 0; j < this->nSdonors[nextnode]; ++j)
				{
					stackhelper.emplace(this->Sdonors[nextnode*8 + j]);
				}

			}

		}
	}
	

	/// this function enforces minimal slope 
	template<class Neighbourer_t, class topo_t>
	void carve_topo(double slope, Neighbourer_t& neighbourer, topo_t& topography)
	{

		std::cout << std::setprecision(8);
		for(int i=this->nnodes-1; i >= 0; --i)
		{
			int node  = this->Sstack[i];
			
			// if(node == 148880)
				// std::cout << "ASSESSED" << std::endl;

			if(neighbourer.can_flow_out_there(node) || neighbourer.can_flow_even_go_there(node) == false)
				continue;
			// if(node == 148880)
				// std::cout << "PASSED" << std::endl;
			int rec = this->Sreceivers[node];
			double dz = topography[node] - topography[rec];
			if(dz <= 0)
			{
				topography[rec] = topography[node] - slope + neighbourer.randu.get()* 1e-7;// * d2rec;
			}

		}
	}

	/// this function enforces minimal slope 
	template<class Neighbourer_t, class topo_t>
	std::vector<int> carve_topo_v2(double slope, Neighbourer_t& neighbourer, topo_t& topography)
	{

		std::cout << std::setprecision(8);
		std::vector<int> to_recompute;
		to_recompute.reserve(1000);
		for(int i=this->nnodes-1; i >= 0; --i)
		{
			int node  = this->Sstack[i];
			
			// if(node == 148880)
				// std::cout << "ASSESSED" << std::endl;

			if(neighbourer.can_flow_out_there(node) || neighbourer.can_flow_even_go_there(node) == false)
				continue;
			// if(node == 148880)
				// std::cout << "PASSED" << std::endl;
			int rec = this->Sreceivers[node];
			double dz = topography[node] - topography[rec];
			if(dz <= 0)
			{
				topography[rec] = topography[node] - slope + neighbourer.randu.get() * 1e-7;// * d2rec;
				to_recompute.emplace_back(rec);
			}

		}
		return to_recompute;
	}
	/// this function enforces minimal slope 
	template<class Neighbourer_t, class topo_t>
	std::vector<int> fill_topo_v2(double slope, Neighbourer_t& neighbourer, topo_t& topography)
	{
		std::vector<int> to_recompute;
		for(int i=0; i < this->nnodes; ++i)
		{
			int node  = this->Sstack[i];
			if(neighbourer.can_flow_out_there(node) || neighbourer.can_flow_even_go_there(node) == false)
				continue;

			int rec = this->Sreceivers[node];
			double dz = topography[node] - topography[rec];

			if(dz <= 0)
			{
				topography[node] = topography[rec] + slope + neighbourer.randu.get() * 1e-6;// * d2rec;
				to_recompute.emplace_back(node);
			}
		}
		return to_recompute;
	}


	/// this function enforces minimal slope 
	template<class Neighbourer_t, class topo_t>
	void fill_topo(double slope, Neighbourer_t& neighbourer, topo_t& topography)
	{

		for(int i=0; i < this->nnodes; ++i)
		{
			int node  = this->Sstack[i];
			if(neighbourer.can_flow_out_there(node) || neighbourer.can_flow_even_go_there(node) == false)
				continue;

			int rec = this->Sreceivers[node];
			double dz = topography[node] - topography[rec];

			if(dz <= 0)
			{
				// double d2rec = this->Sdistance2receivers[node];
				topography[node] = topography[rec] + slope + neighbourer.randu.get()* 1e-7;// * d2rec;
			}
		}
	}

	template<class Neighbourer_t>
	std::vector<int> get_receiver_indices(int i, Neighbourer_t& neighbourer)
	{
		std::vector<int> recs; recs.reserve(8);
		int i_r = neighbourer.get_id_right_SMG(i);
		if(i_r >=0 || i_r < neighbourer.nnodes * 4)
		{
			if(this->isrec[i_r])
				recs.emplace_back(neighbourer.get_right_index(i));
		}
		i_r = neighbourer.get_id_bottomright_SMG(i);
		if(i_r >=0 || i_r < neighbourer.nnodes * 4)
		{
			if(this->isrec[i_r])
				recs.emplace_back(neighbourer.get_bottomright_index(i));
		}
		i_r = neighbourer.get_id_bottom_SMG(i);
		if(i_r >=0 || i_r < neighbourer.nnodes * 4)
		{
			if(this->isrec[i_r])
				recs.emplace_back(neighbourer.get_bottom_index(i));
		}
		i_r = neighbourer.get_id_bottomleft_SMG(i);
		if(i_r >=0 || i_r < neighbourer.nnodes * 4)
		{
			if(this->isrec[i_r])
				recs.emplace_back(neighbourer.get_bottomleft_index(i));
		}
		i_r = neighbourer.get_id_left_SMG(i);
		if(i_r >=0 || i_r < neighbourer.nnodes * 4)
		{
			if(this->isrec[i_r] == false)
				recs.emplace_back(neighbourer.get_left_index(i));
		}
		i_r = neighbourer.get_id_topleft_SMG(i);
		if(i_r >=0 || i_r < neighbourer.nnodes * 4)
		{
			if(this->isrec[i_r] == false)
				recs.emplace_back(neighbourer.get_topleft_index(i));
		}
		i_r = neighbourer.get_id_top_SMG(i);
		if(i_r >=0 || i_r < neighbourer.nnodes * 4)
		{
			if(this->isrec[i_r] == false)
				recs.emplace_back(neighbourer.get_top_index(i));
		}
		i_r = neighbourer.get_id_topright_SMG(i);
		if(i_r >=0 || i_r < neighbourer.nnodes * 4)
		{
			if(this->isrec[i_r] == false)
				recs.emplace_back(neighbourer.get_topright_index(i));
		}
		return recs;

	}

	template<class Neighbourer_t>
	std::vector<std::pair<int,int> > get_receiver_link_indices(int i, Neighbourer_t& neighbourer)
	{
		std::vector<std::pair<int,int>> recs; recs.reserve(8);

		int i_r = neighbourer.get_id_right_SMG(i);
		if(i_r >=0 && i_r < neighbourer.nnodes * 4)
		{
			int ti = neighbourer.get_right_index(i);
			if(this->isrec[i_r] && ti >=0 && ti<this->nnodes)
			{
				recs.emplace_back(std::pair<int,int>{ti,i_r});
			}
		}
		i_r = neighbourer.get_id_bottomright_SMG(i);
		if(i_r >=0 && i_r < neighbourer.nnodes * 4)
		{
			int ti = neighbourer.get_bottomright_index(i);
			if(this->isrec[i_r] && ti >=0 && ti<this->nnodes)
			{
				recs.emplace_back(std::pair<int,int>{ti,i_r});
			}
		}
		i_r = neighbourer.get_id_bottom_SMG(i);
		if(i_r >=0 && i_r < neighbourer.nnodes * 4)
		{
			int ti = neighbourer.get_bottom_index(i);
			if(this->isrec[i_r] && ti >=0 && ti<this->nnodes)
			{
				recs.emplace_back(std::pair<int,int>{ti,i_r});
			}
		}
		i_r = neighbourer.get_id_bottomleft_SMG(i);
		if(i_r >=0 && i_r < neighbourer.nnodes * 4)
		{
			int ti = neighbourer.get_bottomleft_index(i);
			if(this->isrec[i_r] && ti >=0 && ti<this->nnodes)
			{
				recs.emplace_back(std::pair<int,int>{ti,i_r});
			}
		}
		i_r = neighbourer.get_id_left_SMG(i);
		if(i_r >=0 && i_r < neighbourer.nnodes * 4)
		{
			int ti = neighbourer.get_left_index(i);
			if(this->isrec[i_r] == false && ti >=0 && ti<this->nnodes)
			{
				recs.emplace_back(std::pair<int,int>{ti,i_r});
			}
		}
		i_r = neighbourer.get_id_topleft_SMG(i);
		if(i_r >=0 && i_r < neighbourer.nnodes * 4)
		{
			int ti = neighbourer.get_topleft_index(i);
			if(this->isrec[i_r] == false && ti >=0 && ti<this->nnodes)
			{
				recs.emplace_back(std::pair<int,int>{ti,i_r});
			}
		}
		i_r = neighbourer.get_id_top_SMG(i);
		if(i_r >=0 && i_r < neighbourer.nnodes * 4)
		{
			int ti = neighbourer.get_top_index(i);
			if(this->isrec[i_r] == false && ti >=0 && ti<this->nnodes)
			{
				recs.emplace_back(std::pair<int,int>{ti,i_r});
			}
		}
		i_r = neighbourer.get_id_topright_SMG(i);
		if(i_r >=0 && i_r < neighbourer.nnodes * 4)
		{
			int ti = neighbourer.get_topright_index(i);
			if(this->isrec[i_r] == false && ti >=0 && ti<this->nnodes)
			{
				recs.emplace_back(std::pair<int,int>{ti,i_r});
			}
		}
		return recs;

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
				auto receivers = this->get_receiver_indices(node, neighbourer);

				std::vector<double> slopes(receivers.size());
				double sumslopes = 0;
				for(size_t j = 0;j < receivers.size(); ++j)
				{
					int rec = receivers[j];
					slopes[j] = (topography[node] - topography[rec])/neighbourer.dx;
					if(slopes[j] <= 0)
						slopes[j] = 1e-5;
					sumslopes += slopes[j];
				}

				for(size_t j = 0;j < receivers.size(); ++j)
				{
					int rec = receivers[j];
					DA[rec] += DA[node] * slopes[j]/sumslopes;
				}

			}

			// std::cout << std::endl;;

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
			int node = this->Sstack[i];
			DA[node] += neighbourer.cellarea;

			if(neighbourer.can_flow_even_go_there(node) && node != this->Sreceivers[node])
			{
				int Srec = this->Sreceivers[node];
				DA[Srec] += DA[node];
			}
		}

		return format_output(DA);
	}


	template<class Neighbourer_t>
	std::vector<int> get_rowcol_Sreceivers(int row, int col,  Neighbourer_t& neighbourer)
	{
		int node = neighbourer.nodeid_from_row_col(row,col);
		std::vector<int> out_receivers;
		int trow,tcol;
		neighbourer.rowcol_from_node_id(this->Sreceivers[node],trow,tcol);
		out_receivers = std::vector<int>{trow,tcol};
		
		std::cout << "Srec is " << this->Sreceivers[node] << " node was " << node << std::endl;
		return out_receivers;
	}


	template<class Neighbourer_t, class topo_t>
	void print_receivers(int i,Neighbourer_t& neighbourer, topo_t& ttopography)
	{
		std::cout << std::setprecision(12);
		auto topography = format_input(ttopography);
		auto receivers = this->get_receiver_indices(i, neighbourer);

		std::cout << "Topography is " << topography[i] << "# receivers: " << receivers.size() << std::endl;
		for(auto r: receivers)
		{
			int row,col;
			neighbourer.rowcol_from_node_id(r,row,col);
			std::cout << "Rec " << r << " row " << row << " col " << col << " topo " << topography[r] << std::endl;

		}


		auto neighbours = neighbourer.get_neighbours_only_id(i);
		std::cout << "Neighbours are :" << std::endl;

		for(auto r: neighbours)
		{
			int row,col;
			neighbourer.rowcol_from_node_id(r,row,col);
			std::cout << "Neighbour " << r << " row " << row << " col " << col << " topo " << topography[r] << std::endl;
		}




	}


	int get_rec_array_size(){return int(this->isrec.size());}


	/// Takes an array of nnodes size and sum the values at the outlets
	/// This can be useful for checking mass balances for example
	/// if true, include_internal_pits allow the code to add internal unprocessed pits, wether they are on purpose or not
	template<class Neighbourer_t,class array_t, class T>
	T sum_at_outlets(Neighbourer_t& neighbourer, array_t& tarray, bool include_internal_pits = true)
	{
		auto array = format_input(tarray);
		T out = 0;
		for(int i=0; i<this->nnodes; ++i)
		{
			if (this->Sreceivers[i] == i)
			{
				if(include_internal_pits)
				{
					out += array[i];
				}
				else if(neighbourer.can_flow_out_there(i) )
				{
					out += array[i];
				}
			}
		}
		return out;

	}

	/// Takes an array of nnodes size and sum the values at the outlets
	/// This can be useful for checking mass balances for example
	/// if true, include_internal_pits allow the code to add internal unprocessed pits, wether they are on purpose or not
	template<class Neighbourer_t,class array_t, class out_t>
	out_t keep_only_at_outlets(Neighbourer_t& neighbourer, array_t& tarray, bool include_internal_pits = true)
	{
		auto array = format_input(tarray);
		std::vector<double> out = std::vector<double> (this->nnodes,0);
		for(int i=0; i<this->nnodes; ++i)
		{
			if (this->Sreceivers[i] == i)
			{
				if(include_internal_pits)
					out[i] = array[i];
				else if(neighbourer.can_flow_out_there(i) )
					out[i] = array[i];
			}
		}
		return format_output(out);

	}

	bool is_Sstack_full()
	{
		if(int(this->Sstack.size()) != this->nnodes)
		{
			std::cout << "stack size (" << this->Sstack.size() << ") is invalid." << std::endl;
			return false;
		}
		std::vector<int> ntimenodes(this->nnodes,0);

		for(auto v:this->Sstack)
		{
			++ntimenodes[v];
		}

		int n_0 = 0,n_p1 = 0;
		for(int i=0; i<this->nnodes; ++i)
		{
			if(ntimenodes[i] == 0)
				++n_0;
			else if(ntimenodes[i]>1)
				++n_p1;
		}

		if(n_0 > 0 || n_p1 > 0)
		{
			std::cout << "Stack issue: " << n_p1 << " nodes appearing more than once and " << n_0 << " nodes not appearing" << std::endl;
			return false;
		}

		std::vector<bool> isdone(this->nnodes,false);
		for(int i = this->nnodes - 1; i>=0; --i)
		{
			auto v = this->Sstack[i];
			isdone[v] = true;
			if(int(v) !=  this->Sreceivers[v])
			{
				if(isdone[this->Sreceivers[v]])
				{
					std::cout << "Receiver processed before node stack is fucked" << std::endl;
					return false;
				}
			}
		}

		return true;

	}


	std::vector<bool> has_Srecs()
	{
		std::vector<bool> haSrecs(this->nnodes,true);
		for(int i=0;i<this->nnodes;++i)
		{
			if(this->Sreceivers[i] == i)
				haSrecs[i] = false;
		}
		return haSrecs;
	}







};

























#endif
