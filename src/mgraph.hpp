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

// local includes 
// -> General routines and data structures
#include "chonkutils.hpp"
// -> Depression solvers
#include "cordonnier_versatile_2019.hpp"
// -> lightweigth numpy export
#include "npy.hpp"

#include "neighbourer.hpp"

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

	std::vector<int> Sreceivers;
	std::vector<std::vector<int> > Sdonors;
	std::vector<T> Sdistance2receivers;

	std::vector<std::vector<int> > receivers;
	std::vector<std::vector<int> > donors;
	std::vector<std::vector<T> > distance2receivers;


		// #->stack: topological order from downstream to upstream direction
	std::vector<int> stack;




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

	// ------------------------------------------------

	//	                             	              __
	//                                             / _)
	//                                    _.----._/ /
	//   GRAPH CREATION METHODS          /         /
	//                                __/ (  | (  |
	//                               /__.-'|_|--|_|

	// All the methods related to accessing and calculating neighbours
	// ------------------------------------------------

	template<class Neighbourer_t, class topo_t>
	void compute_graph(std::string depression_solver, Neighbourer_t& neighbourer, topo_t& topography)
	{
		// std::cout << "ASDFdsfa->1" << std::endl;
		this->compute_graph_both_v2(neighbourer,topography);
		// std::cout << "ASDFdsfa->2" << std::endl;
		this->compute_TO_SF_stack_version();
		// std::cout << "ASDFdsfa->3" << std::endl;
		if(depression_solver != "none")
		{
			this->solve_depressions(depression_solver, neighbourer, topography);

			topo_t faketopo(topography);
			this->enforce_minimal_slope_SS(1e-3,neighbourer,faketopo);
			this->compute_TO_SF_stack_version();
			// THERE IS A PROBLEM HERE



			this->receivers.clear();
			this->donors.clear();
			this->distance2receivers.clear();
			this->receivers = std::vector<std::vector<int> >(this->nnodes);
			this->donors = std::vector<std::vector<int> >(this->nnodes);
			this->distance2receivers = std::vector<std::vector<T> >(this->nnodes);
			neighbourer.rebuild_mgraph_after_solving_depression(Sdonors, Sreceivers, receivers,donors, Sdistance2receivers, distance2receivers, faketopo);
			this->compute_MF_topological_order();

		// std::cout << "ASDFdsfa->4" << std::endl;
		// std::cout << "ASDFdsfa->5" << std::endl;

		}
	}


	template<class Neighbourer_t, class topo_t>
	void compute_graph_both_v2(Neighbourer_t& neighbourer, topo_t& topography)
	{
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

	template<class Neighbourer_t, class topo_t>
	void solve_depressions(std::string& depression_solver , Neighbourer_t& neighbourer, topo_t& topography)
	{
		
		LMRerouter depsolver(neighbourer, topography, this->Sreceivers,  this->Sdonors,  this->Sdistance2receivers, this->stack);


		if(depsolver.npits > 0)
			depsolver.update_receivers(depression_solver,neighbourer, topography, this->Sreceivers,  this->Sdonors,  this->Sdistance2receivers, this->stack);
	}


	
	/// this function enforces minimal slope 
	template<class Neighbourer_t, class topo_t>
	void enforce_minimal_slope_SS(T slope, Neighbourer_t& neighbourer, topo_t& topography)
	{
		for(auto node:this->stack)
		{
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
		this->stack.clear();
		// reserving the amount of stuff
		this->stack.reserve(this->nnodes_t);

		// The stack container helper
		std::stack<int, std::vector<int> > stackhelper;
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
				this->stack.emplace_back(nextnode);

				// as well as all its donors which will be processed next
				for( int j = 0; j < this->Sdonors[nextnode].size(); ++j)
				{
					stackhelper.emplace(this->Sdonors[nextnode][j]);
				}

			}

		}

		if(this->stack.size() != this->nnodes_t)
		{
			std::cout << "Stack error, "<< this->stack.size() << "|" << this->nnodes << " checking for nans..." << std::endl;
			throw std::runtime_error("Stack error: should be " + std::to_string(this->nnodes_t) + " but is " + std::to_string(this->stack.size()));
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





	template<class Neighbourer_t, class topo_t>
	topo_t get_DA_proposlope(Neighbourer_t& neighbourer, topo_t& topography)
	{
		topo_t DA(neighbourer.nnodes,0.);
		for(int i = neighbourer.nnodes - 1; i>=0; --i)
		{
			int node = this->stack[i];
			DA[node] += neighbourer.cellarea;

			if(neighbourer.is_active(node))
			{
				topo_t slopes(this->receivers[node].size());
				float sumslopes = 0;
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

		return DA;
	}

	template<class Neighbourer_t, class topo_t>
	topo_t get_DA_SS(Neighbourer_t& neighbourer, topo_t& topography)
	{
		topo_t DA(neighbourer.nnodes,0.);
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

		return DA;
	}

};






























#endif