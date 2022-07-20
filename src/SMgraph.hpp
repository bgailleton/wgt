//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef smgraph_HPP
#define smgraph_HPP


/*
SMGraph stands for Static Multiple flow Graph.
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

#include "neighbourer.hpp"
#include "numvec.hpp"





class SMgraph
{
public:


	int nnodes;
	int n_neighbours;
	std::vector<bool> isrec;
	std::vector<int> Sreceivers;
	std::vector<size_t> stack, Sstack;
	std::vector<std::vector<int> > Sdonors;
	std::vector<double> Sdistance2receivers;


	SMgraph(){};
	SMgraph(int nnodes, int n_neighbours){this->nnodes = nnodes; this->n_neighbours = n_neighbours ;}

	template<class Neighbourer_t,class topo_t, class out_t>
	out_t compute_graph(std::string depression_solver, topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		
		ocarina timer;


		timer.tik();

		auto topography = format_input(ttopography);
		add_noise_to_vector(topography,-1e-6,1e-6);

		this->isrec = std::vector<bool>(int(this->nnodes * this->n_neighbours/2), false);
		this->Sreceivers = std::vector<int>(this->nnodes,-1);
		for(int i=0;i<this->nnodes; ++i)
			this->Sreceivers[i] = i;
		this->Sdistance2receivers = std::vector<double >(this->nnodes,-1);
		
		std::vector<double> SS(this->nnodes,0.);
		timer.tok("Initiation");

		timer.tik();
		neighbourer.build_smgraph_SS_only_SS(topography, this->Sreceivers, this->Sdistance2receivers, SS);
		timer.tok("Build_SS");

		timer.tik();
		this->recompute_SF_donors_from_receivers();
		this->compute_TO_SF_stack_version();

		timer.tok("SS_stuff");

		timer.tik();
		this->solve_depressions( depression_solver, neighbourer, topography);
		timer.tok("Depression solver");

		timer.tik();

		std::vector<double> faketopo = to_vec(topography);

		this->compute_TO_SF_stack_version();

		if(depression_solver == "carve")
		{
			this->carve_topo(1e-3,neighbourer,faketopo);

		}
		else if (depression_solver == "fill")
			this->fill_topo(1e-3,neighbourer,faketopo);


		timer.tok("Process topo post depression");


		
		timer.tik();
		neighbourer.build_smgraph_only_MF(faketopo, this->isrec);



		this->compute_MF_topological_order_insort(faketopo);
		
		timer.tok("MF stuffs");


		return format_output(faketopo);
	}

	template<class Neighbourer_t, class topo_t, class out_t>
	out_t compute_graph_multi_filled( topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		
		auto topography = format_input(ttopography);
		ocarina timer;
		timer.tik();

		// auto topography = format_input(ttopography);
		add_noise_to_vector(topography,-1e-6,1e-6);

		this->isrec = std::vector<bool>(int(this->nnodes * this->n_neighbours/2), false);

		timer.tik();
		std::vector<double> faketopo = neighbourer.fill_barne_2014(topography);
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
		py::gil_scoped_release release;
		// std::cout << "made it here -2" << std::endl;
		auto topography = format_input(ttopography);
		// std::cout << "made it here -1,9" << std::endl;
		add_noise_to_vector(topography,-1e-6,1e-6);
		// std::cout << "made it here -1" << std::endl;

		this->isrec = std::vector<bool>(int(this->nnodes * this->n_neighbours/2), false);
		this->Sreceivers = std::vector<int>(this->nnodes,-1);
		for(int i=0;i<this->nnodes; ++i)
			this->Sreceivers[i] = i;
		this->Sdistance2receivers = std::vector<double >(this->nnodes,-1);
		
		std::vector<double> SS(this->nnodes,0.);
		neighbourer.build_smgraph_SS_only_SS_OMP(topography, this->Sreceivers, this->Sdistance2receivers, SS);
		this->recompute_SF_donors_from_receivers();
		this->compute_TO_SF_stack_version();

		this->solve_depressions( depression_solver, neighbourer, topography);
		std::vector<double> faketopo = to_vec(topography);

		this->compute_TO_SF_stack_version();

		if(depression_solver == "carve")
		{
			this->carve_topo(1e-3,neighbourer,faketopo);

		}
		else if (depression_solver == "fill")
			this->fill_topo(1e-3,neighbourer,faketopo);


		neighbourer.build_smgraph_only_MF_OMP(faketopo, this->isrec);
		
		this->compute_MF_topological_order_insort(faketopo);
		py::gil_scoped_acquire acquire;
		
		return format_output(faketopo);
	}

	template<class Neighbourer_t,class topo_t, class out_t>
	out_t update_graph(std::string depression_solver, topo_t& ttopography, Neighbourer_t& neighbourer)
	{
		auto topography = format_input(ttopography);

		for(int i=0;i<this->nnodes; ++i)
			this->Sreceivers[i] = i;

		std::vector<double> SS(this->nnodes,0.);

		neighbourer.build_smgraph_SS_only_SS(topography, this->Sreceivers, this->Sdistance2receivers, SS);
		this->recompute_SF_donors_from_receivers();
		
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


	template<class topo_t>
	void compute_MF_topological_order_insort(topo_t& ttopography)
	{
		auto topography = format_input(ttopography);


		auto yolo = sort_indexes(topography);
		this->stack = std::move(yolo);

	}

	void recompute_SF_donors_from_receivers()
	{
		// Initialising the graph dimesions for the donors
		// All of thenm have the graph dimension
		this->Sdonors = std::vector< std::vector<int> >(this->nnodes);
		for(auto&v:this->Sdonors)
			v.reserve(8);

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
	void solve_depressions(std::string& depression_solver , Neighbourer_t& neighbourer, topo_t& ttopography)
	{
		auto topography = format_input(ttopography);

		
		LMRerouter<int, double,  Neighbourer_t, topo_t> depsolver(neighbourer, topography, this->Sreceivers,  this->Sdistance2receivers, this->Sstack);


		if(depsolver.npits > 0)
		{
			depsolver.update_receivers(depression_solver,neighbourer, topography, this->Sreceivers,  this->Sdistance2receivers, this->Sstack);
			this->recompute_SF_donors_from_receivers();
		}
	}

	void compute_TO_SF_stack_version()
	{
		// Initialising the stack
		this->Sstack.clear();
		// reserving the amount of stuff
		this->Sstack.reserve(this->nnodes);

		// The stack container helper
		std::stack<size_t, std::vector<size_t> > stackhelper;
		// std::vector<bool> isdone(this->nnodes,false);
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
				this->Sstack.emplace_back(nextnode);

				// as well as all its donors which will be processed next
				for( int j = 0; j < this->Sdonors[nextnode].size(); ++j)
				{
					stackhelper.emplace(this->Sdonors[nextnode][j]);
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
				double d2rec = this->Sdistance2receivers[node];
				topography[rec] = topography[node] - slope;// * d2rec;
			}

		}
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
				double d2rec = this->Sdistance2receivers[node];
				topography[node] = topography[rec] + slope;// * d2rec;
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
				for(int j = 0;j < receivers.size(); ++j)
				{
					int rec = receivers[j];
					slopes[j] = (topography[node] - topography[rec])/neighbourer.dx;
					if(slopes[j] <= 0)
						slopes[j] = 1e-5;
					sumslopes += slopes[j];
				}

				for(int j = 0;j < receivers.size(); ++j)
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

			if(neighbourer.is_active(node))
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









};

























#endif
