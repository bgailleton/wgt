//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef graph_HPP
#define graph_HPP

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
class Graph
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
		// #-> number of nodes in x direction.
		int nx = 0;
		// #-> number of nodes in x direction.
		int ny = 0;
		// #-> length in the x direction.
		float dx;
		// #-> length in the y direction.
		float dy;
		float dxy;
		// #-> cell area
		float cellarea;
		float Xmin;
		float Ymin;
		float Xmax;
		float Ymax;

		// Admin variables
		int not_a_node;
		float NDV = -9999;

		std::vector<int> receivers;
		std::vector<std::vector<int> > donors;
		std::vector<float> distance2receivers;

		// #->stack: topological order from downstream to upstream direction
		std::vector<int> stack;



		// Bearer of values
		// #->boundary: index: node index, value: boundary value
		// #-->The possible values are:
		// #---> -1 = no in, no out, not process, internal or external
		// #--->  0 = in, no out, all fluxes leave the system, internal or external
		// #--->  1 = "normal" cell, in, out, internal
		// #--->  2 = external periodic boundary
		std::vector<int> boundary;


		// Helpers for neighbouring operations
		// -> neighbourer holdsthe indices to loop through for each boundary condition
		std::vector<std::vector<int> > neighbourer;
		// -> lengthener is the dx on each directions
		std::vector<float > lengthener;


		// the topography:
		std::vector<float> topography;


		// Coordinate stuff
		// Xs and Ys are vectors of nx and ny size converting row to y and col to X
		// Extents holds the cxmin,xmax,ymin,ymax (extent is an option in matplotlib imshow plots)
		std::vector<float> Xs,Ys, extents;

		// "Hydro" stuff
		// Drainage area
		std::vector<float> area;
		// liste of nodes indices for river sources
		std::vector<int> source_nodes;
		// List of river nodes
		std::vector<int> river_nodes;
		// stack (topological order), for rivers
		std::vector<int> river_stack;
		// converts node index to index in river node
		std::map<int,int> node2rnode;
		// converts node index to junction order
		std::map<int,int> junctions;
		// node index to flow distance
		std::vector<float> flowdistance;
		// I do not remember what it is Oo
		std::vector<int> sourceID;

		// DDB stuff
		// Is the node a drainage boundary
		std::vector<bool> isDD;
		// Label of the drainage basins:
		// In steepest descent, each node belongs to a single watershed
		std::vector<int> basin_labels;
		// number of basin labels
		int nbasins = 0;


		// DEBUG STUFF
		// temporary vector I am using for outputting floating point values
		std::vector<float> debugfloat;


		// Modelling stuff (now moved somewhere else)
		std::vector<Label> labels;


		



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
		Graph(){;};
		Graph(std::vector<float> mock)
		{
			// std::cout << "Yolo yolo beng beng from c++" << std::endl;
			this->topography.clear();
			this->topography = mock;
		};

		Graph(int nx, int ny, int nnodes, float dx, float dy, float xmin, float ymin, std::vector<float>& numtopo)
		{
			this->set_dimensions(nx,ny,nnodes,dx,dy,xmin,ymin);
			this->set_default_boundaries("4edges");
			this->topography = numtopo;
		}
	


		void ingest_topo_from_C(uintptr_t topo, int carray_size)
		{
			const float* ptr = reinterpret_cast<float*>(topo);
			this->topography = std::vector<float>(ptr,ptr+carray_size); 
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


		/// @Description: Initialising the "neighbourer", a data structure managing the iterations though the neighbours of a given node
		/// Function of boundary conditions, one can feed the neighbourer with an indice which returns the index to ADD to the node to get its neighbour.
		/// The neighbour order is (in table referential, be careful if Y axis is inverted) top-left, top, top-right, left, right, bottom-left, bottom, bottom-right
		/// @Authors: B.G.
		/// @Date: 2021
		void initialise_neighbourer()
		{
			float diag = std::sqrt(std::pow(dx,2) + std::pow(dy,2) );
			this->dxy = diag;

			this->lengthener = std::initializer_list<float>{diag,dy,diag,dx,dx,diag,dy,diag};
			this->neighbourer.clear();

			// these vectors are additioned to the node indice to test the neighbors
			this->neighbourer.emplace_back(std::initializer_list<int>{-this->nx - 1, - this->nx, - this->nx + 1, -1,1,this->nx - 1, this->nx, this->nx + 1 }); // internal node 0
			this->neighbourer.emplace_back(std::initializer_list<int>{(this->ny - 1) * this->nx - 1, (this->ny - 1) * this->nx, (this->ny - 1) * this->nx + 1, -1,1,this->nx - 1, this->nx, this->nx + 1 });// periodic_first_row 1
			this->neighbourer.emplace_back(std::initializer_list<int>{-this->nx - 1, - this->nx, - this->nx + 1, -1,1,- (this->ny - 1) * this->nx - 1, - (this->ny - 1) * this->nx, - (this->ny - 1) * this->nx + 1 });// periodic_last_row 2
			this->neighbourer.emplace_back(std::initializer_list<int>{- 1, - this->nx, - this->nx + 1, (this->nx - 1),1, 2 * this->nx - 1, this->nx, this->nx + 1 }); // periodic_first_col 3
			this->neighbourer.emplace_back(std::initializer_list<int>{-this->nx - 1, - this->nx, - 2 * this->nx + 1, -1,-this->nx + 1, this->nx - 1, this->nx, 1 }); // periodic last_col 4
			this->neighbourer.emplace_back(std::initializer_list<int>{this->not_a_node, this->not_a_node, this->not_a_node, -1, 1, this->nx - 1, this->nx, this->nx + 1 }); // normal_first_row 5
			this->neighbourer.emplace_back(std::initializer_list<int>{- this->nx - 1, - this->nx, - this->nx + 1, -1,1, this->not_a_node, this->not_a_node, this->not_a_node}); // normal_last_row 6
			this->neighbourer.emplace_back(std::initializer_list<int>{this->not_a_node, - this->nx, - this->nx + 1, this->not_a_node, 1, this->not_a_node,  this->nx, this->nx + 1 }); // normal_first_col 7
			this->neighbourer.emplace_back(std::initializer_list<int>{-this->nx - 1, - this->nx, this->not_a_node, -1, this->not_a_node, this->nx - 1, this->nx, this->not_a_node }); // normal_last_col 8
			this->neighbourer.emplace_back(std::initializer_list<int>{this->not_a_node, this->not_a_node, this->not_a_node, this->not_a_node, 1, this->not_a_node, this->nx, this->nx + 1 }); // normal_top_left 9
			this->neighbourer.emplace_back(std::initializer_list<int>{ this->not_a_node, this->not_a_node, this->not_a_node, -1, this->not_a_node, this->nx - 1, this->nx, this->not_a_node}); // normal_top_right 10
			this->neighbourer.emplace_back(std::initializer_list<int>{ this->not_a_node,- this->nx, - this->nx + 1,this->not_a_node, 1, this->not_a_node, this->not_a_node, this->not_a_node}); // normal_bottom_left 11
			this->neighbourer.emplace_back(std::initializer_list<int>{-this->nx - 1, - this->nx, this->not_a_node, -1, this->not_a_node, this->not_a_node, this->not_a_node, this->not_a_node}); // normal_bottom_right 12
			this->neighbourer.emplace_back(std::initializer_list<int>{this->ny * this->nx -1, (this->ny - 1) * this->nx, (this->ny - 1) * this->nx + 1, this->nx - 1, 1, 2 * this->nx - 1, this->nx, this->nx+1}); // top_left_periodic 13
			this->neighbourer.emplace_back(std::initializer_list<int>{(this->ny - 1) * this->nx - 1, (this->ny - 1) * this->nx, (this->ny - 2) * this->nx + 1, -1, - this->nx + 1,this->nx - 1, this->nx, 1 }); // top_right_periodic 14
			this->neighbourer.emplace_back(std::initializer_list<int>{- 1, - this->nx, - this->nx + 1, this->nx - 1,1,- (this->ny - 2) * this->nx - 1, - (this->ny - 1) * this->nx, - (this->ny - 1) * this->nx + 1 });// periodic_bottom_left 15
			this->neighbourer.emplace_back(std::initializer_list<int>{-this->nx - 1, - this->nx, - this->nx + 1, -1, 1 - this->nx + 1, - (this->ny - 1) * this->nx - 1, - (this->ny - 1) * this->nx, - (this->ny) * this->nx + 1 });// periodic_bottom_right 16
		
		}




		/// @description: Sets the dimension of the graph and initialise the neighbourer
		void set_dimensions(int nx, int ny, int nnodes, float dx, float dy, float xmin, float ymin)
		{
			this->nx = nx;
			this->ny = ny;
			this->nnodes = nnodes;
			this->nnodes_t = size_t(nnodes);
			this->dx = dx;
			this->dy = dy;
			this->cellarea = this->dx * this->dy;
			this->Xmin = xmin;
			this->Ymin = ymin;

			// Not a node is utilised to detect when a neighbouring operation returns not a node
			this->not_a_node = - nx * ny - 10;

			this->initialise_neighbourer();

			// Initialise coordinate stuff
			this->Xs = std::vector<float>(this->nx);
			this->Ys = std::vector<float>(this->ny);
			// this->extents = std::vector<float>(4,0);

			for(int i=0; i<this->nx; ++i)
				this->Xs[i] = this->Xmin + this->dx/2 + i * this->dx;
			for(int i=this->ny -1; i>=0; --i)
			{
				this->Ys[i] = this->Ymin + this->dy/2 + i * this->dy;
			}

			this->Xmax = this->Xs.back() + this->dx/2;
			this->Ymax = this->Ys.back() + this->dy/2;

			std::reverse(this->Ys.begin(), this->Ys.end());

			this->extents = {this->Xmin, this->Xmin +  (this->nx + 1) * this->dx, this->Ymin, this->Ymin + (this->ny + 1) * this->dy};



			// Initialising A
			this->area = std::vector<float>(this->nnodes);

			// this->save2file("topo", this->topography);
		}

		void set_default_boundaries(std::string bountype)
		{

			this->boundary = std::vector<int>(this->nnodes_t,1);

			if(bountype == "4edges")
			{
				for(size_t i = 0; i < this->nnodes_t; ++i)
				{
					if(this->is_on_dem_edge(i))
						this->boundary[i] = 0;
				}
			}
			else if(bountype == "periodic_EW")
			{
				for(int i = 0; i < this->nnodes; i++)
				{
					if(this->is_on_top_row(i) || this->is_on_bottom_row(i))
						this->boundary[i] = 0;
					else if(this->is_on_leftest_col(i) || this->is_on_rightest_col(i))
						this->boundary[i] = 2;
				}
			}
			else if(bountype == "periodic_NS")
			{
				for(int i = 0; i < this->nnodes; i++)
				{
					if(this->is_on_leftest_col(i) || this->is_on_rightest_col(i))
						this->boundary[i] = 0;
					else if(this->is_on_top_row(i) || this->is_on_bottom_row(i))
						this->boundary[i] = 2;
				}
			}
			else
			{
				throw std::runtime_error("invalid periodic boundaries");
			}
		}

		int get_boundary_at_node(int i){return this->boundary[i];}

		void remove_seas(float sealvl)
		{
			bool is_there_still_0s = false;
			float min_elev_data = std::numeric_limits<float>::max();
			int idf = -1;
			for(int i = 0; i<this->nnodes; ++i)
			{
				if(this->topography[i] < sealvl)
					this->boundary[i] = -1;

				else if(this->topography[i]< min_elev_data)
				{
					idf = i;
					min_elev_data = this->topography[i];
				}

				if(this->boundary[i] == 0)
					is_there_still_0s = true;
			}

			this->compute_graph("none");
			bool has0 = false;
			for(int i =0; i<this->nnodes; ++i)
			{
				if(this->boundary[i] == -1)
					continue;
				if(this->boundary[this->receivers[i]] == 0)
					has0 = true;
			}

			if(has0 == false)
				this->set_lowest_boundary_to_outlet();

			// Done with the check

			if(idf >= 0 && is_there_still_0s == false)
			{
				std::cout << "recasting outlet to " <<idf << std::endl;
				this->boundary[idf] = 0;
			}
			else if (is_there_still_0s == false)
				throw std::runtime_error("no outlet remaining. Is sea level too high?");

		}

		// Draft for a function using connected comonents to solve stuff
		// Check connected components here:
			// this->compute_graph("carve");
			// std::vector<int>  baslab(this->nnodes,-1);
			// int lab = -1;
			// for(int i = 0; i < this->nnodes; ++i)
			// {
				
			// 	int tnode = this->stack[i];
			// 	int trec = this->receivers[tnode];
				
			// 	if(this->boundary[tnode] == -1)
			// 		continue;

			// 	if(tnode == trec)
			// 		lab++;

			// 	baslab[tnode] = lab;	
			// }


			// std::vector<bool> has0(lab,false);

			// for(int i = 0; i < this->nnodes; ++i)
			// {
				
			// 	int tnode = this->stack[i];				
			// 	if(this->boundary[tnode] == -1 || baslab[tnode] == -1)
			// 		continue;

			// 	int lalab = baslab[tnode];
			// 	if(this->boundary[tnode] == 0)
			// 		has0[lalab] = true;
			// }

			// for(int b=0; b<has0.size(); ++b)
			// {
			// 	if(has0[b] == false)
			// 	{
			// 		for(int i = 0; i < this->nnodes; ++i)
			// 		{
			// 			int tnode = this->stack[i];
			// 			if (baslab[tnode] == b)
			// 			{
			// 				this->boundary[tnode] = 0;
			// 				break;
			// 			}
			// 		}
			// 	}
			// }

		void compute(std::string depression_solver)
		{
			this->compute_graph(depression_solver);
			this->calculate_area();
			this->d_sources(1e6);
			this->compute_river_nodes();
		}

		// This function creates the graph of nodes
		// It computes the neighbours of each node in a given direction
		// Direction is a string that can be either "donors", "receivers" or "both"
		void compute_graph_old(std::string depression_solver)
		{
			// const float* ptr = reinterpret_cast<float*>(ttopography);
			// std::vector<float>  topography; topography.reserve(this->nnodes_t);
			// for (int i=0; i< this->nnodes; ++i)
			//	topography.emplace_back(ptr[i]);
			// std::cout << "cpp1" << std::endl;
			// this->compute_graph_SF_both();
			this->compute_graph_SF_both();
			// std::cout << "cpp2" << std::endl;
			this->compute_topological_order();
			// std::cout << "cpp3" << std::endl;
			if(depression_solver != "none")
			{
				// std::cout << "cpp4" << std::endl;
				this->solve_depressions(depression_solver);
				// std::cout << "cpp5" << std::endl;
				this->compute_topological_order();
				// std::cout << "cpp6" << std::endl;
			}
		}
		void compute_graph(std::string depression_solver)
		{
			// const float* ptr = reinterpret_cast<float*>(ttopography);
			// std::vector<float>  topography; topography.reserve(this->nnodes_t);
			// for (int i=0; i< this->nnodes; ++i)
			//	topography.emplace_back(ptr[i]);
			// std::cout << "cpp1" << std::endl;
			// this->compute_graph_SF_both();
			this->compute_graph_SF_both_v2();
			// std::cout << "cpp2" << std::endl;
			this->compute_topological_order();
			// std::cout << "cpp3" << std::endl;
			if(depression_solver != "none")
			{
				// std::cout << "cpp4" << std::endl;
				this->solve_depressions(depression_solver);
				// std::cout << "cpp5" << std::endl;
				this->compute_topological_order();
				// std::cout << "cpp6" << std::endl;
			}
		}

		void compute_graph_SF_both()
		{
			// Initialising the graph dimesions for the receivers
			// All of thenm have the graph dimension
			this->receivers = std::vector<int>(this->nnodes_t, -1);
			this->distance2receivers = std::vector<float>(this->nnodes_t, -1);
			this->donors = std::vector<std::vector<int> >(this->nnodes_t);
			for (auto& v:this->donors)
				v.reserve(8);

			// Iterating through all the nodes and finding the max slope
			int n_donors_iintotal = 0;


			for(int i = 0; i < this->nnodes; ++i)
			{
				this->receivers[i] = i;
				if(this->can_flow_out_there(i) || this->can_flow_even_go_there(i) == false)
				{
					continue;
				}

				auto neighbours = this->get_neighbours(i);

				float this_elev = topography[i];
				float max_slope = -1;
				int id_max_slope = -1;
				float dist_max_slope = -1;
				for (auto& ne:neighbours)
				{
					if(this->can_flow_even_go_there(ne.node) == false)
						continue;

					if(topography[ne.node] < this_elev)
					{
						double this_slope = (this_elev - topography[ne.node])/ne.distance;
						if(this_slope > max_slope)
						{
							max_slope = this_slope;
							id_max_slope = ne.node;
							dist_max_slope = ne.distance;
						}
					}
				}

				if(id_max_slope == -1)
					continue;

				this->receivers[i] = id_max_slope;
				this->distance2receivers[i] = dist_max_slope;
				this->donors[id_max_slope].emplace_back(i);
			}

		}

		void compute_graph_SF_both_v2()
		{
			// Initialising the graph dimesions for the receivers
			// All of thenm have the graph dimension
			this->receivers.clear();
			this->distance2receivers.clear();
			// std::cout << " qwe1" << std::endl;
			this->receivers = std::vector<int>(this->nnodes_t, -1);
			// std::cout << " qwe2" << std::endl;
			this->distance2receivers = std::vector<float>(this->nnodes_t, -1);
			std::vector<float> SS(this->nnodes_t, 0.);
			// std::cout << " qwe3" << std::endl;

			for (int i=0; i < this->nnodes; ++i)
			{
				this->receivers[i] = i;
			}
			// std::cout << " qwe4" << std::endl;

			// Iterating through all the nodes and finding the max slope
			int n_donors_iintotal = 0;
			// bool switc = false;

			for(int row = 0; row < this->ny; ++row)
			{
				// switc = !switc;
				for(int col = 0; col < this->nx; ++col)
				{
					int i = row * this->nx + col;
					// cannot be a neighbour anyway, abort
					if(this->can_flow_even_go_there(i) == false)
					{
						continue;
					}

					if(col > 0 && row > 0 && col < this->nx -1 && row < this->ny - 1 )
					{
						this->check_neighbour_v22(SS,this->get_topleft_index(i), i, this->dxy);
						this->check_neighbour_v22(SS,this->get_top_index(i), i, this->dy);
						this->check_neighbour_v22(SS,this->get_topright_index(i), i, this->dxy);
						this->check_neighbour_v22(SS,this->get_right_index(i), i, this->dx);
					}
					else
					{
						this->check_neighbour_v22(SS,this->get_top_index(i), i, this->dy);
						this->check_neighbour_v22(SS,this->get_left_index(i), i, this->dx);
						this->check_neighbour_v22(SS,this->get_right_index(i), i, this->dx);
						this->check_neighbour_v22(SS,this->get_bottom_index(i), i, this->dy);
						this->check_neighbour_v22(SS,this->get_topright_index(i), i, this->dxy);
						this->check_neighbour_v22(SS,this->get_topleft_index(i), i, this->dxy);
						this->check_neighbour_v22(SS,this->get_bottomright_index(i), i, this->dxy);
						this->check_neighbour_v22(SS,this->get_bottomleft_index(i), i, this->dxy);
					}

				}
			}

			// std::cout << " qwe5" << std::endl;


			this->recompute_SF_donors_from_receivers();

			// std::cout << " qwe6" << std::endl;
		}

		void check_neighbour_v22(std::vector<float>& SS, int ne, int i, float dist)
		{
			if (ne < 0)
				return;

			if(this->can_flow_even_go_there(ne) == false)
				return;
			
			float this_elev = this->topography[i];

			if(topography[ne] < this_elev && this->can_flow_out_there(i) == false)
			{
				double this_slope = (this_elev - topography[ne])/dist;
				if(this_slope > SS[i])
				{
					SS[i] = this_slope;
					this->receivers[i] = ne;
					this->distance2receivers[i] = dist;
				}
			}
			else if (topography[ne] > this_elev && this->can_flow_out_there(ne) == false)
			{
				double this_slope = (topography[ne] - this_elev)/dist;
				if(this_slope > SS[ne])
				{
					SS[ne] = this_slope;
					this->receivers[ne] = i;
					this->distance2receivers[ne] = dist;
				}
			}

		}

		void check_neighbour_v2(std::vector<float>& SS, Neighbour<int,float> ne, int i)
		{
			if (ne.node < 0)
				return;

			if(this->can_flow_even_go_there(ne.node) == false)
				return;
			
			float this_elev = this->topography[i];

			if(topography[ne.node] < this_elev && this->can_flow_out_there(i) == false)
			{
				double this_slope = (this_elev - topography[ne.node])/ne.distance;
				if(this_slope > SS[i])
				{
					SS[i] = this_slope;
					this->receivers[i] = ne.node;
					this->distance2receivers[i] = ne.distance;
				}
			}
			else if (topography[ne.node] > this_elev && this->can_flow_out_there(ne.node) == false)
			{
				double this_slope = (topography[ne.node] - this_elev)/ne.distance;
				if(this_slope > SS[ne.node])
				{
					SS[ne.node] = this_slope;
					this->receivers[ne.node] = i;
					this->distance2receivers[ne.node] = ne.distance;
				}
			}

		}



		void recompute_SF_donors_from_receivers()
		{
			// Initialising the graph dimesions for the donors
			// All of thenm have the graph dimension
			// std::cout << " qwe41" << std::endl;
			this->donors = std::vector<std::vector<int> >(this->nnodes);
			// std::cout << " qwe43" << std::endl;

			for(int i=0; i < this->nnodes; ++i)
			{
				// SF so rid == i cause there is only 1 rec
				int trec = this->receivers[i];
				if(trec == i)
					continue;

				this->donors[trec].emplace_back(i);
			}

			// std::cout << "qwe4end" << std::endl;

		}

		void solve_depressions(std::string& depression_solver)
		{
			// std::cout << "cpp4.1" << std::endl;
			Cordonnier2019_v2<int,float> depsolver(*this);
			// Cordonnier2019<int,float> depsolver(*this, this->topography);
			// std::cout << "cpp4.2" << std::endl;

			depsolver.update_receivers(depression_solver);
			// depsolver.update_receivers(depression_solver, this->topography);
			// std::cout << "cpp4.3" << std::endl;
		}


		void calculate_area()
		{
			this->area = std::vector<float>(this->nnodes,0.);

			for (int i=this->nnodes-1; i>=0;--i)
			{
				int tnode = this->stack[i];
				this->area[tnode] += this->cellarea;
				int trec =this->receivers[tnode];
				if(trec != tnode)
				{
					this->area[trec] += this->area[tnode];
				}
			}
			// std::cout << "Max area is " << this->_max(this->area) << std::endl;
		}


		void diffuse_area(float coeff){this->area = On_gaussian_blur(coeff, this->area, this->nx, this->ny);}

		void calculate_discharge_uniprec(float prec)
		{
			this->area = std::vector<float>(this->nnodes,0.);

			for (int i=this->nnodes-1; i>=0;--i)
			{
				int tnode = this->stack[i];
				this->area[tnode] += this->cellarea * prec;
				int trec =this->receivers[tnode];
				if(trec != tnode)
				{
					this->area[trec] += this->area[tnode];
				}
			}
			// std::cout << "Max area is " << this->_max(this->area) << std::endl;
		}

		void calculate_discharge_prec(std::vector<float>& prec)
		{
			this->area = std::vector<float>(this->nnodes,0.);

			for (int i=this->nnodes-1; i>=0;--i)
			{
				int tnode = this->stack[i];
				this->area[tnode] += this->cellarea * prec[tnode];
				int trec =this->receivers[tnode];
				if(trec != tnode)
				{
					this->area[trec] += this->area[tnode];
				}
			}
			// std::cout << "Max area is " << this->_max(this->area) << std::endl;
		}

		std::vector<float> get_steepest_slope(bool negto0 = true)
		{

			if(this->receivers.size() == 0)
				throw std::runtime_error("Cannot calculate steepest slope if the graph is not computed");

			std::vector<float> out(this->nnodes_t, 0.);

			for(int i = 0; i < this->nnodes_t; ++i)
			{
				
				int rec = this->receivers[i];
				
				if (rec == i)
				{
					out[i] = 0;
					continue;
				}

				out[i] = ( this->topography[i] - this->topography[rec] ) / this->distance2receivers[i];

				if(out[i] < 0 && negto0)
					out[i] = 0;
			}

			return out;
		}


		void set_boundaries_to(float val)
		{
			for(int i=0; i<this->nnodes; ++i)
			{
				if(this->can_flow_out_there(i))
					this->topography[i]= val;
			}
		}


		void set_boundary_at_rowcol(int row, int col, int val)
		{
			int i = this->nodeid_from_row_col(row,col);
			this->set_boundary_at(i,val);
		}

		void set_boundary_at(int i, int val)
		{
			this->boundary[i] = val;
		}

		void set_manual_boundaries(std::vector<int> bounds){this->boundary = bounds;}

		void set_single_outlet_remove_sea(float sealvl)
		{

			for(int i = 0; i<this->nnodes; ++i)
			{
				if(this->topography[i] < sealvl)
					this->boundary[i] = -1;
				else if (this->is_on_dem_edge(i))
					this->boundary[i] = 3;
				else
					this->boundary[i] = 1;
			}

			float minv = std::numeric_limits<float>::max();
			int minid = -9999;
			for(int i=0; i<this->nnodes; ++i)
			{
				if(this->boundary[i] == 1 || this->boundary[i] == 3)
				{
					if(this->topography[i] < minv)
					{
						minid = i;
						minv = this->topography[i];
					}
				}
			}

			if(minid == -9999)
				throw std::runtime_error("bite");

			int row,col;
			this->rowcol_from_node_id(minid,row,col);
			// std::cout << "recasting " << row << "|" << col << std::endl;
			this->boundary[minid] = 0;

		}

		void force_boundary_at_minimal_elevation()
		{
			float min_elev = std::numeric_limits<float>::max();
			int maxid = -9999;
			for(int i=0; i<this->nnodes; ++i)
			{
				if(this->can_flow_even_go_there(i) && this->can_flow_out_there(i) == false && min_elev > this->topography[i])
				{
					min_elev = this->topography[i];
					maxid = i;
				}
			}
			int row,col;
			this->rowcol_from_node_id(maxid,row,col);
			// std::cout << "Recasting " << maxid << " row: " << row << " | col: " << col << " to 0" << std::endl;

			this->boundary[maxid] = 0;
		}


		/// checks all the nodes and sets the lowest one as unique outlet
		void set_lowest_boundary_to_outlet()
		{
			// First trying to get the lowest point neighboured by nodata
			std::vector<bool> has_boundary = this->get_neighbour_to_nodata();
			float max_elev = std::numeric_limits<float>::min();
			int maxid = -9999;
			for(int i=0; i<this->nnodes; ++i)
			{
				if(has_boundary[i] && max_elev< this->topography[i])
				{
					max_elev = this->topography[i];
					maxid = i;
				}
			}


			// if it fails, for example my raster has no nodata, I am chacking the natural boundaries
			if(maxid == -9999)
			{
				for(int i=0; i<this->nnodes; ++i)
				{
					if(this->can_flow_out_there(i) && max_elev< this->topography[i])
					{
						max_elev = this->topography[i];
						maxid = i;
					}
				}
			}

			// if it fails here, there is a problem
			if(maxid == -9999)
			{
				throw std::runtime_error("No availble boundaries for minigraph::set_lowest_boundary_to_outlet. Mke sure there is at least one nodata bordering the");
			}

			
			// finally setting all external nodes to -1 ...
			for(int i=0; i<this->nnodes; ++i)
			{
				if(this->can_flow_even_go_there(i) == false || this->can_flow_out_there(i))
					this->boundary[i] = -1;
			}

			// Except the maxid one which will be the outlet
			this->boundary[maxid] = 0;

		}

		/// Returns a binary vector with node size where true is a data-carrying neighbour of nodata and false is not (nodata is false)
		std::vector<bool> get_neighbour_to_nodata()
		{
			std::vector<bool> has_boundary(this->nnodes,false);
			for(int i=0; i<this->nnodes; ++i)
			{
				if(this->can_flow_even_go_there(i) == false || this->can_flow_out_there(i))
					continue;

				auto neighbours = this->get_neighbours_only_id(i);
				for(auto tn:neighbours)
				{
					if(this->can_flow_even_go_there(tn) == false)
					{
						has_boundary[i] = true;
						break;
					}
				}
			}
			return has_boundary;
		}

		/// this function enforces minimal slope 
		void enforce_minimal_slope(float slope)
		{
			for(auto node:this->stack)
			{
				if(this->can_flow_out_there(node) || this->can_flow_even_go_there(node) == false)
					continue;
				int rec = this->receivers[node];
				float dz = this->topography[node] - this->topography[rec];

				if(dz <= 0)
				{
					float d2rec = this->distance2receivers[node];
					this->topography[node] = this->topography[rec] + slope * d2rec;
				}
			}
		}




		// ------------------------------------------------

		//	                             	            __
		//                                             / _)
		//                                    _.----._/ /
		//   GRAPH Conversion METHODS        /         /
		//                                __/ (  | (  |
		//                               /__.-'|_|--|_|

		// All the methods related to geometrical conversions
		// ------------------------------------------------
		// WILL NEED SOME WORK HERE DEPENDING ON THE GRAPH ORIENTATION AND ALL

		inline void rowcol_from_node_id(int node_index, int& row, int& col)
		{
			col = node_index % this->nx;
			row = (int)std::floor(node_index/this->nx);
		}

		inline int nodeid_from_row_col(int row, int col)
		{
			return row * this->nx + col;
		}

		inline int nodeid_from_XY(float X, float Y)
		{
			int col = floor((X - this->Xmin)/this->dx);
			int row = this->ny - floor((Y - this->Ymin)/this->dy);
			return this->nodeid_from_row_col(row,col);
		}
		// void rowcol_from_XY()


		inline void XY_from_nodeid(int node_index, float& tX, float& tY)
		{
			int row,col;
			this->rowcol_from_node_id(node_index, row, col);
			this->XY_from_rowcol(row, col, tX, tY);
		}

		inline void XY_from_rowcol(int row, int col, float& tX, float& tY)
		{
			tX = this->Xs[col];
			tY = this->Ys[row];
		}

		inline float get_X_from_i(int i)
		{
			int col = i % this->nx;
			return this->Xs[col];
		}

		inline float get_Y_from_i(int i)
		{
			int row = (int)std::floor(i/this->nx);
			return this->Ys[row];
		}



		// ------------------------------------------------

		//	                             	            __
		//                                             / _)
		//                                    _.----._/ /
		//   Neighbouring METHODS            /         /
		//                                __/ (  | (  |
		//                               /__.-'|_|--|_|

		// All the methods related to accessing and calculating neighbours
		// ------------------------------------------------

		std::vector<Neighbour<int,float> > get_neighbours(int i, bool ignore_nodata = false)
		{

			// preformatting the output
			std::vector<Neighbour<int,float> > neighbours;

			// Reserving size depending on the flow topology
			neighbours.reserve(8);

			size_t id_neighbourer = this->_get_neighbourer_id(i);

			// Now I need to determine the index of the neighbourer vector
			// Which provides adder to determine the neighbours f(boundary_conditions)
			// internal node 0
			// periodic_first_row 1
			// periodic_last_row 2
			// periodic_first_col 3
			// periodic last_col 4
			// normal_first_row 5
			// normal_last_row 6
			// normal_first_col 7
			// normal_last_col 8
			// normal_top_left 9
			// normal_top_right 10
			// normal_bottom_left 11
			// normal_bottom_right 12
			// top_left_periodic 13
			// top_right_periodic 14
			// periodic_bottom_left 15
			// periodic_bottom_right 16

			// // internal node, so neighbourer is 0
			// if(this->boundary[i] == 1)
			// 	id_neighbourer = 0;
			// else
			// {
			// 	// Case top left corner
			// 	if(i == 0)
			// 	{
			// 		// Case Periodic
			// 		if(this->boundary[i] == 2)
			// 			id_neighbourer = 13;
			// 		// Case Noraml
			// 		else if ( this->boundary[i] == 0 || this->boundary[i] == -1 )
			// 			id_neighbourer = 9;
			// 	}
			// 	// Case top right corner 
			// 	else if (i == this->nx - 1)
			// 	{
			// 		// Case Periodic
			// 		if(this->boundary[i] == 2)
			// 			id_neighbourer = 14;
			// 		// Case Noraml
			// 		else if ( this->boundary[i] == 0 || this->boundary[i] == -1 )
			// 			id_neighbourer = 10;
			// 	}
			// 	// Case bottom left corner 
			// 	else if (i == (this->nnodes - this->nx) )
			// 	{
			// 		// Case Periodic
			// 		if(this->boundary[i] == 2)
			// 			id_neighbourer = 15;
			// 		// Case Noraml
			// 		else if ( this->boundary[i] == 0 || this->boundary[i] == -1 )
			// 			id_neighbourer = 11;
			// 	}
			// 	// Case bottom right corner 
			// 	else if (i == (this->nnodes - 1) )
			// 	{
			// 		// Case Periodic
			// 		if(this->boundary[i] == 2)
			// 			id_neighbourer = 16;
			// 		// Case Noraml
			// 		else if ( this->boundary[i] == 0 || this->boundary[i] == -1 )
			// 			id_neighbourer = 12;
			// 	}
			// 	// Cases first row (no corner)
			// 	else if(i < this->nx - 1)
			// 	{
			// 		// Case periodic
			// 		if(this->boundary[i] == 2)
			// 			id_neighbourer = 1;
			// 		// Case normal
			// 		else if (this->boundary[i] == 0 || this->boundary[i] == -1)
			// 			id_neighbourer = 5;
			// 	}
			// 	// Cases last row (no corner)
			// 	else if(i > (this->ny - 1) * this->nx)
			// 	{
			// 		// Case periodic
			// 		if(this->boundary[i] == 2)
			// 			id_neighbourer = 2;
			// 		// Case normal
			// 		else if (this->boundary[i] == 0 || this->boundary[i] == -1)
			// 			id_neighbourer = 6;
			// 	}
			// 	// Cases first col (no corner)
			// 	else if(i % this->nx == 0)
			// 	{
			// 		// Case periodic
			// 		if(this->boundary[i] == 2)
			// 			id_neighbourer = 3;
			// 		// Case normal
			// 		else if (this->boundary[i] == 0 || this->boundary[i] == -1)
			// 			id_neighbourer = 7;
			// 	}
			// 	// Cases last col (no corner)
			// 	else if(i % this->nx == this->nx - 1)
			// 	{
			// 		// Case periodic
			// 		if(this->boundary[i] == 2)
			// 			id_neighbourer = 4;
			// 		// Case normal
			// 		else if (this->boundary[i] == 0 || this->boundary[i] == -1)
			// 			id_neighbourer = 8;
			// 	}
			// }

			// Now that I have the id_neighbourer, I can set the neighbours in the vector
			// if(id_neighbourer == 0)
			//	this->set_neighbours_in_vector_test0(neighbours, i, id_neighbourer, ignore_nodata);
			// else
			// std::cout << "id_neighbourer::" << id_neighbourer << std::endl;
				this->set_neighbours_in_vector(neighbours, i, id_neighbourer, ignore_nodata);

			// and return it
			return neighbours;
		}


		// This function is the helper of the neighbouring function: it fills the neighbours vector with the actul values.
		void set_neighbours_in_vector(std::vector<Neighbour<int,float> >& neighbours, int& i, size_t& id_neighbourer, bool ignore_nodata)
		{
			// node index of the current neighbour
			int tn;

			// D4 or D8 I take it all	
			// If you wonder why I am not iterating with a vector and everything here, it is for the small perf gain of doing it this way
			tn = this->neighbourer[id_neighbourer][1];
			if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
				neighbours.emplace_back(Neighbour<int, float>(tn + i, this->lengthener[1]));
			tn = this->neighbourer[id_neighbourer][3];
			if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
				neighbours.emplace_back(Neighbour<int, float>(tn + i, this->lengthener[3]));
			tn = this->neighbourer[id_neighbourer][4];
			if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
				neighbours.emplace_back(Neighbour<int, float>(tn + i, this->lengthener[4]));
			tn = this->neighbourer[id_neighbourer][6];
			if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
				neighbours.emplace_back(Neighbour<int, float>(tn + i, this->lengthener[6]));
			tn = this->neighbourer[id_neighbourer][0];
			if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
				neighbours.emplace_back(Neighbour<int, float>(tn + i, this->lengthener[0]));
			tn = this->neighbourer[id_neighbourer][2];
			if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
				neighbours.emplace_back(Neighbour<int, float>(tn + i, this->lengthener[2]));
			tn = this->neighbourer[id_neighbourer][5];
			if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
				neighbours.emplace_back(Neighbour<int, float>(tn + i, this->lengthener[5]));
			tn = this->neighbourer[id_neighbourer][7];
			if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
				neighbours.emplace_back(Neighbour<int, float>(tn + i, this->lengthener[7]));				
		
		}

		inline size_t _get_neighbourer_id(int i)
		{

			size_t id_neighbourer = -1;
			// internal node, so neighbourer is 0
			if(this->boundary[i] == 1)
				id_neighbourer = 0;
			else
			{
				// Case top left corner
				if(i == 0)
				{
					// Case Periodic
					if(this->boundary[i] == 2)
						id_neighbourer = 13;
					// Case Noraml
					else if ( this->boundary[i] == 0 || this->boundary[i] == -1 || this->boundary[i] == 3 )
						id_neighbourer = 9;
				}
				// Case top right corner 
				else if (i == this->nx - 1)
				{
					// Case Periodic
					if(this->boundary[i] == 2)
						id_neighbourer = 14;
					// Case Noraml
					else if ( this->boundary[i] == 0 || this->boundary[i] == -1 || this->boundary[i] == 3 )
						id_neighbourer = 10;
				}
				// Case bottom left corner 
				else if (i == (this->nnodes - this->nx) )
				{
					// Case Periodic
					if(this->boundary[i] == 2)
						id_neighbourer = 15;
					// Case Noraml
					else if ( this->boundary[i] == 0 || this->boundary[i] == -1 || this->boundary[i] == 3 )
						id_neighbourer = 11;
				}
				// Case bottom right corner 
				else if (i == (this->nnodes - 1) )
				{
					// Case Periodic
					if(this->boundary[i] == 2)
						id_neighbourer = 16;
					// Case Noraml
					else if ( this->boundary[i] == 0 || this->boundary[i] == -1 || this->boundary[i] == 3 )
						id_neighbourer = 12;
				}
				// Cases first row (no corner)
				else if(i < this->nx - 1)
				{
					// Case periodic
					if(this->boundary[i] == 2)
						id_neighbourer = 1;
					// Case normal
					else if (this->boundary[i] == 0 || this->boundary[i] == -1 || this->boundary[i] == 3)
						id_neighbourer = 5;
				}
				// Cases last row (no corner)
				else if(i > (this->ny - 1) * this->nx)
				{
					// Case periodic
					if(this->boundary[i] == 2)
						id_neighbourer = 2;
					// Case normal
					else if (this->boundary[i] == 0 || this->boundary[i] == -1 || this->boundary[i] == 3)
						id_neighbourer = 6;
				}
				// Cases first col (no corner)
				else if(i % this->nx == 0)
				{
					// Case periodic
					if(this->boundary[i] == 2)
						id_neighbourer = 3;
					// Case normal
					else if (this->boundary[i] == 0 || this->boundary[i] == -1 || this->boundary[i] == 3)
						id_neighbourer = 7;
				}
				// Cases last col (no corner)
				else if(i % this->nx == this->nx - 1)
				{
					// Case periodic
					if(this->boundary[i] == 2)
						id_neighbourer = 4;
					// Case normal
					else if (this->boundary[i] == 0 || this->boundary[i] == -1 || this->boundary[i] == 3)
						id_neighbourer = 8;
				}
				else
					id_neighbourer = 0;
			}

			if(id_neighbourer == -1)
			{
				int row,col;
				this->rowcol_from_node_id(i,row,col);
				std::cout << "boundary was " << this->boundary[i] << " i: " << i << "/" << this->nnodes << " row: " << row << " col: " << col << std::endl;
				throw std::runtime_error("neighbouring issue");
			}

			return id_neighbourer;
		}

		// D4 single neighbour routines
		Neighbour<int,float> get_left_neighbour(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return Neighbour<int,float>(this->neighbourer[id_neighbourer][3] + i, this->dx);
		}

		Neighbour<int,float> get_top_neighbour(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return Neighbour<int,float>(this->neighbourer[id_neighbourer][1] + i, this->dy);
		}

		Neighbour<int,float> get_right_neighbour(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return Neighbour<int,float>(this->neighbourer[id_neighbourer][4] + i, this->dx);
		}

		Neighbour<int,float> get_bottom_neighbour(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return Neighbour<int,float>(this->neighbourer[id_neighbourer][6] + i, this->dy);
		}

		// D8 single neighbour extra routines
		Neighbour<int,float> get_topleft_neighbour(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return Neighbour<int,float>(this->neighbourer[id_neighbourer][0] + i, this->dxy);
		}

		Neighbour<int,float> get_topright_neighbour(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return Neighbour<int,float>(this->neighbourer[id_neighbourer][2] + i, this->dxy);
		}

		Neighbour<int,float> get_bottomleft_neighbour(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return Neighbour<int,float>(this->neighbourer[id_neighbourer][5] + i, this->dxy);
		}

		Neighbour<int,float> get_bottomright_neighbour(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return Neighbour<int,float>(this->neighbourer[id_neighbourer][7] + i, this->dxy);
		}

		// D4 single neighbour routines
		int get_left_index(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return this->neighbourer[id_neighbourer][3] + i;
		}

		int get_top_index(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return this->neighbourer[id_neighbourer][1] + i;
		}

		int get_right_index(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return this->neighbourer[id_neighbourer][4] + i;
		}

		int get_bottom_index(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return this->neighbourer[id_neighbourer][6] + i;
		}

		// D8 single neighbour extra routines
		int get_topleft_index(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return this->neighbourer[id_neighbourer][0] + i;
		}

		int get_topright_index(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return this->neighbourer[id_neighbourer][2] + i;
		}

		int get_bottomleft_index(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return this->neighbourer[id_neighbourer][5] + i;
		}

		int get_bottomright_index(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return this->neighbourer[id_neighbourer][7] + i;
		}

		std::vector<int> get_D4_neighbours_only_id(int i)
		{
			std::vector<int> neighs;neighs.reserve(4);
			int tn = this->get_left_index(i);
			if(tn >= 0 )
				neighs.emplace_back(tn);
			tn = this->get_top_index(i);
			if(tn >= 0 )
				neighs.emplace_back(tn);
			tn = this->get_right_index(i);
			if(tn >= 0 )
				neighs.emplace_back(tn);
			tn = this->get_bottom_index(i);
			if(tn >= 0 )
				neighs.emplace_back(tn);
			return neighs;
		}

		std::vector<int> get_neighbours_only_id(int i)
		{
			std::vector<int> neighs;neighs.reserve(8);
			int tn = this->get_left_index(i);
			if(tn >= 0 )
				neighs.emplace_back(tn);
			tn = this->get_top_index(i);
			if(tn >= 0 )
				neighs.emplace_back(tn);
			tn = this->get_right_index(i);
			if(tn >= 0 )
				neighs.emplace_back(tn);
			tn = this->get_bottom_index(i);
			if(tn >= 0 )
				neighs.emplace_back(tn);
			tn = this->get_topleft_index(i);
			if(tn >= 0 )
				neighs.emplace_back(tn);
			tn = this->get_topright_index(i);
			if(tn >= 0 )
				neighs.emplace_back(tn);
			tn = this->get_bottomright_index(i);
			if(tn >= 0 )
				neighs.emplace_back(tn);
			tn = this->get_bottomleft_index(i);
			if(tn >= 0 )
				neighs.emplace_back(tn);


			return neighs;
		}





		// Method to test whether a node can outlet flow OUT of the model
		inline bool can_flow_out_there(int i)
		{
			if(this->boundary[i] == 0)
				return true;
			else
			{
				return false;
			}
		}

		// method to check if a node can even accept flow
		inline bool can_flow_even_go_there(int i)
		{
			if(this->boundary[i] < 0)
				return false;
			else
			{
				return true;
			}

		}

		// Returns true if the node is on the dem edge
		// (i.e. one of the 4 lines)
		bool is_on_dem_edge(int i)
		{
			if( this->is_on_top_row(i) || this->is_on_bottom_row(i) || this->is_on_leftest_col(i) || this->is_on_rightest_col(i) )
				return true;
			return false;
		}

		bool is_on_top_row(int i)
		{
			if(i < this->nx)
				return true;
			return false;
		}

		bool is_on_bottom_row(int i)
		{
			if(i >= this->nnodes - this->nx)
				return true;
			else
				return false;
		}

		bool is_on_leftest_col(int i)
		{
			if(i % this->nx == 0)
				return true;
			else
				return false;
		}

		bool is_on_rightest_col(int i)
		{
			if(i % this->nx == this->nx - 1)
				return true;
			else
				return false;
		}

		// Return an array of receivers of a given node
		// Return array in a DirectedNeighbour structure: node ID receiver ID and distance
		Neighbour<int,float> get_receivers(int i)
		{
			return Neighbour<int,float>(this->receivers[i],this->distance2receivers[i]);
		}

		// Return an array of donors of a given node
		// Return array in a DirectedNeighbour structure: node ID donor ID and distance
		std::vector<Neighbour<int,float> > get_donors(int i)
		{
			// Preprocessing the output to the right size (number of donors for the given node)
			std::vector<Neighbour<int,float> > output;
			output.reserve(this->donors[i].size());

			//Going through all the donor ID 
			for(int j = 0; j < this->donors[i].size(); ++j)
			{
				// And gathering the node ID,donor ID and distance to i of the donor node
				output.emplace_back(this->donors[i][j], this->distance2receivers[this->donors[i][j]]);
			}

			// AANd we are done
			return output;
		}


		std::vector<int> circular_window(int i, float radius, int reserve = 8)
		{
			// Typedef for the neighbour
			typedef Neighbour<int, float> stuff;

			// output
			std::vector<int> output;
			output.reserve(reserve);

			// Queue used to gather stuff
			std::queue<stuff> tQ;

			// The set to avoid doubling thingies
			std::set<int> done;

			// inserting th first element in the Q
			tQ.emplace(stuff(i,0));

			while(tQ.empty() == false)
			{
				// next node in the queue
				stuff next = tQ.front();
				tQ.pop(); // pop

				// std::cout << tQ.size() << std::endl;

				// Putting it in the output and the set
				done.insert(next.node);
				output.emplace_back(next.node);

				auto neighs = this->get_neighbours(next.node);
				for (auto& tn : neighs)
				{
					float tdist = next.distance + tn.distance;
					// std::cout << tdist << std::endl;

					if(tdist > radius)
						continue;

					if( izinset(done, tn.node) )
						continue;

					tQ.emplace(stuff(tn.node, tdist));

				}

			}

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

		void compute_topological_order()
		{
			this->compute_TO_SF_stack_version();
		}


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
				if(this->receivers[i] == i)
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
					for( int j = 0; j < this->donors[nextnode].size(); ++j)
					{
						stackhelper.emplace(this->donors[nextnode][j]);
					}

				}

			}

			if(this->stack.size() != this->nnodes_t)
			{
				std::cout << "Stack error, "<< this->stack.size() << "|" << this->topography.size() << " checking for nans..." << std::endl;
				for(auto v: this->topography)
				{
					if(std::isfinite(v) == false)
					{
						std::cout << "Yes, there are nans..." << std::endl;
					}
				}
				throw std::runtime_error("Stack error: should be " + std::to_string(this->nnodes_t) + " but is " + std::to_string(this->stack.size()));
			}

		}

		

		// ------------------------------------------------

		//	                             	            __
		//                                             / _)
		//                                    _.----._/ /
		//   DATA OUTPUT METHODS.            /         /
		//                                __/ (  | (  |
		//                               /__.-'|_|--|_|

		// All the methods related to accessing and calculating neighbours
		// ------------------------------------------------



		std::vector<int> get_topological_order_copy()
		{
			return this->stack;
		}

		// ------------------------------------------------

		//	                             	            __
		//                                             / _)
		//                                    _.----._/ /
		//   OTHER METHODS.            /         /
		//                                __/ (  | (  |
		//                               /__.-'|_|--|_|

		// Haphazard stuff for testing or really common routines
		// ------------------------------------------------

		std::vector<float> accumulate_downstream(float constant)
		{
			std::vector<float> output(this->nnodes,0.);
			for(int i = this->nnodes - 1; i >= 0; i--)
			{
				int tnode = this->stack[i];
				output[tnode] += constant;
				if(tnode != this->receivers[tnode])
					output[this->receivers[tnode]] += output[tnode];
			}
			return output;
		}

		std::vector<float> accumulate_downstream_var(float constant, std::vector<float>& weights)
		{
			std::vector<float> output(this->nnodes,0.);
			for(int i = this->nnodes - 1; i >= 0; i--)
			{
				int tnode = this->stack[i];
				output[tnode] += constant * weights[tnode];
				if(tnode != this->receivers[tnode])
					output[this->receivers[tnode]] += output[tnode];
			}
			return output;
		}

		std::vector<float> get_topo_array()
		{

			return this->topography;
		}

		void sources_from_Acrit(std::vector<int> labarray)
		{
			this->source_nodes.clear();
			this->river_nodes.clear();

			std::vector<bool> done(this->nnodes,false);
			for(int i = this->nnodes -1; i>=0; --i)
			{
				int tnode = this->stack[i];
				if(done[tnode])
				{
					int trec = this->receivers[tnode];
					if(trec != tnode)
						done[trec] = true;
					continue;
				}

				if(this->area[tnode] >= this->labels[labarray[tnode]].Acrit)
				{
					this->source_nodes.emplace_back(tnode);
					int trec = this->receivers[tnode];
					done[trec] = true;
					done[tnode] = true;
				}
			}
		}

		void d_sources(float Ath)
		{
			// std::cout << "Ath is " << Ath << std::endl;

			this->source_nodes.clear();
			this->river_nodes.clear();

			std::vector<bool> done(this->nnodes,false);
			for(int i = this->nnodes -1; i>=0; --i)
			{
				int tnode = this->stack[i];
				if(done[tnode])
				{
					int trec = this->receivers[tnode];
					if(trec != tnode)
						done[trec] = true;
					continue;
				}

				if(this->area[tnode] >= Ath)
				{
					this->source_nodes.emplace_back(tnode);
					int trec = this->receivers[tnode];
					done[trec] = true;
					done[tnode] = true;
				}
			}
		}


		void compute_river_nodes()
		{
			this->river_nodes.clear();
			this->node2rnode.clear();
			this->sourceID.clear();
			std::vector<bool> done(this->nnodes,false);
			int idex = 0;
			int sourcedex = -1;

			for(auto ts:this->source_nodes)
			{
				++sourcedex;
				int tnode = ts;
				int trec = this->receivers[tnode];
				// std::cout << this->can_flow_even_go_there(tnode) << "::::" << std::endl;
				while(this->can_flow_even_go_there(tnode) && done[tnode] == false)
				{
					// std::cout << "emplacing " << tnode << std::endl;
					this->river_nodes.emplace_back(tnode);
					this->sourceID.emplace_back(sourcedex);
					this->node2rnode[tnode] = idex;
					++idex;
					done[tnode] = true;
					tnode = trec;
					trec = this->receivers[tnode];
					if(this->can_flow_out_there(tnode) && done[tnode] == false)
					{
						this->river_nodes.emplace_back(tnode);
						this->sourceID.emplace_back(sourcedex);
						this->node2rnode[tnode] = idex;
						++idex;
						break;
					}
				}
			}
			// std::cout << "Got " << this->river_nodes.size() << " river pixels, I had " << this->source_nodes.size() << " sources" << std::endl;;


			// creating the river stack

			done = std::vector<bool>(this->nnodes,false);
			for(auto tn : this->river_nodes)
				done[tn] = true;

			this->river_stack.clear();this->river_stack.reserve(this->river_nodes.size());
			for(auto tn:this->stack)
			{
				if(done[tn])
					this->river_stack.emplace_back(tn);
			}

			this->compute_flow_distance_river();
		}


		void compute_flow_distance_river()
		{
			this->flowdistance = std::vector<float>(this->river_nodes.size(),0.);
			for(int i = 0; i < this->river_stack.size(); ++i)
			{
				int tn = this->river_stack[i]; 
				if(this->can_flow_out_there(tn) || this->can_flow_even_go_there(tn) == false)
					continue;
				int rec = this->receivers[tn]; int rrec = this->node2rnode[rec];
				this->flowdistance[this->node2rnode[tn]] = this->flowdistance[rrec] + this->distance2receivers[tn];
			}
		}

		std::vector<float> compute_flow_distance_full()
		{
			std::vector<float>TFD(this->nnodes,0.);
			for(auto tn:this->stack)
			{
				if(this->can_flow_out_there(tn) || this->can_flow_even_go_there(tn) == false)
					continue;

				int rec = this->receivers[tn];

				TFD[tn] = TFD[rec] + this->distance2receivers[tn];
			}
			return TFD;
		}

		std::vector<float> compute_chi_full(float theta, float A0)
		{
			std::vector<float>chi(this->nnodes,0.);
			for(auto tn:this->stack)
			{
				if(this->can_flow_out_there(tn) || this->can_flow_even_go_there(tn) == false)
					continue;

				int rec = this->receivers[tn];

				chi[tn] = chi[rec] + this->distance2receivers[tn] * A0/std::pow(this->area[tn], theta);
			}
			return chi;
		}



		int getNriv(){return int(this->river_nodes.size());}


		std::vector<float> get_normalised_flow_distance_by_basin()
		{
			std::vector<float>TFD(this->nnodes,0.);
			std::vector<int> max_FD_by_basins(this->nbasins,this->dx);
			for(auto tn:this->stack)
			{
				if(this->can_flow_out_there(tn) || this->can_flow_even_go_there(tn) == false)
					continue;

				int rec = this->receivers[tn];

				TFD[tn] = TFD[rec] + this->distance2receivers[tn];
				if(TFD[tn] > max_FD_by_basins[this->basin_labels[tn]]) max_FD_by_basins[this->basin_labels[tn]] = TFD[tn];
			}

			for(auto tn:this->stack)
				TFD[tn] = std::fmin(TFD[tn]/max_FD_by_basins[this->basin_labels[tn]],1.);

			return TFD;
		}

		std::vector<float> get_normalised_chi_by_basin(float theta, float A0)
		{
			std::vector<float>chi(this->nnodes,0.);
			std::vector<int> max_FD_by_basins(this->nbasins,0);

			for(auto tn:this->stack)
			{
				if(this->can_flow_out_there(tn) || this->can_flow_even_go_there(tn) == false)
					continue;

				int rec = this->receivers[tn];

				chi[tn] = chi[rec] + this->distance2receivers[tn] * std::pow(A0/this->area[tn], theta);
				if(chi[tn] > max_FD_by_basins[this->basin_labels[tn]]) max_FD_by_basins[this->basin_labels[tn]] = chi[tn];
			}

			for(auto tn:this->stack)
				chi[tn] = std::fmin(chi[tn]/max_FD_by_basins[this->basin_labels[tn]],1.);

			return chi;
		}

		





// Landscape creation routines
		void generate_white_noise(float intensity)
		{
			this->topography = std::vector<float>(this->nnodes,0);
			std::srand(std::time(0));
			for(int i = 0; i < this->nnodes; ++i)
			{
				if(this->boundary[i] <=0 )
					this->topography[i] = 0;
				else
				{
					this->topography[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/intensity));
				}
			}
		}

		void generate_white_noise_with_centered_divide_PEW(float intensity, float rel)
		{
			this->topography = std::vector<float>(this->nnodes,0);
			for(int i = 0; i < this->nnodes; ++i)
			{
				if(this->boundary[i] <=0 )
					this->topography[i] = 0;
				else
				{
					this->topography[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/intensity));
				}
			}

			int nyd2 = std::round(this->ny/2);
			for(int i=0; i <= nyd2; ++i)
			{
				float tz = intensity * rel / nyd2 * i;
				for(int j = 0; j<this->nx; ++j)
				{
					int node = i * this->nx + j;
					this->topography[node] += tz;
					node = this->nx * (this->nx - i - 1) + j;
					this->topography[node] += tz;
				}

			}
		}

		void add_white_noise(float intensity)
		{
			this->topography = std::vector<float>(this->nnodes,0);
			for(int i = 0; i < this->nnodes; ++i)
			{
				if(this->boundary[i] > 0)
				{
					this->topography[i] += static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/intensity));
				}
				
			}
		}



		std::vector<float> compute_graph_distance_from_rivers()
		{
			std::vector<float> output(this->nnodes, 0);
			std::vector<bool> done(this->nnodes, false);
			for(auto node:this->river_nodes)
			{
				done[node] = true;
			}

			for(auto node:this->stack)
			{
				if(this->can_flow_out_there(node) || this->can_flow_even_go_there(node) == false || done[node])
					continue;
				
				output[node] = output[this->receivers[node]] + this->distance2receivers[node];
			}

			return output;

		}

		template<class T>
		std::vector<T> propagate_from_rivers(std::vector<T>& input)
		{
			std::vector<T> output(input);
			std::vector<bool> done(this->nnodes, false);
			for(auto node:this->river_nodes)
			{
				done[node] = true;
			}

			for(auto node:this->stack)
			{
				if(this->can_flow_out_there(node) || this->can_flow_even_go_there(node) == false || done[node])
					continue;
				
				output[node] = input[this->receivers[node]];
			}
			return output;
		}



		std::vector<float> get_maxA_per_basins()
		{
			std::vector<float> out(this->nbasins,0);
			for(int i = 0; i <this->nnodes; ++i)
			{
				// if(this->can_flow_out_there(i))
				// {
				if(this->area[i] > out[this->basin_labels[i]])
					out[this->basin_labels[i]] = this->area[i];
				// }
			}
			return out;
		}


		// get_normalised_flow_distance_by_basin
		void generate_random_landscape_exp(std::vector<int> theselabels, std::vector<float> Es, float FD_limit, bool simple, float sigma_gblur, bool chi, float min_A, bool floodplains)
		{
			// Graph
			// std::cout << "Yolo>0" << std::endl;
			this->compute_graph("carve");
			// std::cout << "Yolo>0.1" << std::endl;
			// Drainage area
			this->calculate_area();
			// std::cout << "Yolo>0.2" << std::endl;

			this->compute_all_basins();
			// std::cout << "Yolo>0.3" << std::endl;
			auto fd = (chi) ? this->get_normalised_chi_by_basin(0.45,1): this->get_normalised_flow_distance_by_basin();
			// Few temp variables
			float maxZ = std::numeric_limits<float>::min();
			// std::cout << "Yolo>0.4" << std::endl;

			std::vector<float> maxAper_basins = this->get_maxA_per_basins();

			std::vector<float> distfromriv(this->nnodes,0);

			for(auto node:this->stack)
			{
				Label& tlab = this->labels[theselabels[node]];
				if(this->area[node] < tlab.Acrit)
					distfromriv[node] = distfromriv[this->receivers[node]] + this->distance2receivers[node];
			}
			// std::cout << "Yolo>0.5" << std::endl;

			// 5th step: build landscape
			// std::cout << "Yolo>1" << std::endl;
			for(auto node:this->stack)
			{
				// std::cout << "Yolo>" << node << std::endl;

				int rec = this->receivers[node];
				int baslab = this->basin_labels[node];
				if(rec == node || (fd[node] > FD_limit && maxAper_basins[baslab] >= min_A ) )
				{
					this->topography[node] = 0;
					continue;
				}

				Label& tlab = this->labels[theselabels[node]];

				if(simple)
				{
					this->topography[node] = this->topography[rec] + this->distance2receivers[node] * 
																	std::pow(Es[node] / ( tlab.Kref * std::pow(this->area[node],tlab.mref) ), 1/tlab.nref );
				}
				else if (tlab.breaks.size() == 0)
				{
					if(this->area[node] >= tlab.Acrit)
					{
						this->topography[node] = this->topography[rec] + this->distance2receivers[node] * 
																	std::pow(Es[node] / ( tlab.Kref * std::pow(this->area[node],tlab.mref) ), 1/tlab.nref );
					}
					else if (distfromriv[node] < 1000 && floodplains)
					{
						this->topography[node] = this->topography[rec] + 0.01 * this->distance2receivers[node];
					}
					else
					{
						this->topography[node] = this->topography[rec] + tlab.S_c * this->distance2receivers[node];
					}
				}
				else
				{
					float tA = this->area[node];
					if(tA >= tlab.Acrit)
					{
						int ib = -1;
						for(auto tl:tlab.breaks)
						{
							if(tA > tl[0])
								break;
							++ib;

						}

						if(ib == -1)
						{
							this->topography[node] = this->topography[rec] + this->distance2receivers[node] * 
																	std::pow(Es[node] / ( tlab.Kref * std::pow(this->area[node],tlab.mref) ), 1/tlab.nref );
						}
						else
						{
							this->topography[node] = this->topography[rec] + this->distance2receivers[node] * 
																	std::pow(Es[node] / ( tlab.breaks[ib][1] * std::pow(this->area[node],tlab.breaks[ib][2]) ), 1/tlab.breaks[ib][3] );
																	// std::cout << tlab.Kref << "o" << tlab.breaks[ib][1] << "|"<< tlab.mref << "o" << tlab.breaks[ib][2] << "|"<< tlab.nref << "o" << tlab.breaks[ib][3] << std::endl;
						}
						
					}
					else
					{
						this->topography[node] = this->topography[rec] + tlab.S_c * this->distance2receivers[node];
					}

				}

				if(this->topography[node] > maxZ)
					maxZ = 	this->topography[node];																
			}

			// std::cout << "Yolo>2" << std::endl;

			if(sigma_gblur > 0) this->fast_gblur(sigma_gblur);

			for(int i=0; i < this->nnodes; ++i)
			{
				if(fd[i] > FD_limit && maxAper_basins[this->basin_labels[i]] >= min_A )
					this->topography[i] = maxZ + 1 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
			}

			// std::cout << "Yolo>3" << std::endl;


		}


		// get_normalised_flow_distance_by_basin
		void generate_random_landscape_v4(std::vector<int> theselabels, std::vector<float> Es, float FD_limit, bool simple, float sigma_gblur, bool chi, float min_A)
		{
			// Graph
			// std::cout << "Yolo>0" << std::endl;
			this->compute_graph("carve");
			// std::cout << "Yolo>0.1" << std::endl;
			// Drainage area
			this->calculate_area();
			// std::cout << "Yolo>0.2" << std::endl;

			this->compute_all_basins();
			// std::cout << "Yolo>0.3" << std::endl;
			auto fd = (chi) ? this->get_normalised_chi_by_basin(0.45,1): this->get_normalised_flow_distance_by_basin();
			// Few temp variables
			float maxZ = std::numeric_limits<float>::min();
			// std::cout << "Yolo>0.4" << std::endl;

			std::vector<float> maxAper_basins = this->get_maxA_per_basins();
			// std::cout << "Yolo>0.5" << std::endl;

			// 5th step: build landscape
			// std::cout << "Yolo>1" << std::endl;
			for(auto node:this->stack)
			{
				// std::cout << "Yolo>" << node << std::endl;

				int rec = this->receivers[node];
				int baslab = this->basin_labels[node];
				if(rec == node || (fd[node] > FD_limit && maxAper_basins[baslab] >= min_A ) )
				{
					this->topography[node] = 0;
					continue;
				}

				Label& tlab = this->labels[theselabels[node]];

				if(simple)
				{
					this->topography[node] = this->topography[rec] + this->distance2receivers[node] * 
																	std::pow(Es[node] / ( tlab.Kref * std::pow(this->area[node],tlab.mref) ), 1/tlab.nref );
				}
				else if (tlab.breaks.size() == 0)
				{
					if(this->area[node] >= tlab.Acrit)
					{
						this->topography[node] = this->topography[rec] + this->distance2receivers[node] * 
																	std::pow(Es[node] / ( tlab.Kref * std::pow(this->area[node],tlab.mref) ), 1/tlab.nref );
					}
					else
					{
						this->topography[node] = this->topography[rec] + tlab.S_c * this->distance2receivers[node];
					}
				}
				else
				{
					float tA = this->area[node];
					if(tA >= tlab.Acrit)
					{
						int ib = -1;
						for(auto tl:tlab.breaks)
						{
							if(tA > tl[0])
								break;
							++ib;

						}

						if(ib == -1)
						{
							this->topography[node] = this->topography[rec] + this->distance2receivers[node] * 
																	std::pow(Es[node] / ( tlab.Kref * std::pow(this->area[node],tlab.mref) ), 1/tlab.nref );
						}
						else
						{
							this->topography[node] = this->topography[rec] + this->distance2receivers[node] * 
																	std::pow(Es[node] / ( tlab.breaks[ib][1] * std::pow(this->area[node],tlab.breaks[ib][2]) ), 1/tlab.breaks[ib][3] );
																	// std::cout << tlab.Kref << "o" << tlab.breaks[ib][1] << "|"<< tlab.mref << "o" << tlab.breaks[ib][2] << "|"<< tlab.nref << "o" << tlab.breaks[ib][3] << std::endl;
						}
						
					}
					else
					{
						this->topography[node] = this->topography[rec] + tlab.S_c * this->distance2receivers[node];
					}

				}

				if(this->topography[node] > maxZ)
					maxZ = 	this->topography[node];																
			}

			// std::cout << "Yolo>2" << std::endl;

			if(sigma_gblur > 0) this->fast_gblur(sigma_gblur);

			for(int i=0; i < this->nnodes; ++i)
			{
				if(fd[i] > FD_limit && maxAper_basins[this->basin_labels[i]] >= min_A )
					this->topography[i] = maxZ + 1 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
			}

			// std::cout << "Yolo>3" << std::endl;


		}


		void generate_random_landscape_v3(std::vector<int> theselabels, std::vector<float> Es, float FD_limit, bool simple)
		{
			// Graph
			this->compute_graph("fill");
			// Drainage area
			this->calculate_area();
			auto fd = this->compute_flow_distance_full();
			// Few temp variables
			float maxZ = std::numeric_limits<float>::min();
			float maxFD = this->_max(fd);
			FD_limit = maxFD *  FD_limit;

			// 5th step: build landscape
			for(auto node:this->stack)
			{
				int rec = this->receivers[node];
				if(rec == node || fd[node] > FD_limit)
				{
					this->topography[node] = 0;
					continue;
				}

				Label& tlab = this->labels[theselabels[node]];

				if(simple)
				{
					this->topography[node] = this->topography[rec] + this->distance2receivers[node] * 
																	std::pow(Es[node] / ( tlab.Kref * std::pow(this->area[node],tlab.mref) ), 1/tlab.nref );
				}
				else if (tlab.breaks.size() == 0)
				{
					if(this->area[node] >= tlab.Acrit)
					{
						this->topography[node] = this->topography[rec] + this->distance2receivers[node] * 
																	std::pow(Es[node] / ( tlab.Kref * std::pow(this->area[node],tlab.mref) ), 1/tlab.nref );
					}
					else
					{
						this->topography[node] = this->topography[rec] + tlab.S_c * this->distance2receivers[node];
					}
				}
				else
				{
					float tA = this->area[node];
					if(tA >= tlab.Acrit)
					{
						int ib = -1;
						for(auto tl:tlab.breaks)
						{
							if(tA > tl[0])
								break;
							++ib;

						}

						if(ib == -1)
						{
							this->topography[node] = this->topography[rec] + this->distance2receivers[node] * 
																	std::pow(Es[node] / ( tlab.Kref * std::pow(this->area[node],tlab.mref) ), 1/tlab.nref );
						}
						else
						{
							this->topography[node] = this->topography[rec] + this->distance2receivers[node] * 
																	std::pow(Es[node] / ( tlab.breaks[ib][1] * std::pow(this->area[node],tlab.breaks[ib][2]) ), 1/tlab.breaks[ib][3] );
																	// std::cout << tlab.Kref << "o" << tlab.breaks[ib][1] << "|"<< tlab.mref << "o" << tlab.breaks[ib][2] << "|"<< tlab.nref << "o" << tlab.breaks[ib][3] << std::endl;
						}
						
					}
					else
					{
						this->topography[node] = this->topography[rec] + tlab.S_c * this->distance2receivers[node];
					}

				}

				if(this->topography[node] > maxZ)
					maxZ = 	this->topography[node];																
			}

			for(int i=0; i < this->nnodes; ++i)
			{
				if(fd[i] > FD_limit)
					this->topography[i] = maxZ + 1 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
			}
		}



		void generate_random_landscape_v2_spatial(std::vector<float> Ks, float mexp, float nexp, std::vector<float> Es, float FD_limit)
		{
			// Graph
			this->compute_graph("carve");
			// Drainage area
			this->calculate_area();
			auto fd = this->compute_flow_distance_full();
			// Few temp variables
			float maxZ = std::numeric_limits<float>::min();
			float maxFD = this->_max(fd);
			FD_limit = maxFD *  FD_limit;

			// 5th step: build landscape
			for(auto node:this->stack)
			{
				int rec = this->receivers[node];
				if(rec == node || fd[node] > FD_limit)
				{
					this->topography[node] = 0;
					continue;
				}
				this->topography[node] = this->topography[rec] + this->distance2receivers[node] * 
																	std::pow(Es[node] / ( Ks[node] * std::pow(this->area[node],mexp) ), 1/nexp );
				if(this->topography[node] > maxZ)
					maxZ = 	this->topography[node];																
			}
			for(int i=0; i < this->nnodes; ++i)
			{
				if(fd[i] > FD_limit)
					this->topography[i] = maxZ + 1 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
			}
		}

		void generate_random_landscape_v2(float this_K, float mexp, float nexp, float E, float FD_limit)
		{
			// Graph
			this->compute_graph("carve");
			// Drainage area
			this->calculate_area();
			auto fd = this->compute_flow_distance_full();
			// Few temp variables
			float maxZ = std::numeric_limits<float>::min();
			float maxFD = this->_max(fd);
			FD_limit = maxFD *  FD_limit;

			// 5th step: build landscape
			for(auto node:this->stack)
			{
				int rec = this->receivers[node];
				if(rec == node || fd[node] > FD_limit)
				{
					this->topography[node] = 0;
					continue;
				}
				this->topography[node] = this->topography[rec] + this->distance2receivers[node] * 
																	std::pow(E / ( this_K * std::pow(this->area[node],mexp) ), 1/nexp );
				if(this->topography[node] > maxZ)
					maxZ = 	this->topography[node];																
			}
			for(int i=0; i < this->nnodes; ++i)
			{
				if(fd[i] > FD_limit)
					this->topography[i] = maxZ + 1 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
			}
		}

		void generate_random_landscape(float min_K, float max_K, float mexp, float nexp, float E)
		{
			// // 1st step: generate landscape // maybe better to left that separate
			// this->generate_white_noise(1);

			// 2nd step: compute graph, we'll fill the landscape to test it
			this->compute_graph("carve");

			// 3rd step: get A
			std::vector<float> A = this->accumulate_downstream(this->cellarea);

			// 4th step: functionalise K
			float tmax = std::log10(this->_max(A));
			float tmin = std::log10(this->cellarea);

			// 5th step: build landscape
			for(auto node:this->stack)
			{
				int rec = this->receivers[node];
				if(rec == node)
					continue;
				this->topography[node] = this->topography[rec] + this->distance2receivers[node] * 
																	std::pow(E / ( (std::log10(A[node]) * (max_K - min_K)/(tmax - tmin) + min_K) * 
																		std::pow(A[node],mexp) ), 1/nexp );

				if( ( (std::log10(A[node]) * (max_K - min_K)/(tmax - tmin) + min_K) * 
																		std::pow(A[node],mexp) ) < 0)
				{
					std::cout << std::log10(A[node]) << "|" << (max_K - min_K) << "|" << (tmax - tmin) + min_K << std::endl;
					return;					
				}

			}

		}


		void refine_topo(float Ath, float DDth, float FDth)
		{
			// marking junctions with a label
			std::vector<int> tbaslab(this->nnodes, -1);
			int tlab = -1;
			for(auto junc:this->junctions)
			{
				++tlab;
				tbaslab[junc.first] = tlab;
			}
			for(int i =0; i<this->nnodes;++i)
			{
				if(this->can_flow_out_there(i))
				{
					++tlab;
					tbaslab[i] = tlab;
				}
			}

			// Propagating the labels up
			for(auto tn:this->stack)
			{
				if(this->can_flow_out_there(tn) || this->can_flow_even_go_there(tn) == false || tbaslab[tn] != -1 )
				{
					continue;
				}
				tbaslab[tn] = tbaslab[this->receivers[tn]];
			}

			// marking up nodes of interest
			std::vector<bool> marked(this->nnodes,false);
			std::vector<float> FD = this->compute_flow_distance_full();
			std::priority_queue<PQ_helper<int, float>, std::vector<PQ_helper<int, float>>, std::greater<PQ_helper<int, float> > > PQ;

			for(int i=0; i<this->nnodes; ++i)
			{
				if(this->can_flow_out_there(i) || this->can_flow_even_go_there(i) == false || marked[i] || FD[i] < FDth)
					continue;

				int tlab = tbaslab[i];
				auto neighs = this->get_D4_neighbours_only_id(i);
				for(auto tn:neighs)
				{
					if(tn < 0)
						continue;
					if(marked[tn])
						continue;

					if(tlab != tbaslab[tn])
					{
						marked[i] = true;
						marked[tn] = true;
						PQ.emplace(PQ_helper<int,float>(i,0.));
						PQ.emplace(PQ_helper<int,float>(tn,0.));
					}
				}
			}

			std::vector<float> DDD(this->nnodes, 0.);
			while(PQ.empty() == false)
			{
				auto next = PQ.top(); PQ.pop();
				auto neighs = this->get_neighbours(next.node);
				for (auto tn : neighs)
				{
					if(tn.node<0)
						continue;
					if(marked[tn.node] && DDD[tn.node] == 0)
						continue;


					float tg_dist = DDD[next.node] + tn.distance;

					if(tg_dist < DDD[tn.node] || DDD[tn.node] == 0)
					{
						DDD[tn.node] = tg_dist;
						if(this->can_flow_out_there(tn.node) == false && this->can_flow_even_go_there(tn.node))
							PQ.emplace(PQ_helper<int,float>(tn.node,DDD[tn.node]));
					}
				}
				if(PQ.size()> 10 * this->nnodes)
					throw std::runtime_error("ERROR ^%&");
			}



			// std::vector<float> DDD(this->nnodes, 0.);
			// for(int i=this->nnodes-1; i>=0;--i)
			// {
			// 	int tn = this->stack[i];
			// 	int rec = this->receivers[tn];
			// 	if(this->can_flow_out_there(tn) || this->can_flow_even_go_there(tn) == false || marked[rec] || FD[tn] < FDth)
			// 		continue;

			// 	DDD[rec] = (DDD[rec] == 0)? DDD[tn] + this->distance2receivers[tn] : std::fmax(DDD[rec], DDD[tn] + this->distance2receivers[tn]);

			// }

			this->debugfloat = std::vector<float>(FD);

			float maxtopo = this->_max(this->topography);

			for(int i=0; i < this->nnodes; ++i)
			{
				if(DDD[i] < DDth && FD[i] > FDth)
					this->topography[i] = maxtopo + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX)) * 1 + 1;
			}

		}

		void refine_topo_v2(float Ath, float Rth)
		{
			auto itopo = this->_invert_topography(Ath);
			std::vector<bool> isdone(this->nnodes,false);
			std::queue<int> to_replace;
			if(Rth > 0)
				Rth = - Rth;

			for(int i=0; i <this->nnodes; ++i)
			{

				if(itopo[i] <= Rth)
					isdone[i] = true;
				else
					to_replace.emplace(i);
			}

			while(to_replace.empty() == false)
			{
				int next = to_replace.front();to_replace.pop();
				auto neighs = this->get_D4_neighbours_only_id(next);
				int maxN = -1;
				float val = std::numeric_limits<float>::max();
				for(auto tn:neighs)
				{
					if(tn >=0)
					{
						if(isdone[tn] && this->topography[tn] < val)
						{
							val = this->topography[tn];
							maxN = tn;
						}
					}
				}

				if(maxN > -1)
				{
					this->topography[next] = val + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX)) * 1e-3;
					isdone[next] = true;
				}
				else
				{
					to_replace.emplace(next);
				}
			}
		}

		void refine_topo_v1(float area_threshold)
		{
			// 1 recomputing the graph
			this->compute_graph("fill");

			// 2  Getting new drainage area
			std::vector<float> A = this->accumulate_downstream(this->cellarea);

			// 3 setting up a bool array as ref for rivers
			std::vector<bool> isriv(this->nnodes, false);
			for(int i=0; i<this->nnodes; ++i)
			{
				if(A[i]>=area_threshold)
					isriv[i] = true;
			}
			// Freeing up memory
			A.clear();

			// 4 Set up the distance2river
			std::vector<float> distance2river(this->nnodes, 0.);
			std::vector<float> Z0(this->nnodes, 0.);
			for(auto node:this->stack)
			{
				if(isriv[node] || this->can_flow_even_go_there(node) == false || this->can_flow_out_there(node))
				{
					Z0[node] = this->topography[node];
					continue;
				}
				int rec = this->receivers[node];
				distance2river[node] = distance2river[rec] + this->distance2receivers[node];
				Z0[node] = Z0[rec];
			}

			// 5 set up the max elevation
			std::vector<float> targetZ(this->nnodes, -1);
			std::vector<float> targetX(this->nnodes, -1);
			for(int i = this->nnodes - 1; i >= 0; --i)
			{
				int node = this->stack[i];
				int rec = this->receivers[node];
				if(targetZ[node] == -1)
				{
					targetZ[node] = this->topography[node];
					targetX[node] = distance2river[node];
				}
				else if(targetZ[rec] < targetZ[node])
				{
					targetZ[rec] = targetZ[node];
					targetX[rec] = targetZ[node];
				}
			}

			// OK now we have a field of X, X0 is 0, Z and a target Z
			// 6 Calculating topo by increment
			for(auto node:this->stack)
			{
				if(isriv[node] || this->can_flow_even_go_there(node) == false || this->can_flow_out_there(node))
					continue;
				int rec = this->receivers[node];
				this->topography[node] = this->topography[rec] + (targetZ[node] - Z0[node])/targetX[node] * this->distance2receivers[node];
			}


		}

		void replace_top_by_topo(float frac)
		{
			float maxtopo = this->_max(this->topography) * frac;
			for(auto& v:this->topography)
			{
				if(v > maxtopo)
					v = maxtopo + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/1));
			}

		}

		void replace_top_by_A(float frac)
		{
			std::vector<float> A = this->accumulate_downstream(this->cellarea);

			float maxA = this->_max(A) * frac;
			float maxtopo = this->_max(this->topography);
			
			for(int i=0; i<this->nnodes; ++i)
			{
				if(this->can_flow_out_there(i))
					continue;
				if(A[i] < maxA)
					this->topography[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/1)) + maxtopo ;
			}

		}

		void replace_cleverly(float min_K, float max_K, float mexp, float nexp, float E, float factor)
		{
			// std::vector<int> baslab(this->nnodes_t, -1);
			// int tlab = -1;
			// for(auto i:this->stack)
			// {
			// 	if(i == this->receivers[i])
			// 		tlab++;
			// 	baslab[i] = tlab;
			// }

			std::vector<bool> drains_S(this->nnodes_t,false);
			for(int i=0; i < this->nx; ++i)
				drains_S[i] = true;

			for(auto i:this->stack)
			{
				drains_S[i] = drains_S[this->receivers[i]];
			}
			
			std::vector<int> oldrecs(this->receivers);
			std::vector<int> oldstack(this->stack);

			std::priority_queue<PQ_helper<int, float>, std::vector<PQ_helper<int, float>>, std::greater<PQ_helper<int, float> > > Q;

			this->compute_graph("fill");
			float maxS = std::pow(E/(min_K * std::pow(this->cellarea, mexp)), 1/nexp) * factor;

			for(int i=0; i< this->nnodes; ++i)
			{
				int rec = this->receivers[i];

				if(rec == i)
					continue;

				float tS = (this->topography[i] - this->topography[rec])/ this->distance2receivers[i];

				if(tS > maxS && drains_S[rec] !=  drains_S[i])
				{
					// this->topography[i] = -9999;
					Q.emplace(PQ_helper<int, float>(i,this->topography[rec]));
				}

			}

			while(Q.empty() == false)
			{
				PQ_helper<int, float> next = Q.top(); Q.pop();
				if(this->topography[next.node] > next.score)
				{
					this->topography[next.node] = next.score;
					
					auto neighbours = this->get_neighbours(next.node);
					for(auto& tn:neighbours)
					{
						if(this->topography[tn.node] > next.score)
							Q.emplace(PQ_helper<int, float>(tn.node, next.score));
					}
				}
			}

			this->add_white_noise(1e-3);


			// test1
			// labelling S boundary
			// std::vector<bool> drains_S(this->nnodes_t,false);
			// for(int i=0; i < this->nx; ++i)
			// 	drains_S[i] = true;
			// float max_Z_S = 0;
			// float max_Z_N = 0;
			// for(auto i:this->stack)
			// {
			// 	drains_S[i] = drains_S[this->receivers[i]];
			// 	if(drains_S[i] && this->topography[i]>max_Z_S)
			// 		max_Z_S = this->topography[i];
			// 	else if(drains_S[i] == false && this->topography[i]>max_Z_N)
			// 		max_Z_N = this->topography[i];
			// }

			// bool tgb = (max_Z_N > max_Z_S) ? false : true;
			// float tgz = (max_Z_N > max_Z_S) ? max_Z_S : max_Z_N;
			// for(int i=0; i<this->nnodes; ++i)
			// {
			// 	if(this->topography[i] > tgz && drains_S[i] == tgb)
			// 		this->topography[i] = tgz + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/1)) ;
			// }
			
		}


		void calculate_junction_map()
		{
			this->junctions.clear();
			std::vector<int> tempjunc(this->river_nodes.size(),0);
			std::vector<int> junctionnode; junctionnode.reserve(this->river_nodes.size());

			for(int i=this->river_nodes.size() - 1; i>=0; --i)
			{
				int tn = this->river_stack[i];
				int rn = this->node2rnode[tn];
				// std::cout << "tn:" << tn << std::endl;
				if(tempjunc[rn] == 0)
					tempjunc[rn]++;

				int rec = this->receivers[tn];
				if(this->can_flow_out_there(tn) || this->can_flow_even_go_there(tn) == false)
					continue;

				if(this->node2rnode.count(rec) == 0)
					continue;

				int recriv = this->node2rnode[rec];

				if(rn < 0 || recriv < 0 || rn >= this->river_nodes.size() || recriv >= this->river_nodes.size())
					throw std::runtime_error("BIP BIP tn||rn "+ std::to_string(tn) + "||" + std::to_string(rn));
				// std::cout << recriv << "|";

				if(tempjunc[recriv] > 0)
					junctionnode.emplace_back(rec);

				if(tempjunc[recriv] == 0)
					tempjunc[recriv] = tempjunc[rn];
				else if (tempjunc[recriv] < tempjunc[rn])
					tempjunc[recriv] = tempjunc[rn];
				else if (tempjunc[recriv] == tempjunc[rn])
					tempjunc[recriv] = tempjunc[rn] + 1;
			}

			for(auto tn:junctionnode)
			{
				if(this->node2rnode[tn] < 0 || this->node2rnode[tn] >= this->river_nodes.size())
					throw std::runtime_error("YOLO");
				this->junctions[tn] = tempjunc[this->node2rnode[tn]];
			}

		}

		std::vector<int> _get_stream_order()
		{

			std::vector<int> SO(this->river_nodes.size(),0);
			for (int i = this->river_stack.size() - 1 ; i >= 0; --i)
			{
				int tn = this->river_stack[i];
				int rn = this->node2rnode[tn];
				if(this->junctions.count(tn)>0)
					SO[rn] = this->junctions[tn];
				else if(SO[rn] == 0)
				{
					SO[rn] = 1;
				}

				int rec = this->receivers[tn];
				if(this->node2rnode.count(rec) > 0)
				{
					int rrn = this->node2rnode[rec];
					SO[rrn] = SO[rn];
				}

			}

			return SO;

		}


		std::vector<float> _invert_topography(float Ath)
		{
			// marking junctions with a label
			std::vector<int> tbaslab(this->nnodes, -1);
			int tlab = -1;
			for(auto junc:this->junctions)
			{
				++tlab;
				tbaslab[junc.first] = tlab;
			}
			for(int i =0; i<this->nnodes;++i)
			{
				if(this->can_flow_out_there(i))
				{
					++tlab;
					tbaslab[i] = tlab;
				}
			}

			// Propagating the labels up
			for(auto tn:this->stack)
			{
				if(this->can_flow_out_there(tn) || this->can_flow_even_go_there(tn) == false || tbaslab[tn] != -1 )
				{
					continue;
				}
				tbaslab[tn] = tbaslab[this->receivers[tn]];
			}

			// marking up nodes of interest
			std::vector<bool> marked(this->nnodes,false);
			for(int i=0; i<this->nnodes; ++i)
			{
				if(this->can_flow_out_there(i) || this->can_flow_even_go_there(i) == false || marked[i])
					continue;

				int tlab = tbaslab[i];
				auto neighs = this->get_D4_neighbours_only_id(i);
				for(auto tn:neighs)
				{
					if(tn < 0)
						continue;
					if(tlab != tbaslab[tn])
					{
						marked[i] = true;
						marked[tn] = true;
					}
				}
			}

			// inverting topo
			std::vector<float> itopo(this->nnodes, 1);
			//first pass
			for(int i = this->nnodes - 1; i>=0; --i)
			{
				int tn = this->stack[i];
				int rec = this->receivers[tn];
				if(tn == rec)
					continue;

				if(marked[tn])
					itopo[tn] = 0;

				if(itopo[tn] <= 0)
				{
					itopo[rec] = (itopo[rec] == 1) ? itopo[tn] - std::abs(this->topography[tn] - this->topography[rec]) : std::fmax(itopo[tn] - std::abs(this->topography[tn] - this->topography[rec]),itopo[rec]);
				}
			}
			// second pass
			for(auto tn : this->stack)
			{
				if(itopo[tn] == 1)
				{
					int rec = this->receivers[tn];
					if(rec != tn)
						itopo[tn] = std::fmin(itopo[rec] + std::abs(this->topography[tn] - this->topography[rec]), 0);
				}
			}

			return itopo;

		}



		float _max(std::vector<float>& vec)
		{
			float maxv = std::numeric_limits<float>::min();
			for(int i=0; i<this->nnodes; ++i)
			{
				if(this->boundary[i] >= 0)
				{
					if(vec[i] > maxv)
						maxv = vec[i];
				}
			}
			return maxv;
		}






		void compute_all_basins()
		{
			int tlabel = -1;
			this->basin_labels = std::vector<int>(this->nnodes,-1);

			for(int ti=0; ti<this->nnodes; ++ti)
			{
				int i = this->stack[ti];
				if(this->can_flow_even_go_there(i))
				{
					if(this->can_flow_out_there(i))
						++tlabel;
					this->basin_labels[i] = tlabel;
				}
			}
			this->nbasins = tlabel + 1;
		}


		void compute_DD_simple()
		{
			this->isDD = std::vector<bool>(this->nnodes, false);
			for(int i=0; i< this->nnodes ; ++i)
			{
				int tlab = this->basin_labels[i];
				if(tlab == -1)
					continue;

				auto neighs = this->get_D4_neighbours_only_id(i);
				for(auto tn:neighs)
					if(tlab != this->basin_labels[tn])
						isDD[i] = true;
			}

		}


		Graph downscale(float factor)
		{
			if(factor <= 1)
				throw std::runtime_error("Cannot downscale bya factor <= 1");

			int newnx = floor(this->nx/factor);
			float coefx = this->nx/factor;
			int newny = floor(this->ny/factor);
			float coefy = this->ny/factor;
			int newnxy = newnx * newny;
			float newdx = this->dx*factor * newnx/coefx;
			float newdy = this->dy*factor * newny/coefy;

			float coor_factor = this->nnodes/newnxy;
			float coor_factorx = this->nx/newnx;
			float coor_factory = this->ny/newny;
			std::vector<float> newtopo(newnxy, 0.);

			int nx2average = floor(this->nx/newnx);
			int ny2average = floor(this->ny/newny);
			// std::cout << nx2average << "||" << ny2average << std::endl;

			int oro = 0;
			for (int r =0; r < newny; ++r)
			{
				oro = r * ny2average;
				int oco = 0;
				for (int c =0; c < newnx; ++c)
				{
					oco = nx2average * c;
					int newi = r * newnx + c;
					int i = oro * this->nx + oco;

					int N = 1;
					float mean = this->topography[i];
					int right = i;
					for (int iright=i; iright<i + nx2average; ++iright)
					{
						

						int bottom = right;
						for (int ibottom=bottom; ibottom<right + ny2average; ++ibottom)
						{
							bottom = this->get_bottom_index(bottom);
							if(bottom >= this->nnodes || bottom < 0)
								break;
							++N;
							mean += this->topography[bottom];
						}

						right = this->get_right_index(right);
						if(right >= this->nnodes || right < 0)
							break;
					}

					newtopo[newi] = mean/N;

				}
			}


			return Graph(newnx, newny, newnxy, newdx, newdy, this->Xmin, this->Ymin, newtopo);
		}


		void add_label(Label& l){this->labels.emplace_back(l);}





		void fast_gblur(float r)
		{
			std::vector<float>newtopo(this->topography);
			this->gaussBlur_4 (this->topography, newtopo, r);
			this->topography = std::move(newtopo);
		}

		std::vector<int> boxesForGauss(float sigma, int n)  // standard deviation, number of boxes
		{
		    float wIdeal = std::sqrt((12*sigma*sigma/n)+1);  // Ideal averaging filter width 
		    int wl = std::floor(wIdeal);  if(wl%2==0) wl--;
		    int wu = wl+2;

		    float mIdeal = (12*sigma*sigma - n*wl*wl - 4*n*wl - 3*n)/(-4*wl - 4);
		    int m = std::round(mIdeal);
		    // var sigmaActual = Math.sqrt( (m*wl*wl + (n-m)*wu*wu - n)/12 );
						
		    std::vector<int> sizes; sizes.reserve(n);  
		    for(int i=0; i<n; ++i) sizes.emplace_back(i<m?wl:wu);
		    
		    return sizes;
		}


		void gaussBlur_4 (std::vector<float>& scl, std::vector<float>& tcl, float r) 
		{
			auto bxs = this->boxesForGauss(r, 3);
			this->boxBlur_4 (scl, tcl, this->nx, this->ny, (bxs[0]-1)/2);
			this->boxBlur_4 (tcl, scl, this->nx, this->ny, (bxs[1]-1)/2);
	  	this->boxBlur_4 (scl, tcl, this->nx, this->ny, (bxs[2]-1)/2);
		}

		void boxBlur_4 (std::vector<float>& scl, std::vector<float>& tcl, int w, int h, float r) {
		    for(int i=0; i<scl.size(); ++i) tcl[i] = scl[i];
		    this->boxBlurH_4(tcl, scl, w, h, r);
		    this->boxBlurT_4(scl, tcl, w, h, r);
		}

		void boxBlurH_4 (std::vector<float>& scl, std::vector<float>& tcl, int w, int h, float r) {
		    float iarr = 1 / (r+r+1);
		    for(int i=0; i<h; ++i) {
		        int ti = i*w, li = ti, ri = ti+r;
		        float fv = scl[ti], lv = scl[ti+w-1], val = (r+1)*fv;
		        for(int j=0; j<r; ++j) val += scl[ti+j];
		        for(int j=0  ; j<=r ; ++j) { val += scl[ri++] - fv       ;   tcl[ti++] =std::round(val*iarr); }
		        for(int j=r+1; j<w-r; ++j) { val += scl[ri++] - scl[li++];   tcl[ti++] =std::round(val*iarr); }
		        for(int j=w-r; j<w  ; ++j) { val += lv - scl[li++];   tcl[ti++] =std::round(val*iarr); }
		    }
		}
		void boxBlurT_4 (std::vector<float>& scl, std::vector<float>& tcl, int w, int h, float r) {
		    float iarr = 1 / (r+r+1);
		    for(int i=0; i<w; i++) 
		    {
		        int ti = i, li = ti, ri = ti+r*w;
		        float fv = scl[ti], lv = scl[ti+w*(h-1)], val = (r+1)*fv;
		        for(int j=0; j<r; ++j) val += scl[ti+j*w];
		        for(int j=0  ; j<=r ; ++j) { val += scl[ri] - fv     ;  tcl[ti] = std::round(val*iarr);  ri+=w; ti+=w; }
		        for(int j=r+1; j<h-r; ++j) { val += scl[ri] - scl[li];  tcl[ti] = std::round(val*iarr);  li+=w; ri+=w; ti+=w; }
		        for(int j=h-r; j<h  ; ++j) { val += lv      - scl[li];  tcl[ti] = std::round(val*iarr);  li+=w; ti+=w; }
		    }
		}


		void add_label_full(float mref, float nref, float Kref, float S_c, float Acrit)
		{
			Label tlab = Label();
			tlab.mref = mref;
			tlab.nref = nref;
			tlab.Kref = Kref;
			tlab.S_c = S_c;
			tlab.Acrit = Acrit;
			this->labels.emplace_back(tlab);
		}

		void _run_nit_v4(int nit_chi, int nit_gblur, int nit_end, std::vector<int> theselabels, std::vector<float> Es, float final_blur_sigma)
		{

			// std::cout << "A" << std::endl;
			this->generate_white_noise(1);
			for(int i=0;i<this->nnodes; ++i)
			{
				if(this->can_flow_out_there(i))
					continue;
				this->topography[i] += 1;
			}
			float step_chi = 0.4/nit_chi;
			float min_A = this->dx * this->dy * this->nnodes/10;
			// std::cout << "E" << std::endl;

			for (int i=0; i<nit_chi; ++i)
			{
				// std::cout << "E" << i << std::endl;
				this->generate_random_landscape_v4(theselabels, Es, step_chi * i, false, 0, true, min_A);
				for(auto v:this->topography)
					if(std::isfinite(v) == false)
					{
						// std::cout << "gabulon" << std::endl;
						throw std::runtime_error("fdsljk");
					}
			}
			// std::cout << "F" << std::endl;

			for (int i=0; i<nit_gblur; ++i)
			{
				this->generate_random_landscape_v4(theselabels, Es, 1, false, 3, false, min_A);
			}
			// std::cout << "G" << std::endl;

			for (int i=0; i<nit_end; ++i)
			{
				this->generate_random_landscape_v4(theselabels, Es, 1, false, 0, false, min_A);
			}
			// std::cout << "H" << std::endl;

			if(final_blur_sigma>0)
				this->fast_gblur(final_blur_sigma);
			this->set_boundaries_to(0);

		}

		std::vector<float> calculate_diff_from_Sc(float Sc)
		{
			std::vector<float> out(this->nnodes,0.), fake_elev(this->nnodes,0.);

			for(auto tn:this->river_nodes)
				fake_elev[tn] = this->topography[tn];

			for(auto tn:this->stack)
			{
				if(this->can_flow_out_there(tn) || this->can_flow_even_go_there(tn) == false || fake_elev[tn]>0)
					continue;
				int rec = this->receivers[tn];
				fake_elev[tn] = fake_elev[rec] + this->distance2receivers[tn] * Sc;
				out[tn] = this->topography[tn] - fake_elev[tn];
			}
			return out;

		}


		void reset_labels(){this->labels.clear();}

		std::vector<float> compute_distance_to_nodes(std::vector<int>& nodes,float maxv)
		{
			std::vector<float> distfromstuff(this->nnodes,std::numeric_limits<float>::max());
			std::priority_queue<PQ_helper<int, float>, std::vector<PQ_helper<int, float>>, std::greater<PQ_helper<int, float> > > PQ;

			for(auto i:nodes)
				PQ.emplace(PQ_helper<int, float>(i,0));



			while(PQ.empty() == false)
			{
				auto next = PQ.top();PQ.pop();

				if(next.score < distfromstuff[next.node] && next.score < maxv)
				{
					distfromstuff[next.node] = next.score;
					auto neighbours = this->get_neighbours(next.node);
					float td = next.score;
					for(auto& tn:neighbours)
					{
						if( distfromstuff[tn.node] > td + tn.distance )
							PQ.emplace( PQ_helper<int,float>(tn.node,td + tn.distance) );
					}
				}
			}

			for(auto& v:distfromstuff)
			{
				if(v == std::numeric_limits<float>::max())
				{
					v = maxv;
				}
			}

			return distfromstuff;



		}

		std::vector<float> compute_distance_to_boundaries(float maxv)
		{
			std::vector<float> distfromstuff(this->nnodes,std::numeric_limits<float>::max());
			std::priority_queue<PQ_helper<int, float>, std::vector<PQ_helper<int, float>>, std::greater<PQ_helper<int, float> > > PQ;
			for (int i=0; i<this->nnodes;++i)
			{
				if(this->can_flow_out_there(i))
				{
					PQ.emplace(PQ_helper<int, float>(i,0));
				}
			}

			while(PQ.empty() == false)
			{
				auto next = PQ.top();PQ.pop();

				if(next.score < distfromstuff[next.node] && next.score < maxv)
				{
					distfromstuff[next.node] = next.score;
					auto neighbours = this->get_neighbours(next.node);
					float td = next.score;
					for(auto& tn:neighbours)
					{
						if( distfromstuff[tn.node] > td + tn.distance )
							PQ.emplace( PQ_helper<int,float>(tn.node,td + tn.distance) );
					}
				}
			}

			for(auto& v:distfromstuff)
			{
				if(v == std::numeric_limits<float>::max())
				{
					v = maxv;
				}
			}

			return distfromstuff;



		}











//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Emscripten only function
#ifdef __EMSCRIPTEN__

		void run_nit_v4(int nit_chi, int nit_gblur, int nit_end, emscripten::val theselabels, emscripten::val Es, float final_blur_sigma)
		{

			// std::cout << "A" << std::endl;
			this->generate_white_noise(1);
			for(int i=0;i<this->nnodes; ++i)
			{
				if(this->can_flow_out_there(i))
					continue;
				this->topography[i] += 1;
			}
			// std::cout << "B" << std::endl;
			std::vector<float> tEs = emscripten::convertJSArrayToNumberVector<float>(Es);
			// std::cout << "C" << std::endl;
			std::vector<int> ttheselabels = emscripten::convertJSArrayToNumberVector<int>(theselabels);
			// std::cout << "D::" << tEs.size() << "|" << ttheselabels.size() << "|" << this->nnodes << std::endl;

			float step_chi = 0.4/nit_chi;
			float min_A = this->dx * this->dy * this->nnodes/10;
			// std::cout << "E" << std::endl;

			for (int i=0; i<nit_chi; ++i)
			{
				// std::cout << "E" << i << std::endl;
				this->generate_random_landscape_v4(ttheselabels, tEs, step_chi * i, false, 0, true, min_A);
				for(auto v:this->topography)
					if(std::isfinite(v) == false)
					{
						// std::cout << "gabulon" << std::endl;
						throw std::runtime_error("fdsljk");
					}
			}
			// std::cout << "F" << std::endl;

			for (int i=0; i<nit_gblur; ++i)
			{
				this->generate_random_landscape_v4(ttheselabels, tEs, 1, false, 3, false, min_A);
			}
			// std::cout << "G" << std::endl;

			for (int i=0; i<nit_end; ++i)
			{
				this->generate_random_landscape_v4(ttheselabels, tEs, 1, false, 0, false, min_A);
			}
			// std::cout << "H" << std::endl;

			if(final_blur_sigma>0)
				this->fast_gblur(final_blur_sigma);
			this->set_boundaries_to(0);

		}

		emscripten::val gettopo()
		{
			return emscripten::val::array(this->topography);
		}



		void ingest_topo(emscripten::val ttopo)
		{
			this->topography = emscripten::convertJSArrayToNumberVector<float>(ttopo);
			std::cout << "Ingested topo has " << this->topography.size() << " and should be " << this->nnodes << std::endl;
		}

		void save2file(std::string file_prefix, std::vector<float>& data)
		{
				// EM_ASM is a macro to call in-line JavaScript code.
			EM_ASM(
			    // Make a directory other than '/'
			    FS.mkdir('/wgt_data');
			    // Then mount with IDBFS type
			    FS.mount(IDBFS, {}, '/wgt_data');

			    // Then sync
			    FS.syncfs(true, function (err) {
			        // Error
			    });
			);

			long unsigned A = this->ny, B = this->nx;

			std::array<long unsigned, 2> leshape31 = {A, B};

			npy::SaveArrayAsNumpy("/wgt_data/"+ file_prefix +".npy", false, leshape31.size(), leshape31.data(), data);

			// Don't forget to sync to make sure you store it to IndexedDB
			EM_ASM(
			    FS.syncfs(function (err) {
			        // Error
			    });
			);
		}

		emscripten::val getXriv()
		{
			std::vector<float> rXs; rXs.reserve(this->getNriv());
			int i=0;
			for(auto tn:this->river_nodes)
			{
				// ++i;
				// if(i < 20)
				// 	std::cout << this->get_X_from_i(tn) << std::endl;
				rXs.emplace_back(this->get_X_from_i(tn));
			}

			// return emscripten::val{ emscripten::typed_memory_view(rXs.size(), rXs.data()) };
			return emscripten::val::array(rXs);
		}

		emscripten::val getYriv()
		{
			std::vector<float> rYs; rYs.reserve(this->getNriv());
			for(auto tn:this->river_nodes)
				rYs.emplace_back(this->get_Y_from_i(tn));

			// return emscripten::val{ emscripten::typed_memory_view(rYs.size(), rYs.data()) };
			return emscripten::val::array(rYs);
		
		}

		emscripten::val getAriv()
		{
			std::vector<float> rA; rA.reserve(this->getNriv());
			for(auto tn:this->river_nodes)
				rA.emplace_back(this->area[tn]);

			// return emscripten::val{ emscripten::typed_memory_view(rA.size(), rA.data()) };
			return emscripten::val::array(rA);

		}

		

		emscripten::val get_HS()
		{
			// float* hillshade = reinterpret_cast<float*>(yolo);
			// float hillshade[] = ptr;

			float altitude = 45;
			float azimuth = 315;
			float z_factor = 1;
			float pi = 3.1415926;

			// std::vector<float> hillshade(ptr, ptr + this->nnodes);
			std::vector<float> hillshade(this->nnodes_t,0.);


			//convert zenith and azimuth into radians for calculation
			float zenith_rad = (90 - altitude) * pi / 180.0;
			float azimuth_math = 360-azimuth + 90;
			if (azimuth_math >= 360.0) azimuth_math = azimuth_math - 360;
			float azimuth_rad = azimuth_math * pi /180.0;

			for(int i = 0; i<this->nnodes; ++i)
			{
				// Ignoring no data
				if(this->boundary[i] < 0)
					continue;

				float slope_rad = 0;
				float aspect_rad = 0;
				float dzdx = 0;
				float dzdy = 0;

				float ij = this->topography[i];
				float ijp1 = this->topography[this->get_right_neighbour(i).node];
				float ip1j = this->topography[this->get_bottom_neighbour(i).node];
				float ip1jp1 = this->topography[this->get_bottomright_neighbour(i).node];
				float im1jm1 = this->topography[this->get_topleft_neighbour(i).node];
				float im1j = this->topography[this->get_top_neighbour(i).node];
				float im1jp1 = this->topography[this->get_topright_neighbour(i).node];
				float ijm1 = this->topography[this->get_left_neighbour(i).node];
				float ip1jm1 = this->topography[this->get_bottomleft_neighbour(i).node];


				if (ij > 0 )
				{
					dzdx = ((ijp1 + 2*ip1j + ip1jp1) - (im1jm1 + 2*im1j + im1jp1)) / (8 * this->dx);
					dzdy = ((im1jp1 + 2*ijp1 + ip1jp1) - (im1jm1 + 2*ijm1 + ip1jm1)) / (8 * this->dy);
					slope_rad = atan(z_factor * sqrt((dzdx*dzdx) + (dzdy*dzdy)));
					if (dzdx != 0)
					{
						aspect_rad = std::atan2(dzdy, (dzdx*-1));
						if (aspect_rad < 0) aspect_rad = 2*pi + aspect_rad;
					}
					else
					{
						if (dzdy > 0) aspect_rad = pi/2;
						else if (dzdy < 0) aspect_rad = 2 * pi - pi/2;
						else aspect_rad = aspect_rad;
					}

					hillshade[i] = ((std::cos(zenith_rad) * std::cos(slope_rad)) + (std::sin(zenith_rad) * std::sin(slope_rad) * std::cos(azimuth_rad - aspect_rad)));
					// std::cout << hillshade[i] << "|";
					if (hillshade[i] < 0) hillshade[i] = 0;
				}

			}
					// std::cout << "Done" << std::endl;

			return emscripten::val{ emscripten::typed_memory_view(hillshade.size(), hillshade.data()) };
    
		}


// Python only
#else
		Graph(int nx, int ny, int nnodes, float dx, float dy, float xmin, float ymin, py::array_t<float>& numtopo)
		{
			this->set_dimensions(nx,ny,nnodes,dx,dy,xmin,ymin);
			this->set_default_boundaries("4edges");
			this->topography.clear();
			this->topography.reserve(this->nnodes_t);
			py::buffer_info yolo = numtopo.request();
			float *baft = (float *)yolo.ptr;
			for(int i=0; i<this->nnodes; ++i)
				this->topography.emplace_back(baft[i]);
			// std::cout << "Got topo:: " << this->topography.size() << " nnodes:: " << this->nnodes << std::endl;
		}

		py::array get_HS()
		{



			float altitude = 45;
			float azimuth = 315;
			float z_factor = 1;
			float pi = 3.1415926;

			// std::vector<float> hillshade(ptr, ptr + this->nnodes);
			std::vector<float> hillshade(this->nnodes_t,0.);

			float min_value = 0;
			for(int i=0; i<this->nnodes;++i)
			{
				if(this->topography[i] < min_value)
					min_value = this->topography[i];
			}

			// if(min_value)
			//convert zenith and azimuth into radians for calculation
			float zenith_rad = (90 - altitude) * pi / 180.0;
			float azimuth_math = 360-azimuth + 90;
			if (azimuth_math >= 360.0) azimuth_math = azimuth_math - 360;
			float azimuth_rad = azimuth_math * pi /180.0;

			for(int i = 0; i<this->nnodes; ++i)
			{
				// Ignoring no data
				if(this->boundary[i] < 0)
				{
					continue;
				}

				float slope_rad = 0;
				float aspect_rad = 0;
				float dzdx = 0;
				float dzdy = 0;

				float ij = this->topography[i] + std::abs(min_value);
				float ijp1 = this->topography[this->get_right_neighbour(i).node] + std::abs(min_value);
				float ip1j = this->topography[this->get_bottom_neighbour(i).node] + std::abs(min_value);
				float ip1jp1 = this->topography[this->get_bottomright_neighbour(i).node] + std::abs(min_value);
				float im1jm1 = this->topography[this->get_topleft_neighbour(i).node] + std::abs(min_value);
				float im1j = this->topography[this->get_top_neighbour(i).node] + std::abs(min_value);
				float im1jp1 = this->topography[this->get_topright_neighbour(i).node] + std::abs(min_value);
				float ijm1 = this->topography[this->get_left_neighbour(i).node] + std::abs(min_value);
				float ip1jm1 = this->topography[this->get_bottomleft_neighbour(i).node] + std::abs(min_value);


				if (ij > 0 )
				{
					dzdx = ((ijp1 + 2*ip1j + ip1jp1) - (im1jm1 + 2*im1j + im1jp1)) / (8 * this->dx);
					dzdy = ((im1jp1 + 2*ijp1 + ip1jp1) - (im1jm1 + 2*ijm1 + ip1jm1)) / (8 * this->dy);
					slope_rad = atan(z_factor * sqrt((dzdx*dzdx) + (dzdy*dzdy)));
					if (dzdx != 0)
					{
						aspect_rad = std::atan2(dzdy, (dzdx*-1));
						if (aspect_rad < 0) aspect_rad = 2*pi + aspect_rad;
					}
					else
					{
						if (dzdy > 0) aspect_rad = pi/2;
						else if (dzdy < 0) aspect_rad = 2 * pi - pi/2;
						else aspect_rad = aspect_rad;
					}

					hillshade[i] = 255.0 * ((std::cos(zenith_rad) * std::cos(slope_rad)) + (std::sin(zenith_rad) * std::sin(slope_rad) * std::cos(azimuth_rad - aspect_rad)));
					// std::cout << hillshade[i] << "|";
					// if (hillshade[i] < 0) hillshade[i] = 0;
				}

			}
			return py::array(hillshade.size(), hillshade.data());
		}

		void np2topo(py::array_t<float>& ntopo)
		{
			py::buffer_info yolo = ntopo.request();
			float *baft = (float *)yolo.ptr;

			for (int i=0; i<this->nnodes;++i)
				this->topography[i] = baft[i];
		}

		void add2topo(py::array_t<float>& nadd)
		{
			py::buffer_info yolo = nadd.request();
			float *baft = (float *)yolo.ptr;

			for (int i=0; i<this->nnodes;++i)
				this->topography[i] += baft[i];
		}

		py::array get_topo_np(){return py::array(this->topography.size(), this->topography.data());}
		py::array get_boundaries(){return py::array(this->boundary.size(), this->boundary.data());}

		py::array get_receiver_array(){return py::array(this->receivers.size(), this->receivers.data());}
		py::array get_distance_to_receivers_array(){return py::array(this->distance2receivers.size(), this->distance2receivers.data());}

		py::array get_stack(){return py::array(this->stack.size(), this->stack.data());}



		py::array get_faces()
		{
			std::vector<int> faces;
			faces.reserve((this->nx - 1) * (this->ny - 1) * 4);
			for(int row = 0; row < this->ny - 1 ; ++row)
			{
				for(int col = 0; col < this->nx - 1 ; ++col)
				{
					int tid = row * this->nx + col;
					faces.emplace_back(this->get_bottom_index(tid));
					faces.emplace_back(this->get_bottomright_index(tid));
					faces.emplace_back(this->get_right_index(tid));
					faces.emplace_back(tid);
				}

			}

			return py::array(faces.size(), faces.data());

		}

		py::dict get_rivers_dict()
		{
			std::cout << "stuffy1" << std::endl;
			py::dict output;

			std::vector<float> rYs; rYs.reserve(this->getNriv());
			std::vector<float> rXs; rXs.reserve(this->getNriv());
			std::vector<float> rZs; rZs.reserve(this->getNriv());
			std::vector<float> As; As.reserve(this->getNriv());
			for(auto tn:this->river_nodes)
			{
				// std::cout <<"1|" << std::endl;
				rYs.emplace_back(this->get_Y_from_i(tn));
				// std::cout <<"2|" << std::endl;
				rXs.emplace_back(this->get_X_from_i(tn));
				// std::cout <<"3|" << std::endl;
				rZs.emplace_back(this->topography[tn]);
				// std::cout <<"4|" << std::endl;
				As.emplace_back(this->area[tn]);
				// std::cout <<"5|" << std::endl;
			}

			output[py::str("Y")] = py::array(rYs.size(), rYs.data());			
			output[py::str("X")] = py::array(rXs.size(), rXs.data());
			output[py::str("Z")] = py::array(rZs.size(), rZs.data());
			output[py::str("flow_distance")] = py::array(this->flowdistance.size(), this->flowdistance.data());
			output[py::str("A")] = py::array(As.size(), As.data());
			output[py::str("key")] = py::array(sourceID.size(), sourceID.data());
			
			return output;
		}

		py::dict get_rivers_rowcolnode()
		{
			py::dict output;

			std::vector<int> rows; rows.reserve(this->getNriv());
			std::vector<int> cols; cols.reserve(this->getNriv());
			std::vector<int> rivrec; rivrec.reserve(this->getNriv());
			std::vector<int> globrec; globrec.reserve(this->getNriv());
			for(auto tn:this->river_nodes)
			{
				int row,col;
				this->rowcol_from_node_id(tn, row, col);
				rows.emplace_back(row);
				cols.emplace_back(col);
				rivrec.emplace_back((this->node2rnode.count(this->receivers[tn]) > 0) ?this->node2rnode[this->receivers[tn]]:tn );
				globrec.emplace_back(this->receivers[tn]);
			}

			output["row"] = py::array(rows.size(), rows.data());			
			output["col"] = py::array(cols.size(), cols.data());
			output["node"] = py::array(this->river_nodes.size(), this->river_nodes.data());
			output["rivrec"] = py::array(rivrec.size(), rivrec.data());
			output["globrec"] = py::array(globrec.size(), globrec.data());
			return output;
		}

		py::array get_rivers_basin_label()
		{
			std::vector<int> rBz; rBz.reserve(this->getNriv());
			for(auto tn:this->river_nodes)
				rBz.emplace_back(this->basin_labels[tn]);
			return py::array(rBz.size(), rBz.data());
		}

		py::array get_stream_order()
		{
			auto yolo = this->_get_stream_order();
			return py::array(yolo.size(), yolo.data());
		}

		py::array get_basin_array()
		{
			return py::array(this->basin_labels.size(),this->basin_labels.data());
		}

		py::array invert_topography(float Ath)
		{
			auto output = this->_invert_topography(Ath);
			return py::array(output.size(), output.data());
		}

		py::dict get_junction_map()
		{
			py::dict output;
			std::vector<int> OOOOORDEEEEEEEERRRRR;
			std::vector<float> Xs;
			std::vector<float> Ys;

			for(auto junc:this->junctions)
			{
				OOOOORDEEEEEEEERRRRR.emplace_back(junc.second);
				Xs.emplace_back(this->get_X_from_i(junc.first));
				Ys.emplace_back(this->get_Y_from_i(junc.first));
			}

			output["X"] = py::array(Xs.size(),Xs.data());
			output["Y"] = py::array(Ys.size(),Ys.data());
			output["SO"] = py::array(OOOOORDEEEEEEEERRRRR.size(),OOOOORDEEEEEEEERRRRR.data());
			return output;


		}

		// This function returns a python dictionnary with all INDIVIDUAL FULL RIVERS
		// I cannot stress enough how important it is, a lot of nodes will be duplicated.
		// This is a rather specific function intended to select a specific river easily
		py::dict get_full_rivers_of_basin(int bkey)
		{
			py::dict output;

			for(auto sk:this->source_nodes)
			{
				if(this->basin_labels[sk] != bkey)
					continue;

				std::vector<int> nodes, flow_distance, elevation;

				int tnode = sk;
				do
				{
					nodes.emplace_back(tnode);
					flow_distance.emplace_back(this->flowdistance[this->node2rnode[tnode]]);
					elevation.emplace_back(this->topography[tnode]);
					tnode = this->receivers[tnode];
				}while(this->can_flow_out_there(tnode) == false);

				py::dict temp;
				temp["node"] = py::array(nodes.size(), nodes.data());
				temp["flow_distance"] = py::array(flow_distance.size(), flow_distance.data());
				temp["elevation"] = py::array(elevation.size(), elevation.data());
				std::string ttttttttt = std::to_string(sk);
				output[ttttttttt.c_str()] = temp;
			}
			return output;
		}

		py::array get_diff_from_Sc(float Sc)
		{
			auto res = this->calculate_diff_from_Sc(Sc);
			return py::array(res.size(), res.data());

		}


		py::array get_debug_f(){return py::array(this->debugfloat.size(), this->debugfloat.data());}
		py::array get_flow_distance_full(){auto FD = this->compute_flow_distance_full();return py::array(FD.size(), FD.data());}
		py::array get_drainage_area(){this->calculate_area(); return py::array(this->area.size(), this->area.data());}

		py::dict get_dimensions()
		{
			py::dict output;
			output["nx"] = this->nx;
			output["ny"] = this->ny;
			output["dx"] = this->dx;
			output["dy"] = this->dy;
			output["nnodes"] = this->nnodes;
			py::list yolo;
			yolo.append(this->Xmin);
			yolo.append(this->Xmin + this->nx * this->dx + this->dx);
			yolo.append(this->Ymin) ;
			yolo.append(this->Ymin + this->ny * this->dy + this->dy);
			output["extents"] = yolo;
			return output;
		}



#endif



};






























#endif