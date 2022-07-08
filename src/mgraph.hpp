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
template<T>
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

		std::vector<std:vector<int> > receivers;
		std::vector<std::vector<int> > donors;
		std::vector<std:vector<T> > distance2receivers;
ß
		// #->stack: topological order from downstream to upstream direction
		std::vector<int> stack;

		// Coordinate stuff
		// Xs and Ys are vectors of nx and ny size converting row to y and col to X
		// Extents holds the cxmin,xmax,ymin,ymax (extent is an option in matplotlib imshow plots)
		std::vector<T> Xs,Ys, extents;

	


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

		// ------------------------------------------------

		//	                             	              __
		//                                             / _)
		//                                    _.----._/ /
		//   GRAPH CREATION METHODS          /         /
		//                                __/ (  | (  |
		//                               /__.-'|_|--|_|

		// All the methods related to accessing and calculating neighbours
		// ------------------------------------------------

		void compute_graph(std::string depression_solver)
		{
			this->compute_graph_both_v2();
			this->compute_TO_SF_stack_version();
			if(depression_solver != "none")
			{
				this->solve_depressions(depression_solver);
				this->compute_TO_SF_stack_version();
			}
		}


		void compute_graph_both_v2()
		{
			// Initialising the graph dimesions for the receivers
			// All of thenm have the graph dimension
			this->Sreceivers.clear();
			this->Sdistance2receivers.clear();
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
			{
				this->receivers[i] = i;
			}

			// Iterating through all the nodes and finding the max slope
			int n_donors_iintotal = 0;
			// bool switc = false;

			
			}
			this->recompute_SF_donors_from_receivers();
		}

		void recompute_MF_impose_slope_SS()
		{
			std::vector<T> temptopo((*this->topography));
			std::vector<T>* otopo = this->topography;

			this->topography = &temptopo;
			this->enforce_minimal_slope_SS(1e-3);

			this->receivers.clear();
			this->donors.clear();
			this->distance2receivers.clear();
			this->receivers = std::vector<std::vector<int> >(this->nnodes);
			this->donors = std::vector<std::vector<int> >(this->nnodes);
			this->distance2receivers = std::vector<std::vector<T> >(this->nnodes);

			for(int row = 0; row < this->ny; ++row)
			{
				for(int col = 0; col < this->nx; ++col)
				{
					int i = row * this->nx + col;
					// cannot be a neighbour anyway, abort
					if(this->can_flow_even_go_there(i) == false)
					{
						continue;
					}

					bool check = (col > 1 && row > 1 && col < this->nx - 2 && row < this->ny - 2 ) ? false : true;


					if(col > 0 && row > 0 && col < this->nx -1 && row < this->ny - 1 )
					{
						this->check_neighbour_v22_MF(check, this->get_topleft_index(i), i, this->dxy);
						this->check_neighbour_v22_MF(check,this->get_top_index(i), i, this->dy);
						this->check_neighbour_v22_MF(check,this->get_topright_index(i), i, this->dxy);
						this->check_neighbour_v22_MF(check,this->get_right_index(i), i, this->dx);
					}
					else
					{
						this->check_neighbour_v22_MF(check,this->get_top_index(i), i, this->dy);
						this->check_neighbour_v22_MF(check,this->get_left_index(i), i, this->dx);
						this->check_neighbour_v22_MF(check,this->get_right_index(i), i, this->dx);
						this->check_neighbour_v22_MF(check,this->get_bottom_index(i), i, this->dy);
						this->check_neighbour_v22_MF(check,this->get_topright_index(i), i, this->dxy);
						this->check_neighbour_v22_MF(check,this->get_topleft_index(i), i, this->dxy);
						this->check_neighbour_v22_MF(check,this->get_bottomright_index(i), i, this->dxy);
						this->check_neighbour_v22_MF(check,this->get_bottomleft_index(i), i, this->dxy);
					}
				}
			}

			this->topography = otopo;
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

		void solve_depressions(std::string& depression_solver)
		{
			
			Cordonnier2019_v2MF<int,T> depsolver(*this);

			if(depsolver.npits > 0)
				depsolver.update_receivers(depression_solver);
		}


		
		/// this function enforces minimal slope 
		void enforce_minimal_slope_SS(T slope)
		{
			for(auto node:this->stack)
			{
				if(this->can_flow_out_there(node) || this->can_flow_even_go_there(node) == false)
					continue;
				int rec = this->Sreceivers[node];
				T dz = (*this->topography)[node] - (*this->topography)[rec];

				if(dz <= 0)
				{
					T d2rec = this->Sdistance2receivers[node];
					(*this->topography)[node] = (*this->topography)[rec] + slope * d2rec;
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

		inline int nodeid_from_XY(T X, T Y)
		{
			int col = floor((X - this->Xmin)/this->dx);
			int row = this->ny - floor((Y - this->Ymin)/this->dy);
			return this->nodeid_from_row_col(row,col);
		}
		// void rowcol_from_XY()


		inline void XY_from_nodeid(int node_index, T& tX, T& tY)
		{
			int row,col;
			this->rowcol_from_node_id(node_index, row, col);
			this->XY_from_rowcol(row, col, tX, tY);
		}

		inline void XY_from_rowcol(int row, int col, T& tX, T& tY)
		{
			tX = this->Xs[col];
			tY = this->Ys[row];
		}

		inline T get_X_from_i(int i)
		{
			int col = i % this->nx;
			return this->Xs[col];
		}

		inline T get_Y_from_i(int i)
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

		std::vector<Neighbour<int,T> > get_neighbours(int i, bool ignore_nodata = false)
		{

			// preformatting the output
			std::vector<Neighbour<int,T> > neighbours;

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
		void set_neighbours_in_vector(std::vector<Neighbour<int,T> >& neighbours, int& i, size_t& id_neighbourer, bool ignore_nodata)
		{
			// node index of the current neighbour
			int tn;

			// D4 or D8 I take it all	
			// If you wonder why I am not iterating with a vector and everything here, it is for the small perf gain of doing it this way
			tn = this->neighbourer[id_neighbourer][1];
			if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
				neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[1]));
			tn = this->neighbourer[id_neighbourer][3];
			if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
				neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[3]));
			tn = this->neighbourer[id_neighbourer][4];
			if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
				neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[4]));
			tn = this->neighbourer[id_neighbourer][6];
			if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
				neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[6]));
			tn = this->neighbourer[id_neighbourer][0];
			if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
				neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[0]));
			tn = this->neighbourer[id_neighbourer][2];
			if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
				neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[2]));
			tn = this->neighbourer[id_neighbourer][5];
			if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
				neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[5]));
			tn = this->neighbourer[id_neighbourer][7];
			if(tn != this->not_a_node && (ignore_nodata == false || (ignore_nodata && this->can_flow_even_go_there(tn))))
				neighbours.emplace_back(Neighbour<int, T>(tn + i, this->lengthener[7]));				
		
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
		Neighbour<int,T> get_left_neighbour(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return Neighbour<int,T>(this->neighbourer[id_neighbourer][3] + i, this->dx);
		}

		Neighbour<int,T> get_top_neighbour(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return Neighbour<int,T>(this->neighbourer[id_neighbourer][1] + i, this->dy);
		}

		Neighbour<int,T> get_right_neighbour(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return Neighbour<int,T>(this->neighbourer[id_neighbourer][4] + i, this->dx);
		}

		Neighbour<int,T> get_bottom_neighbour(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return Neighbour<int,T>(this->neighbourer[id_neighbourer][6] + i, this->dy);
		}

		// D8 single neighbour extra routines
		Neighbour<int,T> get_topleft_neighbour(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return Neighbour<int,T>(this->neighbourer[id_neighbourer][0] + i, this->dxy);
		}

		Neighbour<int,T> get_topright_neighbour(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return Neighbour<int,T>(this->neighbourer[id_neighbourer][2] + i, this->dxy);
		}

		Neighbour<int,T> get_bottomleft_neighbour(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return Neighbour<int,T>(this->neighbourer[id_neighbourer][5] + i, this->dxy);
		}

		Neighbour<int,T> get_bottomright_neighbour(int i)
		{
			size_t id_neighbourer = this->_get_neighbourer_id(i);	
			return Neighbour<int,T>(this->neighbourer[id_neighbourer][7] + i, this->dxy);
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

		


		void compute_MF_topological_order(bool calculate_gradient = true)
		{

			this->stack.clear();
			this->stack.reserve(this->usize);

		  std::vector<T> vis(this->usize,0), parse(this->usize,-1), ndons(this->usize,0);
		  
		  int nparse = -1;
		  int nstack = -1;
		  for(int i=0; i<this->isize;++i)
		  {
		  	ndons[i] = this->get_donors_ID(i).size();
		  }

		  // we go through the nodes
		  for(int i=0; i<this->isize;++i)
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

		      auto recs = this->get_receivers_ID(ijn);

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
		  
		  // std::cout << "TO full finished with " << this->topological_order.size() << " versus " << this->isize << std::endl;;
		}


#endif



};






























#endif