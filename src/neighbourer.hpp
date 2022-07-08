//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef neighbourer_HPP
#define neighbourer_HPP

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



template<class T>
class D8Neighbourer
{
public:

		// General informations about the graph
	// #-> number of nodes (integer and unsized int for loops or list innit).
	int nnodes = 0;
	size_t nnodes_t = 0;
	// #-> number of nodes in x direction.
	int nx = 0;
	// #-> number of nodes in x direction.
	int ny = 0;
	// #-> length in the x direction.
	T dx;
	// #-> length in the y direction.
	T dy;
	T dxy;
	// #-> cell area
	T cellarea;
	T Xmin;
	T Ymin;
	T Xmax;
	T Ymax;

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
	std::vector<T > lengthener;


	D8Neighbourer(int nx, int ny, int nnodes, T dx, T dy, T xmin, T ymin)
	{
		this->set_dimensions(nx,ny,nnodes,dx,dy,xmin,ymin);
		this->set_default_boundaries("4edges");
	}
	
		/// @Description: Initialising the "neighbourer", a data structure managing the iterations though the neighbours of a given node
	/// Function of boundary conditions, one can feed the neighbourer with an indice which returns the index to ADD to the node to get its neighbour.
	/// The neighbour order is (in table referential, be careful if Y axis is inverted) top-left, top, top-right, left, right, bottom-left, bottom, bottom-right
	/// @Authors: B.G.
	/// @Date: 2021
	void initialise_neighbourer()
	{
		T diag = std::sqrt(std::pow(dx,2) + std::pow(dy,2) );
		this->dxy = diag;

		this->lengthener = std::initializer_list<T>{diag,dy,diag,dx,dx,diag,dy,diag};
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
	void set_dimensions(int nx, int ny, int nnodes, T dx, T dy, T xmin, T ymin)
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
		this->Xs = std::vector<T>(this->nx);
		this->Ys = std::vector<T>(this->ny);
		// this->extents = std::vector<T>(4,0);

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


	template<class U, class V>
	void build_mgraph(U& SS, std::vector<std::vector<int> >& Sdonors, std::vector<int>& Sreceivers, std::vector<std::vector<int> >& receivers,
		std::vector<std::vector<int> >& donors, std::vector<T>& Sdistance2receivers, std::vector<std::vector<T> >& distance2receivers, V& topography )
	{
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
					this->check_neighbour_v22(SS,this->get_topleft_index(i), i, this->dxy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22_MF(check, this->get_topleft_index(i), i, this->dxy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22(SS,this->get_top_index(i), i, this->dy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22_MF(check,this->get_top_index(i), i, this->dy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22(SS,this->get_topright_index(i), i, this->dxy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22_MF(check,this->get_topright_index(i), i, this->dxy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22(SS,this->get_right_index(i), i, this->dx,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22_MF(check,this->get_right_index(i), i, this->dx,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
				}
				else
				{
					this->check_neighbour_v22(SS,this->get_top_index(i), i, this->dy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22_MF(check,this->get_top_index(i), i, this->dy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22(SS,this->get_left_index(i), i, this->dx,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22_MF(check,this->get_left_index(i), i, this->dx,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22(SS,this->get_right_index(i), i, this->dx,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22_MF(check,this->get_right_index(i), i, this->dx,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22(SS,this->get_bottom_index(i), i, this->dy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22_MF(check,this->get_bottom_index(i), i, this->dy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22(SS,this->get_topright_index(i), i, this->dxy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22_MF(check,this->get_topright_index(i), i, this->dxy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22(SS,this->get_topleft_index(i), i, this->dxy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22_MF(check,this->get_topleft_index(i), i, this->dxy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22(SS,this->get_bottomright_index(i), i, this->dxy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22_MF(check,this->get_bottomright_index(i), i, this->dxy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22(SS,this->get_bottomleft_index(i), i, this->dxy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
					this->check_neighbour_v22_MF(check,this->get_bottomleft_index(i), i, this->dxy,Sdonors,Sreceivers,receivers,donors,Sdistance2receivers,distance2receivers, topography);
				}
			}

		}
	}

	template<class U, class V>
	void check_neighbour_v22(std::vector<T>& SS, int ne, int i, T dist, std::vector<std::vector<int> >& Sdonors, std::vector<int>& Sreceivers, std::vector<std::vector<int> >& receivers,
		std::vector<std::vector<int> >& donors, std::vector<T>& Sdistance2receivers, std::vector<std::vector<T> >& distance2receivers, V& topography)
	{
		if (ne < 0)
			return;

		if(this->can_flow_even_go_there(ne) == false)
			return;
		
		T this_elev = topography[i];

		if(topography[ne] < this_elev && this->can_flow_out_there(i) == false)
		{
			double this_slope = (this_elev - topography[ne])/dist;
			if(this_slope > SS[i])
			{
				SS[i] = this_slope;
				Sreceivers[i] = ne;
				Sdistance2receivers[i] = dist;
			}
		}
		else if (topography[ne] > this_elev && this->can_flow_out_there(ne) == false)
		{
			double this_slope = (topography[ne] - this_elev)/dist;
			if(this_slope > SS[ne])
			{
				SS[ne] = this_slope;
				Sreceivers[ne] = i;
				Sdistance2receivers[ne] = dist;
			}
		}
	}

	template<class U, class V>
	void check_neighbour_v22_MF(bool check, int from, int to, T dist, std::vector<std::vector<int> >& Sdonors, std::vector<int>& Sreceivers, std::vector<std::vector<int> >& receivers,
		std::vector<std::vector<int> >& donors, std::vector<T>& Sdistance2receivers, std::vector<std::vector<T> >& distance2receivers, V& topography)
	{
		
		if(topograpphy[from] == topograpphy[to])
			return;

		int temp = from;
		if(topograpphy[from] < topograpphy[to])
		{
			from = to;
			to = temp;
		}

		if(check)
		{
			for(auto rec:receivers[from])
			{
				if(rec == to)
					return;
			}
		}
		receivers[from].emplace_back(to);
		distance2receivers[from].emplace_back(dist);
		donors[to].emplace_back(from);
	}



};























































#endif