//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef minilem_HPP
#define minilem_HPP

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

#include "chonkutils.hpp"
#include "graph.hpp"
// #include "wapart.hpp"
#include "PerlinNoise.hpp"

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








class minilem
{

public:

	// Graph

	// The graph manages all the node connectivities and local minima connections
	Graph graph;

	// Graph options
	// -> method for minima (carve, fill, simple or none)
	std::string minima_solver = "carve";
	std::string boundary_string;


	// Parameters

	// Vertical motion mode
	// -> 0 for uniform value
	// -> 1 for spatial
	// -> other potential modes for future: functional, coupled, ...
	int vmmode = 0;
	float uvertmot = 1e-3;
	std::vector<float> vertmot;

	// QArea mode
	// 0 -> drainage area
	// 1 -> simple precipitation (uniform)
	// 2 -> precipitations weighter
	int qamode = 0;
	float uPrec;
	std::vector<float> prec;

	// Parameter mode: 0 for simple, 1 for labels
	// -> if simple all the params are taken from label 0
	int paramode = 0;

	// Spatial K modifyer
	bool is_K_mod = false;
	std::vector<float> K_mod;

	// -> Labels
	// --> Labels are for any kind of parameter/data belonging to a labelisable group
	// --> For example data function of lithology: the user provide an array of label indices, and each indice has a set of property
	std::vector<Label> labels;
	// --> Array of nnodes size with the label indices correponding to the labels array
	std::vector<int> labelarray;

	// Analytical model options
	// -> fluvial: 0 = deactivated; 1: SPL, + (in the future more options with probably SPACE and the extended SPL)
	int fluvial_mode = 1; 

	// -> hillslope: 0 = deactivated; 1: snap to critical slope, + (in the future more options with probably geometrical stuff)
	int hillslope_mode = 1; 

	// Optimisation
	// -> flowptimode: optimisation level based on flow: stops the solver above a certain flow (see doc or manuscript for the details)
	// -->  0 is no flow based optimisation, 1 is based on flow distance and 2 on chi
	int flowptimode = 0;
	// -> Determine is the threshold for flow is normalised per basin or not
	bool flowptimode_normalised = true;
	// -> minimum area (or discharge) for the basin for the flow opti to be applied
	float flow_opti_min_A = 1e5;
	// -> minimum flow (in flow distance or chi) to be used as threshold
	float flow_opti_threshold = 0.5;

	// float timestep
	float timestep = 1000;

	// TExturing for CGs
	std::vector<int> pixeltype;

	// experimental features
	float Vwidth_coeff = 1e-2;
	float Vwidth_exp = 0.8;
	float Acritbeach = 1e6;


	minilem(){;};
	~minilem(){;};


	// This function construct the graph structure
	void initialise_graph(int nx, int ny, float dx, float dy, std::string bound)
	{
		// Reinitialising the graph
		this->graph = Graph();
		
		// Setting the dimensions, as the name may suggest
		this->graph.set_dimensions(nx, ny, nx * ny, dx, dy, 0., 0.);

		// Boundary condistions
		this->boundary_string = bound;
		this->graph.set_default_boundaries(bound);
	}


	void init_topo(float intensity, std::string method = "white_noise")
	{
		this->graph.topography = std::vector<float>(this->graph.nnodes,0.);
		this->add_noise(intensity);
	}


	// add random noise to the topography
	void add_noise(float intensity)
	{
		std::srand(std::time(0));
		for(int i = 0; i < this->graph.nnodes; ++i)
		{
			if(this->graph.boundary[i] > 0 )
				this->graph.topography[i] += static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/intensity));
		}
	}

	void add_perlin_noise(int dimension, int octaves, float xscale, float yscale, float scale)
	{
		siv::BasicPerlinNoise<float> perlin;

		for(int i = 0; i < this->graph.nnodes; ++i)
		{
			float x,y;
			this->graph.XY_from_nodeid(i,x,y);
			this->graph.topography[i] += perlin.octave2D_01((x * xscale), (y * yscale), octaves) * scale;
		}


	}


	void double_resolution()
	{
		int newnxy = this->graph.nnodes * 4;
		std::vector<float> newtopo(newnxy, 0.);
		for(int col=0; col < this->graph.nx; ++col)
		{
			for(int row=0; row < this->graph.ny; ++row)
			{
				int i = this->graph.nx * row + col;
				int jcol = col * 2;
				int jrow = row * 2;
				int j = (this->graph.nx * 2) * jrow + jcol;
				newtopo[j] = this->graph.topography[i];
				newtopo[j+1] = this->graph.topography[i];
				newtopo[j+this->graph.nx*2] = this->graph.topography[i];
				newtopo[j+1+this->graph.nx*2] = this->graph.topography[i];
			}
		}
		this->initialise_graph(this->graph.nx * 2, this->graph.ny * 2, this->graph.dx/2, this->graph.dy/2, this->boundary_string);
		this->graph.topography = newtopo;
	}

	void blur(float sigma){this->graph.topography = On_gaussian_blur(sigma, this->graph.topography, this->graph.nx, this->graph.ny); this->graph.set_boundaries_to(0.);}

	void node_solve_spl(int tnode, int trec, float tm, float tn, float tK, float tA, float tVM)
	{
		this->graph.topography[tnode] = this->graph.topography[trec] + this->graph.distance2receivers[tnode] * std::pow(tVM / (tK * std::pow(tA,tm) ),1/tn);
	}

	void node_solve_critical_slope(int tnode, int trec, float tSc)
	{
		this->graph.topography[tnode] = this->graph.topography[trec] + this->graph.distance2receivers[tnode] * tSc;
	}

	void node_solve_linear_diffusion(int tnode, int trec, float eu, float alpha, std::vector<float>& distance2river)
	{
		this->graph.topography[tnode] = this->graph.topography[trec] + std::pow(distance2river[tnode],2) * eu/alpha;
	}

	void node_solve_beach(int tnode, int trec, Label& tlab)
	{
		this->graph.topography[tnode] = this->graph.topography[trec] + this->graph.distance2receivers[tnode] * random_number(tlab.beach_slope_min, tlab.beach_slope_max);
	}

	void node_solve_floodplain(int tnode, int trec, float floodplain_slope)
	{
		this->graph.topography[tnode] = this->graph.topography[trec] + this->graph.distance2receivers[tnode] * floodplain_slope;
	}

	void solve_analytically()
	{

		// 1) compute graph
		this->graph.compute_graph(this->minima_solver);

		// 2) calculate area/discharge according to plan
		if(qamode == 0)
			this->graph.calculate_area();
		else if(qamode == 1)
			this->graph.calculate_discharge_uniprec(this->uPrec);
		else if(qamode == 2)
			this->graph.calculate_discharge_prec(this->prec);

		// 3) preparing for metrics for the optimisation
		std::vector<float> fdo;
		std::vector<float> maxAper_basins;
		std::vector<float> distance2river;
		std::vector<float> Scs;
		if(this->hillslope_mode ==2)
		{
			this->graph.d_sources(this->labels[0].Acrit);
			this->graph.compute_river_nodes();
			distance2river = this->graph.compute_graph_distance_from_rivers();
			// distance2river = this->graph.compute_distance_to_nodes(this->graph.river_nodes,1e32);
		}

		else if (this->hillslope_mode == 3)
		{
			Scs = std::vector<float>(this->graph.nnodes,0.2);
			// for(int i=0; i<this->graph.nnodes; ++i)
			// {
			// 	int label = (this->paramode == 0) ? 0:this->labelarray[i];
			// 	Label& tlab = this->labels[label];
			// 	Scs[i] = tlab.S_c;
			// }
		}

		std::vector<bool> needreplace(this->graph.nnodes, false);
		bool need_basin_label = false;
		if(this->flowptimode > 0)
		{
			// case if normalised
			if(this->flowptimode_normalised)
			{
				need_basin_label = true;
				this->graph.compute_all_basins();
				fdo = (flowptimode == 2) ? this->graph.get_normalised_chi_by_basin(0.45,1): this->graph.get_normalised_flow_distance_by_basin();
			}
			// Cas if not normalised
			else
				fdo = (flowptimode == 1) ? this->graph.compute_flow_distance_full(): this->graph.compute_chi_full(0.45,1);

			// min basin size for opti:
			maxAper_basins = this->graph.get_maxA_per_basins();
		}

		// Keeping track of maximum elevation
		float maxZ = std::numeric_limits<float>::min();


		// Let's go for the main iteration
		for(auto node:this->graph.stack)
		{
			// Getting the single flow receiver
			int rec = this->graph.receivers[node];
			// If needed, getting the basin label, for opti testing
			int baslab;
			if(need_basin_label)
				baslab = this->graph.basin_labels[node];

			// testing the flow optimisation
			bool process = true;
			if(this->flowptimode > 0)
			{
				if(fdo[node]>this->flow_opti_threshold)
					if(need_basin_label && maxAper_basins[baslab] >= flow_opti_min_A)
						process = false;
			}

			// Stpping here if flow opti or base level checkers tells so
			if(rec == node || this->graph.can_flow_out_there(node))
			{
				// this->graph.topography[node] = 0;
				continue;
			}

			if(process == false)
			{
				needreplace[node] = true;
				continue;
			}

			// ok let's solve analytically then

			// first getting the current label
			int label = (this->paramode == 0) ? 0:this->labelarray[node];
			Label& tlab = this->labels[label];

			// Getting the local vertical motion to equilibrate to
			float tU = (this->vmmode == 0)? this->uvertmot:this->vertmot[node];

			bool hillslopeSc = false;
			bool hillslopeLin = false;
			bool hillslopeExp = false;
			bool fluvialSPL = false;
			// checking which stuff to process
			if(this->hillslope_mode == 1)
			{
				if(this->graph.area[node] < tlab.Acrit)
					hillslopeSc = true;
			}
			else if (this->hillslope_mode == 2)
			{
				if(this->graph.area[node] < tlab.Acrit)
					hillslopeLin = true;
			}
			else if(this->hillslope_mode == 3)
			{
				if(this->graph.area[node] < tlab.Acrit)
					hillslopeExp = true;
			}

			if(this->fluvial_mode > 0)
			{
				fluvialSPL = true;
			}

			// process hillslope if hillslopes
			if(hillslopeSc)
			{
				this->node_solve_critical_slope(node, rec, tlab.S_c);
			}
			else if(hillslopeLin)
			{
				this->node_solve_linear_diffusion(node, rec, tU, tlab.alpha_D, distance2river);
			}
			else if(hillslopeExp)
			{
				auto tSc = Scs[rec];
				if(tSc <= 0.6)
					tSc += 0.001;
				else
					tSc -= 0.001;
				this->node_solve_critical_slope(node, rec, tSc);
				Scs[node] = tSc;
			}
			// process fluvial if fluvial
			else if (fluvialSPL)
			{
				float this_K = (is_K_mod)? this->K_mod[node] : tlab.Kref ;
				this->node_solve_spl(node, rec, tlab.mref, tlab.nref, this_K, this->graph.area[node], tU);
			}
			
			// Save max topography if needed
			if(this->graph.topography[node] > maxZ)
				maxZ = this->graph.topography[node];
		}

		// Finally, if flow opti mode is the right one
		if(this->flowptimode > 0)
		{
			std::srand(std::time(0));
			for (int i=0; i<this->graph.nnodes; ++i)
			{
				if(needreplace[i])
					this->graph.topography[i] = maxZ + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
			}
		}

	}

	// void experimental_downstuff()

	void experimental_island(float beach_distance, float blurA)
	{

		// 1) compute graph
		this->graph.compute_graph(this->minima_solver);

		// 2) calculate area/discharge according to plan
		if(qamode == 0)
			this->graph.calculate_area();
		else if(qamode == 1)
			this->graph.calculate_discharge_uniprec(this->uPrec);
		else if(qamode == 2)
			this->graph.calculate_discharge_prec(this->prec);

		std::vector<bool> done(this->graph.nnodes,false);

		// this->graph.d_sources(rivAth);
		// this->graph.compute_river_nodes();
		// std::vector<float> dist2riv = this->graph.compute_graph_distance_from_rivers();
		// std::vector<float> dist2riv = this->graph.compute_distance_to_nodes(this->graph.river_nodes,1e32);
		// std::vector<float> area_linked = this->graph.propagate_from_rivers(this->graph.area);
		// auto izFP = this->is_potential_floodplain();
		std::vector<float> dist2sea = this->graph.compute_distance_to_boundaries(1e32);

		// first getting the current label

		if(blurA > 0)
		{
			auto temp = On_gaussian_blur(blurA, this->graph.area,this->graph.nx, this->graph.ny);
			this->graph.area = temp;
		}
		
		// Pixel type: 0 = sea
		// 1 = river
		// 2 = floodplain
		// 3 = beach
		// 4 = colluvial
		// 5 = hillslope
		this->pixeltype = std::vector<int>(this->graph.nnodes,0);
		

		for(int i=0; i<this->graph.nnodes; ++i)
		{
			// Getting the single flow receiver
			int rec = this->graph.receivers[i];

			int label = (this->paramode == 0) ? 0:this->labelarray[i];
			Label& tlab = this->labels[label];


			// Stpping here if flow opti or base level checkers tells so
			if(rec == i || this->graph.can_flow_out_there(i))
			{
				pixeltype[i] = 0;
				continue;
			}

			if(this->graph.area[i] > tlab.Acrit)
			{
				if(dist2sea[i] < beach_distance && this->graph.area[i] < this->Acritbeach)
				{
					pixeltype[i] = 3;
					continue;
				}	
				pixeltype[i] = 1;
				continue;
			}

			// Checking if beach
			if(dist2sea[i] < beach_distance)
			{
				pixeltype[i] = 3;
				continue;
			}

			// checking if floodplain
			// if(izFP[i] )
			// {
			// 	// std::cout << i << std::endl;
			// 	pixeltype[i] = 2;
			// 	continue;
			// }

			

			pixeltype[i] = 5;


		}

		normalise_vector(dist2sea);


		// Keeping track of maximum elevation
		float maxZ = std::numeric_limits<float>::min();


		// Let's go for the main iteration
		for(auto node:this->graph.stack)
		{
			// Getting the single flow receiver
			int rec = this->graph.receivers[node];
			// first getting the current label
			int label = (this->paramode == 0) ? 0:this->labelarray[node];
			Label& tlab = this->labels[label];

			// Getting the local vertical motion to equilibrate to
			float tU = (this->vmmode == 0)? this->uvertmot:this->vertmot[node];

			bool hillslopeSc = false;
			bool fluvialSPL = false;

			// Stpping here if flow opti or base level checkers tells so
			if(pixeltype[node] == 0)
			{
				this->graph.topography[node] = 0;
				continue;
			}
			// ok let's solve analytically then
			else if(pixeltype[node] == 1)
			{
				this->node_solve_spl(node, rec, tlab.mref, tlab.nref, tlab.Kref / dist2sea[node], this->graph.area[node], tU);
				continue;
			}
			// else if(pixeltype[node] == 2)
			// {
			// 	this->node_solve_spl(node, rec, tlab.mref, tlab.nref, tlab.Kref * 1e2, this->graph.area[node], tU);
			// 	continue;
			// }
			else if(pixeltype[node] == 3)
			{
				this->node_solve_beach(node, rec, tlab);
				continue;
			}
			else if(pixeltype[node] == 5)
			{
				this->node_solve_critical_slope(node, rec, tlab.S_c);
				continue;
			}
			
			// Save max topography if needed
			if(this->graph.topography[node] > maxZ)
				maxZ = this->graph.topography[node];
		}



	}

	std::vector<bool> is_potential_floodplain()
	{
		std::vector<bool> izFP(this->graph.nnodes,false);
		std::vector<float> distfromstuff(this->graph.nnodes,0.);
		std::priority_queue<PQ_helper<int, float>, std::vector<PQ_helper<int, float>>, std::less<PQ_helper<int, float> > > PQ;
		// std::vector<float> dist2riv = this->graph.compute_distance_to_nodes(this->graph.river_nodes,1e32);


		for(auto i:this->graph.river_nodes)
			PQ.emplace(PQ_helper<int, float>(i,this->Vwidth_coeff * std::pow(this->graph.area[i], this->Vwidth_coeff)));



		while(PQ.empty() == false)
		{
			auto next = PQ.top();PQ.pop();

			if(next.score > distfromstuff[next.node])
			{
				izFP[next.node] = true;
				distfromstuff[next.node] = next.score;
				auto neighbours = this->graph.get_neighbours(next.node);
				float td = next.score;
				for(auto& tn:neighbours)
				{
					if( distfromstuff[tn.node] < td - tn.distance )
						PQ.emplace( PQ_helper<int,float>(tn.node,td - tn.distance) );
				}
			}
		}

		// for (int i=0' i< this->graph.nnodes;++i)
		// 	if()
		return izFP;

	}



	void set_flow_options(std::string mode = "none", bool normalised = true, float min_A = 1e6, float threshold = 0.5)
	{
		this->flowptimode = (mode == "chi") ? 2: ((mode == "flow") ? 1:0 );
		this->flowptimode_normalised = normalised;
		this->flow_opti_min_A = min_A;
		this->flow_opti_threshold = threshold;
	}

	void set_vertical_motions_uniform(float val)
	{
		this->uvertmot = val;
		this->vmmode = 0;
	}

	void set_K_mod(std::vector<float> tK_mod)
	{
		this->is_K_mod = true;
		this->K_mod = tK_mod;
	}

	void set_vertical_motions_spatial(std::vector<float>& val)
	{
		this->vertmot = val;
		this->vmmode = 1;
	}

	void set_local_minima_option(std::string mode)
	{
		if(mode == "carve" || mode == "fill" || mode == "none")
			this->minima_solver = mode;
		else
		{
			std::cout << "Mode needs to be on of carve, fill or none" << std::endl;
			throw std::runtime_error("Mode needs to be on of carve, fill or none");
		}
	}


	void set_precipitation_uniform(float val){this->qamode = 1; this->uPrec = val;}
	void set_precipitation_spatial(std::vector<float> val){this->qamode = 1; this->prec = val;}
	void use_drainage_area(){this->qamode = 0;}


	void set_uniform_parameters(float tm, float tn, float tKref, float tAcrit, float tS_c, float alpha_D )
	{
		this->paramode = 0;
		Label lab = Label();
		lab.mref = tm;
		lab.nref = tn;
		lab.Kref = tKref;
		lab.Acrit = tAcrit;
		lab.S_c = tS_c;
		lab.alpha_D = alpha_D;
		this->labels.clear();
		this->labels.emplace_back(lab);
	}

	void use_spatial_labels(){this->paramode =1; this->labels.clear();}

	void add_label(float tm, float tn, float tKref, float tAcrit, float tS_c, float alpha_D )
	{
		Label lab = Label();
		lab.mref = tm;
		lab.nref = tn;
		lab.Kref = tKref;
		lab.Acrit = tAcrit;
		lab.S_c = tS_c;
		lab.alpha_D = alpha_D;
		this->labels.emplace_back(lab);
	}

	void set_hillslope_mode(std::string mode)
	{
		if(mode == "none")
			this->hillslope_mode = 0;
		else if (mode == "critical_slope")
			this->hillslope_mode = 1;
		else if (mode == "linear")
			this->hillslope_mode = 2;
		else if (mode == "experimental")
			this->hillslope_mode = 3;
		else
			throw std::runtime_error("Unknow hillslope mode, needs to be one of 'none' or 'critical_slope' .");

	}

	void set_fluvial_mode(std::string mode)
	{
		if(mode == "none")
			this->fluvial_mode = 0;
		else if (mode == "SPL")
			this->fluvial_mode = 1;
		else
			throw std::runtime_error("Unknow fluvial mode, needs to be one of 'none' or 'SPL' .");

	}

	void set_timestep(float tt){this->timestep = tt;}

	void run_SPL_basic()
	{

		// 1) compute graph
		this->graph.compute_graph(this->minima_solver);

		// 2) calculate area/discharge according to plan
		if(qamode == 0)
			this->graph.calculate_area();
		else if(qamode == 1)
			this->graph.calculate_discharge_uniprec(this->uPrec);
		else if(qamode == 2)
			this->graph.calculate_discharge_prec(this->prec);

		Label& tlab = this->labels[0];
		for (auto tnode:this->graph.stack)
		{
			if(this->graph.can_flow_out_there(tnode))
				continue;

			// Newton-rahpson method to solve non linear equation, see Braun and Willett 2013
			float epsilon;

			// A bit of factorisation to clarify the equation
			float this_K = (is_K_mod)? this->K_mod[tnode] : tlab.Kref ;
			float streamPowerFactor = this_K * pow(this->graph.area[tnode], tlab.mref) * this->timestep;
			float slope; // Slope.
			float new_zeta = this->graph.topography[tnode]; float old_zeta = this->graph.topography[tnode]; // zeta = z = elevation
			float dx = this->graph.distance2receivers[tnode];
			// iterate until you converge on a solution. Uses Newton's method.
			int nit = 0;
			do
			{
				// Get the slope
				slope = (new_zeta - this->graph.topography[this->graph.receivers[tnode]]) / dx;
				// Check backslope or no slope ie no erosion
				if(slope <= 0)
				{
					epsilon = 0;
				}
				else
				{
					// Applying the newton's method
					epsilon = (new_zeta - old_zeta + streamPowerFactor * std::pow(slope, tlab.nref)) /
					(1 + streamPowerFactor * (tlab.nref/dx) * std::pow(slope, tlab.nref - 1));
				}
				// iterate the result
				old_zeta = new_zeta;
				new_zeta -= epsilon;
				++nit;
				if(nit > 100)
				{
					std::cout << dx << "||" << slope << "||" << streamPowerFactor << "||" << std::pow(slope, tlab.nref - 1) << std::endl;
					epsilon = 0;	
				}
				// I want it to run while it can still have an effect on the elevation
			} while (abs(epsilon) > 1e-3);

			this->graph.topography[tnode] = new_zeta;

		}
	
	}

	void run_cidre_hillslope_basic(float kappa, float maxD, float oSc)
	{
		std::vector<float> Qs(this->graph.nnodes,0);

		for (int i = this->graph.nnodes - 1 ; i >= 0 ; --i)
		{
			int tnode = this->graph.stack[i];
			if(this->graph.can_flow_out_there(tnode) || this->graph.can_flow_even_go_there(tnode) == false)
				continue;
			int rec = this->graph.receivers[tnode];
			float tS = std::fmax(1e-6,this->graph.gradient(tnode));
			float dx = this->graph.distance2receivers[tnode];
			
			bool isbellowSc = (tS < oSc - 1e-6);

			 

			float L = dx / (1 - std::pow(tS/oSc,2) );

			float D = std::fmin(isbellowSc ? ((L > dx) ?  Qs[tnode]/L : Qs[tnode] / dx) : 0, maxD * Qs[tnode]/dx);
			float E = isbellowSc ? kappa * tS: (std::abs(tS - oSc) * dx)/this->timestep;

			if(E > 1 || D > 1)
				std::cout << E << "||" << D << "||" << D * dx << " vs " << Qs[tnode] << "||" << Qs[tnode] << "||" << L  << " = " <<dx << "/" << (1 - std::pow(tS/oSc,2) ) << std::endl;

			Qs[rec] += Qs[tnode] + (E - D) * dx;
			this->graph.topography[tnode] += (D - E) * this->timestep;
		}

	}

	void run_SPL_basic_for_hydraulic_erosion(float kcof, std::vector<bool>& mask)
	{

		// 1) compute graph
		this->graph.compute_graph(this->minima_solver);

		// 2) calculate area/discharge according to plan
		if(qamode == 0)
			this->graph.calculate_area();
		else if(qamode == 1)
			this->graph.calculate_discharge_uniprec(this->uPrec);
		else if(qamode == 2)
			this->graph.calculate_discharge_prec(this->prec);

		Label& tlab = this->labels[0];
		for (auto tnode:this->graph.stack)
		{
			if(this->graph.can_flow_out_there(tnode) || this->graph.can_flow_even_go_there(tnode) == false || mask[tnode] == false)
				continue;

			// Newton-rahpson method to solve non linear equation, see Braun and Willett 2013
			float epsilon;

			// A bit of factorisation to clarify the equation
			float this_K = (is_K_mod)? this->K_mod[tnode] : tlab.Kref ;
			this_K *= kcof;
			float streamPowerFactor = this_K * pow(this->graph.area[tnode], tlab.mref) * this->timestep;
			float slope; // Slope.
			float new_zeta = this->graph.topography[tnode]; float old_zeta = this->graph.topography[tnode]; // zeta = z = elevation
			float dx = this->graph.distance2receivers[tnode];
			// iterate until you converge on a solution. Uses Newton's method.
			int nit = 0;
			do
			{
				// Get the slope
				slope = (new_zeta - this->graph.topography[this->graph.receivers[tnode]]) / dx;
				// Check backslope or no slope ie no erosion
				if(slope <= 0)
				{
					epsilon = 0;
				}
				else
				{
					// Applying the newton's method
					epsilon = (new_zeta - old_zeta + streamPowerFactor * std::pow(slope, tlab.nref)) /
					(1 + streamPowerFactor * (tlab.nref/dx) * std::pow(slope, tlab.nref - 1));
				}
				// iterate the result
				old_zeta = new_zeta;
				new_zeta -= epsilon;
				++nit;
				if(nit > 100)
				{
					std::cout << dx << "||" << slope << "||" << streamPowerFactor << "||" << std::pow(slope, tlab.nref - 1) << std::endl;
					epsilon = 0;	
				}
				// I want it to run while it can still have an effect on the elevation
			} while (abs(epsilon) > 1e-3);

			this->graph.topography[tnode] = new_zeta;

		}
	
	}

	void apply_vertical_motions()
	{
		for(int i=0; i<this->graph.nnodes; ++i)
		{
			if(this->graph.can_flow_out_there(i))
				continue;
			float tU = (this->vmmode == 0)? this->uvertmot:this->vertmot[i];
			this->graph.topography[i] += tU * this->timestep;
		}
	}


	void init_circle_boundaries(float x_centre, float y_centre, float radius, int boundtype)
	{
		for(int row = 0; row < this->graph.ny; ++row)
		{
			for(int col = 0 ; col < this->graph.nx; ++col)
			{
				float tX,tY;
				this->graph.XY_from_rowcol(row, col,  tX,  tY);
				if(std::pow(tX - x_centre,2) + std::pow(tY - y_centre,2) > std::pow(radius,2))
					this->graph.boundary[this->graph.nodeid_from_row_col(row,col)] = boundtype;
				if(row ==0 ||row == this->graph.ny-1 || col == 0 || col == this->graph.nx-1)
					this->graph.boundary[this->graph.nodeid_from_row_col(row,col)] = boundtype;
			}
		}

		if (boundtype < 0)
		{
			for(int i=0; i<this->graph.nnodes; ++i)
			{
				if(this->graph.boundary[i] == 1)
				{
					auto neighbours = this->graph.get_D4_neighbours_only_id(i);
					for (auto tn:neighbours)
					{
						if(this->graph.boundary[tn] == -1)
						{
							this->graph.boundary[i] = 0;
							break;
						}
					}
				}
			}

		}

	}


	void init_boundaries_from_binary_array(std::vector<bool> barr, bool edges_base_level = true)
	{
		// Error checking
		if(barr.size() != this->graph.nnodes_t)
			throw std::runtime_error("The input boundaries need to be the same dimention as the graph");

		// Init all the boundaries to -1 (no data, no flow)
		this->graph.boundary = std::vector<int>(this->graph.nnodes,-1);
		
		// Cheecking every nodes
		for(int i=0; i<this->graph.nnodes; ++i)
		{
			if(barr[i] == false)
			{
				// if mask is false: no flow
				this->graph.boundary[i] = -1;
				this->graph.topography[i] = 0;
				continue;
			}
			// Else normal flow
			this->graph.boundary[i] = 1;

			// If I need to check dem edges, I assign base level if the node is on dem edge
			if(edges_base_level)
			{
				if(this->graph.is_on_dem_edge(i))
				{
					this->graph.boundary[i] = 0;
					continue;
				}
			}

			// Otherwise, if any of the neighbour is a no data
			auto neighbours = this->graph.get_neighbours_only_id(i);

			for(auto tn:neighbours)
			{
				if(barr[tn] == false)
				{
					// Then this node is a base level
					this->graph.boundary[i] = 0;
					// I do not need to check other neighbours so I break the loop
					break;
				}
			}
		}
	}

	void burn_data_to_base_levels(std::vector<float> data)
	{
		for(int i=0; i< this->graph.nnodes; ++i)
		{
			if(this->graph.can_flow_out_there(i))
				this->graph.topography[i] = data[i];
		}
	}

	std::vector<int> get_pixeltype(){return this->pixeltype;}

	void set_Vwidth_params(float A, float B){this->Vwidth_coeff = A; this->Vwidth_exp = B;};
	void set_Acritbeach(float A){this->Acritbeach = A;}


	// void hydraulic_erosion_graph_sources(float coefK, int n_sources)
	// {
	// 	std::vector<int> sources; sources.reserve(this->graph.nnodes);
	// 	for(int i=0; i<this->graph.nnodes; ++i)
	// 	{
	// 		if(this->graph.can_flow_even_go_there(i) == false)
	// 			continue;
	// 		if(this->graph.donors[i].size() == 0)
	// 			sources.emplace_back(i);
	// 	}

	// 	sources.shrink_to_fit();
	// 	std::vector<bool> source_done(sources.size(),false), need_proc(this->graph.nnodes, false);

	// 	for(int nsamples = 0; nsamples < n_sources; ++nsamples)
	// 	{
	// 		while(true)
	// 		{
	// 			int randpos = floor(static_cast <float> (rand()) / (static_cast <float> (RAND_MAX))  * int( sources.size() ) );
	// 			if(source_done[randpos])
	// 				continue;
				
	// 			source_done[randpos] = true;
	// 			need_proc[sources[randpos]] = true;
	// 			break;
	// 		}
	// 	}

	// 	for(int i =  this->graph.nnodes - 1; i >= 0; --i)
	// 	{
	// 		int tnode = this->graph.stack[i];
	// 		if(need_proc[tnode] == false)
	// 			continue;
	// 		int trec = this->graph.receivers[tnode];
	// 		need_proc[trec] = true;
	// 	}

	// 	this->run_SPL_basic_for_hydraulic_erosion(coefK, need_proc);

	// }

	// void hydraulic_erosion_graph_all(float coefK, int n_sources)
	// {

	// 	std::vector<bool> need_proc(this->graph.nnodes, false);

	// 	for(int nsamples = 0; nsamples < n_sources; ++nsamples)
	// 	{
	// 		while(true)
	// 		{
	// 			int randpos = floor(static_cast <float> (rand()) / (static_cast <float> (RAND_MAX))  * int( this->graph.nnodes ) );
	// 			if(need_proc[randpos])
	// 				continue;
				
	// 			need_proc[randpos] = true;
	// 			break;
	// 		}
	// 	}

	// 	for(int i =  this->graph.nnodes - 1; i >= 0; --i)
	// 	{
	// 		int tnode = this->graph.stack[i];
	// 		if(need_proc[tnode] == false)
	// 			continue;
	// 		int trec = this->graph.receivers[tnode];
	// 		need_proc[trec] = true;
	// 	}

	// 	this->run_SPL_basic_for_hydraulic_erosion(coefK, need_proc);

	// }



	// void simple_hydraulic_erosion(int n_drops, float tdt, float friction, float scale, float depositionRate, float evapRate, float k_err)
	// {
	// 	int n_break = 0;
	// 	for (int i=0; i < n_drops; ++i)
	// 	{
	// 		// Determining random position
	// 		// std::srand(std::time(0));
	// 		int randpos = floor(static_cast <float> (rand()) / (static_cast <float> (RAND_MAX))  * this->graph.nnodes);

	// 		// std::cout << randpos << std::endl;
	// 		// creating particles
	// 		wapart<float,int> tparticle(this->graph, randpos);
	// 		// std::cout << tparticle.check_if_in_dem() << std::endl;

	// 		// running the while loop
	// 		int n_speed_0 = 0;

	// 		while(tparticle.check_if_in_dem())
	// 		{
	// 			int temp_i = tparticle.i;
	// 			// std::cout << "1" << std::endl;
	// 			tparticle.compute_speed(tdt, friction, scale);

	// 			if(temp_i ==tparticle.i )
	// 				n_speed_0++;

	// 			if(this->graph.is_on_dem_edge(tparticle.i))
	// 				continue;
	// 			// std::cout << tparticle.i << "|";
	// 			// std::cout << "2" << std::endl;
	// 			// float dz = this->graph.topography[temp_i] - this->graph.topography[tparticle.i];
	// 			tparticle.compute_c_eq(this->graph.topography[temp_i] - this->graph.topography[tparticle.i]);
	// 			// std::cout << dz << " vs " << tparticle.c_eq << "||";
	// 			// std::cout << "3" << std::endl;
	// 			tparticle.compute_ED(tdt, depositionRate, k_err);
	// 			// std::cout << "4" << std::endl;
	// 			tparticle.evaporate(tdt, evapRate);
	// 			// std::cout << "5" << std::endl;
	// 			// if(tparticle.speed.x + tparticle.speed.y < 1e-3)
	// 			// 	n_speed_0 ++;
	// 			if(n_speed_0 > 10)
	// 			{
	// 				n_break++;
	// 				break;
	// 			}
	// 			if(tparticle.volume < 0.05)
	// 				break;
	// 		}
			
	// 	}

	// 	std::cout << n_break << " were stopped in the middle " << std::endl;
	// }

	// void simple_hydraulic_erosion_v2(int n_drops, float k_err, float exp_err, float k_dep, float exp_dep, float factor)
	// {
	// 	int n_break = 0;
	// 	std::vector<float> topo(this->graph.topography);
	// 	for (int i=0; i < n_drops; ++i)
	// 	{
	// 		// Determining random position
	// 		// std::srand(std::time(0));
	// 		int randpos = floor(static_cast <float> (rand()) / (static_cast <float> (RAND_MAX))  * this->graph.nnodes);

	// 		// std::cout << randpos << std::endl;
	// 		// creating particles
	// 		wapart2<float,int> tparticle(this->graph, randpos);
	// 		// std::cout << tparticle.check_if_in_dem() << std::endl;
	// 		int n_loop = 0;
	// 		while(true)
	// 		{
	// 			n_loop++;
	// 			// std::cout << tparticle.sediments << "|";
	// 			bool test = tparticle.compute_sloperec(topo);
	// 			if(this->graph.is_on_dem_edge(tparticle.i) || test == false)
	// 			{
	// 				// if(test == false && this->graph.is_on_dem_edge(tparticle.i) == false)
	// 				// 	tparticle.fill();

	// 				break;
	// 			}
	// 			tparticle.mass_transfer(k_err ,exp_err,k_dep,exp_dep,factor);
	// 			tparticle.move();
	// 		}
	// 		// std::cout << std::endl;
	// 	}
	// }


	// void simple_hydraulic_erosion_v3(int n_drops, float K, float tm, float tn, float L, float factor)
	// {
	// 	int n_break = 0;
	// 	std::vector<float> topo(this->graph.topography);
	// 	for (int i=0; i < n_drops; ++i)
	// 	{
	// 		// Determining random position
	// 		// std::srand(std::time(0));
	// 		int randpos = floor(static_cast <float> (rand()) / (static_cast <float> (RAND_MAX))  * this->graph.nnodes);

	// 		// std::cout << randpos << std::endl;
	// 		// creating particles
	// 		wapart2<float,int> tparticle(this->graph, randpos);
	// 		// std::cout << tparticle.check_if_in_dem() << std::endl;
	// 		int n_loop = 0;
	// 		while(true)
	// 		{
	// 			n_loop++;
	// 			// std::cout << tparticle.sediments << "|";
	// 			bool test = tparticle.compute_sloperec(topo);
	// 			if(this->graph.is_on_dem_edge(tparticle.i) || test == false)
	// 			{
	// 				// if(test == false && this->graph.is_on_dem_edge(tparticle.i) == false)
	// 				// 	tparticle.fill();

	// 				break;
	// 			}
	// 			// tparticle.mass_transfer(k_err ,exp_err,k_dep,exp_dep,factor);
	// 			tparticle.mass_transfer_SPL( K,  tm,  tn,  L, factor);
	// 			tparticle.move();
	// 		}
	// 		// std::cout << std::endl;
	// 	}
	// }





	
};


























#endif