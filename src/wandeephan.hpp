//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef wandeephan_HPP
#define wandeephan_HPP

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



template<class T, class U>
class wandeephan
{

public:

U nnodes = 0;

T dx = 0;
T Qw = 0;
T Ks_gravel = 0;
T Ks_sand = 0;
T deplength_gravel = 0;
T deplength_sand = 0;
T dt = 1000;
T m_exp = 0.45;
T n_exp = 1;
T width = 50;
T kp = 0.6;
T pp = 2;
T dz = 1;

std::vector<T> Z, Qs_gravel, Qs_sand, Vm, hs;
std::vector<std::vector<T> > strat;



wandeephan(){};
~wandeephan(){};


void init(U nnodes, 
	std::vector<T>& topo, 
	std::vector<T>& Vm, 
	T Qs_in_gravel, 
	T Qs_in_sand, 
	T Qw_in, 
	T dx, 
	T Ks_gravel, 
	T Ks_sand, 
	T m, 
	T n, 
	T deplength_gravel, 
	T deplength_sand, 
	T width
	)
{
	this->Z = topo;
	this->Vm = Vm;
	this->nnodes = nnodes;
	this->Qs_gravel = std::vector<T>(nnodes,0.);
	this->Qs_sand = std::vector<T>(nnodes,0.);
	this->Qs_gravel[0] = Qs_in_gravel;
	this->Qs_sand[0] = Qs_in_sand;
	this->Qw = Qw_in;
	this->Ks_gravel = Ks_gravel;
	this->Ks_sand = Ks_sand;
	this->dx = dx;
	this->m_exp = m;
	this->n_exp = n;
	this->deplength_gravel = deplength_gravel;
	this->deplength_sand = deplength_sand;
	this->width = width;
	this->hs = std::vector<T>(nnodes, 0.);
	this->strat = std::vector<std::vector<T> >(nnodes);
	// for(int i=0; i<this->nnodes; ++i)
		// this->strat[i].emplace_back(0.5);
}

void update_Qs_in(T Qs_gravel_in, T Qs_sand_in){this->Qs_gravel[0] = Qs_gravel_in; this->Qs_sand[0] = Qs_sand_in;};
void update_Qw_in(T Qw_in){this->Qw = Qw_in;};
void update_dt(T dt){this->dt = dt;};;
std::vector<T> get_topo(){return this->Z;};
std::vector<T> get_Qs_gravel(){return this->Qs_gravel;};
std::vector<T> get_Qs_sand(){return this->Qs_sand;};
std::vector< std::vector<T> > get_full_prop(){return this->strat;};
std::vector<T> get_surface_prop()
{
	std::vector<T> out(this->nnodes,0.);
	for(int i=0; i < this->nnodes; ++i)
	{
		if(this->strat[i].size() > 0)
			out[i] = this->strat[i].back();
	}
	return out;
};

void run()
{
	T Qw2use = this->Qw;
	T increments = this->kp * std::pow(this->dx, this->pp);

	for (size_t i=0; i < this->nnodes - 1; ++i)
	{
		// this->width = 9.5 * std::pow(Qw2use,0.35);
		// Getting hte topographic gradient and recasting it to 0 if negative
		T S = (this->Z[i] - this->Z[i+1])/this->dx;
		if(S<0) S = 0;

		// Calculating the stream power
		T SP = std::pow(Qw2use, this->m_exp) * std::pow(S, this->n_exp);
		T used = 0;

		// Sediment entrainment
		T Es_gravel = this->Ks_gravel * SP;
		T Es_sand = this->Ks_sand * SP;

		T Etot = (Es_gravel + Es_sand)/2;
		T prop_gravel_tot = 0;

		// while(true)
		// {
		// 	T prop_gravel = 0;
		// 	if(this->strat[i].size() == 0)
		// 	{
		// 		prop_gravel = 0;
		// 		T this_E = (Es_sand * (1 - prop_gravel) + Es_gravel * prop_gravel) * (1 -  used);
		// 		prop_gravel_tot = (this_E + Etot > 0) ? (this_E * prop_gravel + Etot * prop_gravel_tot)/(this_E + Etot) : 0;

		// 		Etot += this_E;
		// 		this->hs[i] = 0;
		// 		break;
		// 	}

		// 	prop_gravel = this->strat[i].back();
		// 	T rest = remainderf(this->hs[i], this->dz);
		// 	T this_E = (Es_sand * (1 - prop_gravel) + Es_gravel * prop_gravel) * (1 - used);
		// 	if(this_E > rest && this_E > 0)
		// 	{
		// 		used += rest/this_E;
		// 		this_E = rest;
		// 		this->strat[i].pop_back();
		// 	}
		// 	else
		// 	{
		// 		used = 1;
		// 	}


		// 	prop_gravel_tot = (this_E + Etot > 0)? (this_E * prop_gravel + Etot * prop_gravel_tot)/(this_E + Etot):0;
		// 	Etot += this_E;
		// 	this->hs[i] -= this_E * this->dt;

		// 	if(used >= 1)
		// 		break;

		// }

		// and deposition
		T Ds_gravel = this->Qs_gravel[i]/(this->Qw * this->deplength_gravel);
		T Ds_sand = this->Qs_sand[i]/(this->Qw * this->deplength_sand);

		if(Ds_gravel * this->width * this->dx > this->Qs_gravel[i])
			Ds_gravel = this->Qs_gravel[i]/(this->width* this->dx);

		if(Ds_sand * (this->width * this->dx) > this->Qs_sand[i])
			Ds_sand = this->Qs_sand[i]/(this->width * this->dx);

		T prop_gravel_depo = Ds_gravel/(Ds_sand + Ds_gravel);

		// T rest2fill = this->dz - remainderf(this->hs[i] ,this->dz);
		// T n_other = (Ds_gravel + Ds_sand - this->dz)/this->dz;
		// if(n_other < 0)
		// 	n_other = 0;

		// if ( (Ds_gravel + Ds_sand) > rest2fill  )
		// {
		// 	if(this->strat[i].size() > 0)
		// 		this->strat[i].back() = rest2fill/this->dz * prop_gravel_depo + (1 - rest2fill)/this->dz * this->strat[i].back();
		// 	else
		// 		this->strat[i].emplace_back(prop_gravel_depo);
		// }


		// for(int i = 0; i<n_other; ++i)
		// {
		// 	this->strat[i].emplace_back(prop_gravel_depo);
		// }

		// this->hs[i] += (Ds_gravel + Ds_sand) * this->dt;


		if(i > 0)
		{
			this->Qs_gravel[i] += Etot * this->width * this->dx * prop_gravel_tot;
			this->Qs_sand[i] -=  Ds_sand * this->width * this->dx;
			this->Qs_sand[i] += Etot * this->width * this->dx * (1 - prop_gravel_tot);
			this->Qs_gravel[i] -= Ds_gravel * this->width * this->dx;

		}

		this->Qs_gravel[i + 1] = this->Qs_gravel[i];
		this->Qs_sand[i + 1] = this->Qs_sand[i];

		this->Z[i] += (-1 * Etot + (Ds_sand + Ds_gravel) ) * this->dt;

		Qw2use += increments;
	}

	for (size_t i=0; i < this->nnodes - 1; ++i)
		this->Z[i] += this->Vm[i] * this->dt;

}










};































#endif