//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef chonkutils_HPP
#define chonkutils_HPP

// STL imports
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <ctime>
#include <fstream>
#include <random>
#include <queue>
#include <iostream>
#include <numeric>
#include <cmath>
#include <initializer_list>
#include <stdio.h>
#include <string.h>



// This enum contains the different graph types
// SF -> Single Flow direction
// DInf -> Multiple flow, following the Dinf logic (only 2 recs)
// MF -> Multiple flow
enum class graphtype {SF, Dinf, MF};

// This enum contains the topology of flow directions
// D4 -> deliver to neighbors using the Queen direction
// D8 -> deliver to neighbors using the King direction
// Dangle -> not sure I need this
enum class flowdirtype {D4, D8, Dangle};

// This enum contains the type of direction to take account when calling a graph
// donors -> upstream direction
// receivers -> downstream direction
// both -> well, both
enum class graphdir {none,donors, receivers, both};

enum class when2update {dynamic, timestep, phase};

enum class CordonnierMethod {CARVE, FILL, SIMPLE, NONE};


enum class when2reinitialise {never, timestep, phase};

enum class BouType {OUT, NORM, PER, NONODE};
enum class NeiType {DON, REC, NOPE, FLAT, ALL};


// Quick function checking if an element is in a set or not
template<class T>
inline bool izinset(std::set<T>& tset, T& element )
{
	 return std::binary_search(tset.begin(), tset.end(), element);
};



// Small container class for neighbourers
// It helps the navigation when accessing a neighbouring node
// node is the node ID of the neightbours
// distance is the distance to the reference node
// Note that I keep the number of members low on purpose (no mention of original node, azimuth, elevation, ...)
// to remain efficient.
template<class n_t, class dist_t>
class Neighbour
{
	public:
		Neighbour(){;}
		Neighbour(n_t node, dist_t distance){this->node = node; this->distance = distance;}
		n_t node;
		dist_t distance;
};


// Small container class for directed neighbour (i.e receiver or donor)
// It helps the navigation when accessing a neighbouring node
// node is the node ID of the neightbours
// sid is the special ID of the node -> rid or did depending on neighbour type
// distance is the distance to the reference node
// Note that I keep the number of members low on purpose (no mention of original node, azimuth, elevation, ...)
// to remain efficient.
template<class n_t, class dist_t>
class DirectedNeighbour : public Neighbour<n_t,dist_t>
{
	public:
		DirectedNeighbour(){;}
		DirectedNeighbour(n_t node, n_t sid,  dist_t distance){this->node = node; this->sid = sid; this->distance = distance;}
		n_t node;
		n_t sid;
		dist_t distance;
};

template<class T, class U>
class PQ_helper
{
  public:
    // empty constructor
    PQ_helper(){};
    // Constructor by default
    PQ_helper(T node,U score){this->node = node; this->score = score;};
    // Elevation data
    U score;
    // Node index
    T node;
};


template<class T, class U>
inline bool operator>( const PQ_helper<T,U>& lhs, const PQ_helper<T,U>& rhs )
{
  return lhs.score > rhs.score;
}
template<class T, class U>
inline bool operator<( const PQ_helper<T,U>& lhs, const PQ_helper<T,U>& rhs )
{
  return lhs.score < rhs.score;
}

// Hack the container behind
template <class T, class S, class C>
    S& Container(std::priority_queue<T, S, C>& q) {
        struct HackedQueue : private std::priority_queue<T, S, C> {
            static S& Container(std::priority_queue<T, S, C>& q) {
                return q.*&HackedQueue::c;
            }
        };
    return HackedQueue::Container(q);
};

// Calss utilised in the marine topology (or else if needed)
// It sorts with 2 prioritise: normally by score_1 and inversely by score 2
template<class T, class U>
class DoubleSorterMarine
{
  public:
    // empty constructor
    DoubleSorterMarine(){};
    // Constructor by default
    DoubleSorterMarine(T node,U score_1, U score_2){this->node = node; this->score_1 = score_1; this->score_2 = score_2;};
    // Elevation data
    U score_1;
    U score_2;
    // Node index
    T node;
};

template<class T, class U>
inline bool operator>( const DoubleSorterMarine<T,U>& lhs, const DoubleSorterMarine<T,U>& rhs )
{
	if(lhs.score_1 != rhs.score_1)
	  return lhs.score_1 > rhs.score_1;
	else
		return lhs.score_2 < rhs.score_2;
}

template<class T, class U>
inline bool operator<( const DoubleSorterMarine<T,U>& lhs, const DoubleSorterMarine<T,U>& rhs )
{
	if(lhs.score_1 != rhs.score_1)
	  return lhs.score_1 < rhs.score_1;
	else
		return lhs.score_2 > rhs.score_2;
}




// Convert a graphdir to a string
inline std::string graphdir2string(graphdir this_dir)
{
	if(this_dir == graphdir::none)
		return "none";
	if(this_dir == graphdir::donors)
		return "donors";
	if(this_dir == graphdir::receivers)
		return "receivers";
	if(this_dir == graphdir::both)
		return "both";
	return "UNKNOWN_TYPE";
}

// Convert a graphtype to string
inline std::string graphtype2string(graphtype this_type)
{
	if(this_type == graphtype::SF)
		return "SF";
	if(this_type == graphtype::Dinf)
		return "Dinf";
	if(this_type == graphtype::MF)
		return "MF";
	return "UNKNOWN_TYPE";
}

inline graphdir string2graphdir(std::string strin)
{
    if(strin == "none")
        return graphdir::none;
    if(strin == "donors")
        return graphdir::donors;
    if(strin == "receivers")
        return graphdir::receivers;
    if(strin == "both")
        return graphdir::both;
    else
        return graphdir::none;
}


template<class T>
T random_number(T min, T max)
{ 
    T random = ((T) rand()) / (T) RAND_MAX;
    T range = max - min;  
    return (random*range) + min;
}

// Return a vector of sorted indices
// Equivalent of argsort in numpy or other platforms
// Shamelessly ripped from https://stackoverflow.com/a/12399290/7114716 - credits to Lukasz Wiklendt
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  std::stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}



inline std::vector<std::string> split_mah_string(std::string str,std::string sep)
{
    char* cstr=const_cast<char*>(str.c_str());
    char* current;
    std::vector<std::string> arr;
    current=strtok(cstr,sep.c_str());
    while(current!=NULL){
        arr.push_back(current);
        current=strtok(NULL,sep.c_str());
    }
    return arr;
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






void boxBlurH_4 (std::vector<float>& scl, std::vector<float>& tcl, int w, int h, float r) 
{
	float iarr = 1 / (r+r+1);
	for(int i=0; i<h; ++i) 
	{
		int ti = i*w, li = ti, ri = ti+r;
		float fv = scl[ti], lv = scl[ti+w-1], val = (r+1)*fv;
		for(int j=0; j<r; ++j) val += scl[ti+j];
		for(int j=0  ; j<=r ; ++j) { val += scl[ri++] - fv       ;   tcl[ti++] =std::round(val*iarr); }
		for(int j=r+1; j<w-r; ++j) { val += scl[ri++] - scl[li++];   tcl[ti++] =std::round(val*iarr); }
		for(int j=w-r; j<w  ; ++j) { val += lv - scl[li++];   tcl[ti++] =std::round(val*iarr); }
	}
}


void boxBlurT_4 (std::vector<float>& scl, std::vector<float>& tcl, int w, int h, float r) 
{
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


void boxBlur_4 (std::vector<float>& scl, std::vector<float>& tcl, int w, int h, float r) 
{
	for(size_t i=0; i<scl.size(); ++i) tcl[i] = scl[i];
	boxBlurH_4(tcl, scl, w, h, r);
	boxBlurT_4(scl, tcl, w, h, r);
}

void gaussBlur_4 (std::vector<float>& scl, std::vector<float>& tcl, float r, int nx, int ny) 
{
	auto bxs = boxesForGauss(r, 3);
	boxBlur_4 (scl, tcl, nx, ny, (bxs[0]-1)/2);
	boxBlur_4 (tcl, scl, nx, ny, (bxs[1]-1)/2);
	boxBlur_4 (scl, tcl, nx, ny, (bxs[2]-1)/2);
}

std::vector<float> On_gaussian_blur(float r, std::vector<float> topo, int nx, int ny)
{
	std::vector<float>newtopo(topo);
	gaussBlur_4 (topo, newtopo, r, nx, ny);
	return newtopo;
}



class Label
{
	public:
		float mref = 0.45;
		float nref = 1;
		float Kref = 1e-4;
		float S_c = 0.4;
		float Acrit = 1e6;
		float floodplain_slope = 0.001;
		float beach_slope_min = 0.001;
		float beach_slope_max = 0.005;
		float alpha_D = 1e-2;
		std::vector<std::vector<float> > breaks;

		Label(){};
		void add_break(float tA, float tK, float tm, float tn)
		{
			this->breaks.emplace_back(std::vector<float>({tA,tK,tm,tn}));
		}

		void set_m(float t){this->mref = t;}
		void set_n(float t){this->nref = t;}
		void set_Kref(float t){this->Kref = t;}
		void set_S_c(float t){this->S_c = t;}
		void set_Acrit(float t){this->Acrit = t;}

};



void normalise_vector(std::vector<float>& vec)
{
	float min = std::numeric_limits<float>::max();
	float max = std::numeric_limits<float>::min();

	for(auto v:vec)
	{
		if(v<min)
			min = v;
		else if(v>max)
			max = v;
	}

	for(auto&v:vec)
		v = (v-min)/(max-min);
}




// Small class emulating some aspect of a vector vector from its pointer
template<class T>
class pvector
{
public:
	std::shared_ptr<std::vector<T> > data;

	pvector(){};
	pvector(std::vector<T>& dat){this->data = std::make_shared<std::vector<T> >(dat);}
	pvector(std::vector<T>* dat){this->data = std::make_shared<std::vector<T> >(dat);}
	size_t size(){return this->data->size();}
  T & operator [](int i) {return (*this->data)[i];}
  void emplace_back(T& tin){this->data->emplace_back(tin);}
  void push_back(T& tin){this->data->push_back(tin);}
  void shrink_to_fit(){this->data->shrink_to_fit();}
  void clear(){this->data->clear();}

};


// template< class T>
// std::vector<T> to_vec(pvector<T>& pv)
// {
// 	return *pv.data;
// }


// template<class Neighbourer, class Holder>
// class RasterData
// {
// public:

// 	Holder data;

// 	RasterData(){}
// 	RasterData(Holder& ho){this->data = ho;}

// 	T operator [](int i) const    {return this->data[i];}
//   T & operator [](int i) {return this->data[i];}

// };









template <typename T>
std::vector<size_t> sort_indexes(T &v) {

  // initialize original index locations
  int size = int(v.size());
  std::vector<size_t> idx(size);
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  std::sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}


template<class T, class U>
void add_noise_to_vector(T& vec, U min, U max)
{

	std::random_device rd; // obtain a random number from hardware
  std::mt19937 gen(rd()); // seed the generator
  std::uniform_real_distribution<> distr(min, max); // define the range

  for(size_t i=0; i<vec.size(); ++i)
  {
  	auto cd = distr(gen);
  	vec[i] += cd;
  }

}

template<class T>
pvector<T> format_input(std::vector<T>& in){return pvector<T>(in);}

template<class T>
pvector<T> format_input(pvector<T>& in){return in;}

template<class T>
std::vector<T> to_vec(pvector<T>& in)
{
	return std::vector<T>(*in.data);
}


template<class T>
std::vector<T> to_vec(std::vector<T>& in)
{
	return std::vector<T>(in);
}



// template<class T>
// std::vector<T> format_output(std::vector<T>& in){return in;}

// int sort(std::vector<int>& array, int l, int h, int k)
// {
//   int mid = l + (h - l) / 2; //Chose middle element as pivot
//   int i = max(l, mid - k), j = i, end = min(mid + k, h); // Set appropriate range
//   swap(array[mid], array[end]); //Swap middle and last element to avoid extra complications
//   while (j < end) {
//       if (array[j] < array[end]) {
//           swap(array[i++], array[j]);
//       }
//       j = j + 1;
//   }
//   swap(array[end], array[i]);
//   return i;
// }
 
// void ksorter(std::vector<int>& array, int l, int h, int k)
// {
//   if (l < h) {
//       int q = sort(array, l, h, k);
//       ksorter(array, l, q - 1, k);
//       ksorter(array, q + 1, h, k);
//   }
// }


class easyRand
{
public:
	std::random_device rd; // obtain a random number from hardware
  std::mt19937 gen; // seed the generator
  std::uniform_real_distribution<> distr; // define the range

  easyRand()
  {
		std::random_device rd;
	  this->gen = std::mt19937(rd()); // seed the generator
	  this->distr = std::uniform_real_distribution<>(0, 1); // define the range
  }

  easyRand(double min, double max)
  {
		std::random_device rd;
	  this->gen = std::mt19937(rd()); // seed the generator
	  this->distr = std::uniform_real_distribution<>(min, max); // define the range
  }

  double get(){return this->distr(this->gen);}
};




class ocarina
{
public:

	std::chrono::high_resolution_clock::time_point start,end;

	ocarina(){};

	void tik(){this->start = std::chrono::high_resolution_clock::now();}
	void tok(std::string message)
	{
		this->end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> ms_double = this->end - this->start;
		double timing = ms_double.count();
		std::cout << message << " took " << timing << " ms" << std::endl;
	}


	
};




























#endif