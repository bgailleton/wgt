/*
This header file extends the graph to provide routines to process depressions using Cordonnier et al, 2019
*/

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifndef cordonnier_versatile_2019_HPP
#define cordonnier_versatile_2019_HPP

// STL imports
#include <limits>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <ctime>
#include <fstream>
#include <queue>
#include <iostream>
#include <numeric>
#include <cmath>
#include <initializer_list>
#include <chrono>
#include <unordered_map>

#include "chonkutils.hpp"
// #include "graph.hpp"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;


class Graph;
template<class T>
class MGraph;

template<class n_t, class dist_t, class Neighbourer_t, class topo_t>
class LMRerouter;

template<class n_t, class dist_t>
class Cordonnier2019_v2;

template<class n_t, class dist_t>
class Cordonnier2019_v2MF;

template<class n_t, class dist_t>
class Link
{
public:

	n_t from,to, node_to, node_from;
	
	dist_t score;
	
	Link(){;}
	
	Link(n_t from, n_t to, n_t node_from, n_t node_to, dist_t score)
	{
		this->from = from;
		this->to = to;
		this->score = score;
		this->node_from = node_from;
		this->node_to = node_to;
	}

	void inverse()
	{
		auto temp = Link(this->to,this->from,this->node_to,this->node_from,this->score); 
		this->to = temp.to;
		this->from = temp.from;
		this->node_to = temp.node_to;
		this->node_from = temp.node_from;
		this->score = temp.score;
	}

	bool is_bas_in(int bas){return (bas == this->from || bas == this->to) ? true : false;}

};

template<class T, class U>
inline bool operator>( const Link<T,U>& lhs, const Link<T,U>& rhs )
{
	return lhs.score > rhs.score;
}
template<class T, class U>
inline bool operator<( const Link<T,U>& lhs, const Link<T,U>& rhs )
{
	return lhs.score < rhs.score;
}

class UnionFind
{
	public:
		UnionFind(int size)
		{
			this->_parent =std::vector<int>(size);
			for(int i=0; i<size; i++)
				this->_parent[i] = i;
			this->_rank = std::vector<int>(size,0);
		};

		void Union(int& x, int& y)
		{
			int xroot = this->Find(x);
			int yroot = this->Find(y);

			if (xroot != yroot)
			{
				if(this->_rank[xroot] < this->_rank[yroot])
						this->_parent[xroot] = yroot;
				else
				{
					this->_parent[yroot] = xroot;
					if(this->_rank[xroot] == this->_rank[yroot])
						this->_rank[xroot] ++;
				}
			}
		}

		int Find(int& x)
		{
			int xp = x,xc;
			while (true)
			{
				xc = xp;
				xp = this->_parent[xc];
				if (xp == xc)
					break;
			}
			this->_parent[x] = xc;
			return xc;
		}

		std::vector<int> _parent;
		std::vector<int> _rank;

};


template<class n_t, class dist_t>
class UnionFindv2
{
	public:
		UnionFindv2(int size,Cordonnier2019_v2<n_t,dist_t> &cod )
		{
			this->cod = &cod;
			this->_parent = std::vector<int>(size);
			this->_open = std::vector<bool>(size);
			for(int i=0; i<size; i++)
			{
				this->_parent[i] = i;
				this->_open[i] = this->cod->is_open_basin[i];
			}
			this->_rank = std::vector<int>(size,0);
		};

		void Union(int& x, int& y)
		{
			int xroot = this->Find(x);
			int yroot = this->Find(y);

			if (xroot != yroot)
			{
				if(this->_rank[xroot] < this->_rank[yroot])
						this->_parent[xroot] = yroot;
				else
				{
					this->_parent[yroot] = xroot;
					if(this->_rank[xroot] == this->_rank[yroot])
						this->_rank[xroot] ++;
				}

				if(this->_open[xroot] || this->_open[yroot])
				{
					this->_open[xroot] = true;
					this->_open[yroot] = true;
				}
			}
		}

		int Find(int& x)
		{
			int xp = x,xc;
			while (true)
			{
				xc = xp;
				xp = this->_parent[xc];
				if (xp == xc)
					break;
			}
			this->_parent[x] = xc;
			return xc;
		}

		std::vector<int> _parent;
		std::vector<int> _rank;
		std::vector<bool> _open;
		Cordonnier2019_v2<n_t,dist_t> *cod;


};


template<class n_t, class dist_t, class Neighbourer_t, class topo_t, class LM_t>
class UnionFindv3
{
	public:
		UnionFindv3(int size,LM_t &cod )
		{
			this->cod = &cod;
			this->_parent = std::vector<int>(size);
			this->_open = std::vector<bool>(size);
			for(int i=0; i<size; i++)
			{
				this->_parent[i] = i;
				this->_open[i] = this->cod->is_open_basin[i];
			}
			this->_rank = std::vector<int>(size,0);
		};

		void Union(int& x, int& y)
		{
			int xroot = this->Find(x);
			int yroot = this->Find(y);

			if (xroot != yroot)
			{
				if(this->_rank[xroot] < this->_rank[yroot])
						this->_parent[xroot] = yroot;
				else
				{
					this->_parent[yroot] = xroot;
					if(this->_rank[xroot] == this->_rank[yroot])
						this->_rank[xroot] ++;
				}

				if(this->_open[xroot] || this->_open[yroot])
				{
					this->_open[xroot] = true;
					this->_open[yroot] = true;
				}
			}
		}

		int Find(int& x)
		{
			int xp = x,xc;
			while (true)
			{
				xc = xp;
				xp = this->_parent[xc];
				if (xp == xc)
					break;
			}
			this->_parent[x] = xc;
			return xc;
		}

		std::vector<int> _parent;
		std::vector<int> _rank;
		std::vector<bool> _open;
		LM_t *cod;


};


template<class n_t, class dist_t>
class Cordonnier2019
{

	public:

		// Pointer to the mother graph; NNEDS TO BE SINGLE FLOW
		Graph* graph;

		// Labelling the watersheds
		std::vector<n_t> basin_labels;
		std::vector<n_t> basin_to_outlets;
		std::vector<n_t> pits_to_reroute;
		std::vector<n_t> mstree;
		std::vector< std::vector<int> > conn_basins;
		std::vector< std::vector<int> > conn_nodes;

		int nbasins;
		int npits = 0;


		Cordonnier2019(){;};
		Cordonnier2019(Graph& graph, std::vector<dist_t>& elevation)
		{
			this->graph = &graph;
			this->compute_basins_and_pits();
			this->preprocess_flowrouting(elevation);

		}

		void update_receivers(std::string& method, std::vector<dist_t>& elevation)
		{
			if(method == "simple" || method == "Simple")
				this->_update_pits_receivers_sompli(this->conn_basins, this->conn_nodes, mstree, elevation);
			else if (method == "carve")
			{
				this->_update_pits_receivers_carve(this->conn_basins, this->conn_nodes, mstree, elevation);
			}
			else if (method == "fill")
			{
				// std::cout << "gwamoulg" << std::endl;
				this->_update_pits_receivers_fill(this->conn_basins, this->conn_nodes, mstree, elevation);
			}
			this->graph->recompute_SF_donors_from_receivers();
		}

		// Uses the stack structure to build a quick basin array
		void compute_basins_and_pits()
		{
			this->basin_labels = std::vector<n_t>(this->graph->nnodes_t, -1);
			this->basin_to_outlets.reserve(200);
			this->pits_to_reroute.reserve(200);
			n_t lab = -1;
			for(auto tnode: this->graph->stack)
			{
				if(this->graph->can_flow_even_go_there(tnode) == false)
				{
					// std::cout << "bagul?" << std::endl;
					continue;
				}
				
				if(this->graph->receivers[tnode] == tnode)
				{
					lab++;
					this->basin_to_outlets.emplace_back(tnode);
					if(this->graph->can_flow_out_there(tnode) == false)
					{
						this->pits_to_reroute.emplace_back(tnode);
						this->npits++;
					}

				}

				if(lab == -1)
					throw std::runtime_error("flonflon?");

				this->basin_labels[tnode] = lab;
			
			}
			
			// for(auto tl : this->basin_labels)
			// 	if(tl == -1)
			// 		throw std::runtime_error("flonflonww?");

			this->nbasins = lab + 1;
		}


	void preprocess_flowrouting(std::vector<dist_t>& elevation)
	{
		// """Ensure that no flow is captured in sinks.

		// If needed, update `receivers`, `dist2receivers`, `ndonors`,
		// `donors` and `stack`.

		// """
		// int nnodes = this->graph->nnodes;

		// # theory of planar graph -> max nb. of connections known
		int nconn_max = this->nbasins * 8;
		this->conn_basins.clear(); this->conn_basins.reserve(nconn_max); // np.empty((nconn_max, 2), dtype=np.intp)
		this->conn_nodes.clear(); this->conn_nodes.reserve(nconn_max); // np.empty((nconn_max, 2), dtype=np.intp)
	 	std::vector<dist_t> conn_weights; conn_weights.reserve(nconn_max);
	 	for(int i=0; i<nconn_max; i++)
	 	{
	 		this->conn_basins.emplace_back(std::vector<n_t>({-2,-2}));
	 		this->conn_nodes.emplace_back(std::vector<n_t>({-2,-2}));
	 		conn_weights.emplace_back(-2);
	 	}

		int nconn, basin0;
		auto t1 = high_resolution_clock::now();

		this->_connect_basins(this->conn_basins, this->conn_nodes, conn_weights, elevation, nconn, basin0);
		auto t2 = high_resolution_clock::now();
		duration<double, std::milli> ms_double = t2 - t1;
		std::cout << "_compute_links --> " << ms_double.count() << " milliseconds" << std::endl;;
		int scb = nconn, scn = nconn, scw = nconn;

		this->conn_basins = std::vector<std::vector<n_t> >(this->conn_basins.begin(), this->conn_basins.begin() + scb - 1);
		this->conn_nodes = std::vector<std::vector<n_t> >(this->conn_nodes.begin(), this->conn_nodes.begin() + scn - 1);
		conn_weights = std::vector<dist_t>(conn_weights.begin(), conn_weights.begin() + scw - 1);

		mstree = this->_compute_mst_kruskal(this->conn_basins, conn_weights);

		this->_orient_basin_tree(this->conn_basins, this->conn_nodes, basin0, mstree);
		// this->_update_pits_receivers(this->conn_basins, conn_nodes, mstree, elevation);		
	}


	void _connect_basins(std::vector<std::vector<n_t> >& conn_basins,
				 std::vector<std::vector<n_t> >& conn_nodes, std::vector<dist_t>& conn_weights,					
				 std::vector<dist_t>& elevation, int& nconn, int& basin0)
	{
		int iconn = 0;

		basin0 = -1; //intp? intp. intp!??? intp!
		int ibasin = 0;

		std::vector<int> conn_pos = std::vector<int>(this->nbasins, -1); //np.full(nbasins, -1, dtype=np.intp)
		std::vector<int> conn_pos_used = std::vector<int>(this->nbasins, -1); //np.full(nbasins, -1, dtype=np.intp)
		
		int conn_pos_used_size = 0;
		
		bool iactive = false;

		for(auto istack : this->graph->stack)
		{
			if(this->graph->can_flow_even_go_there(istack) == false)
				continue;

			int irec = this->graph->receivers[istack];

			// # new basin
			if(irec == istack)
			{
				ibasin = this->basin_labels[istack];
				iactive = !this->graph->can_flow_out_there(istack);

				// std::cout << "Active::" << iactive << std::endl;
				// for iused in conn_pos_used[:conn_pos_used_size]
				for(int iused = 0; iused < conn_pos_used_size; iused++)
				{
					conn_pos[conn_pos_used[iused]] = -1;
				}

				conn_pos_used_size = 0;

				if (iactive == false)
				{
					if(basin0 == -1)
					{
						int row,col;
						this->graph->rowcol_from_node_id(istack,row,col);
						std::cout << "SHOULD HPPEN ONLY ONCE AND THE LOCATION IS " << row << "||" << col << std::endl;
						basin0 = ibasin;
					}
					else
					{
						this->conn_basins[iconn][0] = basin0;
						this->conn_basins[iconn][1] = ibasin;
						this->conn_nodes[iconn][0] = -1;
						this->conn_nodes[iconn][1] = -1;

						conn_weights[iconn] = - std::numeric_limits<double>::max();
						iconn ++;
					}
				}
			}

			if(iactive)
			{
				// std::cout << "iactive" << std::endl;

				auto neighbours = this->graph->get_neighbours(istack, false);

				for(auto& ineighbor:neighbours)
				{
					// std::cout << "1::" << ineighbor.node << "|" << this->basin_labels.size() << std::endl;
					if(this->graph->can_flow_even_go_there(ineighbor.node) == false)
					{
						// std::cout << "flowout" << std::endl;
						continue;
					}

					// std::cout << "sg1" << std::endl;
					if(ineighbor.node >= this->basin_labels.size() || ineighbor.node < 0)
					{
						std::cout << "Node::" << istack << "/" << this->graph->nnodes << std::endl;;
						std::cout << "neighbour::" << ineighbor.node << "/" << this->graph->nnodes << std::endl;;
						std::cout << "NEIGHBOURS::";
						for(auto nn:neighbours)
							std::cout << nn.node << "::";
						std::cout << std::endl;
						std::cout << "ineighb::" << this->graph->_get_neighbourer_id(istack) << std::endl;;
						std::cout << this->graph->boundary[istack] << std::endl;
						throw std::runtime_error("bite2");
					}

					int ineighbor_basin = this->basin_labels[ineighbor.node];
					// std::cout << "fault" << std::endl;
					
					if(ineighbor_basin == -1)
					{
						// std::cout << "gaboune" << std::endl;
						continue;
						throw std::runtime_error("wuuuut??? " + std::to_string(this->graph->can_flow_out_there(ineighbor.node)) );
					}

					// std::cout << "sg2" << std::endl;
					int ineighbor_outlet = this->basin_to_outlets[ineighbor_basin];
					// std::cout << "fault" << std::endl;

					// # skip same basin or already connected adjacent basin
					// # don't skip adjacent basin if it's an open basin
					if(ibasin >= ineighbor_basin && this->graph->can_flow_out_there(ineighbor_outlet) == false)
						continue;

					// std::cout << "sg3" << std::endl;
					double weight = std::max(elevation[istack], elevation[ineighbor.node]);
					// std::cout << "fault" << std::endl;
					// std::cout << "sg4" << std::endl;
					int conn_idx = conn_pos[ineighbor_basin];
					// std::cout << "fault" << std::endl;

					// # add new connection
					if(conn_idx == -1)
					{
						// conn_basins[iconn] = (ibasin, ineighbor_basin);

						this->conn_basins[iconn][0] = ibasin;
						this->conn_basins[iconn][1] = ineighbor_basin;
						// conn_nodes[iconn] = (istack, ineighbor);
						this->conn_nodes[iconn][0] = istack;
						this->conn_nodes[iconn][1] = ineighbor.node;

						conn_weights[iconn] = weight;

						conn_pos[ineighbor_basin] = iconn;
						iconn ++;

						conn_pos_used[conn_pos_used_size] = ineighbor_basin;
						conn_pos_used_size ++;

					}

					// # update existing connection
					else if (weight < conn_weights[conn_idx])
					{
						// std::cout << "sg1" << std::endl;
						// conn_nodes[conn_idx] = (istack, ineighbor);
						this->conn_nodes[conn_idx][0] = istack;
						this->conn_nodes[conn_idx][1] = ineighbor.node;
						conn_weights[conn_idx] = weight;
						// std::cout << "fault" << std::endl;
					}

					// std::cout << "2::" << ineighbor.node << std::endl;

				}
				// std::cout << "iactive done" << std::endl;

			}

		}

		nconn = iconn;

		return;
	}


	std::vector<int> _compute_mst_kruskal(std::vector<std::vector<n_t> >& conn_basins, std::vector<dist_t>& conn_weights)
	{
		std::vector<int> mstree = std::vector<int>(this->nbasins - 1);
		int mstree_size = 0;

		// # sort edges by elevation
		auto sort_id = sort_indexes<dist_t>(conn_weights);


		UnionFind uf(nbasins);

		for (auto eid : sort_id)
		{
			// if(eid > 1000000 || eid <0)
				// std::cout << "FLURB::" << eid << std::endl; 
			int b0 = this->conn_basins[eid][0];
			int b1 = this->conn_basins[eid][1];

			if (uf.Find(b0) != uf.Find(b1))
			{
				mstree[mstree_size] = eid;
				// mstree_translated.emplace_back(std::initializer_list<int>{SBasinOutlets[b0],SBasinOutlets[b1]});
				mstree_size ++;
				uf.Union(b0, b1);
			}
		}
		return mstree;
	}


	void _orient_basin_tree(std::vector<std::vector<n_t> >& conn_basins, std::vector<std::vector<n_t> >& conn_nodes, int& basin0, std::vector<int>& tree)
	{
		// # nodes connections
		std::vector<int> nodes_connects_size(this->nbasins,0);
		std::vector<int> nodes_connects_ptr(this->nbasins);

		// # parse the edges to compute the number of edges per node
		for (auto i : tree)
		{
			nodes_connects_size[ this->conn_basins[i][0] ]++;
			nodes_connects_size[ this->conn_basins[i][1] ]++;
		}

		// # compute the id of first edge in adjacency table
		nodes_connects_ptr[0] = 0;

		//lsdkfjlsjd
		// gruklb.
		for (int i = 1; i < nbasins; i++)
		{
			nodes_connects_ptr[i] = (nodes_connects_ptr[i - 1] + nodes_connects_size[i - 1]);
			nodes_connects_size[i - 1] = 0;
		}

		// # create the adjacency table
		int nodes_adjacency_size = nodes_connects_ptr[nbasins - 1] + nodes_connects_size[nbasins - 1];
		nodes_connects_size[this->nbasins -1] = 0;
		std::vector<int> nodes_adjacency(nodes_adjacency_size,0);

		// # parse the edges to update the adjacency
		for (auto i : tree)
		{

			int n1 = this->conn_basins[i][0];
			int n2 = this->conn_basins[i][1];
			nodes_adjacency[nodes_connects_ptr[n1] + nodes_connects_size[n1]] = i;
			nodes_adjacency[nodes_connects_ptr[n2] + nodes_connects_size[n2]] = i;
			nodes_connects_size[n1] ++;
			nodes_connects_size[n2] ++;
		}


		// # depth-first parse of the tree, starting from basin0
		// # stack of node, parent
		std::vector<std::vector<n_t> > stack; stack.reserve(nbasins);
		for(int i=0; i < nbasins; i++)
			stack.emplace_back(std::vector<n_t>({-2,-2}));


		int stack_size = 1;
		stack[0][0] = basin0;// (basin0, basin0)
		stack[0][1] = basin0;

		int n_turn = 0;
		while (stack_size > 0)
		{
			n_turn++;
			// # get parsed node
			stack_size = stack_size - 1;
			int node = stack[stack_size][0];
			int parent = stack[stack_size][1];

			// # for each edge of the graph
			// for i in range(nodes_connects_ptr[node], nodes_connects_ptr[node] + nodes_connects_size[node])
			for( int i = nodes_connects_ptr[node]; i < (nodes_connects_ptr[node] + nodes_connects_size[node]); i++) 
			{
				int edge_id = nodes_adjacency[i];

				// # the edge comming from the parent node has already been updated.
				// # in this case, the edge is (parent, node)
				if (this->conn_basins[edge_id][0] == parent && node != parent)
				{
						continue;
				}
				// # we want the edge to be (node, next)
				// # we check if the first node of the edge is not "node"
				if(node != this->conn_basins[edge_id][0])
				{
					// # swap n1 and n2
					int cb1 = this->conn_basins[edge_id][1];
					int cb0 = this->conn_basins[edge_id][0];
					this->conn_basins[edge_id][0] = cb1;
					this->conn_basins[edge_id][1] = cb0;

					cb1 = this->conn_nodes[edge_id][1];
					cb0 = this->conn_nodes[edge_id][0];
					// # swap p1 and p2
					this->conn_nodes[edge_id][0] = cb1;
					this->conn_nodes[edge_id][1] = cb0;
				}
				// # add the opposite node to the stack
				stack[stack_size][0] = this->conn_basins[edge_id][1];
				stack[stack_size][1] = node;
				stack_size ++;
			}
		}

	}

	void _update_pits_receivers_sompli(std::vector<std::vector<n_t> >& conn_basins,std::vector<std::vector<n_t> >& conn_nodes, 
		std::vector<n_t>& mstree, std::vector<dist_t>& elevation)
	{
		// for i in mstree:
		for(auto i : mstree)
		{

			int node_to = this->conn_nodes[i][0];
			int node_from = this->conn_nodes[i][1];

			// # skip open basins
			if (node_from == -1)
			{
					continue;
			}

			int outlet_from = this->basin_to_outlets[conn_basins[i][1] ];
			this->graph->receivers[outlet_from] = node_to;
			this->graph->distance2receivers[outlet_from] = this->graph->dx; // just to have a length but it should not actually be used
		}

	 
	}


	void _update_pits_receivers_carve(std::vector<std::vector<n_t> >& conn_basins,std::vector<std::vector<n_t> >& conn_nodes, 
		std::vector<n_t>& mstree, std::vector<dist_t>& elevation)
	{
		// for i in mstree:
		// std::cout <<"ID,rfrom,cfrom,rto,cto,rout,cout" << std::endl;

		std::vector<bool> pitdone(this->nbasins,false);

		for(auto i : mstree)
		{
			// int i = mstree[ti];
			int node_from = this->conn_nodes[i][1];
			int onode_to = this->conn_nodes[i][0];
			int onode_from = this->conn_nodes[i][1];
			int outlet_from = this->basin_to_outlets[conn_basins[i][1] ];

			if(pitdone[this->basin_labels[outlet_from]])
				continue;
			pitdone[this->basin_labels[outlet_from]] = true;

			// int rowf,colf,rowt,colt, rowo, colo;
			// this->graph->rowcol_from_node_id(node_from,rowf,colf);
			// this->graph->rowcol_from_node_id(node_to,rowt,colt);
			// this->graph->rowcol_from_node_id(outlet_from,rowo,colo);
			// std::cout << lab <<"," <<rowf <<"," <<colf <<"," << rowt<<"," <<colt <<"," << rowo<<"," << colo << std::endl;
			// lab++;
			// std::cout << "Basin:" << i << std::endl << 
			// " can outlet_from escape? " << this->graph->can_flow_out_there(outlet_from) << "|" << rowo << "|" << colo << std::endl <<
			// "node_from::" << node_from << "|" << rowf << "|" << colf << " node_to::" << node_to << "|" << rowt << "|" << colt << std::endl;;


			// # skip open basins
			if (node_from == -1)
			{
					continue;
			}

			int next_node = this->graph->receivers[node_from];
			// if(next_node == node_from)
			// 	std::cout << "WABUL WABUL YOLO" << std::endl;
			int temp = node_from;
			bool keep_on = true;
			do
			{
				if(next_node == outlet_from)
					keep_on = false;

				temp = this->graph->receivers[next_node];
				this->graph->receivers[next_node] = node_from;
				this->graph->distance2receivers[next_node] = this->graph->dx; // just to have a length but it should not actually be used
				node_from = next_node;
				next_node = temp;
				// std::cout << next_node << "|";
			} while(keep_on);



			this->graph->receivers[onode_from] = onode_to;
			this->graph->distance2receivers[onode_from] = this->graph->dx; // just to have a length but it should not actually be used
		}
		
	}

	void _update_pits_receivers_fill(std::vector<std::vector<n_t> >& conn_basins,std::vector<std::vector<n_t> >& conn_nodes, 
		std::vector<n_t>& mstree, std::vector<dist_t>& elevation)
	{
		// for i in mstree:
		// std::cout <<"yolo";
		std::vector<bool> isdone(this->graph->nnodes_t,false);
		for(auto i : mstree)
		{
			// std::cout << i << "|" << this->graph->can_flow_even_go_there(this->conn_nodes[i][0]) << "|" << this->graph->can_flow_even_go_there(this->conn_nodes[i][1]) << "|" << conn_basins[i][1] << "|" << this->nbasins << std::endl;;
			int node_to = this->conn_nodes[i][0];
			int node_from = this->conn_nodes[i][1];
			int tg_bas = this->basin_labels[node_from];


			// # skip open basins
			if (node_from == -1)
			{
					continue;
			}
			// std::cout << "rec nodde_from beef " << this->graph->receivers[node_from];;
			float otopo = std::fmax(elevation[node_to],elevation[node_from]);
			std::queue<int> yonode;
			yonode.push(node_to);
			while(yonode.size() > 0)
			{
				int tnode = yonode.front(); yonode.pop();
	 			auto neighbours = this->graph->get_neighbours(tnode, false);
				for(auto& ineighbor:neighbours)
				{
					if(this->graph->can_flow_even_go_there(ineighbor.node) == false || tg_bas != this->basin_labels[ineighbor.node] || isdone[ineighbor.node] == true)
						continue;

					if(elevation[ineighbor.node] <= otopo)
					{
						// std::cout << ineighbor.node << "<--" << tnode << "|||";
						this->graph->receivers[ineighbor.node] = tnode;
						this->graph->distance2receivers[ineighbor.node] = ineighbor.distance;
						yonode.push(ineighbor.node);
						isdone[ineighbor.node] = true;
					}
				}

			}
			// std::cout << " after " << this->graph->receivers[node_from] << std::endl;;;

		}
		
	}

};





template<class n_t, class dist_t>
class Cordonnier2019_v2
{
public:

	// Pointer to the mother graph; NNEDS TO BE SINGLE FLOW
	Graph* graph;

	// Labelling the watersheds
	std::vector<n_t> basin_labels;
	std::vector<bool> is_open_basin;
	std::vector<n_t> basin_to_outlets;
	std::vector<n_t> pits_to_reroute;
	std::vector<n_t> mstree;
	std::vector<n_t> receivers, nodes_to, nodes_from, outlets_from;
	std::vector<std::vector<n_t> > donors;

	std::vector< Link<n_t, dist_t>* > stack;


	std::vector<Link<n_t, dist_t> > active_links;
	std::vector< std::vector< Link<n_t, dist_t>* > > bas2links;
	std::priority_queue<Link<n_t, dist_t>, std::vector<Link<n_t, dist_t>>, std::greater<Link<n_t, dist_t> > > pqlinks;




	int nbasins;
	int npits = 0;


	Cordonnier2019_v2(){;};
	Cordonnier2019_v2(Graph& graph)
	{
		this->graph = &graph;
		auto t1 = high_resolution_clock::now();
		this->compute_basins_and_pits();
		auto t2 = high_resolution_clock::now();
		duration<double, std::milli> ms_double = t2 - t1;
		// std::cout << "Computing basins and pits -> " << ms_double.count() << " milliseconds" << std::endl;;
		if(this->npits > 0)
		{
			t1 = high_resolution_clock::now();
			this->preprocess_flowrouting();
			t2 = high_resolution_clock::now();
			ms_double = t2 - t1;
			// std::cout << "Preprocess_flowrouting -> " << ms_double.count() << " milliseconds" << std::endl;;
			
		}
	}

	void preprocess_flowrouting()
	{
		auto t1 = high_resolution_clock::now();

		this->_compute_links();
		auto t2 = high_resolution_clock::now();
		duration<double, std::milli> ms_double = t2 - t1;
		// std::cout << "_compute_links --> " << ms_double.count() << " milliseconds" << std::endl;;
		t1 = high_resolution_clock::now();
		this->_compute_mst_kruskal();
		t2 = high_resolution_clock::now();
		ms_double = t2 - t1;
		// std::cout << "_compute_mst_kruskal --> " << ms_double.count() << " milliseconds" << std::endl;;
		t1 = high_resolution_clock::now();
		this->_orient_basin_tree();
		t2 = high_resolution_clock::now();
		ms_double = t2 - t1;
		// std::cout << "_orient_basin_tree --> " << ms_double.count() << " milliseconds" << std::endl;;
	}

	// Uses the stack structure to build a quick basin array
	void compute_basins_and_pits()
	{
		this->basin_labels = std::vector<n_t>(this->graph->nnodes_t, -1);
		this->basin_to_outlets.reserve(200);
		this->pits_to_reroute.reserve(200);
		n_t lab = -1;
		for(auto tnode: this->graph->stack)
		{
			if(this->graph->can_flow_even_go_there(tnode) == false)
			{
				// std::cout << "bagul?" << std::endl;
				continue;
			}
			
			if(this->graph->receivers[tnode] == tnode)
			{
				lab++;
				this->basin_to_outlets.emplace_back(tnode);
				if(this->graph->can_flow_out_there(tnode) == false)
				{
					this->pits_to_reroute.emplace_back(tnode);
					this->is_open_basin.emplace_back(false);
					this->npits++;
				}
				else
					this->is_open_basin.emplace_back(true);


			}

			// if(lab == -1)
			// 	throw std::runtime_error("flonflon?");

			this->basin_labels[tnode] = lab;
		
		}
		
		this->nbasins = lab + 1;
	}

	void _compute_links()
	{

		// Initialising a matrix of links for each basins in order to stor the minimum elevation links between each pair of basins
		// std::vector<std::vector<dist_t> > mat_of_scores(this->nbasins,std::vector<dist_t>(this->nbasins,std::numeric_limits<dist_t>::max()));
		// std::vector<std::vector<std::vector<n_t> > > mat_of_nodes(this->nbasins,std::vector<std::vector<n_t> >(this->nbasins));

		std::unordered_map<int, std::unordered_map<int,dist_t > > mat_of_scores;
		std::unordered_map<int, std::unordered_map<int, std::vector<n_t> > > mat_of_nodes;
		auto t1 = high_resolution_clock::now();
		

		for(int i=0; i<this->graph->nnodes; ++i)
		{
			if(this->graph->can_flow_even_go_there(i) == false)
				continue;

			auto neighbours = this->graph->get_neighbours_only_id(i);
			int tbas = this->basin_labels[i];
			dist_t telev = this->graph->topography[i];
			for(auto tn:neighbours)
			{
				if(this->graph->can_flow_even_go_there(tn) == false)
					continue;

				int obas = this->basin_labels[tn];
				if(tbas != obas)
				{

					if(this->is_open_basin[tbas] && this->is_open_basin[obas])
						continue;

					dist_t score = std::min(telev, this->graph->topography[tn]);
					int basA = (tbas > obas) ? obas : tbas;
					int basB = (tbas > obas) ? tbas : obas;

					bool isinmap = mat_of_scores.count(basA) > 0;
					if(isinmap)
						isinmap = mat_of_scores[basA].count(basB) > 0;

					bool need = false;
					if(isinmap)
					{
						if(mat_of_scores[basA][basB] > score)
						{
							need = true;
						}
					}
					else
						need = true;


					if(need)
					{
						if(basA == tbas)
						{
							std::vector<n_t> temp = std::initializer_list<n_t>{i,tn};
							mat_of_nodes[basA][basB] = temp;
						}
						else
						{
							std::vector<n_t> temp = std::initializer_list<n_t>{tn,i};
							mat_of_nodes[basA][basB] = temp;
						}

						mat_of_scores[basA][basB] = score;
					}
				}
			}
		}
		auto t2 = high_resolution_clock::now();
		duration<double, std::milli> ms_double = t2 - t1;
		// std::cout << "_compute_links::building ---> " << ms_double.count() << " milliseconds" << std::endl;;


		t1 = high_resolution_clock::now();
		for (auto it=mat_of_scores.begin(); it!=mat_of_scores.end(); ++it)
		{
			int basA = it->first;
			for (auto it2=mat_of_scores[basA].begin(); it2!=mat_of_scores[basA].end(); ++it2)
			{
				int basB = it2->first;
				dist_t score = it2->second;
				this->pqlinks.emplace(Link<n_t,dist_t>(basA,basB,mat_of_nodes[basA][basB][0],mat_of_nodes[basA][basB][1],score));
			}


			

		}

		t2 = high_resolution_clock::now();
		ms_double = t2 - t1;
		// std::cout << "_compute_links::sorting ---> " << ms_double.count() << " milliseconds" << std::endl;;

		// dist_t maxval = std::numeric_limits<dist_t>::max();
		// for(int i = 0; i < this->nbasins;++i)
		// {
		// 	for(int j=0; j<this->nbasins;++j)
		// 	{

		// 		if(j <= i)
		// 			continue;

		// 		if(mat_of_scores[i][j] ==	maxval)
		// 			continue;

		// 		this->pqlinks.emplace(Link<n_t,dist_t>(i,j,mat_of_nodes[i][j][0],mat_of_nodes[i][j][1],mat_of_scores[i][j]));
		// 		// std::cout << "LIIINK::"<<i << "|" << j << std::endl;
		// 	}
		// }

	}

	void _compute_mst_kruskal()
	{
		this->bas2links = std::vector< std::vector< Link<n_t, dist_t>* > >( this->nbasins, std::vector< Link<n_t, dist_t>* >() );

		UnionFindv2<n_t,dist_t> uf(nbasins, (*this) );
		int j = 0;

		while(this->pqlinks.empty() == false)
		{
			auto next = this->pqlinks.top();
			this->pqlinks.pop();

			int b0 = next.from;
			int b1 = next.to;
			int fb0 = uf.Find(b0);
			int fb1 = uf.Find(b1) ;
			if (fb0 != fb1)
			{

				if(uf._open[fb0] && uf._open[fb1])
					continue;
				// std::cout << "GOUGN::" << b0 << "|" << b1 << "|" << this->nbasins << std::endl;
				
				uf.Union(b0, b1);
				this->active_links.emplace_back(next);
				++j;
			}
		}

		for(size_t i = 0; i < this->active_links.size(); ++i)
		{
			this->bas2links[this->active_links[i].from].emplace_back(&this->active_links[i]);
			this->bas2links[this->active_links[i].to].emplace_back(&this->active_links[i]);
		}

		// for (int tb = 0 ; tb < this->nbasins; ++tb)
		// for(size_t j=0; j<this->bas2links[tb].size(); ++j)
		// {
		// 	if(this->bas2links[tb][j]->from >= this->nbasins || this->bas2links[tb][j]->to >= this->nbasins)
		// 			throw std::runtime_error("waffon2");
		// }

	}

	void _orient_basin_tree()
	{

		this->receivers = std::vector<n_t>(this->nbasins, -1);
		this->nodes_to = std::vector<n_t>(this->nbasins, -1);
		this->nodes_from = std::vector<n_t>(this->nbasins, -1);
		this->outlets_from = std::vector<n_t>(this->nbasins, -1);
		this->stack.clear();this->stack.reserve(this->nbasins);

		std::vector<bool> isdone(this->nbasins,false);



		// std::cout << "here::" << this->active_links.size() << std::endl;		
		for(int i = this->active_links.size() - 1; i >= 0; --i)
		{
			// getting the link
			// auto& tlink = this->active_links[i]; 

			if(this->is_open_basin[this->active_links[i].from] == false && this->is_open_basin[this->active_links[i].to] == false)
				continue;

			if(this->is_open_basin[this->active_links[i].from])
			{
				// std::cout << "BOULOUF" << std::endl;
				this->active_links[i].inverse();
			}

			std::stack<int, std::vector<int> > stackhelper;
			stackhelper.emplace(this->active_links[i].to);

			while(stackhelper.empty() == false)
			{
				int next = stackhelper.top();
				isdone[next] = true;

				stackhelper.pop();
				for(auto ptr: this->bas2links[next])
				{
					// std::cout << "Wafulb1::" << ptr->from << "|" << ptr->to << std::endl;
					if(ptr->from == next)
					{
						if(isdone[ptr->to] == false)
						{
							stackhelper.emplace(ptr->to);
							ptr->inverse();
							// if(isdone[ptr->from] == false)
							// {
								// std::cout << "STACKADD1" << std::endl;
								this->stack.emplace_back(ptr);
								// std::cout << this->is_open_basin[ptr->from] << "||" << this->is_open_basin[ptr->to] << std::endl;
								isdone[ptr->from] = true;
							// }
						}
					}
					else if(isdone[ptr->from] == false)
					{
						stackhelper.emplace(ptr->from);
						if(isdone[ptr->from] == false)
						{
							// std::cout << "STACKADD2" << std::endl;
							this->stack.emplace_back(ptr);
							isdone[ptr->from] = true;
						}
					}
				}
			}
		}

	}

	void compute_TO_SF_stack_version()
	{
		// Initialising the stack
		this->stack.clear();
		// reserving the amount of stuff
		this->stack.reserve(this->nbasins);

		// The stack container helper
		std::stack<int, std::vector<int> > stackhelper;
		// std::vector<bool> isdone(this->nbasins,false);
		// going through all the nodes
		for(int i=0; i<this->nbasins; ++i)
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



	}

	void _update_pits_receivers_carve()
	{

		for(int i=this->stack.size() - 1; i >=0 ; --i)
		{

			auto tlink = this->stack[i];
			int node_from = tlink->node_from;
			int onode_to = tlink->node_to;
			int onode_from = tlink->node_from;
			int outlet_from = this->basin_to_outlets[tlink->from];

			int next_node = this->graph->receivers[node_from];
			int temp = node_from;
			bool keep_on = true;
			do
			{


				if(next_node == outlet_from)
				{
					keep_on = false;
				}


				temp = this->graph->receivers[next_node];
				this->graph->receivers[next_node] = node_from;
				this->graph->distance2receivers[next_node] = this->graph->distance2receivers[node_from]; // just to have a length but it should not actually be used
				node_from = next_node;
				next_node = temp;
			} while(keep_on);



			this->graph->receivers[onode_from] = onode_to;
			this->graph->distance2receivers[onode_from] = this->graph->dx; // just to have a length but it should not actually be used

		}
		
	}



	void _update_pits_receivers_sompli()
	{
		// for i in mstree:
		for(int i=this->stack.size() - 1; i >=0 ; --i)
		{

			auto tlink = this->stack[i];


			int node_to = tlink->node_to;
			int node_from = tlink->node_from;

			// # skip open basins
			if (node_from == -1)
			{
					continue;
			}

			// int outlet_from = this->basin_to_outlets[conn_basins[i][1] ];
			int outlet_from = this->basin_to_outlets[tlink->from];

			this->graph->receivers[outlet_from] = node_to;
			this->graph->distance2receivers[outlet_from] = this->graph->dx; // just to have a length but it should not actually be used
		}

	 
	}

	void _update_pits_receivers_fill()
	{
		// for i in mstree:
		// std::cout <<"yolo";
		std::vector<bool> isdone(this->graph->nnodes_t,false);
		for(int i=this->stack.size() - 1; i >=0 ; --i)
		// for(int i=0;  i<this->stack.size(); ++i)
		{

			auto tlink = this->stack[i];
			// int i = mstree[ti];
			int node_to = tlink->node_to;
			int node_from = tlink->node_from;

			int tg_bas = this->basin_labels[node_from];


			// # skip open basins
			if (node_from == -1)
			{
					continue;
			}
			// std::cout << "rec nodde_from beef " << this->graph->receivers[node_from];;
			float otopo = std::fmax(this->graph->topography[node_to],this->graph->topography[node_from]);
			std::queue<int> yonode;
			yonode.push(node_to);
			while(yonode.size() > 0)
			{
				int tnode = yonode.front(); yonode.pop();
	 			auto neighbours = this->graph->get_neighbours(tnode, false);
				for(auto& ineighbor:neighbours)
				{
					if(this->graph->can_flow_even_go_there(ineighbor.node) == false || tg_bas != this->basin_labels[ineighbor.node] || isdone[ineighbor.node] == true)
						continue;

					if(this->graph->topography[ineighbor.node] <= otopo)
					{
						// std::cout << ineighbor.node << "<--" << tnode << "|||";
						this->graph->receivers[ineighbor.node] = tnode;
						this->graph->distance2receivers[ineighbor.node] = ineighbor.distance;
						yonode.push(ineighbor.node);
						isdone[ineighbor.node] = true;
					}
				}

			}
		}
		
	}



	// bool is_open(int i)
	// {
	// 	return this->graph->can_flow_out_there(this->outlets_from[i]);
	// }


	void update_receivers(std::string& method)
	{
		// auto t1 = high_resolution_clock::now();

		if(method == "simple" || method == "Simple")
			this->_update_pits_receivers_sompli();
		else if (method == "carve")
		{
			this->_update_pits_receivers_carve();
		}
		else if (method == "fill")
		{
			// std::cout << "gwamoulg" << std::endl;
			this->_update_pits_receivers_fill();
		}
		// std::cout << "fabul" << std::endl;
		this->graph->recompute_SF_donors_from_receivers();

		// auto t2 = high_resolution_clock::now();
		// duration<double, std::milli> ms_double = t2 - t1;
		// std::cout << "Update recs -> " << ms_double.count() << " milliseconds" << std::endl;;

	}


};



template<class n_t, class dist_t>
class Cordonnier2019_v2MF
{
public:

	// Pointer to the mother graph; NNEDS TO BE SINGLE FLOW
	MGraph<dist_t>* graph;

	// Labelling the watersheds
	std::vector<n_t> basin_labels;
	std::vector<bool> is_open_basin;
	std::vector<n_t> basin_to_outlets;
	std::vector<n_t> pits_to_reroute;
	std::vector<n_t> mstree;
	std::vector<n_t> receivers, nodes_to, nodes_from, outlets_from;
	std::vector<std::vector<n_t> > donors;

	std::vector< Link<n_t, dist_t>* > stack;


	std::vector<Link<n_t, dist_t> > active_links;
	std::vector< std::vector< Link<n_t, dist_t>* > > bas2links;
	std::priority_queue<Link<n_t, dist_t>, std::vector<Link<n_t, dist_t>>, std::greater<Link<n_t, dist_t> > > pqlinks;




	int nbasins;
	int npits = 0;


	Cordonnier2019_v2MF(){;};
	Cordonnier2019_v2MF(MGraph<dist_t>& graph)
	{
		this->graph = &graph;
		auto t1 = high_resolution_clock::now();
		this->compute_basins_and_pits();
		auto t2 = high_resolution_clock::now();
		duration<double, std::milli> ms_double = t2 - t1;
		// std::cout << "Computing basins and pits -> " << ms_double.count() << " milliseconds" << std::endl;;
		if(this->npits > 0)
		{
			t1 = high_resolution_clock::now();
			this->preprocess_flowrouting();
			t2 = high_resolution_clock::now();
			ms_double = t2 - t1;
			// std::cout << "Preprocess_flowrouting -> " << ms_double.count() << " milliseconds" << std::endl;;
			
		}
	}

	void preprocess_flowrouting()
	{
		auto t1 = high_resolution_clock::now();

		this->_compute_links();
		auto t2 = high_resolution_clock::now();
		duration<double, std::milli> ms_double = t2 - t1;
		// std::cout << "_compute_links --> " << ms_double.count() << " milliseconds" << std::endl;;
		t1 = high_resolution_clock::now();
		this->_compute_mst_kruskal();
		t2 = high_resolution_clock::now();
		ms_double = t2 - t1;
		// std::cout << "_compute_mst_kruskal --> " << ms_double.count() << " milliseconds" << std::endl;;
		t1 = high_resolution_clock::now();
		this->_orient_basin_tree();
		t2 = high_resolution_clock::now();
		ms_double = t2 - t1;
		// std::cout << "_orient_basin_tree --> " << ms_double.count() << " milliseconds" << std::endl;;
	}

	// Uses the stack structure to build a quick basin array
	void compute_basins_and_pits()
	{
		this->basin_labels = std::vector<n_t>(this->graph->nnodes_t, -1);
		this->basin_to_outlets.reserve(200);
		this->pits_to_reroute.reserve(200);
		n_t lab = -1;
		for(auto tnode: this->graph->stack)
		{
			if(this->graph->can_flow_even_go_there(tnode) == false)
			{
				// std::cout << "bagul?" << std::endl;
				continue;
			}
			
			if(this->graph->Sreceivers[tnode] == tnode)
			{
				lab++;
				this->basin_to_outlets.emplace_back(tnode);
				if(this->graph->can_flow_out_there(tnode) == false)
				{
					this->pits_to_reroute.emplace_back(tnode);
					this->is_open_basin.emplace_back(false);
					this->npits++;
				}
				else
					this->is_open_basin.emplace_back(true);


			}

			// if(lab == -1)
			// 	throw std::runtime_error("flonflon?");

			this->basin_labels[tnode] = lab;
		
		}
		
		this->nbasins = lab + 1;
	}

	void _compute_links()
	{

		// Initialising a matrix of links for each basins in order to stor the minimum elevation links between each pair of basins
		// std::vector<std::vector<dist_t> > mat_of_scores(this->nbasins,std::vector<dist_t>(this->nbasins,std::numeric_limits<dist_t>::max()));
		// std::vector<std::vector<std::vector<n_t> > > mat_of_nodes(this->nbasins,std::vector<std::vector<n_t> >(this->nbasins));

		std::unordered_map<int, std::unordered_map<int,dist_t > > mat_of_scores;
		std::unordered_map<int, std::unordered_map<int, std::vector<n_t> > > mat_of_nodes;
		auto t1 = high_resolution_clock::now();
		

		for(int i=0; i<this->graph->nnodes; ++i)
		{
			if(this->graph->can_flow_even_go_there(i) == false)
				continue;

			auto neighbours = this->graph->get_neighbours_only_id(i);
			int tbas = this->basin_labels[i];
			dist_t telev = (*this->graph->topography)[i];
			for(auto tn:neighbours)
			{
				if(this->graph->can_flow_even_go_there(tn) == false)
					continue;

				int obas = this->basin_labels[tn];
				if(tbas != obas)
				{

					if(this->is_open_basin[tbas] && this->is_open_basin[obas])
						continue;

					dist_t score = std::min(telev, (*this->graph->topography)[tn]);
					int basA = (tbas > obas) ? obas : tbas;
					int basB = (tbas > obas) ? tbas : obas;

					bool isinmap = mat_of_scores.count(basA) > 0;
					if(isinmap)
						isinmap = mat_of_scores[basA].count(basB) > 0;

					bool need = false;
					if(isinmap)
					{
						if(mat_of_scores[basA][basB] > score)
						{
							need = true;
						}
					}
					else
						need = true;


					if(need)
					{
						if(basA == tbas)
						{
							std::vector<n_t> temp = std::initializer_list<n_t>{i,tn};
							mat_of_nodes[basA][basB] = temp;
						}
						else
						{
							std::vector<n_t> temp = std::initializer_list<n_t>{tn,i};
							mat_of_nodes[basA][basB] = temp;
						}

						mat_of_scores[basA][basB] = score;
					}
				}
			}
		}
		auto t2 = high_resolution_clock::now();
		duration<double, std::milli> ms_double = t2 - t1;
		// std::cout << "_compute_links::building ---> " << ms_double.count() << " milliseconds" << std::endl;;


		t1 = high_resolution_clock::now();
		for (auto it=mat_of_scores.begin(); it!=mat_of_scores.end(); ++it)
		{
			int basA = it->first;
			for (auto it2=mat_of_scores[basA].begin(); it2!=mat_of_scores[basA].end(); ++it2)
			{
				int basB = it2->first;
				dist_t score = it2->second;
				this->pqlinks.emplace(Link<n_t,dist_t>(basA,basB,mat_of_nodes[basA][basB][0],mat_of_nodes[basA][basB][1],score));
			}


			

		}

		t2 = high_resolution_clock::now();
		ms_double = t2 - t1;
		// std::cout << "_compute_links::sorting ---> " << ms_double.count() << " milliseconds" << std::endl;;

		// dist_t maxval = std::numeric_limits<dist_t>::max();
		// for(int i = 0; i < this->nbasins;++i)
		// {
		// 	for(int j=0; j<this->nbasins;++j)
		// 	{

		// 		if(j <= i)
		// 			continue;

		// 		if(mat_of_scores[i][j] ==	maxval)
		// 			continue;

		// 		this->pqlinks.emplace(Link<n_t,dist_t>(i,j,mat_of_nodes[i][j][0],mat_of_nodes[i][j][1],mat_of_scores[i][j]));
		// 		// std::cout << "LIIINK::"<<i << "|" << j << std::endl;
		// 	}
		// }

	}

	void _compute_mst_kruskal()
	{
		// this->mstree = std::vector<n_t>(this->nbasins - 1);
		this->bas2links = std::vector< std::vector< Link<n_t, dist_t>* > >( this->nbasins, std::vector< Link<n_t, dist_t>* >() );
		UnionFindv2<n_t,dist_t> uf(nbasins, (*this) );
		int j = 0;

		while(this->pqlinks.empty() == false)
		{
			auto next = this->pqlinks.top();
			this->pqlinks.pop();

			int b0 = next.from;
			int b1 = next.to;
			int fb0 = uf.Find(b0);
			int fb1 = uf.Find(b1) ;
			if (fb0 != fb1)
			{

				if(uf._open[fb0] && uf._open[fb1])
					continue;
				// std::cout << "GOUGN::" << b0 << "|" << b1 << "|" << this->nbasins << std::endl;
				
				uf.Union(b0, b1);
				this->active_links.emplace_back(next);
				++j;
			}
		}

		for(size_t i = 0; i < this->active_links.size(); ++i)
		{
			this->bas2links[this->active_links[i].from].emplace_back(&this->active_links[i]);
			this->bas2links[this->active_links[i].to].emplace_back(&this->active_links[i]);
		}

		// for (int tb = 0 ; tb < this->nbasins; ++tb)
		// for(size_t j=0; j<this->bas2links[tb].size(); ++j)
		// {
		// 	if(this->bas2links[tb][j]->from >= this->nbasins || this->bas2links[tb][j]->to >= this->nbasins)
		// 			throw std::runtime_error("waffon2");
		// }

	}

	void _orient_basin_tree()
	{

		this->receivers = std::vector<n_t>(this->nbasins, -1);
		this->nodes_to = std::vector<n_t>(this->nbasins, -1);
		this->nodes_from = std::vector<n_t>(this->nbasins, -1);
		this->outlets_from = std::vector<n_t>(this->nbasins, -1);
		this->stack.clear();this->stack.reserve(this->nbasins);

		std::vector<bool> isdone(this->nbasins,false);



		// std::cout << "here::" << this->active_links.size() << std::endl;		
		for(int i = this->active_links.size() - 1; i >= 0; --i)
		{
			// getting the link
			// auto& tlink = this->active_links[i]; 

			if(this->is_open_basin[this->active_links[i].from] == false && this->is_open_basin[this->active_links[i].to] == false)
				continue;

			if(this->is_open_basin[this->active_links[i].from])
			{
				// std::cout << "BOULOUF" << std::endl;
				this->active_links[i].inverse();
			}

			std::stack<int, std::vector<int> > stackhelper;
			stackhelper.emplace(this->active_links[i].to);

			while(stackhelper.empty() == false)
			{
				int next = stackhelper.top();
				isdone[next] = true;

				stackhelper.pop();
				for(auto ptr: this->bas2links[next])
				{
					// std::cout << "Wafulb1::" << ptr->from << "|" << ptr->to << std::endl;
					if(ptr->from == next)
					{
						if(isdone[ptr->to] == false)
						{
							stackhelper.emplace(ptr->to);
							ptr->inverse();
							// if(isdone[ptr->from] == false)
							// {
								// std::cout << "STACKADD1" << std::endl;
								this->stack.emplace_back(ptr);
								// std::cout << this->is_open_basin[ptr->from] << "||" << this->is_open_basin[ptr->to] << std::endl;
								isdone[ptr->from] = true;
							// }
						}
					}
					else if(isdone[ptr->from] == false)
					{
						stackhelper.emplace(ptr->from);
						if(isdone[ptr->from] == false)
						{
							// std::cout << "STACKADD2" << std::endl;
							this->stack.emplace_back(ptr);
							isdone[ptr->from] = true;
						}
					}
				}
			}
		}

	}

	void compute_TO_SF_stack_version()
	{
		// Initialising the stack
		this->stack.clear();
		// reserving the amount of stuff
		this->stack.reserve(this->nbasins);

		// The stack container helper
		std::stack<int, std::vector<int> > stackhelper;
		// std::vector<bool> isdone(this->nbasins,false);
		// going through all the nodes
		for(int i=0; i<this->nbasins; ++i)
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



	}

	void _update_pits_receivers_carve()
	{

		// std::cout <<"ID,rfrom,cfrom,rto,cto,rout,cout" << std::endl;
			
		// int lab=0;
		for(int i=this->stack.size() - 1; i >=0 ; --i)
		// for(int i=0;  i<this->stack.size(); ++i)
		{

			auto tlink = this->stack[i];
			// int i = mstree[ti];
			// int node_to = tlink->node_to;
			int node_from = tlink->node_from;
			int onode_to = tlink->node_to;
			int onode_from = tlink->node_from;
			int outlet_from = this->basin_to_outlets[tlink->from];


			int next_node = this->graph->Sreceivers[node_from];
			int temp = node_from;
			bool keep_on = true;
			do
			{


				if(next_node == outlet_from)
				{
					keep_on = false;
				}


				temp = this->graph->Sreceivers[next_node];

				this->graph->Sreceivers[next_node] = node_from;
				this->graph->Sdistance2receivers[next_node] = this->graph->Sdistance2receivers[node_from]; // just to have a length but it should not actually be used
				node_from = next_node;
				next_node = temp;
			} while(keep_on);



			this->graph->Sreceivers[onode_from] = onode_to;
			this->graph->Sdistance2receivers[onode_from] = this->graph->dx; // just to have a length but it should not actually be used

		}
		
	}



	void _update_pits_receivers_sompli()
	{
		// for i in mstree:
		for(int i=this->stack.size() - 1; i >=0 ; --i)
		{

			auto tlink = this->stack[i];


			int node_to = tlink->node_to;
			int node_from = tlink->node_from;

			// # skip open basins
			if (node_from == -1)
			{
					continue;
			}

			// int outlet_from = this->basin_to_outlets[conn_basins[i][1] ];
			int outlet_from = this->basin_to_outlets[tlink->from];

			this->graph->Sreceivers[outlet_from] = node_to;
			this->graph->Sdistance2receivers[outlet_from] = this->graph->dx; // just to have a length but it should not actually be used
		}

	 
	}

	void _update_pits_receivers_fill()
	{
		// for i in mstree:
		// std::cout <<"yolo";
		std::vector<bool> isdone(this->graph->nnodes_t,false);
		for(int i=this->stack.size() - 1; i >=0 ; --i)
		// for(int i=0;  i<this->stack.size(); ++i)
		{

			auto tlink = this->stack[i];
			// int i = mstree[ti];
			int node_to = tlink->node_to;
			int node_from = tlink->node_from;
			// int onode_to = tlink->node_to;
			// int onode_from = tlink->node_from;
			// int outlet_from = this->basin_to_outlets[tlink->from];

			int tg_bas = this->basin_labels[node_from];


			// # skip open basins
			if (node_from == -1)
			{
					continue;
			}
			// std::cout << "rec nodde_from beef " << this->graph->receivers[node_from];;
			float otopo = std::fmax((*this->graph->topography)[node_to],(*this->graph->topography)[node_from]);
			std::queue<int> yonode;
			yonode.push(node_to);
			while(yonode.size() > 0)
			{
				int tnode = yonode.front(); yonode.pop();
	 			auto neighbours = this->graph->get_neighbours(tnode, false);
				for(auto& ineighbor:neighbours)
				{
					if(this->graph->can_flow_even_go_there(ineighbor.node) == false || tg_bas != this->basin_labels[ineighbor.node] || isdone[ineighbor.node] == true)
						continue;

					if((*this->graph->topography)[ineighbor.node] <= otopo)
					{
						// std::cout << ineighbor.node << "<--" << tnode << "|||";
						this->graph->Sreceivers[ineighbor.node] = tnode;
						this->graph->Sdistance2receivers[ineighbor.node] = ineighbor.distance;
						yonode.push(ineighbor.node);
						isdone[ineighbor.node] = true;
					}
				}

			}
		}
		
	}



	// bool is_open(int i)
	// {
	// 	return this->graph->can_flow_out_there(this->outlets_from[i]);
	// }


	void update_receivers(std::string& method)
	{

		if(method == "simple" || method == "Simple")
			this->_update_pits_receivers_sompli();
		else if (method == "carve")
		{
			this->_update_pits_receivers_carve();
		}
		else if (method == "fill")
		{
			this->_update_pits_receivers_fill();
		}

		this->graph->recompute_SF_donors_from_receivers();

		this->graph->recompute_MF_impose_slope_SS();


	}


};



















template<class n_t, class dist_t, class Neighbourer_t, class topo_t>
class LMRerouter
{
public:

	// Labelling the watersheds
	std::vector<n_t> basin_labels;
	std::vector<bool> is_open_basin;
	std::vector<n_t> basin_to_outlets;
	std::vector<n_t> pits_to_reroute;
	std::vector<n_t> mstree;
	std::vector<n_t> receivers, nodes_to, nodes_from, outlets_from;
	std::vector<std::vector<n_t> > donors;
	std::vector< Link<n_t, dist_t>* > stack;
	std::vector<Link<n_t, dist_t> > active_links;
	std::vector< std::vector< Link<n_t, dist_t>* > > bas2links;
	std::priority_queue<Link<n_t, dist_t>, std::vector<Link<n_t, dist_t>>, std::greater<Link<n_t, dist_t> > > pqlinks;




	int nbasins;
	int npits = 0;


	LMRerouter(){;};

	LMRerouter(Neighbourer_t& neighbourer,topo_t& topography, std::vector<n_t>& Sreceivers,  std::vector<dist_t>& Sdistance2receivers, std::vector<size_t>& stack)
	{
		auto t1 = high_resolution_clock::now();
		this->compute_basins_and_pits(neighbourer,topography,Sreceivers,Sdistance2receivers, stack);
		auto t2 = high_resolution_clock::now();
		duration<double, std::milli> ms_double = t2 - t1;
		// std::cout << "Computing basins and pits -> " << ms_double.count() << " milliseconds" << std::endl;;
		if(this->npits > 0)
		{
			t1 = high_resolution_clock::now();
			this->preprocess_flowrouting(neighbourer,topography,Sreceivers,Sdistance2receivers);
			t2 = high_resolution_clock::now();
			ms_double = t2 - t1;
			// std::cout << "Preprocess_flowrouting -> " << ms_double.count() << " milliseconds" << std::endl;;
			
		}
	}

	void preprocess_flowrouting(Neighbourer_t& neighbourer,topo_t& topography, std::vector<n_t>& Sreceivers,  std::vector<dist_t>& Sdistance2receivers)
	{
		auto t1 = high_resolution_clock::now();

		this->_compute_links(neighbourer,topography,Sreceivers,Sdistance2receivers);
		auto t2 = high_resolution_clock::now();
		duration<double, std::milli> ms_double = t2 - t1;
		// std::cout << "_compute_links --> " << ms_double.count() << " milliseconds" << std::endl;;
		t1 = high_resolution_clock::now();
		this->_compute_mst_kruskal(neighbourer,topography,Sreceivers,Sdistance2receivers);
		t2 = high_resolution_clock::now();
		ms_double = t2 - t1;
		// std::cout << "_compute_mst_kruskal --> " << ms_double.count() << " milliseconds" << std::endl;;
		t1 = high_resolution_clock::now();
		this->_orient_basin_tree(neighbourer,topography,Sreceivers,Sdistance2receivers);
		t2 = high_resolution_clock::now();
		ms_double = t2 - t1;
		// std::cout << "_orient_basin_tree --> " << ms_double.count() << " milliseconds" << std::endl;;
	}

	// Uses the stack structure to build a quick basin array
	void compute_basins_and_pits(Neighbourer_t& neighbourer,topo_t& topography, std::vector<n_t>& Sreceivers,  std::vector<dist_t>& Sdistance2receivers, std::vector<size_t>& stack)
	{
		this->basin_labels = std::vector<n_t>(neighbourer.nnodes_t, -1);
		this->basin_to_outlets.reserve(200);
		this->pits_to_reroute.reserve(200);
		n_t lab = -1;
		for(auto tnode: stack)
		{
			if(neighbourer.can_flow_even_go_there(tnode) == false)
			{
				// std::cout << "bagul?" << std::endl;
				continue;
			}
			
			if(Sreceivers[tnode] == int(tnode))
			{
				lab++;
				this->basin_to_outlets.emplace_back(tnode);
				if(neighbourer.can_flow_out_there(tnode) == false)
				{
					this->pits_to_reroute.emplace_back(tnode);
					this->is_open_basin.emplace_back(false);
					this->npits++;
				}
				else
					this->is_open_basin.emplace_back(true);


			}

			// if(lab == -1)
			// 	throw std::runtime_error("flonflon?");

			this->basin_labels[tnode] = lab;
		
		}
		
		this->nbasins = lab + 1;
	}

	void _compute_links(Neighbourer_t& neighbourer,topo_t& topography, std::vector<n_t>& Sreceivers,  std::vector<dist_t>& Sdistance2receivers)
	{

		// Initialising a matrix of links for each basins in order to stor the minimum elevation links between each pair of basins
		// std::vector<std::vector<dist_t> > mat_of_scores(this->nbasins,std::vector<dist_t>(this->nbasins,std::numeric_limits<dist_t>::max()));
		// std::vector<std::vector<std::vector<n_t> > > mat_of_nodes(this->nbasins,std::vector<std::vector<n_t> >(this->nbasins));

		std::unordered_map<int, std::unordered_map<int,dist_t > > mat_of_scores;
		std::unordered_map<int, std::unordered_map<int, std::vector<n_t> > > mat_of_nodes;
		auto t1 = high_resolution_clock::now();
		
		// #pragma omp parallel for
		for(int i=0; i<neighbourer.nnodes; ++i)
		{
			if(neighbourer.can_flow_even_go_there(i) == false)
				continue;

			auto neighbours = neighbourer.get_neighbours_only_id(i);
			int tbas = this->basin_labels[i];
			dist_t telev = topography[i];
			for(auto tn:neighbours)
			{
				if(neighbourer.can_flow_even_go_there(tn) == false)
					continue;

				int obas = this->basin_labels[tn];
				if(tbas != obas)
				{
					if(this->is_open_basin[tbas] && this->is_open_basin[obas])
						continue;

					dist_t score = std::min(telev, topography[tn]);
					int basA = (tbas > obas) ? obas : tbas;
					int basB = (tbas > obas) ? tbas : obas;

					bool isinmap = mat_of_scores.count(basA) > 0;
					if(isinmap)
						isinmap = mat_of_scores[basA].count(basB) > 0;

					bool need = false;
					if(isinmap)
					{
						if(mat_of_scores[basA][basB] > score)
						{
							need = true;
						}
					}
					else
						need = true;


					if(need)
					{
						if(basA == tbas)
						{
							std::vector<n_t> temp = std::initializer_list<n_t>{i,tn};
							mat_of_nodes[basA][basB] = temp;
						}
						else
						{
							std::vector<n_t> temp = std::initializer_list<n_t>{tn,i};
							mat_of_nodes[basA][basB] = temp;
						}

						mat_of_scores[basA][basB] = score;
					}
				}
			}
		}
		auto t2 = high_resolution_clock::now();
		duration<double, std::milli> ms_double = t2 - t1;
		// std::cout << "_compute_links::building ---> " << ms_double.count() << " milliseconds" << std::endl;;


		t1 = high_resolution_clock::now();
		for (auto it=mat_of_scores.begin(); it!=mat_of_scores.end(); ++it)
		{
			int basA = it->first;
			for (auto it2=mat_of_scores[basA].begin(); it2!=mat_of_scores[basA].end(); ++it2)
			{
				int basB = it2->first;
				dist_t score = it2->second;
				this->pqlinks.emplace(Link<n_t,dist_t>(basA,basB,mat_of_nodes[basA][basB][0],mat_of_nodes[basA][basB][1],score));
			}


			

		}

		t2 = high_resolution_clock::now();
		ms_double = t2 - t1;


	}

	void _compute_mst_kruskal(Neighbourer_t& neighbourer,topo_t& topography, std::vector<n_t>& Sreceivers,  std::vector<dist_t>& Sdistance2receivers)
	{
		// this->mstree = std::vector<n_t>(this->nbasins - 1);
		this->bas2links = std::vector< std::vector< Link<n_t, dist_t>* > >( this->nbasins, std::vector< Link<n_t, dist_t>* >() );
		// int mstree_size = 0;

		UnionFindv3<n_t,dist_t,Neighbourer_t,topo_t, LMRerouter<n_t,dist_t,Neighbourer_t,topo_t> > uf(nbasins, (*this) );
		int j = 0;

		while(this->pqlinks.empty() == false)
		{
			auto next = this->pqlinks.top();
			this->pqlinks.pop();

			int b0 = next.from;
			int b1 = next.to;
			int fb0 = uf.Find(b0);
			int fb1 = uf.Find(b1) ;
			if (fb0 != fb1)
			{

				if(uf._open[fb0] && uf._open[fb1])
					continue;
				// std::cout << "GOUGN::" << b0 << "|" << b1 << "|" << this->nbasins << std::endl;
				
				uf.Union(b0, b1);
				this->active_links.emplace_back(next);
				++j;
			}
		}

		for(size_t i = 0; i < this->active_links.size(); ++i)
		{
			this->bas2links[this->active_links[i].from].emplace_back(&this->active_links[i]);
			this->bas2links[this->active_links[i].to].emplace_back(&this->active_links[i]);
		}

		// for (int tb = 0 ; tb < this->nbasins; ++tb)
		// for(size_t j=0; j<this->bas2links[tb].size(); ++j)
		// {
		// 	if(this->bas2links[tb][j]->from >= this->nbasins || this->bas2links[tb][j]->to >= this->nbasins)
		// 			throw std::runtime_error("waffon2");
		// }

	}

	void _orient_basin_tree(Neighbourer_t& neighbourer,topo_t& topography, std::vector<n_t>& Sreceivers,  std::vector<dist_t>& Sdistance2receivers)
	{

		this->receivers = std::vector<n_t>(this->nbasins, -1);
		this->nodes_to = std::vector<n_t>(this->nbasins, -1);
		this->nodes_from = std::vector<n_t>(this->nbasins, -1);
		this->outlets_from = std::vector<n_t>(this->nbasins, -1);
		this->stack.clear();this->stack.reserve(this->nbasins);

		std::vector<bool> isdone(this->nbasins,false);



		// std::cout << "here::" << this->active_links.size() << std::endl;		
		for(int i = this->active_links.size() - 1; i >= 0; --i)
		{
			// getting the link
			// auto& tlink = this->active_links[i]; 

			if(this->is_open_basin[this->active_links[i].from] == false && this->is_open_basin[this->active_links[i].to] == false)
				continue;

			if(this->is_open_basin[this->active_links[i].from])
			{
				// std::cout << "BOULOUF" << std::endl;
				this->active_links[i].inverse();
			}

			std::stack<int, std::vector<int> > stackhelper;
			stackhelper.emplace(this->active_links[i].to);

			while(stackhelper.empty() == false)
			{
				int next = stackhelper.top();
				isdone[next] = true;

				stackhelper.pop();
				for(auto ptr: this->bas2links[next])
				{
					// std::cout << "Wafulb1::" << ptr->from << "|" << ptr->to << std::endl;
					if(ptr->from == next)
					{
						if(isdone[ptr->to] == false)
						{
							stackhelper.emplace(ptr->to);
							ptr->inverse();
							// if(isdone[ptr->from] == false)
							// {
								// std::cout << "STACKADD1" << std::endl;
								this->stack.emplace_back(ptr);
								// std::cout << this->is_open_basin[ptr->from] << "||" << this->is_open_basin[ptr->to] << std::endl;
								isdone[ptr->from] = true;
							// }
						}
					}
					else if(isdone[ptr->from] == false)
					{
						stackhelper.emplace(ptr->from);
						if(isdone[ptr->from] == false)
						{
							// std::cout << "STACKADD2" << std::endl;
							this->stack.emplace_back(ptr);
							isdone[ptr->from] = true;
						}
					}
				}
			}
		}

	}

	void compute_TO_SF_stack_version(Neighbourer_t& neighbourer,topo_t& topography, std::vector<n_t>& Sreceivers,  std::vector<dist_t>& Sdistance2receivers)
	{
		// Initialising the stack
		this->stack.clear();
		// reserving the amount of stuff
		this->stack.reserve(this->nbasins);

		// The stack container helper
		std::stack<int, std::vector<int> > stackhelper;
		// std::vector<bool> isdone(this->nbasins,false);
		// going through all the nodes
		for(int i=0; i<this->nbasins; ++i)
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



	}

	void _update_pits_receivers_carve(Neighbourer_t& neighbourer,topo_t& topography, std::vector<n_t>& Sreceivers, 
		 std::vector<dist_t>& Sdistance2receivers, std::vector<size_t>& stack)
	{

		// std::cout <<"ID,rfrom,cfrom,rto,cto,rout,cout" << std::endl;
			
		// int lab=0;
		for(int i=this->stack.size() - 1; i >=0 ; --i)
		// for(int i=0;  i<this->stack.size(); ++i)
		{

			auto tlink = this->stack[i];
			// int i = mstree[ti];
			// int node_to = tlink->node_to;
			int node_from = tlink->node_from;
			int onode_to = tlink->node_to;
			int onode_from = tlink->node_from;
			int outlet_from = this->basin_to_outlets[tlink->from];

			int next_node = Sreceivers[node_from];
			int temp = node_from;
			bool keep_on = true;
			do
			{

				// int this_row,this_col, this_rfrom, this_cfrom;
				// this->graph->rowcol_from_node_id(next_node, this_row, this_col);
				// this->graph->rowcol_from_node_id(node_from, this_rfrom, this_cfrom);

				// std::cout << "This Rerouting " << this_row << "/" <<this_col << " to " << this_rfrom << "/" << this_cfrom << std::endl;


				if(next_node == outlet_from)
				{
					keep_on = false;
				}


				temp = Sreceivers[next_node];

				// if(viz[temp])
				// {
				// 	std::cout << "CYCLICITY ON LINK:" << std::endl;
				// 	std::cout << tlink->from << " to " << tlink->to << std::endl;;
				// 	std::cout << rowf << "|" << colf << " to " << rowt << "|" << colt  << std::endl;;
				// 	throw std::runtime_error("asdfasdfsadfasdf");
				// }
				// viz[next_node] = true;

				// if(temp == next_node || this->graph->can_flow_out_there(next_node))
				// 	keep_on = false;

				Sreceivers[next_node] = node_from;
				Sdistance2receivers[next_node] = Sdistance2receivers[node_from]; // just to have a length but it should not actually be used
				node_from = next_node;
				next_node = temp;
			} while(keep_on);



			Sreceivers[onode_from] = onode_to;
			Sdistance2receivers[onode_from] = neighbourer.dx; // just to have a length but it should not actually be used
		}
		
	}



	void _update_pits_receivers_sompli(Neighbourer_t& neighbourer,topo_t& topography, std::vector<n_t>& Sreceivers, 
		 std::vector<dist_t>& Sdistance2receivers, std::vector<size_t>& stack)
	{
		// for i in mstree:
		for(int i=this->stack.size() - 1; i >=0 ; --i)
		{

			auto tlink = this->stack[i];


			int node_to = tlink->node_to;
			int node_from = tlink->node_from;

			// # skip open basins
			if (node_from == -1)
			{
					continue;
			}

			// int outlet_from = this->basin_to_outlets[conn_basins[i][1] ];
			int outlet_from = this->basin_to_outlets[tlink->from];

			Sreceivers[outlet_from] = node_to;
			Sdistance2receivers[outlet_from] = neighbourer.dx; // just to have a length but it should not actually be used
		}

	 
	}

	void _update_pits_receivers_fill(Neighbourer_t& neighbourer,topo_t& topography, std::vector<n_t>& Sreceivers, 
		 std::vector<dist_t>& Sdistance2receivers, std::vector<size_t>& stack)
	{
		// for i in mstree:
		// std::cout <<"yolo";
		std::vector<bool> isdone(neighbourer.nnodes_t,false);
		// for(int i=this->stack.size() - 1; i >=0 ; --i)
		for(size_t i=0;  i<this->stack.size(); ++i)
		{

			auto tlink = this->stack[i];
			// int i = mstree[ti];
			int node_to = tlink->node_to;
			int node_from = tlink->node_from;
			// int onode_to = tlink->node_to;
			// int onode_from = tlink->node_from;
			// int outlet_from = this->basin_to_outlets[tlink->from];

			int tg_bas = this->basin_labels[node_from];


			// # skip open basins
			if (this->is_open_basin[this->basin_labels[node_to]])
			{
				continue;
			}

			// std::cout << "rec nodde_from beef " << this->graph->receivers[node_from];;
			float otopo = std::fmax(topography[node_to],topography[node_from]);
			std::queue<int> yonode;
			yonode.push(node_to);
			while(yonode.size() > 0)
			{
				int tnode = yonode.front(); yonode.pop();
	 			auto neighbours = neighbourer.get_neighbours(tnode, false);
				for(auto& ineighbor:neighbours)
				{
					if(neighbourer.can_flow_even_go_there(ineighbor.node) == false || tg_bas != this->basin_labels[ineighbor.node] || isdone[ineighbor.node] == true)
						continue;

					if(topography[ineighbor.node] <= otopo)
					{
						// std::cout << ineighbor.node << "<--" << tnode << "|||";
						Sreceivers[ineighbor.node] = tnode;
						Sdistance2receivers[ineighbor.node] = ineighbor.distance;
						yonode.push(ineighbor.node);
						isdone[ineighbor.node] = true;
					}
				}

			}
		}
		
	}



	// bool is_open(int i)
	// {
	// 	return this->graph->can_flow_out_there(this->outlets_from[i]);
	// }


	void update_receivers(std::string& method,Neighbourer_t& neighbourer,topo_t& topography, std::vector<n_t>& Sreceivers,  std::vector<dist_t>& Sdistance2receivers, std::vector<size_t>& stack)
	{

		if(method == "simple" || method == "Simple")
			this->_update_pits_receivers_sompli(neighbourer,topography,Sreceivers,Sdistance2receivers,stack);
		else if (method == "carve")
		{
			this->_update_pits_receivers_carve(neighbourer,topography,Sreceivers,Sdistance2receivers,stack);
		}
		else if (method == "fill")
		{
			this->_update_pits_receivers_fill(neighbourer,topography,Sreceivers,Sdistance2receivers,stack);
		}
	}


};

// Only for pairs of std::hash-able types for simplicity.
// You can of course template this struct to allow other hash functions
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        // Mainly for demonstration purposes, i.e. works but is overly simple
        // In the real world, use sth. like boost.hash_combine
        return h1 ^ h2;  
    }
};



class LMRerouter_II
{

public:

	int nbas;

	// node size
	std::vector<int> basins;

	// nbasins size
	std::vector<bool> is_open_basin;
	std::vector<int> receivers;
	std::vector<int> pitnode;
	std::vector<std::pair<int,int> > receivers_node;
	std::vector<int> stack;
	std::vector<std::vector<int> > donors;


	//
	std::unordered_map<std::pair<int,int> , double, pair_hash> edges;
	std::unordered_map<std::pair<int,int> , std::pair<int,int>, pair_hash > edges_nodes;

	LMRerouter_II(){};

	template<class topo_t, class Neighbourer_t>
	bool run(std::string method, topo_t& topography, Neighbourer_t& neighbourer, std::vector<int>& Sreceivers, std::vector<double>& Sdistance2receivers, std::vector<size_t>& Sstack, std::vector<int>& links)
	{
		// std::cout << "DEBUGLM_II::1" <<std::endl;
		// tracking the number of basins
		this->nbas = -1;
		// tracking to which basin each node belongs to
		this->basins = std::vector<int>(neighbourer.nnodes,0);
		// number of internal basins to reroute
		int nbas2solve = 0;

		// First comuting the basin array
		for(int i=0; i<neighbourer.nnodes; ++i)
		{
			// getting the current nodes
			size_t node = Sstack[i];

			// If it is its own receiver, then
			if(Sreceivers[node] == int(node))
			{
				// incrementing the basin label
				++this->nbas;
				// is it a base level
				if(neighbourer.can_flow_out_there(node))
				{
					// then it is an open basin
					this->is_open_basin.emplace_back(true);
					// saving its pit node
					this->pitnode.emplace_back(node);
				}
				else
				{
					// otherwise it is a basin to solve
					++nbas2solve;
					// not open
					this->is_open_basin.emplace_back(false);
					// saving its pit node
					this->pitnode.emplace_back(node);
				}
			}

			// labelling the node
			this->basins[node] = this->nbas;
		}
		
		// need a last increment as first label is 0
		++this->nbas;


		// std::cout << "DEBUGLM_II::2" <<std::endl;

		// Relabelling 0 all the open basins to gain time
		for(int i=0; i<neighbourer.nnodes; ++i)
		{
			if(this->is_open_basin[this->basins[i]]) this->basins[i] = 0;
		}

		// std::cout << "DEBUGLM_II::3" <<std::endl;
		// if there is literally no basins to solve, then I am done
		if(nbas2solve == 0)
			return false;

		// tracking the number of links between a basin to another
		int nlinks = 0;

		// std::cout << "DEBUGLM_II::4" <<std::endl;
		
		// going through each and every link
		for(int i=0; i < int(links.size()); ++i)
		{
			// if the link is not valid, I continue
			if(links[i] < 0)
			{
				// REALLY IMPORTANT: incrementing i as i is a noe and i+1 its counterpart
				++i;
				continue;
			}

			// j is first index and k the next
			int j = i;
			int k = j+1;

			// REALLY IMPORTANT: incrementing i as i is a noe and i+1 its counterpart
			++i;
			
			// translating the nodes to basin IDs 
			int bj = this->basins[links[j]];
			int bk = this->basins[links[k]];

			// if in same basin or both open -> I skip
			if(bj == bk || (this->is_open_basin[bj] && this->is_open_basin[bk]) )
				continue;

			// The score is the minimum elevation of the pass
			double score = std::min(topography[links[j]],topography[links[k]]);
			// is bj < bk (the std::pair storing the pass always starts from the lowes to the highest by convention to keep the std::pair map keys unique)
			bool bjmin = bj<bk;
			// if (bj<bk) pair is {bj,bk} else {bk,bj} (I love ternary operators)
			std::pair<int,int> tp = {(bjmin)?bj:bk, (bjmin)?bk:bj};
			// is the pair already in the map-e
			auto it_e = this->edges.find(tp);
			if(it_e == this->edges.end())
			{
				// Nope.
				// counting that link
				++nlinks;
				// registering the elev of the pass
				this->edges[tp] = score;
				// registering nodes of the pass
				this->edges_nodes[tp] = std::pair<int,int>{(bjmin)?links[j]:links[k], (bjmin)?links[k]:links[j]};
			}
			else
			{
				// The basins are already connected
				// Checking if the current connection is lower in Z
				if(score < it_e->second)
				{
					// It is!
					// registering the new score ...
					it_e->second = score;
					// ... and nodes
					this->edges_nodes[it_e->first] = std::pair<int,int>{(bjmin)?links[j]:links[k], (bjmin)?links[k]:links[j]};
				}
			}
		}
		// Done with the link construction

		// std::cout << "DEBUGLM_II::5" <<std::endl;

		// Gathering all the links in a vector
		std::vector<PQ_helper<std::pair<int,int>, double > > basinlinks;basinlinks.reserve(nlinks);
		for(auto it: this->edges)
		{
			basinlinks.emplace_back(PQ_helper<std::pair<int,int>, double >(it.first,it.second));
		}

		// And sorting it
		std::sort(basinlinks.begin(), basinlinks.end());

		// This will track which links are active or not
		std::vector<bool> isactive(basinlinks.size(), false);


		// std::cout << "DEBUGLM_II::6" <<std::endl;


		// This part is applying the kruskal algorithm (I think)
		UnionFindv3<int, double, Neighbourer_t,topo_t, LMRerouter_II> uf(this->nbas, (*this) );

		// trackng the receiver of all basins
		this->receivers = std::vector<int>(this->nbas);
		// and the subsequent node pair
		this->receivers_node = std::vector<std::pair<int,int> >(this->nbas);
		// init recs to themselves (base level)
		for(int i =0; i<this->nbas; ++i)
			this->receivers[i] = i;

		// ok going through all the links from lowest pass to the highest
		for(size_t i = 0; i < basinlinks.size(); ++i)
		{
			// getting next link
			auto& next = basinlinks[i];

			// getting basin IDs
			int b1 = next.node.first;
			int b2 = next.node.second;

			// getting basin IDs unionised ☭
			int fb1 = uf.Find(b1);
			int fb2 = uf.Find(b2);

			// If they are united, I skip (they already merged)
			if (fb1 != fb2)
			{

				// if both are open, I skip
				if(uf._open[fb1] && uf._open[fb2])
					continue;
				
				// Unification of both basin				
				uf.Union(b1, b2);
				// this link is active then
				isactive[i] = true;

				// // If basin one is open
				// if(this->is_open_basin[b1])
				// {
				// 	std::cout << "gulg::1" << this->is_open_basin[b2]  << std::endl;
				// 	// rec of b2 is b1
				// 	this->receivers[b2] = b1;
				// 	// connecting node are node b2 to node b1
				// 	this->receivers_node[b2] = std::pair<int,int>{this->edges_nodes[next.node].second ,this->edges_nodes[next.node].first};
				// 	// b2 is now open
				// 	this->is_open_basin[b2] = true;
				// }
				// else if(this->is_open_basin[b2])
				// {
				// 	std::cout << "gulg::2" << this->is_open_basin[b1] << std::endl;
				// 	this->receivers[b1] = b2;
				// 	this->receivers_node[b1] = std::pair<int,int>{this->edges_nodes[next.node].first ,this->edges_nodes[next.node].second};
				// 	this->is_open_basin[b1] = true;
				// }
			}
		}

		// std::cout << "DEBUGLM_II::7" <<std::endl;

		while(true)
		{
			bool alltrue = true;
			for(size_t i=0; i<basinlinks.size();++i)
			{
				if(isactive[i] == false)
					continue;

				int b1 = basinlinks[i].node.first, b2 = basinlinks[i].node.second;
				if(this->is_open_basin[b1] && this->is_open_basin[b2])
					continue;
				auto& next =  basinlinks[i];

				// std::cout << "bulf";

				if(this->is_open_basin[b1])
				{
					// std::cout << "pluf" << std::endl;
					this->receivers[b2] = b1;
					this->receivers_node[b2] = std::pair<int,int>{this->edges_nodes[next.node].second ,this->edges_nodes[next.node].first};
					this->is_open_basin[b2] = true;
				}
				else if(this->is_open_basin[b2])
				{
					// std::cout << "pluf" << std::endl;
					this->receivers[b1] = b2;
					this->receivers_node[b1] = std::pair<int,int>{this->edges_nodes[next.node].first ,this->edges_nodes[next.node].second};
					this->is_open_basin[b1] = true;
				}
				else
					alltrue = false;
			}

			if(alltrue)
				break;
		}

		// std::cout << "DEBUGLM_II::8" <<std::endl;


		this->donors = std::vector<std::vector<int> >(this->nbas, std::vector<int>());
		for(int i=0; i<this->nbas; ++i)
		{
			if(this->receivers[i] != i)
			{
				this->donors[this->receivers[i]].emplace_back(i);
			}
		}

		// std::cout << "DEBUGLM_II::9" <<std::endl;
		this->compute_TO_SF_stack_version();

		// std::cout << "DEBUGLM_II::10::" << this->stack.size() <<std::endl;

		if(method == "carve")
		{
			for(int i =  this->nbas-1; i>=0; --i)
			{
				int bas = this->stack[i];
				// std::cout << bas << "/" << this->nbas << std::endl;;
				if(neighbourer.can_flow_out_there(this->pitnode[bas]))
					continue;
				// std::cout << "A" << std::endl;
				int from = this->receivers_node[bas].first; 
				int to = this->receivers_node[bas].second;
				// std::cout << "B" << std::endl;
				// std::cout << Sreceivers[this->pitnode[bas]] << "|";

				int A = from;
				int B = Sreceivers[A];
				int C = B;
				// std::cout << "C" << std::endl;

				while(A != this->pitnode[bas])
				{
					// std::cout << B << std::endl;
					C = Sreceivers[B];
					Sreceivers[B] = A;

					A = B;
					B = C;
				}
				// std::cout << "D" << std::endl;

				Sreceivers[from] = to;

				// std::cout << Sreceivers[this->pitnode[bas]] << std::endl;
			}

		}
		else if (method == "fill")
		{
			std::vector<char> isinQ(neighbourer.nnodes,false);
			std::vector<char> isfilled(neighbourer.nnodes,false);
			std::vector<bool> basinDone(this->nbas,false);
			std::vector<int> basfam(this->nbas,-1);
			for(int i = 0; i< this->nbas; ++i)
			{
				if(neighbourer.can_flow_out_there(this->pitnode[i]))
					basinDone[i] = true;

				int node = this->stack[i];
				if(this->receivers[node] != node)
					basfam[node] = basfam[this->receivers[node]];
				else
					basfam[node] = node;
			}


			for(int i = 0; i < this->nbas; ++i)
			{
				int bas = this->stack[i];
				if(neighbourer.can_flow_out_there(this->pitnode[bas]))
					continue;
				int from = this->receivers_node[bas].first; 
				int to = this->receivers_node[bas].second;
				double zref = std::max(topography[from], topography[to]);
				Sreceivers[from] = to;
				Sdistance2receivers[from] = neighbourer.dx;
				isinQ[from] = true;
				std::queue<int> Q;Q.emplace(from);
				while(Q.empty() == false)
				{
					int next = Q.front();Q.pop();
					isfilled[next] = true;
					auto neighs = neighbourer.get_neighbours_only_id(next);
					double lowest_z = std::max(topography[Sreceivers[next]],topography[next]);
					int nznodeext = Sreceivers[next];
					for(auto n : neighs )
					{
						int basn = this->basins[n];

						if(isfilled[n] || basinDone[basn] || basfam[basn] != basfam[bas])
						{
							if(lowest_z > topography[n])
							{
								lowest_z = topography[n];
								nznodeext = n;
							}
						}

						if(isinQ[n])
							continue;
						// if(basfam[basn] != basfam[bas])
						if(basn != bas)
							continue;
						if(basinDone[basn])
							continue;
						if(topography[n] <= zref)
						{
							Sreceivers[n] = next;
							isinQ[n] = true;
							Q.emplace(n);
						}
					}

					topography[next] = std::max(lowest_z + 1e-4 + neighbourer.randu.get() * 1e-6, topography[next]);
					zref = std::max(topography[next],zref);
					Sreceivers[next] = nznodeext;
					Sdistance2receivers[next] = neighbourer.dx;
				}
				basinDone[bas] = true;
			}
		}

		// std::cout << "DEBUGLM_II::11" <<std::endl;

		return true;




	}

	void compute_TO_SF_stack_version()
	{
		// Initialising the stack
		this->stack.clear();
		// reserving the amount of stuff
		this->stack.reserve(this->nbas);

		// The stack container helper
		std::stack<int, std::vector<int> > stackhelper;
		// std::vector<bool> isdone(this->nbas,false);
		// going through all the nodes
		for(int i=0; i<this->nbas; ++i)
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
				for(size_t j = 0; j < this->donors[nextnode].size(); ++j)
				{
					stackhelper.emplace(this->donors[nextnode][j]);
				}

			}

		}

		// if(this->nbas != this->stack.size())
		// 	throw std::runtime_error("stacksize issue in LMRerouter_II::" + std::to_string(this->stack.size()) + " vs " + std::to_string(this->stack.size()));

		// for(auto v:this->stack)
		// 	std::cout << v << "|";
		// std::cout << std::endl;

	}


};





























































#endif