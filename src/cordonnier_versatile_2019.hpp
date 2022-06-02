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

#include "chonkutils.hpp"
// #include "graph.hpp"

class Graph;

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

	void inverse(){*this = Link(this->to,this->from,this->node_to,this->node_from,this->score);}

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
			std::cout << "wololo-2" << std::endl;

			sleep(1);
			std::cout << "wololo-1" << std::endl;
			this->graph = &graph;
			std::cout << "wololo0" << std::endl;
			this->compute_basins_and_pits();
			std::cout << "wololo1" << std::endl;
			this->preprocess_flowrouting(elevation);
			std::cout << "wololo2" << std::endl;

		}

		void update_receivers(std::string& method, std::vector<dist_t>& elevation)
		{
			if(method == "simple" || method == "Simple")
				this->_update_pits_receivers_sompli(this->conn_basins, this->conn_nodes, mstree, elevation);
			else if (method == "carve")
			{
				std::cout << "karper" << std::endl;
				this->_update_pits_receivers_carve(this->conn_basins, this->conn_nodes, mstree, elevation);
				std::cout << "ski" << std::endl;
			}
			else if (method == "fill")
			{
				// std::cout << "gwamoulg" << std::endl;
				this->_update_pits_receivers_fill(this->conn_basins, this->conn_nodes, mstree, elevation);
			}
			std::cout << "fabul" << std::endl;
			this->graph->recompute_SF_donors_from_receivers();
			std::cout << "fabudl" << std::endl;
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
		std::cout << "connect basins" << std::endl;
		this->_connect_basins(this->conn_basins, this->conn_nodes, conn_weights, elevation, nconn, basin0);
		std::cout << "connect bagul" << std::endl;
		
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

		int lab = 1;
		for(auto i : mstree)
		// for(int ti = mstree.size() - 1 ; ti >=0; --ti)
		{
			// int i = mstree[ti];
			int node_to = this->conn_nodes[i][0];
			int node_from = this->conn_nodes[i][1];
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



			this->graph->receivers[this->conn_nodes[i][1]] = node_to;
			this->graph->distance2receivers[this->conn_nodes[i][1]] = this->graph->dx; // just to have a length but it should not actually be used
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
			int outlet_from = this->basin_to_outlets[conn_basins[i][1] ];
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
	std::vector<n_t> basin_to_outlets;
	std::vector<n_t> pits_to_reroute;
	std::vector<n_t> mstree;
	std::vector<n_t> receivers, stack, nodes_to, nodes_from, outlets_from;
	std::vector<std::vector<n_t> > donors;


	std::vector<Link<n_t, dist_t> > active_links;
	std::vector< std::vector< *Link<n_t, dist_t> > > bas2links;
	std::priority_queue<Link<n_t, dist_t>, std::vector<Link<n_t, dist_t>>, std::greater<Link<n_t, dist_t> > > pqlinks;




	int nbasins;
	int npits = 0;


	Cordonnier2019_v2(){;};
	Cordonnier2019_v2(Graph& graph)
	{
		this->graph = &graph;
		std::cout << "Compute basepit" << std::endl;
		this->compute_basins_and_pits();
		this->preprocess_flowrouting();
	}

	void preprocess_flowrouting()
	{
		std::cout << "Compute links" << std::endl;
		this->_compute_links();
		std::cout << "Kurskal" << std::endl;
		this->_compute_mst_kruskal();
		std::cout << "btree links" << std::endl;
		this->_orient_basin_tree();
		std::cout << "dolone" << std::endl;
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
		
		this->nbasins = lab + 1;
	}

	void _compute_links()
	{

		// Initialising a matrix of links for each basins in order to stor the minimum elevation links between each pair of basins
		std::vector<std::vector<dist_t> > mat_of_scores(this->nbasins,std::vector<dist_t>(this->nbasins,std::numeric_limits<dist_t>::max()));
		std::vector<std::vector<std::vector<n_t> > > mat_of_nodes(this->nbasins,std::vector<std::vector<n_t> >(this->nbasins));

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
					dist_t score = std::min(telev, this->graph->topography[tn]);
					int basA = (tbas > obas) ? obas : tbas;
					int basB = (tbas > obas) ? tbas : obas;
					if(mat_of_scores[basA][basB] > score)
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
						// std::cout << "happens??" << std::endl;

					}
				}
			}
		}

		dist_t maxval = std::numeric_limits<dist_t>::max();
		for(int i = 0; i < this->nbasins;++i)
		{
			for(int j=0; j<this->nbasins;++j)
			{

				if(j <= i)
					continue;

				if(mat_of_scores[i][j] ==	maxval)
					continue;

				this->pqlinks.emplace(Link<n_t,dist_t>(i,j,mat_of_nodes[i][j][0],mat_of_nodes[i][j][1],mat_of_scores[i][j]));
			}
		}

	}

	void _compute_mst_kruskal()
	{
		// this->mstree = std::vector<n_t>(this->nbasins - 1);
		this->bas2links = std::vector< std::vector< *Link<n_t, dist_t> > >( this->nbasins, std::vector< *Link<n_t, dist_t> >() );
		int mstree_size = 0;

		UnionFind uf(nbasins);

		while(this->pqlinks.empty() == false)
		{
			auto next = this->pqlinks.top();
			this->pqlinks.pop();


			int b0 = next.from;
			int b1 = next.to;
			std::cout << "b0->" << b0 << std::endl;;

			if (uf.Find(b0) != uf.Find(b1))
			{
				// mstree[mstree_size] = eid;
				// mstree_size ++;
				uf.Union(b0, b1);
				this->active_links.emplace_back(next);
				this->bas2links[b0].emplace_back(&this->active_links.back());
				this->bas2links[b1].emplace_back(&this->active_links.back());
			}
		}
	}

	void _orient_basin_tree()
	{

		this->receivers = std::vector<n_t>(this->nbasins, -1);
		this->nodes_to = std::vector<n_t>(this->nbasins, -1);
		this->nodes_from = std::vector<n_t>(this->nbasins, -1);
		this->outlets_from = std::vector<n_t>(this->nbasins, -1);
		this->stack.clear();this->stack.reserve(this->nbasins);

		std::vector<bool> isdone(this->nbasins,false);

		std::cout << "here::" << this->active_links.size() << std::endl;		
		for(int i = this->active_links.size() - 1; i >= 0; --i)
		{
			// getting the link
			auto& tlink = this->active_links[i];

			if(this->graph->can_flow_out_there(this->basin_to_outlets[tlink.from]) == false && this->graph->can_flow_out_there(this->basin_to_outlets[tlink.to]) == false)
				continue;

			if(this->graph->can_flow_out_there(this->basin_to_outlets[tlink.to]) == false)
				tlink.inverse();

			isdone[tlink.to] = true;

			std::stack<int, std::vector<int> > stackhelper;
			stackhelper.emplace(tlink.from);

			while(stackhelper.empty() == false)
			{
				int next = stackhelper.top();
				isdone[next] = true;

				stackhelper.pop();
				for(auto ptr: this->bas2links[next])
				{
					if(ptr->from == next)
					{
						if(isdone[ptr->to] == false)
						{
							stackhelper.emplace(ptr->to);
							ptr->inverse();
						}
					}
					else if(isdone[ptr->from] == false)
						stackhelper.emplace(ptr->from);
				}
			}
		}





		// 	std::cout << i << "::" << tlink.node_from << "|" << tlink.node_to << std::endl;
		// 	if(this->graph->can_flow_out_there(tlink.node_from) == false && this->graph->can_flow_out_there(tlink.node_to) == false)
		// 	{
		// 		if(this->receivers[tlink.from] != -1)
		// 		{
		// 			this->receivers[tlink.to] = tlink.from;
		// 			this->nodes_to[tlink.to] = tlink.node_from;
		// 			this->nodes_from[tlink.to] = tlink.node_to;
		// 			this->outlets_from[tlink.to] = this->basin_to_outlets[tlink.to];
		// 		}
		// 		else if (this->receivers[tlink.to] != -1)
		// 		{
		// 			this->receivers[tlink.from] = tlink.to;
		// 			this->nodes_to[tlink.from] = tlink.node_to;
		// 			this->nodes_from[tlink.from] = tlink.node_from;
		// 			this->outlets_from[tlink.from] = this->basin_to_outlets[tlink.from];
		// 		}
		// 		else
		// 			throw std::runtime_error("Erm, did not work...");
		// 	}

		// 	else
		// 	{
		// 		if(this->graph->can_flow_out_there(tlink.from))
		// 		{
		// 			this->receivers[tlink.from] = tlink.from;
		// 			this->receivers[tlink.to] = tlink.from;

		// 			this->outlets_from[tlink.from] = this->basin_to_outlets[tlink.from];
		// 			this->outlets_from[tlink.to] = this->basin_to_outlets[tlink.to];

		// 		}
		// 		else
		// 		{
		// 			this->receivers[tlink.to] = tlink.to;
		// 			this->receivers[tlink.from] = tlink.to;
		// 			this->receivers[tlink.from] = tlink.to;
		// 			this->nodes_to[tlink.from] = tlink.node_to;
		// 			this->nodes_from[tlink.from] = tlink.node_from;
		// 			this->outlets_from[tlink.from] = this->basin_to_outlets[tlink.from];
		// 			this->outlets_from[tlink.to] = this->basin_to_outlets[tlink.to];
		// 		}

		// 	}
		// }

		// for(int i=0; i< this->nbasins; ++i)
		// {
		// 	if(this->receivers[i] == i)
		// 		continue;
		// 	this->donors[this->receivers[i]].emplace_back(i);
		// }


		// this->compute_TO_SF_stack_version();



		// // # nodes connections
		// std::vector<int> nodes_connects_size(this->nbasins,0);
		// std::vector<int> nodes_connects_ptr(this->nbasins);

		// // # parse the edges to compute the number of edges per node
		// for (auto i : tree)
		// {
		// 	nodes_connects_size[ this->conn_basins[i][0] ]++;
		// 	nodes_connects_size[ this->conn_basins[i][1] ]++;
		// }

		// // # compute the id of first edge in adjacency table
		// nodes_connects_ptr[0] = 0;

		// //lsdkfjlsjd
		// // gruklb.
		// for (int i = 1; i < nbasins; i++)
		// {
		// 	nodes_connects_ptr[i] = (nodes_connects_ptr[i - 1] + nodes_connects_size[i - 1]);
		// 	nodes_connects_size[i - 1] = 0;
		// }

		// // # create the adjacency table
		// int nodes_adjacency_size = nodes_connects_ptr[nbasins - 1] + nodes_connects_size[nbasins - 1];
		// nodes_connects_size[this->nbasins -1] = 0;
		// std::vector<int> nodes_adjacency(nodes_adjacency_size,0);

		// // # parse the edges to update the adjacency
		// for (auto i : tree)
		// {

		// 	int n1 = this->conn_basins[i][0];
		// 	int n2 = this->conn_basins[i][1];
		// 	nodes_adjacency[nodes_connects_ptr[n1] + nodes_connects_size[n1]] = i;
		// 	nodes_adjacency[nodes_connects_ptr[n2] + nodes_connects_size[n2]] = i;
		// 	nodes_connects_size[n1] ++;
		// 	nodes_connects_size[n2] ++;
		// }


		// // # depth-first parse of the tree, starting from basin0
		// // # stack of node, parent
		// std::vector<std::vector<n_t> > stack; stack.reserve(nbasins);
		// for(int i=0; i < nbasins; i++)
		// 	stack.emplace_back(std::vector<n_t>({-2,-2}));


		// int stack_size = 1;
		// stack[0][0] = basin0;// (basin0, basin0)
		// stack[0][1] = basin0;

		// int n_turn = 0;
		// while (stack_size > 0)
		// {
		// 	n_turn++;
		// 	// # get parsed node
		// 	stack_size = stack_size - 1;
		// 	int node = stack[stack_size][0];
		// 	int parent = stack[stack_size][1];

		// 	// # for each edge of the graph
		// 	// for i in range(nodes_connects_ptr[node], nodes_connects_ptr[node] + nodes_connects_size[node])
		// 	for( int i = nodes_connects_ptr[node]; i < (nodes_connects_ptr[node] + nodes_connects_size[node]); i++) 
		// 	{
		// 		int edge_id = nodes_adjacency[i];

		// 		// # the edge comming from the parent node has already been updated.
		// 		// # in this case, the edge is (parent, node)
		// 		if (this->conn_basins[edge_id][0] == parent && node != parent)
		// 		{
		// 				continue;
		// 		}
		// 		// # we want the edge to be (node, next)
		// 		// # we check if the first node of the edge is not "node"
		// 		if(node != this->conn_basins[edge_id][0])
		// 		{
		// 			// # swap n1 and n2
		// 			int cb1 = this->conn_basins[edge_id][1];
		// 			int cb0 = this->conn_basins[edge_id][0];
		// 			this->conn_basins[edge_id][0] = cb1;
		// 			this->conn_basins[edge_id][1] = cb0;

		// 			cb1 = this->conn_nodes[edge_id][1];
		// 			cb0 = this->conn_nodes[edge_id][0];
		// 			// # swap p1 and p2
		// 			this->conn_nodes[edge_id][0] = cb1;
		// 			this->conn_nodes[edge_id][1] = cb0;
		// 		}
		// 		// # add the opposite node to the stack
		// 		stack[stack_size][0] = this->conn_basins[edge_id][1];
		// 		stack[stack_size][1] = node;
		// 		stack_size ++;
		// 	}
		// }

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

		for(auto i : this->active_links)
		{
			// int i = mstree[ti];
			int node_to = i.node_to;
			int node_from = i.node_from;
			int outlet_from = this->basin_to_outlets[i.from];

			// // # skip open basins
			// if (this->is_open(i))
			// {
			// 	continue;
			// }

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
				this->graph->distance2receivers[next_node] = this->graph->distance2receivers[node_from]; // just to have a length but it should not actually be used
				node_from = next_node;
				next_node = temp;
				// std::cout << next_node << "|";
			} while(keep_on);



			this->graph->receivers[node_from] = node_to;
			this->graph->distance2receivers[node_from] = this->graph->dx; // just to have a length but it should not actually be used
		}
		
	}


	// bool is_open(int i)
	// {
	// 	return this->graph->can_flow_out_there(this->outlets_from[i]);
	// }


	void update_receivers(std::string& method)
	{
		// if(method == "simple" || method == "Simple")
		// 	this->_update_pits_receivers_sompli(this->conn_basins, this->conn_nodes, mstree, elevation);
		// else if (method == "carve")
		// {
		// 	std::cout << "karper" << std::endl;
			this->_update_pits_receivers_carve();
		// 	std::cout << "ski" << std::endl;
		// }
		// else if (method == "fill")
		// {
		// 	// std::cout << "gwamoulg" << std::endl;
		// 	this->_update_pits_receivers_fill(this->conn_basins, this->conn_nodes, mstree, elevation);
		// }
		// std::cout << "fabul" << std::endl;
		this->graph->recompute_SF_donors_from_receivers();
		// std::cout << "fabudl" << std::endl;
	}


};





















































































#endif