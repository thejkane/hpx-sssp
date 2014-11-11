/**
 * This is the compressed sparse matrix representation of the graph for HPX-3.
 * The row index array, column array, vertex property array and edge property
 * array data structures are allocated in AGAS.
 * 
 * Authors : Thejaka Kaanewala
 *           Marcin Zaleswski
 *           Andrew Lumsdaine
 **/

#ifndef HPX_CSR_GRAPH
#define HPX_CSR_GRAPH

#include <limits>
#include <stdint.h>
#include <cstddef>
#include <assert.h>
#include <iostream>
#include <iterator>
#include <set>
#include <list>
#include <map>
#include <vector>

#include "distributed_control.hpp"

#define INVALID_VERTEX -1
#define INVALID_EDGE -1
#define INVALID_EDGE_WEIGHT -1

#ifdef GRAPH_64
typedef int64_t edge_t;
typedef int64_t vertex_t;
typedef int64_t index_t;
typedef int64_t edge_property_t;
typedef int64_t vertex_property_t;

#else
typedef int32_t edge_t;
typedef int32_t vertex_t;
typedef int32_t index_t;
typedef int32_t edge_property_t;
typedef int32_t vertex_property_t;
#endif

typedef std::iterator<std::forward_iterator_tag, int> EdgeIterator_t;

//===========================================
// This is needed to sort edges
//===========================================
struct edge_comparator {
  bool operator() (const EdgeType_t& e1, const EdgeType_t& e2) const {
    if (e1.first < e2.first)
      return true;
    else if (e1.first == e2.first) {
      return e1.second < e2.second;
    } else {
      return false;
    }
  }
};

//===========================================
// HPX CSR Implementation. 
//===========================================
class hpx_csr_graph {

public:  
  //===========================================
  // To iterate through all vertices. 
  //===========================================
  class vertex_iterator : public std::iterator<std::input_iterator_tag, vertex_t> {
  private:
    vertex_t vertex_array;

  public:
    vertex_iterator(vertex_t vertices) :vertex_array(vertices) {}
    vertex_iterator(const vertex_iterator& vit) : vertex_array(vit.vertex_array) {}
    vertex_iterator& operator++() {++vertex_array;return *this;}
    vertex_iterator operator++(int) {vertex_iterator tmp(*this); operator++(); return tmp;}
    bool operator==(const vertex_iterator& rhs) {return vertex_array == rhs.vertex_array;}
    bool operator!=(const vertex_iterator& rhs) {return vertex_array != rhs.vertex_array;}
    vertex_t& operator*() {return vertex_array;}
  };
  //===========================================

  //===========================================
  // To iterate through all edges
  //===========================================
  class edge_iterator : 
    public std::iterator<std::input_iterator_tag, EdgeType_t> {
  private:
    edge_t* row_index;
    edge_t* columns;
    edge_t ri; // the progressing row index
    edge_t ci; // the progressing colomn index

  public:
    edge_iterator(edge_t* row_i, edge_t* r, int rowind, int colind) : 
      row_index(row_i), 
      columns(r), 
      ri(rowind), 
      ci(colind) 
    {}

    edge_iterator(const edge_iterator& eit) : 
      row_index(eit.row_index), 
      columns(eit.columns), 
      ri(eit.ri), 
      ci(eit.ci) 
    {}

    void print() {
      std::cout << "Row index : " << ri 
		<< " Column index : " 
		<< ci 
		<< std::endl;
    }

    edge_iterator& operator++() {
      edge_t ci_index_end = *(row_index+ri+1);
      if (ci == (ci_index_end-1)) {
	++ri;
	++ci;
      } else {
	++ci;
      }
      return *this;
    }

    edge_iterator operator++(int) {
      edge_iterator tmp(*this); operator++(); return tmp;
    }

    bool operator==(const edge_iterator& rhs) {
      return ((row_index == rhs.row_index) && 
	      (columns == rhs.columns) && (ci == rhs.ci));
    }

    bool operator!=(const edge_iterator& rhs) {
      return ((row_index != rhs.row_index) ||
	      (columns != rhs.columns) || (ci != rhs.ci)); 
    }

    EdgeType_t operator*() {
      while(row_index[ri] == -1) {
	assert(false);
      }
      edge_t val = *(columns+ci);
      return EdgeType_t(ri, val, ci); 
    }
  };
  //===========================================

  //===========================================
  // The constructor.
  //===========================================
  hpx_csr_graph(std::size_t num_vertices,
		std::size_t num_edges, 
		bool is_undirected): vertices(num_vertices),
				     edges(num_edges),
				     undirected(is_undirected){
    if (undirected)
      edges = num_edges * 2;

    // for row indices we need vertices+1
    row_indices.resize(vertices+1); // use hpx-alloc
    columns.resize(edges);

    // TODO we are wasting space. need to figure out a better approach for undirected
    // graphs
    edge_weight_map.resize(edges); 
    vertex_distance_map.resize(vertices);

    init();
  }

  ~hpx_csr_graph(){
  }

  void buildHistogram(HistogramMap_t& histogram_map, 
		      vertex_t source,
		      vertex_t target, 
		      edge_property_t weight, 
		      bool flipped) {

    HistogramMap_t::iterator iteFind = histogram_map.find(source);
    if (iteFind == histogram_map.end()) { // new source
      EdgeList_t target_list;
      target_list.insert(target_weight(target, weight));
      histogram_map.insert(std::make_pair(source, target_list));
    } else { // source already exists
      (*iteFind).second.insert(target_weight(target, weight));
    }
    
    if (undirected) {
      if (!flipped) {
	buildHistogram(histogram_map, target, source, weight, true);
      }
    }

  }

  // Assumes iterator also includes reversed edges for undirected graphs.
  template <typename RandomAccessIterator, typename EdgePropertyIterator>
  void addEdges(RandomAccessIterator begin, RandomAccessIterator end,
		EdgePropertyIterator epiter) {

    HistogramMap_t histogram_map;
    for (; begin != end; ++begin, ++epiter) {
      //      std::cout << "Before histogram : (" 
      //<< (*begin).first << ", " << (*begin).second
      //		<< ")" << std::endl;
      buildHistogram(histogram_map, (*begin).first, 
		     (*begin).second, *epiter, false); 
      // flipped=false => edge need to flip for undirected
    }

    edge_t row_ind = 1; // starts with 1
    edge_t col_ind = 0;
    row_indices[0] = 0;
    // check whether edges are sorted
    HistogramMap_t::iterator ite = histogram_map.begin();
    for(; ite != histogram_map.end(); ++ite) {
      EdgeList_t list = (*ite).second;
      row_indices[row_ind] = row_indices[row_ind-1] + list.size();
      ++row_ind;

      EdgeList_t::iterator iteList = list.begin();
      for (; iteList != list.end(); ++iteList) {
	columns[col_ind] = (*iteList).target;
	edge_weight_map[col_ind] = (*iteList).weight;
	++col_ind;
      }
    }
    assert(edges == col_ind);

    // Need to handle vertices that doesnt have
    // edges. For rest of the edges put last updated
    // value.
    vertex_t last_seen_id = row_indices[row_ind-1];
    for (int k=row_ind; k<(vertices+1); ++k) {
      row_indices[k] = last_seen_id;
    }
  }

  vertex_iterator vertices_begin() {
    return vertex_iterator(0);
  }

  vertex_iterator vertices_end() {
    return vertex_iterator(vertices);
  }

  edge_iterator edges_begin() {
    return edge_iterator(&row_indices[0], &columns[0], 0, 0);
  }

  edge_iterator edges_end() {
    return edge_iterator(&row_indices[0], &columns[0], vertices, edges);
  }

  // Get a start iterator to edges going out from vertex v
  std::pair<edge_iterator, edge_iterator> out_going_edges(vertex_t v) {
    edge_iterator starte = edge_iterator(&row_indices[0], 
					 &columns[0], v, row_indices[v]);
    edge_iterator ende = edge_iterator(&row_indices[0], 
				       &columns[0], v, row_indices[v+1]);
    return std::make_pair(starte, ende);
  }

  edge_property_t get_edge_weight(EdgeType_t e) {
    assert(e.eid != -1 && e.eid < edges);
    return edge_weight_map[e.eid];
  }

  void set_vertex_property(vertex_property_t val, vertex_t v) {
    vertex_distance_map[v] = val;
  }

  vertex_property_t get_vertex_property(vertex_t v) {
    return vertex_distance_map[v];
  }

  void partition_graph(partition_client_map_t& partitions, int num_qs,
		       int& num_vert_per_local) {

    std::vector<hpx::naming::id_type> localities =
      hpx::find_all_localities();
    std::size_t num_locs = localities.size();
    std::cout << "Number of localities : " << num_locs << std::endl;

     // equally distribute vertices among localities
    num_vert_per_local = vertices / num_locs;
    int vert_balance = vertices % num_locs;

    //std::cout << "partitions X : " << pdatas.size() << std::endl;

    for(std::size_t i=0; i<num_locs; ++i) {
      int startv = i*num_vert_per_local;

      // if this is last locality add balance vertices to last
      int endv;
      if (i == num_locs-1) {
	endv = num_vert_per_local+i*num_vert_per_local + 1 + vert_balance;
      } else {
	endv = num_vert_per_local+i*num_vert_per_local + 1;
      }
     
      index_t starte = row_indices[startv];
      index_t ende = row_indices[endv-1];

      std::cout << "startv : " << startv << " endv : " << endv
		<< " starte : " << starte << " ende : " << ende 
		<< std::endl;
      graph_partition_data pd(startv,
			      endv,
			      num_vert_per_local,
			      num_qs,
			      undirected);

      pd.vertex_distances.resize(endv-startv);
      pd.vertex_distances.assign((endv-startv), 
				 std::numeric_limits<vertex_t>::max());
           
      // assign row indices
      for (int k=startv; k < endv; ++k) {
	pd.row_indices.push_back(row_indices[k]);
      }

      // assign columns and weights
      for (int k=starte; k < ende; ++k) {
	pd.columns.push_back(columns[k]);
	pd.weights.push_back(edge_weight_map[k]);
      }

      // Distribute data
      // To distribute we invoke component client, i.e. partition
      // and give locality and graph_partition. This operation will 
      // distribute graph_partition to remote locality
      std::cout << "Pushing to locality : " << 
	hpx::naming::get_locality_id_from_id(localities[i]) << std::endl;
      //pd.print();
      partition p(localities[i], pd);
      partitions.insert(std::make_pair(hpx::naming::get_locality_id_from_id(localities[i]), p));
    }
  }

  void validate_partitions(const std::vector<graph_partition_data>&
			   partitions) {

    int loc = 0;

    // checking vertices are equal
    for(std::size_t i=0; i<vertices+1; ++i) {
      if (i == partitions[loc].vertex_end) {
	++loc;
      }
      std::cout << "index : " << i
		<< " main indices : " << row_indices[i]
		<< " partition indices : " 
		<< partitions[loc].row_indices[(i - 
						partitions[loc].vertex_start)]
	        << " locality : " << loc <<  std::endl;
      int row_val = (i - partitions[loc].vertex_start);
      assert(row_indices[i] == partitions[loc].row_indices[row_val]);
    }

    loc = 0;
    std::size_t col_val = 0;
    // checking edges are equals
    for(std::size_t k=0; k < edges; ++k, ++col_val) {
      if (col_val == partitions[loc].columns.size()) {
	++loc;
	col_val = 0;
      }

      assert(columns[k] == partitions[loc].columns[col_val]);
      assert(edge_weight_map[k] == partitions[loc].weights[col_val]);
    }
  }

  void print() {
    std::cout << "Vertices - " << vertices << ", Edges - " << edges << std::endl;
    std::cout << "Printing row index array ...." << std::endl;
    std::cout << "[";
    // printing raw index array
    for (std::size_t i=0; i<(vertices+1); ++i) {
      std::cout << row_indices[i];
      if (i != vertices) {
        std::cout << ", ";
      }
    }

    std::cout << "]\n" << std::endl;
    
    std::cout << "Printing row array ...." << std::endl;
    std::cout << "[";
    // printing raw index array
    for (std::size_t i=0; i<edges; ++i) {
      std::cout << columns[i];
      if (i != (edges-1)) {
        std::cout << ", ";
      }
    }
    std::cout << "]\n" << std::endl;
  }
  
private:
  std::size_t vertices;
  std::size_t edges;
  std::vector<index_t> row_indices;
  std::vector<index_t> columns;
  bool undirected;
  std::vector<edge_property_t> edge_weight_map;
  std::vector<vertex_property_t> vertex_distance_map;
  std::size_t last_updated_ri = 0; // last updated row index; only used during graph construction time

  void init() {
    for(std::size_t j=0; j < (vertices+1); ++j) {
      row_indices[j] = INVALID_VERTEX;
    }
    
    for(std::size_t k = 0; k < edges; ++k) {
      columns[k] = INVALID_EDGE;
    }
  }
};

#endif


