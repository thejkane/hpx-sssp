//=============================================================================
// Copyright (C) 2014 The Trustees of Indiana University.
//
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
//  Authors: Thejaka Kanewala
//           Marcin Zalewski
//           Andrew Lumsdaine
//
// Description : The Distributed Control based SSSP implementation on HPX 
//               platform.
//=============================================================================

#ifndef HPX_DC_SSSP
#define HPX_DC_SSSP

//#define WORK_STATS 1

// This is needed to increase the number of 
// parameters to new operator
#define HPX_LIMIT 6

#include <atomic>
#include <queue>
#include <algorithm>

#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/lcos/reduce.hpp>

#include <boost/format.hpp>
#include <boost/cstdint.hpp>
#include <boost/shared_array.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

// For graph generation
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/graph/graph500_generator.hpp>


#include "boost/graph/parallel/thread_support.hpp"
#include "common_types.hpp"

boost::uint32_t total_q_count = 0;

struct partition;
typedef std::map<boost::uint32_t, partition> partition_client_map_t;

partition_client_map_t global_partitions;

// stores all component ids
typedef std::vector<hpx::naming::id_type> component_ids_t;

typedef std::vector<int> graph_array_t;

//==========================================//
// Work stats
#ifdef WORK_STATS
std::atomic_int_fast64_t useful(0);
std::atomic_int_fast64_t invalidated(0);
std::atomic_int_fast64_t useless(0);
std::atomic_int_fast64_t rejected(0);
std::atomic_int_fast64_t partial_buffers(0);
std::atomic_int_fast64_t full_buffers(0);
#endif
//==========================================//


//===========================================
// This use to build the histogram.
// contains target vertex and edge weight
//===========================================
struct target_weight {
  vertex_t target;
  edge_property_t weight;
      
  target_weight(vertex_t t, edge_property_t w) : target(t), weight(w){}
};

  // Edges are stored as a multiset in the histogram
  // sort comparer to maintain edges in order
struct tw_comparator {
  bool operator() (const target_weight& tw1, 
		   const target_weight& tw2) const {
    if (tw1.target < tw2.target)
      return true;
    else if (tw1.target == tw2.target) {
      return tw1.weight < tw2.weight;
    } else {
      return false;
    }
  }
};

// target vertex and edge weight
typedef std::multiset<target_weight, tw_comparator> EdgeList_t;
// source vertex and list of targets
// for (1,3)-w1, (1,5)-w2 we have
// 1 - 3-w1, 5-w2 etc ...
typedef std::map<vertex_t, EdgeList_t> HistogramMap_t;



//====================================================================//
// Vertex distance type
struct vertex_distance {
  vertex_t vertex;
  vertex_property_t distance;

  vertex_distance(const vertex_distance& other):
    vertex(other.vertex), distance(other.distance)
  {}

  vertex_distance(vertex_t v, vertex_property_t d) : vertex(v),
						     distance(d)
  {}

  vertex_distance():vertex(0),distance(0)
  {}

private:
  // Serialization support: even if all of the code below runs on one
  // locality only, we need to provide an (empty) implementation for the
  // serialization as all arguments passed to actions have to support this.
  friend class boost::serialization::access;

  template <typename Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar & vertex & distance;
  }
};

// Comparer for priority queue
struct default_comparer {
  bool operator()(const vertex_distance& vd1, const vertex_distance& vd2) {
    return vd1.distance > vd2.distance;
  }
};

// Comparer for sorting coalesced messages
struct sort_comparer {
  bool operator()(const vertex_distance& vd1, const vertex_distance& vd2) {
    return vd1.distance < vd2.distance;
  }
} sc;


struct dc_priority_queue;
// Priority queue type
typedef std::priority_queue<vertex_distance, 
			    std::vector<vertex_distance>, default_comparer > priority_q_t;
// All priority queues
typedef std::vector<dc_priority_queue> all_q_t;

//==========================================//
// Stores all priority queues
// This is ugly. But cant help; cos - HPX
// reduction support at component level is not
// working as expected.
//==========================================//
all_q_t buckets;


// Coalesced message type
typedef std::vector<vertex_distance> coalesced_message_t;
//====================================================================//


//=========================================
// Used to transfer partition information
// across different localities. This class represent
// a portion of graph which resides in a single
// locality.
//=========================================
struct graph_partition_data {
  int vertex_start;
  int vertex_end;

  // This is needed to calculate
  // the index of partition particulat vertex belongs
  int number_vertices_per_locality; 
  
  // Number of queues to create
  int num_queues;

  // graph is undirected or not
  bool undirected;

  graph_array_t row_indices;
  graph_array_t columns;
  graph_array_t weights;
  graph_array_t vertex_distances;

private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & vertex_start;
    ar & vertex_end;
    ar & number_vertices_per_locality;
    ar & num_queues;
    ar & undirected;
    ar & row_indices;
    ar & columns;
    ar & weights;
    ar & vertex_distances;
  }

public:

  //=====================================================
  // An iterator class to iterate 
  // graph vertices.
  //=====================================================
  class vertex_iterator : 
    public std::iterator<std::input_iterator_tag, vertex_t> {
  private:
    vertex_t vertex_array;

  public:
    vertex_iterator(vertex_t vertices) :
      vertex_array(vertices) {}

    vertex_iterator(const vertex_iterator& vit) : 
      vertex_array(vit.vertex_array) {}

    vertex_iterator& operator++() {
      ++vertex_array;return *this;
    }

    vertex_iterator operator++(int) {
      vertex_iterator tmp(*this); operator++(); return tmp;
    }

    bool operator==(const vertex_iterator& rhs) {
      return vertex_array == rhs.vertex_array;
    }

    bool operator!=(const vertex_iterator& rhs) {
      return vertex_array != rhs.vertex_array;
    }

    vertex_t& operator*() {return vertex_array;}
  };

  //=====================================================
  // An iterator class to iterate 
  // graph edges.
  //=====================================================
  // To iterate through all edges
  class edge_iterator : 
    public std::iterator<std::input_iterator_tag, EdgeType_t> {
  private:
    edge_t* row_index;
    edge_t* rows;
    edge_t ri; // the progressing row index
    edge_t ci; // the progressing colomn index
    vertex_t offset;

  public:
    edge_iterator(edge_t* row_i, edge_t* r, int rowind, int colind, int off) : 
      row_index(row_i), 
      rows(r), 
      ri(rowind), 
      ci(colind),
      offset(off)
    {}

    edge_iterator(const edge_iterator& eit) : 
      row_index(eit.row_index), 
      rows(eit.rows), 
      ri(eit.ri), 
      ci(eit.ci),
      offset(eit.offset)
    {}

    void print() {
      std::cout << "Row index : " << ri 
		<< " Column index : " 
		<< ci 
		<< " Offset : "
		<< offset
		<< std::endl;
    }

    edge_iterator& operator++() {
      edge_t ci_index_end = row_index[ri+1] - row_index[0];
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
	      (offset == rhs.offset) &&
	      (rows == rhs.rows) && (ci == rhs.ci));
    }

    bool operator!=(const edge_iterator& rhs) {
      return ((row_index != rhs.row_index) ||
	      (offset != rhs.offset) ||
	      (rows != rhs.rows) || (ci != rhs.ci)); 
    }

    EdgeType_t operator*() {
      while(row_index[ri] == -1) {
	assert(false);
      }
      edge_t val = *(rows+ci);
      return EdgeType_t(ri+offset, val, ci+row_index[0]); 
    }
  };

  typedef std::pair<graph_partition_data::edge_iterator, 
		    graph_partition_data::edge_iterator> OutgoingEdgePair_t;

  graph_partition_data() {}

  graph_partition_data(int _vstart, 
		       int _vend, 
		       int vert_loc,
		       int num_q,
		       bool und) : 
    vertex_start(_vstart), 
    vertex_end(_vend),
    number_vertices_per_locality(vert_loc),
    num_queues(num_q),
    undirected(und) {

    row_indices.resize((vertex_end - vertex_start) + 1);
    vertex_distances.resize(vertex_end - vertex_start);
  }

  graph_partition_data(int _vstart, 
		       int _vend,
		       int vert_loc,
		       int num_q,
		       bool und,
		       const graph_array_t& _ri,
		       const graph_array_t& _cl,
		       const graph_array_t& _wt,
		       const graph_array_t& _vd) : 
    vertex_start(_vstart), 
    vertex_end(_vend),
    number_vertices_per_locality(vert_loc),
    num_queues(num_q),
    undirected(und),
    row_indices(_ri), 
    columns(_cl),
    weights(_wt), 
    vertex_distances(_vd)
  {}

  // copy constructor
  graph_partition_data(const graph_partition_data& other):
    vertex_start(other.vertex_start),
    vertex_end(other.vertex_end),
    number_vertices_per_locality(other.number_vertices_per_locality),
    num_queues(other.num_queues),
    undirected(other.undirected),
    row_indices(other.row_indices),
    columns(other.columns),
    weights(other.weights),
    vertex_distances(other.vertex_distances)
  {}

  vertex_iterator vertices_begin() {
    return vertex_iterator(vertex_start);
  }

  vertex_iterator vertices_end() {
    return vertex_iterator(vertex_end-1);
  }

  edge_iterator edges_begin() {
    return edge_iterator(&row_indices[0], &columns[0], 0, 0, vertex_start);
  }

  edge_iterator edges_end() {
    int local_vert_index = vertex_end-vertex_start;
    return edge_iterator(&row_indices[0], &columns[0], 
			 local_vert_index, 
			 (row_indices[local_vert_index-1] - row_indices[0]),
			 vertex_start);
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
      //std::cout << "Before histogram : (" 
      //		<< (*begin).first << ", " << (*begin).second
      //		<< ")" << std::endl;
      // check whether source vertex belong to current
      // partition. If not continue.
      if (vertex_start <= (*begin).first &&
	  (*begin).first < vertex_end) {
	buildHistogram(histogram_map, (*begin).first, 
		       (*begin).second, *epiter, false); 
      }
      // flipped=false => edge need to flip for undirected
    }

    //    std::cout << "coming here : " << vertex_start << " - "
    //	      << vertex_end << std::endl;

    edge_t row_ind = 1; // starts with 1
    edge_t col_ind = 0;
    row_indices[0] = 0;

    for(vertex_t k=vertex_start; k < vertex_end; ++k) {
      HistogramMap_t::iterator iteFind = histogram_map.find(k);
      if (iteFind != histogram_map.end()) {
	EdgeList_t list = (*iteFind).second;
	row_indices[row_ind] = row_indices[row_ind-1] + list.size();
	++row_ind;
	EdgeList_t::iterator iteList = list.begin();
	for (; iteList != list.end(); ++iteList) {
	  columns.push_back((*iteList).target);
	  weights.push_back((*iteList).weight);
	  ++col_ind;
	}
      } else {
	// no edge
	row_indices[row_ind] = row_indices[row_ind-1];
	++row_ind;
      }
    }

    HPX_ASSERT(row_ind == row_indices.size());
    HPX_ASSERT(col_ind == columns.size());
    HPX_ASSERT(col_ind == weights.size());
  }


  void generate_local_graph_partition(boost::uint32_t scale,
				      edge_t n, /* number of edges */
				      boost::uint32_t max_weight,
				      uint64_t a,
				      uint64_t b) {
    // The modified graph 500 iterator
    typedef boost::graph500_iterator<vertex_t, edge_t> Graph500Iter;

    // Random weight generation
    typedef boost::uniform_int<edge_property_t> distribution_type;

    boost::minstd_rand edge_weight_gen;

    typedef boost::variate_generator<boost::minstd_rand&, distribution_type> gen_type;

    gen_type die_gen(edge_weight_gen, distribution_type(1, max_weight));
    boost::generator_iterator<gen_type> die(&die_gen);

    // Randome edge generation
    boost::uniform_int<uint64_t> 
      rand_64(0, std::numeric_limits<uint64_t>::max());

    addEdges(Graph500Iter(scale, 0, a, b), 
	     Graph500Iter(scale, n, a, b), die);

    // Finally initialize vertex distance map
    vertex_distances.resize(vertex_end-vertex_start);
    vertex_distances.assign((vertex_end-vertex_start), 
			       std::numeric_limits<vertex_t>::max());

  }

  // Get a start iterator to edges going out from vertex v
  OutgoingEdgePair_t out_going_edges(vertex_t v) {
#ifdef PRINT_DEBUG
    std::cout << "v-" << v << " start-" << vertex_start << " end-"
    	      << vertex_end << std::endl;
#endif

    assert(vertex_start <= v && v < vertex_end);
    vertex_t local_v = v - vertex_start;

    edge_iterator starte = edge_iterator(&row_indices[0], 
					 &columns[0], local_v, (row_indices[local_v] - row_indices[0]),
					 vertex_start);
    edge_iterator ende = edge_iterator(&row_indices[0], 
				       &columns[0], local_v, (row_indices[local_v+1] - row_indices[0]),
				       vertex_start);
    return std::make_pair(starte, ende);
  }

  edge_property_t get_edge_weight(EdgeType_t e) {
    //    std::cout << "e.eid : " << e.eid << " row_indices[0]  : " << row_indices[0]
    //	      << " row_indices[(vertex_end-vertex_start)-1] : "
    //	      << row_indices[(vertex_end-vertex_start)]
    //	      << std::endl;
    assert(e.eid != -1 && 
	   (row_indices[0] <= e.eid) &&  
	   (e.eid < row_indices[(vertex_end-vertex_start)]));

    return weights[(e.eid - row_indices[0])];
  }

  vertex_property_t get_vertex_distance(vertex_t vid) {
#ifdef PRINT_DEBUG
    std::cout << "vertex_start : " << vertex_start
	      << " vid : " << vid
	      << " vertex_end : " << vertex_end << std::endl;
#endif

    HPX_ASSERT(vertex_start <= vid && vid < vertex_end);
    return vertex_distances[vid-vertex_start];
  }

  inline bool set_vertex_distance_atomic(vertex_t vid, 
					 vertex_property_t new_distance) {
#ifdef PRINT_DEBUG
    std::cout << "vid - " << vid << " start - " << vertex_start
	      << " end - " << vertex_end << std::endl;
#endif
    HPX_ASSERT(vertex_start <= vid && vid < vertex_end);

    vertex_t localvid = vid - vertex_start;

    int old_dist = vertex_distances[localvid], last_old_dist;
    while (new_distance < old_dist) {
      last_old_dist = old_dist;
      old_dist = boost::parallel::val_compare_and_swap
	(&vertex_distances[localvid], old_dist, new_distance);
      if (last_old_dist == old_dist) {
#ifdef WORK_STATS
	if(old_dist < std::numeric_limits<vertex_t>::max()) invalidated++;
#endif
	return true;
      }
    }

    return false;
  }

  boost::uint32_t find_locality_id(vertex_t v) {
    HPX_ASSERT(number_vertices_per_locality != 0);

    std::vector<hpx::naming::id_type> localities =
      hpx::find_all_localities();
    std::size_t num_locs = localities.size();

    vertex_t max_partition_vertex 
      = (num_locs * number_vertices_per_locality) - 1;

    if (v > max_partition_vertex)
      return (num_locs-1);
    else 
      return (v / number_vertices_per_locality);
  }

  hpx::naming::id_type find_component_id(vertex_t v) {

    std::vector<hpx::naming::id_type> localities =
      hpx::find_all_localities();
    std::size_t num_locs = localities.size();

    boost::uint32_t loc_id = find_locality_id(v);
    
    for (std::size_t k=0; k<num_locs; ++k) {
      if (hpx::naming::get_locality_id_from_id(localities[k])
	  == loc_id)
	return localities[k];
    }
    
    HPX_ASSERT(false); // should not come here
  }


  void print() {
    std::cout << "Vertex start : " 
	      << vertex_start << " vertex end : " 
	      << vertex_end << std::endl;

    std::cout << "Row indices : {";
    // copy raw indices
    for(int i=0; i < (vertex_end-vertex_start); ++i) {
      std::cout << row_indices[i] << ", ";
    }
    std::cout << "}" << std::endl;

    int num_edges = row_indices[vertex_end-vertex_start-1] - row_indices[0];
    std::cout << "Num Edges : " << num_edges << std::endl;

    std::cout << "Columns : {";
    // copy columns & weights
    for(int i=0; i<num_edges; ++i) {
      std::cout << columns[i] << ", ";
    }
    std::cout << "}" << std::endl;

    std::cout << "Weights : {";
    // copy columns & weights
    for(int i=0; i<num_edges; ++i) {
      std::cout << weights[i] << ", ";
    }
    std::cout << "}" << std::endl;
  }
};

typedef std::vector< hpx::future <void> > future_collection_t;
typedef hpx::lcos::local::spinlock mutex_type;
typedef std::map<boost::uint32_t, coalesced_message_t> coalsced_message_map_t;
//==========================================//
// completed_count - stores the number of messages
// sent through flush_task
// receive_count - stores the number of messages
// received through relax
// these are useful for termination.
//==========================================//  
std::atomic_int_fast64_t active_count(0);
std::atomic_int_fast64_t completed_count(0); // Source does not send a message
boost::uint32_t empty_q_count = 0; 

hpx::lcos::local::condition_variable q_count_cv;
boost::mutex q_count_mutex;

void increase_empty_q_count() {
  boost::mutex::scoped_lock scopedLock(q_count_mutex); 
  empty_q_count++;

#ifdef PRINT_DEBUG
  std::cout << "[inc] Rank " << hpx::naming::get_locality_id_from_id(hpx::find_here())
	    << " The total q count " << total_q_count 
	    << " The empty q count : " << empty_q_count << std::endl;
#endif

  if (empty_q_count == total_q_count) {
    // notify all threads waiting on this condition variable
    q_count_cv.notify_all();
  }
}

void wait_till_all_qs_empty() {
  boost::mutex::scoped_lock scopedLock(q_count_mutex);

#ifdef PRINT_DEBUG
  std::cout << "[wait] Rank " << hpx::naming::get_locality_id_from_id(hpx::find_here())
	    << " The total q count " << total_q_count 
	    << " The empty q count : " << empty_q_count << std::endl;
#endif

  if (empty_q_count != total_q_count) {
    q_count_cv.wait(scopedLock);
  }
}

void decrease_empty_q_count() {
  boost::mutex::scoped_lock scopedLock(q_count_mutex);
  empty_q_count--;

#ifdef PRINT_DEBUG
  std::cout << "[dec] Rank " << hpx::naming::get_locality_id_from_id(hpx::find_here())
	    << " The total q count " << total_q_count 
	    << " The empty q count : " << empty_q_count << std::endl;
#endif
}

///////////////////////////////////////////////////////////////////////////////
// Represents a single priority queue
///////////////////////////////////////////////////////////////////////////////
struct dc_priority_queue {
  
  dc_priority_queue():
    termination(false)
  {
  }

  dc_priority_queue(const dc_priority_queue& other):
    termination(other.termination),
    pq(other.pq)
  {}

  void push(const vertex_distance& vd) {
    // lock the queue and insert element
    {
      boost::mutex::scoped_lock scopedLock(mutex);      
      pq.push(vd);

#ifdef WORK_STATS
      useful++;
#endif
    }

    // notify waiting threads
    cv.notify_all();
  }

  void handle_queue(graph_partition_data& graph_partition,
		    const boost::uint32_t yield_count);

  void static init() {
    // for each locality initialize a coalesced buffer
    std::vector<hpx::naming::id_type> localities =
      hpx::find_all_localities();

    std::vector<hpx::naming::id_type>::iterator iteLoc = localities.begin();
    for (; iteLoc != localities.end(); ++iteLoc) {
      boost::uint32_t locId = hpx::naming::get_locality_id_from_id(*iteLoc);
      cmap.insert(std::make_pair(locId, coalesced_message_t()));
      std::cout << "pushing to mutexes ... " << std::endl;
      cmap_mutexes.push_back(new boost::mutex);
    }
    
    //    cmap_mutexes.resize(localities.size());
  }

  void static send_all();

  void static send(const vertex_distance vd,
	    boost::uint32_t target_locality,
	    const partition& partition_client);

  // Terminates the algorithms
  void terminate() {
    termination = true;
    cv.notify_all();
  }

  void reset() {
    termination = false;
    
    {
      //boost::mutex::scoped_lock scopedLock(cmap_mutex);      
      // No residues from previous runs
      // Just make sure
      coalsced_message_map_t::iterator ite = cmap.begin();
      for (; ite != cmap.end(); ++ite) {
	HPX_ASSERT((*ite).second.empty());
      }
    }

  }

private:
  bool termination;
  priority_q_t pq;
  hpx::lcos::local::condition_variable cv;
  boost::mutex mutex;
  static coalsced_message_map_t cmap;
  static boost::ptr_vector<boost::mutex> cmap_mutexes;
};


//===================================
// Sends all remaining messages in 
// buffers.
//===================================
void send_all_remaining() {
  // all qs are empty
  // send all remaining messages: but need to lock
  all_q_t::iterator ite = buckets.begin();
  for (; ite != buckets.end(); ++ite) {
    (*ite).send_all();
  }
}


// This is the server side representation of the data. We expose this as a HPX
// component which allows for it to be created and accessed remotely through
// a global address (hpx::id_type).
struct partition_server
  : hpx::components::simple_component_base<partition_server> {

  // construct new instances
  partition_server() {}

  partition_server(graph_partition_data const& data)
    : graph_partition(data) {
    
    init();
  }

  partition_server(partition_server const& ps)
    : graph_partition(ps.graph_partition) {

    init();
  }


  // Access data. The parameter specifies what part of the data should be
  // accessed. As long as the result is used locally, no data is copied,
  // however as soon as the result is requested from another locality only
  // the minimally required amount of data will go over the wire.
  graph_partition_data get_data() const
  {
    return graph_partition;
  }

  // Every member function which has to be invoked remotely needs to be
  // wrapped into a component action. The macro below defines a new type
  // 'get_data_action' which represents the (possibly remote) member function
  // partition::get_data().
  HPX_DEFINE_COMPONENT_CONST_DIRECT_ACTION(partition_server, get_data, get_data_action);

  //==============================================================
  // Creates remote partition clients.
  //==============================================================
  void create_partition_clients();


  HPX_DEFINE_COMPONENT_ACTION(partition_server, create_partition_clients,
			      dc_create_partition_clients_action);


  //==============================================================
  // Experimenting... First with chaotic algorithm.
  // In this we will relax each vertex parallely
  //==============================================================
  void relax(const vertex_distance& vd);

  HPX_DEFINE_COMPONENT_ACTION(partition_server, relax,
			      dc_relax_action);

  //==============================================================
  // Count locally visited vertices
  //==============================================================
  boost::uint32_t count_visited_vertices() {
    boost::uint32_t visited = 0;
    int num_local_verts = (graph_partition.vertex_end - graph_partition.vertex_start) - 1;
    for(int i=0; i < num_local_verts; ++i) {
      if (graph_partition.vertex_distances[i] < std::numeric_limits<vertex_t>::max()) {
	++visited;
      }
    }

    return visited;
  }

  HPX_DEFINE_COMPONENT_ACTION(partition_server, count_visited_vertices,
			      dc_count_visited_vertices_action);

  //==============================================================
  // Validate calculated distances are correct.
  //==============================================================
  void verify_partition_results();

  HPX_DEFINE_COMPONENT_ACTION(partition_server, verify_partition_results,
			      dc_verify_partition_results_action);

  //==============================================================
  // Send coalesced messages
  // In this we will relax each vertex parallely
  //==============================================================
  void coalesced_relax(const coalesced_message_t& vds);

  HPX_DEFINE_COMPONENT_ACTION(partition_server, coalesced_relax,
			      dc_coalesced_relax_action);


  //==============================================================
  // Gets the stored distance for a given vertex
  //==============================================================
  vertex_property_t get_vertex_distance(const vertex_t v) {
    return graph_partition.get_vertex_distance(v);
  }

  HPX_DEFINE_COMPONENT_ACTION(partition_server, get_vertex_distance,
			      dc_get_vd_action);


  //==============================================================
  // wait till all futures complete their work
  // idx is the queue index to work on
  //==============================================================
  void flush_tasks(int idx,
		   const boost::uint32_t yield_count);

  HPX_DEFINE_COMPONENT_ACTION(partition_server, flush_tasks,
			      dc_flush_action);

  //==============================================================
  // Graph generation is distributed. Generate the local part.
  //==============================================================
  void generate_local_graph(boost::uint32_t scale,
			    edge_t n, /* number of edges */
			    boost::uint32_t max_weight,
			    uint64_t a, // seeds
			    uint64_t b) {

    graph_partition.generate_local_graph_partition(scale,
						   n,
						   max_weight,
						   a,
						   b);
  }

  HPX_DEFINE_COMPONENT_ACTION(partition_server, generate_local_graph,
			      dc_local_graph_gen_action);


  //==============================================================
  // Reset counters for a new run
  //==============================================================
  void reset_counters() {
    // reset counters
    completed_count = 0;
    active_count = 0;
    
    // reset termination variable in each q
    for(int i=0; i<graph_partition.num_queues; ++i) {
      buckets[i].reset();
    }
    
#ifdef PRINT_DEBUG
    std::cout << "@@@@@@@@@@@@@@@@ Inside Reset - Before resetting @@@@@@@@@@@@@@" << std::endl;
    std::cout << count_visited_vertices() << std::endl;
#endif

    // reset vertex distance
    graph_partition.
      vertex_distances.
      assign(graph_partition.vertex_distances.size(),
	     std::numeric_limits<vertex_t>::max());

#ifdef PRINT_DEBUG
    std::cout << count_visited_vertices() << std::endl;
    std::cout << "@@@@@@@@@@@@@@@@ Inside Reset - After resetting @@@@@@@@@@@@@@" << std::endl;
#endif

    // reset work stats
#ifdef WORK_STATS
    useful = 0;
    useless = 0;
    rejected = 0;
    invalidated = 0;
    partial_buffers = 0;
    full_buffers = 0;
#endif
  }

  HPX_DEFINE_COMPONENT_ACTION(partition_server, reset_counters,
			      dc_reset_counters_action);

  //==============================================================
  // Handles termination.
  //==============================================================
  void terminate() {
    for(int i=0; i<graph_partition.num_queues; ++i) {
      buckets[i].terminate();
    }
  }

  HPX_DEFINE_COMPONENT_ACTION(partition_server, terminate,
			      dc_terminate_action);


  // reduction action for total_msg_difference
  //  HPX_DEFINE_COMPONENT_ACTION(partition_server, total_msg_difference,
  //			      tot_msg_diff_action);

private:
  graph_partition_data graph_partition;
  
  //==========================================//
  // To select a pq randomly
  //==========================================//
  boost::random::mt19937 gen;
  
  // Initializes thread queues
  void init() {
    
    active_count = 0;

    // create queues
    buckets.resize(graph_partition.num_queues);

    // initialize total q count
    total_q_count = graph_partition.num_queues;
    
    dc_priority_queue::init();
  }

  // Finds a random index value to find a queue
  int select_random_q_idx() {
    boost::random::uniform_int_distribution<> dist(0, (graph_partition.num_queues-1));
    return dist(gen);
  }

};

//============== Non Reduction Action Definitions===================//
HPX_REGISTER_ACTION_DECLARATION(partition_server::dc_relax_action, partition_relax_action);
HPX_REGISTER_ACTION_DECLARATION(partition_server::dc_create_partition_clients_action, partition_create_partition_clients_action);
HPX_REGISTER_ACTION_DECLARATION(partition_server::dc_verify_partition_results_action, partition_verify_partition_results_action);
HPX_REGISTER_ACTION_DECLARATION(partition_server::dc_count_visited_vertices_action, partition_count_visited_vertices_action);
HPX_REGISTER_ACTION_DECLARATION(partition_server::dc_coalesced_relax_action, partition_coalesced_relax_action);
HPX_REGISTER_ACTION_DECLARATION(partition_server::dc_get_vd_action, partition_get_vd_action);
HPX_REGISTER_ACTION_DECLARATION(partition_server::dc_flush_action, partition_flush_action);
HPX_REGISTER_ACTION_DECLARATION(partition_server::dc_terminate_action, partition_terminate_action);
HPX_REGISTER_ACTION_DECLARATION(partition_server::dc_reset_counters_action, partition_counter_reset_action);
HPX_REGISTER_ACTION_DECLARATION(partition_server::dc_local_graph_gen_action, partition_local_graph_gen_action);



// The macros below are necessary to generate the code required for exposing
// our partition type remotely.
//
// HPX_REGISTER_MINIMAL_COMPONENT_FACTORY() exposes the component creation
// through hpx::new_<>().
typedef hpx::components::simple_component<partition_server> partition_server_type;
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(partition_server_type, partition_server);

// HPX_REGISTER_ACTION() exposes the component member function for remote
// invocation.
typedef partition_server::get_data_action get_data_action;
HPX_REGISTER_ACTION(get_data_action);

typedef partition_server::dc_create_partition_clients_action partition_create_partition_clients_action;
HPX_REGISTER_ACTION(partition_create_partition_clients_action);

typedef partition_server::dc_relax_action partition_relax_action;
HPX_REGISTER_ACTION(partition_relax_action);

typedef partition_server::dc_verify_partition_results_action partition_verify_partition_results_action;
HPX_REGISTER_ACTION(partition_verify_partition_results_action);

typedef partition_server::dc_count_visited_vertices_action partition_count_visited_vertices_action;
HPX_REGISTER_ACTION(partition_count_visited_vertices_action);

typedef partition_server::dc_coalesced_relax_action partition_coalesced_relax_action;
HPX_REGISTER_ACTION(partition_coalesced_relax_action);

typedef partition_server::dc_get_vd_action partition_get_vd_action;
HPX_REGISTER_ACTION(partition_get_vd_action);

typedef partition_server::dc_flush_action partition_flush_action;
HPX_REGISTER_ACTION(partition_flush_action);

typedef partition_server::dc_terminate_action partition_terminate_action;
HPX_REGISTER_ACTION(partition_terminate_action);

typedef partition_server::dc_reset_counters_action partition_counter_reset_action;
HPX_REGISTER_ACTION(partition_counter_reset_action);

typedef partition_server::dc_local_graph_gen_action partition_local_graph_gen_action;
HPX_REGISTER_ACTION(partition_local_graph_gen_action);


// TODO encapsulate parameters

char const* base_name = "/partition";
///////////////////////////////////////////////////////////////////////////////
// This is a client side helper class allowing to hide some of the tedious
// boilerplate.
struct partition : hpx::components::client_base<partition, partition_server> {
  typedef hpx::components::client_base<partition, partition_server> base_type;

  partition() {}

  // Create a new component on the locality co-located to the id 'where'. The
  // new instance will be initialized from the given partition_data.
  partition(hpx::id_type where, graph_partition_data const& data)
    : base_type(hpx::new_colocated<partition_server>(where, data))
  {
    hpx::register_id_with_basename(base_name, 
				   get_gid(),
				   hpx::naming::get_locality_id_from_id(where));
  }

  // Attach a future representing a (possibly remote) partition.
  partition(hpx::future<hpx::id_type> && id)
    : base_type(std::move(id))
  {}

  // Attach a future representing a (possibly remote) partition.
  partition(hpx::id_type where)
    : base_type(hpx::find_id_from_basename(base_name, 
					   hpx::naming::get_locality_id_from_id(where)))
  {}


  // Unwrap a future<partition> (a partition already holds a future to the
  // id of the referenced object, thus unwrapping accesses this inner future).
  partition(hpx::future<partition> && c)
    : base_type(std::move(c))
  {}

  ///////////////////////////////////////////////////////////////////////////
  // Invoke the (remote) member function which gives us access to the data.
  // This is a pure helper function hiding the async.
  hpx::future<graph_partition_data> get_data() const {
    partition_server::get_data_action act;
    return hpx::async(act, get_gid());
  }

  // NOT USED YET
  /* void send(vertex_distance const& vd) {
    messages.push_back(vd);

    if (messages.size() == coalesced_message_size) {
      // do sorting if enabled
      if(sort_coalesced_buffer) {
	std::sort(messages.begin(), messages.end(), sc);
      }
      coalesced_relax(messages);
      messages.clear();
    }
    }*/

  ///////////////////////////////////////////////////////////////////////////
  // Invoke remote relax
  ///////////////////////////////////////////////////////////////////////////
  void relax(vertex_distance const& vd) const {
    partition_server::dc_relax_action act;
    hpx::apply(act, get_gid(), vd);
  }

  ///////////////////////////////////////////////////////////////////////////
  // Invoke remote coalesced relax
  ///////////////////////////////////////////////////////////////////////////
  void coalesced_relax(const coalesced_message_t vds) const {
    partition_server::dc_coalesced_relax_action act;
    hpx::apply(act, get_gid(), vds);
  }


  ///////////////////////////////////////////////////////////////////////////
  // Gets the currently stored vertex distance
  ///////////////////////////////////////////////////////////////////////////
  vertex_property_t get_vertex_distance(const vertex_t v) const {
    partition_server::dc_get_vd_action act;
    return act(get_gid(), v);
  }


  ///////////////////////////////////////////////////////////////////////////
  // Invoke remote or local flush
  ///////////////////////////////////////////////////////////////////////////
  void start_flush_tasks(int num_qs,
			 const boost::uint32_t yield_count) const {
    partition_server::dc_flush_action act;
    for (int i=0; i<num_qs; ++i) {
      hpx::apply(act, get_gid(), i, yield_count);
    }
  }


  ///////////////////////////////////////////////////////////////////////////
  // Verify final distances calculated
  ///////////////////////////////////////////////////////////////////////////
  hpx::future<void> verify_distance_results() {
    partition_server::dc_verify_partition_results_action act;
    return hpx::async(act, get_gid());
  }


  ///////////////////////////////////////////////////////////////////////////
  // Calculate locally visited vertices
  ///////////////////////////////////////////////////////////////////////////
  hpx::future<boost::uint32_t> calculate_visited_vertices() {
    partition_server::dc_count_visited_vertices_action act;
    return hpx::async(act, get_gid());
  }

  ///////////////////////////////////////////////////////////////////////////
  // Invokes termination
  ///////////////////////////////////////////////////////////////////////////
  hpx::future<void> terminate() {
    partition_server::dc_terminate_action act;
    return hpx::async(act, get_gid());
  }

  ///////////////////////////////////////////////////////////////////////////
  // Generates the local graph
  ///////////////////////////////////////////////////////////////////////////
  hpx::future<void> generate_local_graph(boost::uint32_t scale,
				      edge_t n, /* number of edges */
				      boost::uint32_t max_weight,
				      uint64_t a,
				      uint64_t b) {

    partition_server::dc_local_graph_gen_action act;
    return hpx::async(act, get_gid(), scale, n,
		      max_weight, a, b);
  }


  ///////////////////////////////////////////////////////////////////////////
  // Create partition client map
  ///////////////////////////////////////////////////////////////////////////
  void create_partition_clients() {
    partition_server::dc_create_partition_clients_action act;
    act(get_gid());
  }


  ///////////////////////////////////////////////////////////////////////////
  // Reset counters for next run
  ///////////////////////////////////////////////////////////////////////////
  void reset_counters() {
    partition_server::dc_reset_counters_action act;
    act(get_gid());
  }

};

//====================================================================//
namespace hpx { namespace traits {
    template <>
    struct serialize_as_future<partition_client_map_t>
      : boost::mpl::true_
    {
      static void call(partition_client_map_t& m)
      {
	auto r = boost::adaptors::values(m);
	//std::vector<partition> ps(boost::begin(r), boost::end(r)); 
	//auto r = boost::adaptors::values(m);
	hpx::lcos::wait_all(boost::begin(r), boost::end(r));
	//	hpx::lcos::wait_all(ps);
      }
    };
    }}
//====================================================================//

//==============================================================
// wait till all futures complete their work
// idx is the queue index to work on
//==============================================================
// Total completed count reduction
boost::int64_t total_completed_count() {

  wait_till_all_qs_empty();
  send_all_remaining();
  boost::uint64_t cc = completed_count.load(std::memory_order_relaxed);
  return cc;
}


HPX_PLAIN_ACTION(total_completed_count);
typedef std::plus<boost::int64_t> std_plus_type;
HPX_REGISTER_REDUCE_ACTION_DECLARATION(total_completed_count_action, std_plus_type)
HPX_REGISTER_REDUCE_ACTION(total_completed_count_action, std_plus_type)

// Total active count reduction
boost::int64_t total_active_count() {

  wait_till_all_qs_empty();
  send_all_remaining();
  boost::uint64_t ac = active_count.load(std::memory_order_relaxed);
  return ac;
}

HPX_PLAIN_ACTION(total_active_count);
typedef std::plus<boost::int64_t> std_plus_type;
HPX_REGISTER_REDUCE_ACTION_DECLARATION(total_active_count_action, std_plus_type)
HPX_REGISTER_REDUCE_ACTION(total_active_count_action, std_plus_type)

#ifdef WORK_STATS
// Total useful work reduction
boost::int64_t total_useful_work() {

  boost::uint64_t use = useful.load(std::memory_order_relaxed);
  return use;
}

HPX_PLAIN_ACTION(total_useful_work);
typedef std::plus<boost::int64_t> std_plus_type;
HPX_REGISTER_REDUCE_ACTION_DECLARATION(total_useful_work_action, std_plus_type)
HPX_REGISTER_REDUCE_ACTION(total_useful_work_action, std_plus_type)

// Total useless work reduction
boost::int64_t total_useless_work() {

  boost::uint64_t uless = useless.load(std::memory_order_relaxed);
  return uless;
}

HPX_PLAIN_ACTION(total_useless_work);
typedef std::plus<boost::int64_t> std_plus_type;
HPX_REGISTER_REDUCE_ACTION_DECLARATION(total_useless_work_action, std_plus_type)
HPX_REGISTER_REDUCE_ACTION(total_useless_work_action, std_plus_type)

// Total invalidated work reduction
boost::int64_t total_invalidated_work() {

  boost::uint64_t invalid = invalidated.load(std::memory_order_relaxed);
  return invalid;
}

HPX_PLAIN_ACTION(total_invalidated_work);
typedef std::plus<boost::int64_t> std_plus_type;
HPX_REGISTER_REDUCE_ACTION_DECLARATION(total_invalidated_work_action, std_plus_type)
HPX_REGISTER_REDUCE_ACTION(total_invalidated_work_action, std_plus_type)

// Total rejected work reduction
boost::int64_t total_rejected_work() {

  boost::uint64_t reject = rejected.load(std::memory_order_relaxed);
  return reject;
}

HPX_PLAIN_ACTION(total_rejected_work);
typedef std::plus<boost::int64_t> std_plus_type;
HPX_REGISTER_REDUCE_ACTION_DECLARATION(total_rejected_work_action, std_plus_type)
HPX_REGISTER_REDUCE_ACTION(total_rejected_work_action, std_plus_type)

// Total full buffer
boost::int64_t total_full_buffer() {

  boost::uint64_t full = full_buffers.load(std::memory_order_relaxed);
  return full;
}

HPX_PLAIN_ACTION(total_full_buffer);
typedef std::plus<boost::int64_t> std_plus_type;
HPX_REGISTER_REDUCE_ACTION_DECLARATION(total_full_buffer_action, std_plus_type)
HPX_REGISTER_REDUCE_ACTION(total_full_buffer_action, std_plus_type)

// Total partial buffer
boost::int64_t total_partial_buffer() {

  boost::uint64_t partial = partial_buffers.load(std::memory_order_relaxed);
  return partial;
}

HPX_PLAIN_ACTION(total_partial_buffer);
typedef std::plus<boost::int64_t> std_plus_type;
HPX_REGISTER_REDUCE_ACTION_DECLARATION(total_partial_buffer_action, std_plus_type)
HPX_REGISTER_REDUCE_ACTION(total_partial_buffer_action, std_plus_type)


#endif


#endif
