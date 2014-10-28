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

// This is needed to increase the number of 
// parameters to new operator
#define HPX_LIMIT 6

#include <atomic>

#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include <boost/format.hpp>
#include <boost/cstdint.hpp>
#include <boost/shared_array.hpp>
#include <boost/serialization/vector.hpp>

#include "boost/graph/parallel/thread_support.hpp"
#include "common_types.hpp"

struct partition;
typedef std::map<boost::uint32_t, partition> partition_client_map_t;

typedef std::vector<int> graph_array_t;
//=========================================
// Used to transfer partition information
// across different localities. This class represent
// a portion of graph which resides in a single
// locality.
//=========================================
struct graph_partition_data {
  int vertex_start;
  int vertex_end;
  int number_vertices_per_locality; // This is needed to calculate
  // the index of partition particulat vertex belongs
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
		       int vert_loc) : 
    vertex_start(_vstart), 
    vertex_end(_vend),
    number_vertices_per_locality(vert_loc){
  }

  graph_partition_data(int _vstart, 
		       int _vend,
		       int vert_loc,
		       const graph_array_t& _ri,
		       const graph_array_t& _cl,
		       const graph_array_t& _wt,
		       const graph_array_t& _vd) : 
    vertex_start(_vstart), 
    vertex_end(_vend),
    number_vertices_per_locality(vert_loc),
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

  // Get a start iterator to edges going out from vertex v
  OutgoingEdgePair_t out_going_edges(vertex_t v) {
    //    std::cout << "v-" << v << " start-" << vertex_start << " end-"
    //	      << vertex_end << std::endl;
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
    assert(e.eid != -1 && 
	   (row_indices[0] <= e.eid) &&  
	   (e.eid < row_indices[(vertex_end-vertex_start)-1]));

    return weights[(e.eid - row_indices[0])];
  }

  vertex_property_t get_vertex_distance(vertex_t vid) {
    return vertex_distances[vid];
  }

  inline bool set_vertex_distance_atomic(vertex_t vid, 
					 vertex_property_t new_distance) {
    int old_dist = vertex_distances[vid], last_old_dist;
    while (new_distance < old_dist) {
      last_old_dist = old_dist;
      old_dist = boost::parallel::val_compare_and_swap
	(&vertex_distances[vid], old_dist, new_distance);
      if (last_old_dist == old_dist) {
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
    
    for (int k=0; k<num_locs; ++k) {
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

struct vertex_distance {
  vertex_t vertex;
  vertex_property_t distance;

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


typedef std::vector< hpx::future <void> > future_collection_t;
typedef hpx::lcos::local::spinlock mutex_type;
///////////////////////////////////////////////////////////////////////////////
// This is the server side representation of the data. We expose this as a HPX
// component which allows for it to be created and accessed remotely through
// a global address (hpx::id_type).
struct partition_server
  : hpx::components::simple_component_base<partition_server> {

  // construct new instances
  partition_server() {}

  partition_server(graph_partition_data const& data)
    : graph_partition(data)
  {}

  partition_server(partition_server const& ps)
    : graph_partition(ps.graph_partition)
  {}


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
  // Experimenting... First with chaotic algorithm.
  // In this we will relax each vertex parallely
  //==============================================================
  void relax(const vertex_distance& vd, const partition_client_map_t& pmap);

  HPX_DEFINE_COMPONENT_ACTION(partition_server, relax,
			      dc_relax_action);

  // wait till all futures complete their work
  std::size_t flush_tasks() {
    //    mutex_type::scoped_lock l(mtx);

    std::size_t qsize = futures.size();
    std::cout << "flush task : " << qsize << std::endl;
    hpx::wait_all(futures);
    futures.clear(); // TODO check whether we need to do any de-allocation

    return qsize;
  }

  HPX_DEFINE_COMPONENT_ACTION(partition_server, flush_tasks,
			      dc_flush_action);

private:
  graph_partition_data graph_partition;
  future_collection_t futures; // collect all locally spawned future tasks
  mutex_type mtx;
};

HPX_REGISTER_ACTION_DECLARATION(partition_server::dc_relax_action, partition_relax_action);
HPX_REGISTER_ACTION_DECLARATION(partition_server::dc_flush_action, partition_flush_action);

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

typedef partition_server::dc_relax_action partition_relax_action;
HPX_REGISTER_ACTION(partition_relax_action);

typedef partition_server::dc_flush_action partition_flush_action;
HPX_REGISTER_ACTION(partition_flush_action);

// TODO encapsulate parameters

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
  {}

  // Attach a future representing a (possibly remote) partition.
  partition(hpx::future<hpx::id_type> && id)
    : base_type(std::move(id))
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

  ///////////////////////////////////////////////////////////////////////////
  // Invoke remote relax
  ///////////////////////////////////////////////////////////////////////////
  hpx::future<void> relax(vertex_distance const& vd,
			  const partition_client_map_t& pmap) const {
    partition_server::dc_relax_action act;
    return hpx::async(act, get_gid(), vd, pmap);
  }

  std::size_t flush() const {
    partition_server::dc_flush_action act;
    return hpx::async(act, get_gid()).get();
  }
};

#endif
