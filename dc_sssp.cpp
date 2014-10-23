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

// This is needed to increase the number of 
// parameters to new operator
#define HPX_LIMIT 6


#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include <boost/format.hpp>
#include <boost/cstdint.hpp>
#include <boost/shared_array.hpp>
#include <boost/serialization/vector.hpp>

#include "distributed_control.hpp"

bool verify = true; // Whether to verify results (i.e. shortest path or not)

typedef std::vector<int> graph_array_t;
//=========================================
// Used to transfer partition information
// across different localities.
//=========================================
struct graph_partition_data {
  int vertex_start;
  int vertex_end;
  graph_array_t row_indices;
  graph_array_t columns;
  graph_array_t weights;

private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & vertex_start;
    ar & vertex_end;
    ar & row_indices;
    ar & columns;
    ar & weights;
  }

public:
  graph_partition_data() {}

  graph_partition_data(int _vstart, int _vend) : 
    vertex_start(_vstart), vertex_end(_vend) {}

  graph_partition_data(int _vstart, int _vend,
		       graph_array_t _ri,
		       graph_array_t _cl,
		       graph_array_t _wt) : 
    vertex_start(_vstart), vertex_end(_vend),
    row_indices(_ri), columns(_cl),
    weights(_wt)
  {}

  // copy constructor
  graph_partition_data(const graph_partition_data& other):
    vertex_start(other.vertex_start),
    vertex_end(other.vertex_end),
    row_indices(other.row_indices),
    columns(other.columns),
    weights(other.weights)
  {}

 void print() {
    std::cout << "Vertex start : " << vertex_start << " vertex end : " << vertex_end << std::endl;

    std::cout << "Row indices : {";
    // copy raw indices
    for(int i=0; i < ((vertex_end-vertex_start)+1); ++i) {
      std::cout << row_indices[i] << ", ";
    }
    std::cout << "}" << std::endl;

    std::cout << "Columns : {";
    // copy columns & weights
    for(int i=0; i<(row_indices[vertex_end] - row_indices[vertex_start]); ++i) {
      std::cout << columns[i] << ", ";
    }
    std::cout << "}" << std::endl;

    std::cout << "Weights : {";
    // copy columns & weights
    for(int i=0; i<(row_indices[vertex_end] - row_indices[vertex_start]); ++i) {
      std::cout << weights[i] << ", ";
    }
    std::cout << "}" << std::endl;
  }
};

///////////////////////////////////////////////////////////////////////////////
// This is the server side representation of the data. We expose this as a HPX
// component which allows for it to be created and accessed remotely through
// a global address (hpx::id_type).
struct partition_server
  : hpx::components::simple_component_base<partition_server>
{
  // construct new instances
  partition_server() {}

  partition_server(graph_partition_data const& data)
    : graph_partition(data)
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

private:
  graph_partition_data graph_partition;
};

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
  hpx::future<graph_partition_data> get_data() const
  {
    partition_server::get_data_action act;
    return hpx::async(act, get_gid());
  }
};


struct distributed_control {
  void partition_graph() {

    // TODO : Graph generation
    // Let's create 2 arrays 2 represent row_indices and columns
    // Then lets partition those 2 arrays - These 2 arrays represent the graph
    int numvertices = 7;
    int numcolumns = 22;

    int rowindices[] = {0, 3, 6, 10, 14, 17, 20, 22}; // size of row indices array = vertices + 1
    int columns[] = {1, 2, 4, 0, 2, 3, 0, 1, 3, 4, 1, 2, 5, 6, 0, 2, 5, 3, 4, 6, 3, 5};
    int weights[] = {5, 10, 8, 20, 12, 3, 10, 15, 3, 6, 10, 22, 35, 16, 20, 32, 25, 23, 34, 26, 33, 15};
    // TODO : Permute graph

    std::vector<hpx::naming::id_type> localities =
      hpx::find_all_localities();
    std::size_t num_locs = localities.size();
    std::cout << "Number of localities : " << num_locs << std::endl;


    // equally distribute vertices among localities
    int num_vert_per_local = numvertices / num_locs;

    for(int i=0; i<num_locs; ++i) {
      int startv = i*num_vert_per_local;
      int endv = num_vert_per_local+i*num_vert_per_local + 1;

      int starte = rowindices[startv];
      int ende = rowindices[endv];

      graph_partition_data pd(startv,
			      endv);

      // assign row indices
      for (int k=startv; k < endv; ++k) {
	pd.row_indices.push_back(rowindices[k]);
      }

      // assign columns and weights
      for (int k=starte; k < ende; ++k) {
	pd.columns.push_back(columns[k]);
	pd.weights.push_back(weights[k]);
      }

      pd.print();
    }
  }
};

///////////////////////////////////////////////////////////////////////////////
// Note : Irrespective of number of localities we are running, hpx_main will only
//        get called for locality 0.
int hpx_main(boost::program_options::variables_map& vm) {

  boost::uint64_t scale = vm["scale"].as<boost::uint64_t>();   // Scale of the graph to be generated

  if (!vm.count("verify"))
    verify = false;
  
  distributed_control dc;
  dc.partition_graph(); 
  
  return hpx::finalize();
}

int main(int argc, char* argv[])
{
  using namespace boost::program_options;

  options_description desc_commandline;
  desc_commandline.add_options()
    ("scale", value<boost::uint64_t>()->default_value(10),
     "The scale of the graph to be generated.")
    ("verify", "Verify results")
    ;

  std::cout << "Inside main ...." << std::endl;

  // Initialize and run HPX
  return hpx::init(desc_commandline, argc, argv);
}




