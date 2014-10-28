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
#include "distributed_control.hpp"
#include <iostream>       // std::cout, std::endl
#include <limits>
#include <map>

// Whether to verify results (i.e. shortest path or not)
bool verify = true; 



struct distributed_control {
  
private:
  // create a partition client array
  // creates a partition for each locality
  // map key - locality id, value - remote partition
  partition_client_map_t partitions;

  int num_vert_per_local = 0;

  inline boost::uint32_t find_locality_id(vertex_t v,
					  int num_vertices_per_local) {
    HPX_ASSERT(num_vertices_per_local != 0);

    std::vector<hpx::naming::id_type> localities =
      hpx::find_all_localities();
    std::size_t num_locs = localities.size();

    vertex_t max_partition_vertex 
      = (num_locs * num_vertices_per_local) - 1;

    if (v > max_partition_vertex)
      return (num_locs-1);
    else 
      return (v / num_vertices_per_local);
  }

public:
  void run_chaotice_sssp(vertex_t source);

  void partition_graph(vertex_t source) {
    // TODO : Graph generation
    // Let's create 2 arrays 2 represent row_indices and columns
    // Then lets partition those 2 arrays - These 2 arrays represent the graph
    int numvertices = 7;
    int numcolumns = 22;

    int rowindices[] = {0, 3, 6, 10, 14, 17, 20, 22}; // size of row indices array = vertices + 1
    int columns[] = {1, 2, 4, 0, 2, 3, 0, 1, 3, 4, 1, 2, 5, 6, 0, 2, 5, 3, 4, 6, 3, 5};
    //    int weights[] = {5, 10, 8, 20, 12, 3, 10, 15, 3, 6, 10, 22, 35, 16, 20, 32, 25, 23, 34, 26, 33, 15};
    int weights[] = {20, 10, 8, 20, 12, 10, 10, 12, 3, 32, 10, 3, 35, 16, 8, 32, 34, 35, 34, 15, 16, 15};
    // TODO : Permute graph

    std::vector<hpx::naming::id_type> localities =
      hpx::find_all_localities();
    std::size_t num_locs = localities.size();
    std::cout << "Number of localities : " << num_locs << std::endl;

    // equally distribute vertices among localities
    num_vert_per_local = numvertices / num_locs;
    int vert_balance = numvertices % num_locs;

    for(int i=0; i<num_locs; ++i) {
      int startv = i*num_vert_per_local;

      // if this is last locality add balance vertices to last
      int endv;
      if (i == num_locs-1) {
	endv = num_vert_per_local+i*num_vert_per_local + 1 + vert_balance;
      } else {
	endv = num_vert_per_local+i*num_vert_per_local + 1;
      }
     
      int starte = rowindices[startv];
      int ende = rowindices[endv-1];
      
      std::cout << "startv : " << startv << " endv : " << endv
		<< " starte : " << starte << " ende : " << ende 
		<< std::endl;
      graph_partition_data pd(startv,
			      endv,
			      num_vert_per_local);

      pd.vertex_distances.resize(endv-startv);
      pd.vertex_distances.assign((endv-startv), 
				 std::numeric_limits<vertex_t>::max());
      pd.vertex_distances[source] = 0;
      
      // assign row indices
      for (int k=startv; k < endv; ++k) {
	pd.row_indices.push_back(rowindices[k]);
      }

      // assign columns and weights
      for (int k=starte; k < ende; ++k) {
	pd.columns.push_back(columns[k]);
	pd.weights.push_back(weights[k]);
      }

      // For debugging
      //pd.print();

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

  void print_results() {
    std::cout << "===================== Printing Results ==============================" << std::endl;
    partition_client_map_t::iterator ite = partitions.begin();
    for (; ite != partitions.end(); ++ite) {
      std::cout << "Partition locality : " << 
	hpx::naming::get_locality_id_from_id((*ite).second.get_gid()) << std::endl;
      // What we get from get_data is a future.
      // We have to call get to get the actual graph_partition_data
      // and then call print on it
      graph_partition_data pd = (*ite).second.get_data().get();

      int num_local_verts = (pd.vertex_end - pd.vertex_start) - 1;
      for(int i=0; i < num_local_verts; ++i) {
	std::cout << "vertex - " << (pd.vertex_start + i)
		  << " distance - " << pd.vertex_distances[i]
		  << std::endl;
      }
    }
  }

  // This function iterates all partitions (remote & local)
  // prints partition data
  void print_all_partitions() {
    partition_client_map_t::iterator ite = partitions.begin();
    for (; ite != partitions.end(); ++ite) {
      HPX_ASSERT(hpx::naming::get_locality_id_from_id((*ite).second.get_gid()) == (*ite).first);

      std::cout << "Partition locality : " << 
	hpx::naming::get_locality_id_from_id((*ite).second.get_gid()) << std::endl;
      // What we get from get_data is a future.
      // We have to call get to get the actual graph_partition_data
      // and then call print on it
      graph_partition_data pd = (*ite).second.get_data().get();
      //pd.print();

      { 
	// iterate vertices
	std::cout << "Vertices {";
	graph_partition_data::vertex_iterator vbegin = pd.vertices_begin();
	graph_partition_data::vertex_iterator vend = pd.vertices_end();
	for (; vbegin != vend; ++vbegin) {
	  std::cout << *vbegin << ", ";
	}
	std::cout << "}" << std::endl;
      }

      {
	// iterate edges
	graph_partition_data::edge_iterator ebegin = pd.edges_begin();
	graph_partition_data::edge_iterator eend = pd.edges_end();

	std::cout << "Edges {";
	for(; ebegin != eend; ++ebegin) {
	  std::cout << "[(" << (*ebegin).first << ", "
		    << (*ebegin).second << ")-" << (*ebegin).eid 
		    << "-" << pd.get_edge_weight(*ebegin)
		    << "], ";
	}
	std::cout << "}" << std::endl;
	std::cout << std::endl;
      }

      // traversing outgoing edges of vertices
      { 
	//typedef std::pair<graph_partition_data::edge_iterator, 
	//		  graph_partition_data::edge_iterator> OutgoingEdgePair_t;

	// iterate vertices
	graph_partition_data::vertex_iterator vbegin = pd.vertices_begin();
	graph_partition_data::vertex_iterator vend = pd.vertices_end();
	for (; vbegin != vend; ++vbegin) {
	  std::cout << "Vertex : " << *vbegin << " - Edges : {";

	  graph_partition_data::OutgoingEdgePair_t pair = pd.out_going_edges(*vbegin);
	  graph_partition_data::edge_iterator vebegin = pair.first;
	  graph_partition_data::edge_iterator veend = pair.second;
	  for(; vebegin != veend; ++vebegin) {
	    std::cout << "(" << (*vebegin).first << ", " << (*vebegin).second << "), ";
	  }
	  std::cout << "}" << std::endl;
	}
      }
    }
  }
};

//HPX_PLAIN_ACTION(distributed_control::relax, dc_relax_action);

void distributed_control::run_chaotice_sssp(vertex_t source) {
  // Find the locality of the source
  boost::uint32_t locality = find_locality_id(source, 
					      num_vert_per_local);
  // set distance to 0
  vertex_distance vd(source, 0);

  // The locality of the source might be different from
  // root locality. Therefore we need to bind and invoke the
  // function.
  // Get the partition the belongs to locality
  partition_client_map_t::iterator iteFind = partitions.find(locality);
  HPX_ASSERT(iteFind != partitions.end());

  //  future_collection_t futures;
  // relax source vertex
  // future_collection_t = vector <future <void> >
  hpx::future<future_collection_t> f = (*iteFind).second.relax<future_collection_t> (vd, partitions);
  hpx::future<future_collection_t> f = (*iteFind).second.relax (vd, partitions);

  hpx::wait_all(f.get());

  print_results();

}

//======================================================
// The actual relax code. For each adjacent edge we
// calculate the new distance and compare against
// existing distance in distance map. If distance is 
// better we update distance map with new distance
// and relax all adjacent edges with new distance.
// Updates to distance map are atomic.
//======================================================
// Here T = vector < future <void> >
template <typename T>
T partition_server::relax(const vertex_distance& vd) {

  std::cout << "Invoking relax in locality : "
	    << hpx::naming::get_locality_id_from_id(hpx::find_here())
	    << " for vertex : " << vd.vertex
	    << " and distance : " << vd.distance
	    << std::endl;

  // graph_partition.print();
  // Populate all future to a vector
  T futures;

  graph_partition_data::OutgoingEdgePair_t pair 
    = graph_partition.out_going_edges(vd.vertex);
  graph_partition_data::edge_iterator vebegin = pair.first;
  graph_partition_data::edge_iterator veend = pair.second;

  for(; vebegin != veend; ++vebegin) {
    std::cout << "Relaxing - (" << (*vebegin).first << ", " << (*vebegin).second << "), ";
    HPX_ASSERT(vd.vertex == (*vebegin).first);
    vertex_t target = (*vebegin).second;
      
    // calculate new distance
    int new_distance = vd.distance + 
      graph_partition.get_edge_weight(*vebegin);

    // check whether new distance is better than existing
    if (graph_partition.get_vertex_distance(target) > 
	new_distance){
      // update distance atomically
      if (graph_partition.set_vertex_distance_atomic
	  (target, new_distance)){
	// returned true. i.e. successfully updated.
	// time to relax. First find the appropriate
	// partition id (i.ee component id)
	boost::uint32_t target_locality = graph_partition.find_locality_id(target);

	partition_client_map_t::const_iterator iteClient =
	  pmap.find(target_locality);
	HPX_ASSERT(iteClient != pmap.end());

	hpx::future<T> f = (*iteClient).second.relax(vertex_distance(target, 
								     new_distance),
						     pmap);

	futures.insert(futures.end(), std::make_move_iterator(f.get().begin()),
		       std::make_move_iterator(f.get().end()));

      }
    }
  }

  return futures;
}


///////////////////////////////////////////////////////////////////////////////
// Note : Irrespective of number of localities we are running, hpx_main will only
//        get called for locality 0.
int hpx_main(boost::program_options::variables_map& vm) {

  boost::uint64_t scale = vm["scale"].as<boost::uint64_t>();   // Scale of the graph to be generated

  if (!vm.count("verify"))
    verify = false;

  vertex_t source = 0;
  distributed_control dc;
  dc.partition_graph(source);

  std::cout << "=============================================" << std::endl;
  std::cout << "=============================================" << std::endl;
  
  //  dc.print_all_partitions();
  dc.run_chaotice_sssp(source);
  
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




