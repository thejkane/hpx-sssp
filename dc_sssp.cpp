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

#include <map>

// Whether to verify results (i.e. shortest path or not)
bool verify = true; 



struct distributed_control {
  
private:
  // create a partition client array
  // creates a partition for each locality
  // map key - locality id, value - remote partition
  typedef std::map<boost::uint32_t, partition> partition_client_map_t;
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
  vertex_distance vd(source, 0);

  // The locality of the source might be different from
  // root locality. Therefore we need to bind and invoke the
  // function.
  // Get the partition the belongs to locality
  partition_client_map_t::iterator iteFind = partitions.find(locality);
  HPX_ASSERT(iteFind != partitions.end());

  // Time to invoke relax for source
  partition_relax_action act;
  hpx::future<void> f = hpx::async(act, (*iteFind).second.get_gid(), vd);
  //				   vd, (*iteFind).second);

  f.get(); // wait till relax is done
}


///////////////////////////////////////////////////////////////////////////////
// Note : Irrespective of number of localities we are running, hpx_main will only
//        get called for locality 0.
int hpx_main(boost::program_options::variables_map& vm) {

  boost::uint64_t scale = vm["scale"].as<boost::uint64_t>();   // Scale of the graph to be generated

  if (!vm.count("verify"))
    verify = false;
  
  distributed_control dc;
  dc.partition_graph();

  std::cout << "=============================================" << std::endl;
  std::cout << "=============================================" << std::endl;
  
  //  dc.print_all_partitions();
  dc.run_chaotice_sssp(5);
  
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




