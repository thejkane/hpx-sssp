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
#include "hpx_csr.hpp"

#include <iostream>
#include <limits>
#include <map>

#include <boost/random/linear_congruential.hpp>
#include <boost/graph/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/graph/graph500_generator.hpp>

// Whether to verify results (i.e. shortest path or not)
bool verify = true; 

struct distributed_control {
  
private:
  // create a partition client array
  // creates a partition for each locality
  // map key - locality id, value - remote partition
  partition_client_map_t partitions;

  // The graph scale
  boost::int32_t scale;  

  // Represents number of queues to be run per locality
  int num_qs;

  // Number of vertitions per locality
  int num_vert_per_local = 0;

  // Maximum weight per edge
  edge_property_t max_weight;

  // Random number generator is needed to select a random vertex
  boost::random::mt19937 gen;

  // average number of outgoing edges
  double edgefactor = 16;

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
  // Number of vertices
  vertex_t m;

  // Number of edges
  edge_t n;

  //========================================================
  // Constructor
  // Args : sc - scale, qs - number of queues to spawn
  //        max_w - maximum weight
  //========================================================
  distributed_control(boost::uint32_t sc,
		      boost::uint32_t qs,
		      boost::uint32_t max_w): 
    scale(sc),
    num_qs(qs),
    max_weight(max_w),
    m((unsigned long long)(1) << scale),
    n(static_cast<edge_t>(floor(m * edgefactor)))
  {}

  //========================================
  // Select a random source vertex
  //========================================
  vertex_t select_random_source() {
    boost::random::uniform_int_distribution<> dist(1, m);
    return (dist(gen) - 1);
  }

  //========================================================
  // Generates edges for the given graph
  // Mainly graph generation happens here.
  //========================================================
  void generate_graph(hpx_csr_graph& g) {

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

    boost::minstd_rand gen;
    uint64_t a = rand_64(gen);
    uint64_t b = rand_64(gen);

    g.addEdges(Graph500Iter(scale, 0, a, b), 
	       Graph500Iter(scale, n, a, b), die);
  }

  //========================================================
  // Invokes actual distributed control algorithm.
  //========================================================
  void run_dc(vertex_t source);

  //========================================================
  // Small scale test case for debugging
  //========================================================
  void partition_graph_test(vertex_t source) {
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

      num_qs = 5;

      std::cout << "startv : " << startv << " endv : " << endv
		<< " starte : " << starte << " ende : " << ende 
		<< std::endl;
      graph_partition_data pd(startv,
			      endv,
			      num_vert_per_local,
			      num_qs);

      pd.vertex_distances.resize(endv-startv);
      pd.vertex_distances.assign((endv-startv), 
				 std::numeric_limits<vertex_t>::max());
      //      pd.vertex_distances[source] = 0;
     
      
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

  //========================================================
  // Partition the graph for the running localities
  //========================================================
  void partition_graph(hpx_csr_graph& g) {
    g.partition_graph(partitions, num_qs, num_vert_per_local);
    HPX_ASSERT(num_vert_per_local != 0);
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

  //========================================================
  // This function iterates all partitions (remote & local)
  // prints partition data
  //========================================================
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


//======================================================
// Start invoking Distributed Control algorithm.
// source - The vertex we need to calculate distance
// from.
//======================================================
void distributed_control::run_dc(vertex_t source) {
  
  // start flush tasks
  partition_client_map_t::iterator ite = 
    partitions.begin();
  for(; ite != partitions.end(); ++ite) {
    (*ite).second.start_flush_tasks(num_qs, partitions);
  }

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
  (*iteFind).second.relax(vd, partitions);

  // termination
  std::vector<hpx::id_type> localities = hpx::find_all_localities();
  boost::int64_t result = -1;
  int attempt = 0;
  while(result != 0 || attempt != 2) {
    result = hpx::lcos::reduce<total_msg_difference_action>
      (localities, std::plus<boost::int64_t>()).get();
    result = result - 1;

    std::cout << "result - " << result
	      << " attempt - " << attempt << std::endl;

    if (result == 0)
      ++attempt;
    else
      attempt = 0;
  }

  std::vector< hpx::future<void> > terminating_ps;
  // invoke termination parallely
  for(ite = partitions.begin(); ite != partitions.end(); ++ite) {
    terminating_ps.push_back((*ite).second.terminate());
  }

  // wait for everyone to terminate
  hpx::when_all(terminating_ps).get();

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
void partition_server::relax(const vertex_distance& vd,
			     const partition_client_map_t& pmap) {

  std::cout << "Invoking relax in locality : "
	    << hpx::naming::get_locality_id_from_id(hpx::find_here())
	    << " for vertex : " << vd.vertex
	    << " and distance : " << vd.distance
	    << std::endl;

  // We got a message so increase receive count
  receive_count++;

  // graph_partition.print();

  std::cout << "v - " << vd.vertex << " dis - " << vd.distance << " stored distance - "
	    << graph_partition.get_vertex_distance(vd.vertex) << std::endl;
  // check whether new distance is better than existing
  if (graph_partition.get_vertex_distance(vd.vertex) > 
      vd.distance){
    // update distance atomically
    if (graph_partition.set_vertex_distance_atomic
	(vd.vertex, vd.distance)) {

      // Get a random priority queue
      int idx = select_random_q_idx();
      buckets[idx].push(vd);
    }
  }
}

//======================================================
// Starts flush tasks
//======================================================
void partition_server::flush_tasks(int idx,
				   const partition_client_map_t& pmap) {  
  HPX_ASSERT(0 <= idx && idx < graph_partition.num_queues);
  buckets[idx].handle_queue(pmap, graph_partition);
}

//======================================================
// This function will go through priority queue
// found by idx and will relax all edges found in 
// the queue.
//======================================================
void dc_priority_queue::handle_queue(const partition_client_map_t& pmap,
				     graph_partition_data& graph_partition) {

  while(!termination) {
    // If queue is empty wait till an element is pushed
    {
      boost::mutex::scoped_lock scopedLock(mutex);
      if (pq.empty()) {
	q_empty = true; // This helps termination
	cv.wait(scopedLock);
      }
    }
    
    if (!pq.empty())
      q_empty = false;

    // TODO we need to lock this ? ...
    while(!pq.empty()) {
      vertex_distance vd;

      // lock and pop the element
      {
	boost::mutex::scoped_lock scopedLock(mutex);
	vd = pq.top();
	pq.pop();
      }
    
      // relax vd
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

	boost::uint32_t target_locality = graph_partition.find_locality_id(target);

	partition_client_map_t::const_iterator iteClient =
	  pmap.find(target_locality);
	HPX_ASSERT(iteClient != pmap.end());

	(*iteClient).second.relax(vertex_distance(target, 
						  new_distance),
				  pmap);
	// sending messages. increase send count
	send_count++;
    }
  }
}
}



//========================================================
// Note : Irrespective of number of localities we are running, 
// hpx_main will only
//        get called for locality 0.
//========================================================
int hpx_main(boost::program_options::variables_map& vm) {

  // Scale of the graph to be generated
  boost::uint32_t scale = vm["scale"].as<boost::uint32_t>();
  boost::uint32_t queues = vm["queues"].as<boost::uint32_t>();
  boost::uint32_t max_weight = vm["max_weight"].as<boost::uint32_t>();

  std::cout << "[Input Read] Graph scale : " << scale 
  << ", queues : " << queues 
  << ", max_weight : " << max_weight << std::endl;

  if (!vm.count("verify"))
    verify = false;

  vertex_t source = 0;
  distributed_control dc(scale, queues, max_weight);

  std::cout << "Generating the graph ..." << std::endl;
  std::cout << "Vertices : " << dc.m << " Edges : " << dc.n << std::endl;

  {
    // Create the graph
    hpx_csr_graph g(dc.m, dc.n, true);
    // Generate edges
    dc.generate_graph(g);

    std::cout << "Graph generation ....... done" << std::endl;

#ifdef DC_TEST_CASE
    dc.partition_graph_test(source);
#else
    dc.partition_graph(g);
#endif
  } // HPX CSR graph is only needed for graph generation
  // Once graph is partitioned & distributed we dont need to
  // holde on to HPX CSR graph. The memory can be released.
  // Therefore we need these additional braces.

  std::cout << "=============================================" << std::endl;
  std::cout << "=============================================" << std::endl;
  
#ifndef DC_TEST_CASE
  // Select a random source vertex
  source = dc.select_random_source();
#endif
  
  std::cout << "Running distributed control for source : " 
            << source << std::endl;

  dc.run_dc(source);
  
  return hpx::finalize();
}


int main(int argc, char* argv[])
{
  using namespace boost::program_options;

  options_description desc_commandline;
  desc_commandline.add_options()
    ("scale", value<boost::uint32_t>()->default_value(10),
     "The scale of the graph to be generated.")
    ("queues", value<boost::uint32_t>()->default_value(5),
     "The number of queues per each locality.")
    ("max_weight", value<boost::uint32_t>()->default_value(100),
     "The number of queues per each locality.")
    ("verify", "Verify results")
    ;

  std::cout << "Inside main ...." << std::endl;

  // Initialize and run HPX
  return hpx::init(desc_commandline, argc, argv);
}




