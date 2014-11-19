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
#include <iostream>
#include <stdexcept>
#include <string>

#include <boost/format.hpp>
#include <boost/cstdint.hpp>
#include <boost/graph/relax.hpp>

// Whether to verify results (i.e. shortest path or not)
bool verify = true; 
// If the completed count is much less than active count
// we can notify thread to run till it empties the queue
// for this we are using following boolean
//bool run_till_empty = false;
// This is to compare active count and completed
// count. 
//int active_to_complete_factor = 2;

// Max number of messages per coalesced buffer
std::size_t coalesced_message_size;

// whether to sort coalesced buffers
// default - yes
bool sort_coalesced_buffer = true;

// undirected graph
bool undirected = true;

// work stats
#ifdef WORK_STATS
boost::int64_t total_useful = 0;
boost::int64_t total_useless = 0;
boost::int64_t total_rejected = 0;
boost::int64_t total_invalidated = 0;
#endif

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

  // Count that flushing thread should yield
  boost::uint32_t yield_count;  

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
		      boost::uint32_t max_w,
		      boost::uint32_t yield): 
    scale(sc),
    num_qs(qs),
    max_weight(max_w),
    yield_count(yield),
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

  //========================================
  // Generates separate partitions
  //========================================
  void generate_distributed_graph() {

    std::vector<hpx::naming::id_type> localities =
      hpx::find_all_localities();
    std::size_t num_locs = localities.size();
    std::cout << "Number of localities : " << num_locs << std::endl;

     // equally distribute vertices among localities
    num_vert_per_local = m / num_locs;
    int vert_balance = m % num_locs;

    //std::cout << "partitions X : " << pdatas.size() << std::endl;

    std::vector< hpx::future<void> > futures;

    //seeds
    // Randome edge generation
    boost::uniform_int<uint64_t> 
      rand_64(0, std::numeric_limits<uint64_t>::max());

    boost::minstd_rand gen;
    uint64_t a = rand_64(gen);
    uint64_t b = rand_64(gen);

    for(std::size_t i=0; i<num_locs; ++i) {
      vertex_t startv = i*num_vert_per_local;

      // if this is last locality add balance vertices to last
      vertex_t endv;
      if (i == num_locs-1) {
	endv = num_vert_per_local+i*num_vert_per_local + vert_balance;
      } else {
	endv = num_vert_per_local+i*num_vert_per_local;
      }
     
      graph_partition_data pd(startv,
			      endv,
			      num_vert_per_local,
			      num_qs,
			      undirected);
           
      // Distribute data
      // To distribute we invoke component client, i.e. partition
      // and give locality and graph_partition. This operation will 
      // distribute graph_partition to remote locality
      std::cout << "Pushing to locality : " << 
	hpx::naming::get_locality_id_from_id(localities[i]) << std::endl;
      //pd.print();
      partition p(localities[i], pd);

      // generate local graphs
      futures.push_back(p.generate_local_graph(scale,
					       n,
					       max_weight,
					       a,
					       b));

      partitions.insert(std::make_pair(hpx::naming::get_locality_id_from_id(localities[i]), p));
    }

    hpx::wait_all(futures);
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
  // Reset remote and local counters for a next run
  //========================================================
  void reset_counters() {
    partition_client_map_t::iterator ite = partitions.begin();
    for (; ite != partitions.end(); ++ite) {
      (*ite).second.reset_counters();
    }
  }

  //========================================================
  // Small scale test case for debugging
  //========================================================
  void partition_graph_test(vertex_t source) {
    // TODO : Graph generation
    // Let's create 2 arrays 2 represent row_indices and columns
    // Then lets partition those 2 arrays - These 2 arrays represent the graph
    int numvertices = 7;

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

    for(std::size_t i=0; i<num_locs; ++i) {
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
			      num_qs,
			      undirected);

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

  //========================================================
  // Counts total number of vertices visited
  //========================================================
  boost::uint32_t count_total_visited_vertices() {
    
    boost::uint32_t visited = 0;

    std::vector< hpx::future<boost::uint32_t> > visited_counts;

    partition_client_map_t::iterator ite = partitions.begin();
    for (; ite != partitions.end(); ++ite) {
      // What we get from get_data is a future.
      // We have to call get to get the actual graph_partition_data
      // and then call print on it
      visited_counts.push_back((*ite).second.calculate_visited_vertices());
    }

    std::vector< hpx::future<boost::uint32_t> >::iterator iteFutures 
      = visited_counts.begin();

    for(; iteFutures != visited_counts.end(); ++iteFutures) {
      visited += (*iteFutures).get();
    }
    
    std::cout << "Total visited vertices : " << visited << std::endl;
    return visited;
  }


  //========================================================
  // Verify generated distances
  //========================================================
  void verify_results() {
    std::cout << "===================== Verifying Results ==============================" << std::endl;
    std::vector< hpx::future<void> > futures;
    partition_client_map_t::iterator ite = partitions.begin();
    for (; ite != partitions.end(); ++ite) {
      futures.push_back((*ite).second.verify_distance_results(partitions));
    }

    hpx::wait_all(futures);
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


//============================================
// Verifies whether calculated distances 
// are correct.
//============================================
void partition_server::verify_partition_results(const partition_client_map_t& partitions) {

  std::cout << "Verifying partition : " << hpx::naming::get_locality_id_from_id(hpx::find_here())
	    << std::endl;

  graph_partition_data::vertex_iterator ite =
    graph_partition.vertices_begin();
  for(; ite != graph_partition.vertices_end(); ++ite) {
    graph_partition_data::OutgoingEdgePair_t p = 
      graph_partition.out_going_edges(*ite);
    graph_partition_data::edge_iterator eite =
      p.first;
	
    for (; eite != p.second; ++eite) {
      EdgeType_t e = *eite;

      // edge weight starting from vertex *ite
      edge_property_t w = graph_partition.get_edge_weight(e);
	  
      // distance for source
      vertex_property_t vs = graph_partition.get_vertex_distance(*ite);

      // get the target vertex; but target vertex might be in a 
      // different locality. Therefore we need to find the appropriate
      // locality and find graph data

      // Find the locality of the target
      boost::uint32_t locality = graph_partition.find_locality_id(e.second);
      // Get the partition the belongs to locality
      partition_client_map_t::const_iterator iteFind = partitions.find(locality);
      HPX_ASSERT(iteFind != partitions.end());

      // distance for target
      vertex_property_t vt = (*iteFind).second.get_vertex_distance(e.second);
	  
      if (vt > boost::closed_plus<vertex_property_t>()(vs, w)) {
	// failure
	std::cout << "The target vertex : " << e.second
		  << " distance : " << vt
		  << " and source vertex : " << e.first 
		  << " distance : " << vs << " + "
		  << " weight : " << w
		  << " does not match." 
		  << std::endl;
	HPX_ASSERT(false);
      }
    }
	  
  }
}


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
    (*ite).second.start_flush_tasks(num_qs, partitions, yield_count);
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
  // source message is active, increase active count
  active_count++;
  (*iteFind).second.relax(vd);

  // termination
  std::vector<hpx::id_type> localities = hpx::find_all_localities();

  int phase = 1; // There are 2 phases
  boost::int64_t prev_g_completed;

  bool terminate = false;

  while(!terminate) {
    boost::int64_t g_completed = hpx::lcos::reduce<total_completed_count_action>
      (localities, std::plus<boost::int64_t>(), partitions).get();
    boost::int64_t g_active = hpx::lcos::reduce<total_active_count_action>
      (localities, std::plus<boost::int64_t>(), partitions).get();

    // We need to add 1 to global count
    // Why ? We match every send with a receive. But for source
    // we do not have a send. Therefore we need add 1 to g_completed
    //g_completed = g_completed + 1;

    //#ifdef PRINT_DEBUG
    std::cout << "g_completed : " << g_completed
	      << " g_active : " << g_active
	      << std::endl;
    //#endif
    
    if (g_completed != g_active) {
      phase = 1;
      prev_g_completed = g_completed;
      hpx::this_thread::yield();
    } else if (phase == 1 ||
	       g_completed != prev_g_completed) {
      phase = 2;
      prev_g_completed = g_completed; 
      hpx::this_thread::yield();
    } else {
      terminate = true;
    }
  }

#ifdef PRINT_DEBUG
  std::cout << "TERMINATION DETECTION DONE" << std::endl;
#endif

  std::vector< hpx::future<void> > terminating_ps;
  // invoke termination parallely
  for(ite = partitions.begin(); ite != partitions.end(); ++ite) {
    terminating_ps.push_back((*ite).second.terminate());
  }

  // wait for everyone to terminate
  hpx::when_all(terminating_ps).get();

#ifdef PRINT_DEBUG
  std::cout << "TERMINATION DONE" << std::endl;
  print_results();
#endif

  if (verify) {
    verify_results();
  }

  // collect work stats if enabled
#ifdef WORK_STATS
  total_useful = hpx::lcos::reduce<total_useful_work_action>
      (localities, std::plus<boost::int64_t>()).get();
  total_useless = hpx::lcos::reduce<total_useless_work_action>
      (localities, std::plus<boost::int64_t>()).get(); 
  total_invalidated = hpx::lcos::reduce<total_invalidated_work_action>
      (localities, std::plus<boost::int64_t>()).get();
  total_rejected = hpx::lcos::reduce<total_rejected_work_action>
      (localities, std::plus<boost::int64_t>()).get();  
#endif

}

//======================================================
// Receives a vector of messages. The vector is sorted
// This will inturn call relax.
//======================================================
void partition_server::coalesced_relax(const coalesced_message_t& vds) {
  coalesced_message_t::const_iterator ite = vds.begin();
  for(; ite != vds.end(); ++ite) {
    relax(*ite);
  }
}

//======================================================
// The actual relax code. For each adjacent edge we
// calculate the new distance and compare against
// existing distance in distance map. If distance is 
// better we update distance map with new distance
// and relax all adjacent edges with new distance.
// Updates to distance map are atomic.
//======================================================
void partition_server::relax(const vertex_distance& vd) {

#ifdef PRINT_DEBUG
  std::cout << "Invoking relax in locality : "
	    << hpx::naming::get_locality_id_from_id(hpx::find_here())
	    << " for vertex : " << vd.vertex
	    << " and distance : " << vd.distance
	    << " thread id : " << hpx::get_worker_thread_num()
	    << std::endl;
#endif

  // graph_partition.print();
#ifdef PRINT_DEBUG
  std::cout << "v - " << vd.vertex << " dis - " << vd.distance << " stored distance - "
	    << graph_partition.get_vertex_distance(vd.vertex) << std::endl;
#endif

  // check whether new distance is better than existing
  if (graph_partition.get_vertex_distance(vd.vertex) > 
      vd.distance){
    // update distance atomically
    if (graph_partition.set_vertex_distance_atomic
	(vd.vertex, vd.distance)) {
#ifdef PRINT_DEBUG
      std::cout << "Inserting to queue vertex : " << vd.vertex
		<< " distance : " << vd.distance << std::endl;
#endif
      // Get a random priority queue
      int idx = select_random_q_idx();
      buckets[idx].push(vd);
    } else {
      // We did not put vertex to queue
      // Therefore the vertex-distance should be marked as
      // completed.
      // completed processing a message
      // increase completed count
      completed_count++;

#ifdef WORK_STATS
      useless++;
#endif
    }
  } else {
    // We did not put vertex to queue
    // Therefore the vertex-distance should be marked as
    // completed.
    // completed processing a message
    // increase completed count
    completed_count++;

    //TODO not sure
#ifdef WORK_STATS
      useless++;
#endif
  }
}

//======================================================
// Starts flush tasks
//======================================================
void partition_server::flush_tasks(int idx,
				   const partition_client_map_t& pmap,
				   const boost::uint32_t yield_count) {  
  HPX_ASSERT(0 <= idx && idx < graph_partition.num_queues);
  buckets[idx].handle_queue(pmap, graph_partition, yield_count);
}

//======================================================
// Send coalesced messages. The vd will be put into 
// appropriate target buffer. If buffer has reached the
// coalescing size limit messages will be sent.
//======================================================
void dc_priority_queue::send(const vertex_distance vd,
			     boost::uint32_t target_locality,
			     const partition& partition_client) {

  // We need to interleave send and send_all. i.e. when we 
  // are doing a send_all we shouldnt be doing a send
  // also this doesnt create too much contention as we have destination
  // buffers per each priority queue. So when sending data we only have
  // single thread working inside this function.
  boost::mutex::scoped_lock scopedLock(cmap_mutex);      
  coalsced_message_map_t::iterator iteFind = cmap.find(target_locality);
  HPX_ASSERT(iteFind != cmap.end());

  // put message to buffer
  (*iteFind).second.push_back(vd);

  // check whether we have reached coalescing limit
  if ((*iteFind).second.size() == coalesced_message_size) {
    // do sorting if sorting enabled
    if (sort_coalesced_buffer) {
      std::sort((*iteFind).second.begin(), (*iteFind).second.end(), sc);
    }

    active_count += (*iteFind).second.size(); // atomic addition
    // send coalesced message
    partition_client.coalesced_relax((*iteFind).second);
    (*iteFind).second.clear();
  }
}

//======================================================
// When queues are empty we need to send the data
// that is not send in buffers. This function will
// go through all the buffers for all destinations and
// clear them up.
//======================================================
void dc_priority_queue::send_all(const partition_client_map_t& pmap) {

  // See for explanation in send function
  boost::mutex::scoped_lock scopedLock(cmap_mutex);      

  coalsced_message_map_t::iterator ite = cmap.begin();
  for (; ite != cmap.end(); ++ite) {
    // do sorting if sorting enabled
    if (sort_coalesced_buffer) {
      std::sort((*ite).second.begin(), (*ite).second.end(), sc);
    }

#ifdef PRINT_DEBUG
    std::cout << "Sending to locality id : " << (*ite).first << std::endl;
#endif
    partition_client_map_t::const_iterator iteClient = pmap.find((*ite).first);
    HPX_ASSERT(iteClient != pmap.end());

    if (!(*ite).second.empty()) {
      active_count += (*ite).second.size(); // atomic addition
      (*iteClient).second.coalesced_relax((*ite).second);
      (*ite).second.clear();
    }
  }
}

//======================================================
// This function will go through priority queue
// found by idx and will relax all edges found in 
// the queue.
//======================================================
void dc_priority_queue::handle_queue(const partition_client_map_t& pmap,
				     graph_partition_data& graph_partition,
				     const boost::uint32_t yield_count) {

  boost::uint32_t ite_count = 0;

  while(!termination) {
    // If queue is empty wait till an element is pushed
    {
      boost::mutex::scoped_lock scopedLock(mutex);
      if (pq.empty()) {
	// increase empty q count
	increase_empty_q_count();
	cv.wait(scopedLock);
	// decrease empty q count
	decrease_empty_q_count();
      }
    }
    
    // TODO we need to lock this ? ...
    while(!pq.empty()) {
      ++ite_count;

      // if iteration count is equal to yield count
      // then yield current thread
      if (ite_count == yield_count) {
#ifdef PRINT_DEBUG
	std::cout << "Invoking thread yield .......... Reached thread count " << yield_count << std::endl;
#endif
	// reset counter
	ite_count = 0;
	hpx::this_thread::yield();
      }

      bool was_empty = false;
      vertex_distance vd;
      // lock and pop the element
      {
	boost::mutex::scoped_lock scopedLock(mutex);
	vd = pq.top();
	pq.pop();

	was_empty = pq.empty();
      }
    
      // relax vd
      graph_partition_data::OutgoingEdgePair_t pair 
	= graph_partition.out_going_edges(vd.vertex);
      graph_partition_data::edge_iterator vebegin = pair.first;
      graph_partition_data::edge_iterator veend = pair.second;
      
      for(; vebegin != veend; ++vebegin) {

#ifdef PRINT_DEBUG
	std::cout << "Relaxing - (" << (*vebegin).first << ", " 
		  << (*vebegin).second << "), Termination - " 
		  << termination << std::endl;
#endif

	HPX_ASSERT(vd.vertex == (*vebegin).first);
	vertex_t target = (*vebegin).second;
      
	// calculate new distance
	int new_distance = vd.distance + 
	  graph_partition.get_edge_weight(*vebegin);

	boost::uint32_t target_locality = graph_partition.find_locality_id(target);

	// check whether distance already in target
	// is better than what we calculated (new_distance)
	partition_client_map_t::const_iterator iteClient =
	  pmap.find(target_locality);
	HPX_ASSERT(iteClient != pmap.end());
	
	if ((*iteClient).second.get_vertex_distance(target) >
	    new_distance) {
#ifdef PRINT_DEBUG
	  std::cout << "Target locality : " << target_locality
		    << " for vertex : " << target
		    << " new distance : " << new_distance
		    << " thread id : " << hpx::get_worker_thread_num()
		    << std::endl;
#endif
	  // We are spawning a new message. Therefore increase
	  // active count
	  send(vertex_distance(target, new_distance),
	       target_locality,
	       (*iteClient).second);
        } else {
#ifdef WORK_STATS
          rejected++;
#endif
        } 
      }

      // if(was_empty) {
      completed_count++;
	//}
    }
  }

  // Termination took place;                                                                                                                                                             
  // All buffers must be empty                                                                                                                                                           
  coalsced_message_map_t::iterator ite = cmap.begin();
  for (; ite != cmap.end(); ++ite) {
    if (!(*ite).second.empty()) {
      std::cout << "All buffers must be empty " << std::endl;
      HPX_ASSERT(false);
    }
  }

  // All queues must be empty                                                                                                                                                            
  HPX_ASSERT(pq.empty());

}

typedef std::vector<boost::uint64_t> all_timing_t;

//======================================================
// Prints total timing results.
//======================================================
void print_summary_results(
			   boost::uint32_t num_localities,
			   boost::uint64_t num_os_threads,
			   const all_timing_t& all_timings,
			   boost::uint32_t sc,
			   boost::uint32_t mw,
			   boost::uint32_t n_qs,
			   boost::uint32_t num_sources,
			   boost::int64_t cum_useful = 0,
			   boost::int64_t cum_useless = 0,
			   boost::int64_t cum_invalidated = 0,
			   boost::int64_t cum_rejected = 0) {

  // calculate avg timing
  all_timing_t::const_iterator ite = all_timings.begin();
  boost::uint64_t tot_readings = all_timings.size();
  boost::uint64_t total_time = 0;
  for (; ite != all_timings.end(); ++ite) {
    total_time += (*ite);
  } 

  boost::uint64_t avg_time = total_time / tot_readings;

  std::string const locs_str = boost::str(boost::format("%u,") % num_localities);
  std::string const threads_str = boost::str(boost::format("%lu,") % num_os_threads);
  std::string const scale = boost::str(boost::format("%u,") % sc);
  std::string const max_weight = boost::str(boost::format("%u,") % mw);
  std::string const num_qs = boost::str(boost::format("%u ") % n_qs);
  std::string const num_srces = boost::str(boost::format("%u ") % num_sources);
  std::string const coal_sz = boost::str(boost::format("%u ") % coalesced_message_size);

  std::string sorting = "true";
  if (!sort_coalesced_buffer)
    sorting = "false";

  std::cout << "Total time to execute DC : " << (total_time / 1e9)
	    << " (Per Source Avg : " << ((total_time / tot_readings) / 1e9) << ")"
	    << ", Localities : " << num_localities
	    << ", OS Threads : " << num_os_threads
	    << ", Scale : " << sc
	    << ", Max Weight : " << mw
	    << ", Number of Queues : " << n_qs
	    << ", Number of Sources : " << num_sources
	    << ", Coalescing Size : " << coalesced_message_size
	    << ", Outbuffer Sorting : " << sorting
	    << std::endl;

#ifdef WORK_STATS
  std::cout << "Work Stats (per source) - Useful : " << (cum_useful / tot_readings)
	    << ", Useless : " << (cum_useless / tot_readings)
	    << ", Invalidated :" << (cum_invalidated / tot_readings)
	    << ", Rejected :" << (cum_rejected / tot_readings)
	    << ", Rejected/Useful Ratio :" << ((double)cum_rejected / (double)cum_useful)
	    << ", Invalidated/Useful Ratio :" << ((double)cum_invalidated / (double)cum_useful)
	    << std::endl;
#endif
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
  boost::uint32_t num_sources = vm["num_sources"].as<boost::uint32_t>();
  boost::uint32_t yield_count = vm["yield_count"].as<boost::uint32_t>();
  coalesced_message_size = vm["coalescing_size"].as<boost::uint32_t>();

  if (!vm.count("verify"))
    verify = false;

  if (vm.count("disable_sort"))
    sort_coalesced_buffer = false;


  std::cout << "[Input Read] Graph scale : " << scale 
	    << ", queues : " << queues 
	    << ", max_weight : " << max_weight 
	    << ", num_sources : " << num_sources 
	    << ", yield_count : " << yield_count 
	    << ", coalescing_size : " << coalesced_message_size
	    << ", disable_sort : " << sort_coalesced_buffer 
	    << ", verify results : " << verify << std::endl;

  vertex_t source = 0;
  distributed_control dc(scale, queues, max_weight, yield_count);

  std::cout << "Generating the graph ..." << std::endl;
  std::cout << "Vertices : " << dc.m << " Edges : " << dc.n << std::endl;

  {
    // Create the graph
    //hpx_csr_graph g(dc.m, dc.n, true);
    // Generate edges
    //dc.generate_graph(g);

    dc.generate_distributed_graph();
    std::cout << "Graph generation ....... done" << std::endl;
    //dc.print_all_partitions();
#ifdef DC_TEST_CASE
    //dc.partition_graph_test(source);
#else
    //dc.partition_graph(g);
#endif
  } // HPX CSR graph is only needed for graph generation
  // Once graph is partitioned & distributed we dont need to
  // holde on to HPX CSR graph. The memory can be released.
  // Therefore we need these additional braces.

  std::cout << "=============================================" << std::endl;
  std::cout << "=============================================" << std::endl;
  
  all_timing_t all_timings;

#ifdef WORK_STATS
  boost::int64_t cum_total_useful = 0;
  boost::int64_t cum_total_useless = 0;
  boost::int64_t cum_total_rejected = 0;
  boost::int64_t cum_total_invalidated = 0;
#endif

#ifndef DC_TEST_CASE
  for (boost::uint32_t i=0; i<num_sources; ++i) {
    // Select a random source vertex
    source = dc.select_random_source();
#endif
    dc.reset_counters();

    std::cout << "Running distributed control for scale : " << scale
	      << ", num_queues : " << queues
	      << ", yield count : " << yield_count
	      << ", max_weight : " << max_weight
	      << ", coalescing_size : " << coalesced_message_size
	      << ", disable_sort : " << sort_coalesced_buffer 
	      << ", source : " << source 
	      << ", iteration : " << i << std::endl;
    
    boost::uint64_t before = hpx::util::high_resolution_clock::now();
    dc.run_dc(source);
    boost::uint64_t after = hpx::util::high_resolution_clock::now();    

    // check the total visited count
    boost::uint32_t tot_visited = dc.count_total_visited_vertices();
    // if scale > 10 and total visited is less than 100 ignore the run
    /*    if (scale > 10 && tot_visited < 100) {
      std::cout << "decrementing iteration ... " << std::endl;
      --i;
      continue;
      }*/
    
    boost::uint64_t elapsed = after - before;
    std::cout << "Time for distributed control run with scale : " << scale
	      << ", num_queues : " << queues
	      << ", yield count : " << yield_count
	      << ", max_weight : " << max_weight
	      << ", coalescing_size : " << coalesced_message_size
	      << ", disable_sort : " << sort_coalesced_buffer 
	      << ", source : " << source 
	      << ", total visited count : " << tot_visited
	      << ", iteration : " << i 
	      << " is : " << (elapsed / 1e9) << std::endl;

// work stats
#ifdef WORK_STATS
    std::cout << "Work statistics - Useful Work : " << total_useful
	      << ", Useless Work : " << total_useless
	      << ", Rejected Work : " << total_rejected
	      << ", Invalidated Work : " << total_invalidated
	      << std::endl;

    cum_total_useful += total_useful;
    cum_total_useless += total_useless;
    cum_total_rejected += total_rejected;
    cum_total_invalidated += total_invalidated;
#endif

    all_timings.push_back(elapsed);

#ifndef DC_TEST_CASE
    }
#endif

    boost::uint64_t const num_worker_threads = hpx::get_num_worker_threads();
    hpx::future<boost::uint32_t> locs = hpx::get_num_localities();

#ifdef WORK_STATS
    print_summary_results(locs.get(),
			  num_worker_threads,
			  all_timings,
			  scale,
			  max_weight,
			  queues,
			  num_sources,
			  cum_total_useful,
			  cum_total_useless,
			  cum_total_invalidated,
			  cum_total_rejected);
#else
    print_summary_results(locs.get(),
			  num_worker_threads,
			  all_timings,
			  scale,
			  max_weight,
			  queues,
			  num_sources);
#endif

  
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
    ("num_sources", value<boost::uint32_t>()->default_value(18),
     "The number of sources to run.")
    ("coalescing_size", value<boost::uint32_t>()->default_value(1000),
     "The maximum messages to be coalesced.")
    ("yield_count", value<boost::uint32_t>()->default_value(100),
     "After how many counts the flushing thread should yield, so that it give space to run other threads.")
    ("verify", "Verify results")
    ("disable_sort", "Disable sorting sending message buffer")
    ;

  std::cout << "Inside main ...." << std::endl;

  // Initialize and run HPX
  return hpx::init(desc_commandline, argc, argv);
}




