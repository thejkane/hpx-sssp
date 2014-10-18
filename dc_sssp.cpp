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

#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include <boost/shared_array.hpp>

#include "distributed_control.hpp"

bool verify = true; // Whether to verify results (i.e. shortest path or not)

///////////////////////////////////////////////////////////////////////////////
int hpx_main(boost::program_options::variables_map& vm) {

  boost::uint64_t scale = vm["scale"].as<boost::uint64_t>();   // Scale of the graph to be generated

  if (!vm.count("verify"))
    verify = false;

  // get the current locality and get the 
  // current thread id and print
  hpx::naming::id_type const here = hpx::find_here();
  std::cout << "Locality : " << hpx::get_locality_id() 
	    << " Here : " << here << std::endl;
  

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

  // Initialize and run HPX
  return hpx::init(desc_commandline, argc, argv);
}




