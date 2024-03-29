// Copyright (C) 2018 Thejaka Amila Kanewala, Marcin Zalewski, Andrew Lumsdaine.

// Boost Software License - Version 1.0 - August 17th, 2003

// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:

// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

//  Authors: Thejaka Kanewala
//           Andrew Lumsdaine

#ifndef SSSP_DC_COMMON
#define SSSP_DC_COMMON

#define INVALID_VERTEX -1
#define INVALID_EDGE -1
#define INVALID_EDGE_WEIGHT -1

#ifdef GRAPH_64
typedef int64_t edge_t;
typedef int64_t vertex_t;
typedef int64_t array_t;
typedef int64_t edge_property_t;
typedef int64_t vertex_property_t;
#else
typedef int32_t edge_t;
typedef int32_t vertex_t;
typedef int32_t array_t;
typedef int32_t edge_property_t;
typedef int32_t vertex_property_t;
#endif

// Edge iterator type. Pair <source, target>
struct EdgeType_t {
  vertex_t first;
  vertex_t second;
  edge_t eid;

  EdgeType_t(vertex_t f, vertex_t s):first(f), second(s), eid(-1){}
  // Constructor with edge id
  EdgeType_t(vertex_t f, vertex_t s, edge_t id):first(f), second(s), eid(id){}
};
#endif
