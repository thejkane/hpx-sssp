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
