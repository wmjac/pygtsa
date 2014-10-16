/* Copyright (C) 2014 William M. Jacobs
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GRAPH_H__
#define __GRAPH_H__

#include <stdlib.h>

typedef struct
{
  unsigned short v1, v2;
}
edge_t;

typedef struct
{
  size_t max_nvertices, nvertices;
  size_t max_nedges, nedges;
  size_t * nvedges;
  unsigned short ** edges;
  size_t nbridges;
  edge_t * bridges;
  size_t nremovable;
  edge_t * removable;
  size_t nadjacents;
  edge_t * adjacents;

  /* Workspace for dfs bridge detection */
  size_t _dfs_iter;
  char * _dfs_visited;
  unsigned short * _dfs_discovery;
  unsigned short * _dfs_low;
  unsigned short * _dfs_parent;
}
graph_t;

graph_t * graph_alloc (size_t nv);
void graph_free (graph_t * graph);
void graph_copy (const graph_t * orig, graph_t * dest);

void graph_add_edge (graph_t * graph, unsigned short v1, unsigned short v2);
void graph_remove_edge (graph_t * graph, unsigned short v1, unsigned short v2);

#define set_edge(edge, va, vb)			\
  do {						\
    (edge).v1 = va;				\
    (edge).v2 = vb;				\
  } while (0)

void graph_find_bridges (graph_t * graph);
void graph_find_adjacents (const graph_t * graph, graph_t * subgraph);

#endif /* __GRAPH_H__ */
