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

#include "config.h"

#include "graph.h"

static void find_bridges_within (graph_t * graph, unsigned short u)
{

  graph->_dfs_visited[u] = 1;
  graph->_dfs_discovery[u] = graph->_dfs_low[u] = ++(graph->_dfs_iter);

  size_t i;
  for (i = 0; i < graph->nvedges[u]; i++)
    {
      unsigned short v = graph->edges[u][i];

      //      printf ("touch %zu %zu, parent %zu\n", u, v, graph->_dfs_parent[u]);

      if (!graph->_dfs_visited[v])
        {
	  graph->_dfs_parent[v] = u;
	  find_bridges_within (graph, v);

	  graph->_dfs_low[u] = (graph->_dfs_low[u] < graph->_dfs_low[v]) ? graph->_dfs_low[u] : graph->_dfs_low[v];

	  if (graph->_dfs_low[v] > graph->_dfs_discovery[u])
	    {
	      set_edge (graph->bridges[graph->nbridges], u, v);
	      graph->nbridges++;
	      //	      printf ("found bridge: %zu -- %zu; nvedges %zu %zu\n", u, v, graph->nvedges[u], graph->nvedges[v]);
	      if (graph->nvedges[u] > 1 && graph->nvedges[v] > 1)
		{
		  continue;
		}
	    }
        }
      else if (v != graph->_dfs_parent[u])
	{
	  graph->_dfs_low[u] = (graph->_dfs_low[u] < graph->_dfs_discovery[v]) ? graph->_dfs_low[u] : graph->_dfs_discovery[v];
	}

      if (graph->_dfs_discovery[u] < graph->_dfs_discovery[v])
	{
	  //	  printf ("removable: %zu -- %zu\n", u, v);
	  if (u < v)
	    set_edge (graph->removable[graph->nremovable], u, v);
	  else
	    set_edge (graph->removable[graph->nremovable], v, u);
	  graph->nremovable++;
	}
    }
}
 
void graph_find_bridges (graph_t * graph)
{
  graph->nbridges = 0;
  graph->nremovable = 0;
  graph->_dfs_iter = 0;

  //  printf ("bridges calculation...\n");
  //  print_graph (graph);

  size_t i;
  for (i = 0; i < graph->max_nvertices; i++)
    {
      graph->_dfs_visited[i] = 0;
      graph->_dfs_parent[i] = -1;
    }
 
  /* Assume only one connected component. */

  for (i = 0; i < graph->max_nvertices; i++)
    {
      if (graph->nvedges[i] > 0)
	{
	  find_bridges_within (graph, i);
	  return;
	}
    }
}
