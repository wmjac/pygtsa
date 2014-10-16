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
 
void graph_find_adjacents (const graph_t * graph, graph_t * subgraph)
{
  subgraph->nadjacents = 0;

  /* printf ("--> adj call\n"); */

  size_t i, j, k;
  for (i = 0; i < subgraph->max_nvertices; i++)
    {
      if (subgraph->nvedges[i] > 0 && subgraph->nvedges[i] != graph->nvedges[i])
	{
	  for (j = 0, k = 0; j < graph->nvedges[i]; j++) // && graph->edges[i][j] < i
	    {
	      /* if (k >= subgraph->nvedges[i]) */
	      /* 	{ */
	      /* 	  printf ("adj (%zu -- %u) [%zu]\n", i, graph->edges[i][j], k); */
	      /* 	} */
	      /* else if (subgraph->edges[i][k] > graph->edges[i][j]) */
	      /* 	{ */
	      /* 	  printf ("adj (%zu -- %u / %u) [%zu]\n", i, graph->edges[i][j], subgraph->edges[i][k], k); */
	      /* 	} */
	      if (k >= subgraph->nvedges[i] || subgraph->edges[i][k] > graph->edges[i][j])
		{
		  /* printf ("adj yes\n"); */
		  if (graph->edges[i][j] > i)
		    {
		      set_edge (subgraph->adjacents[subgraph->nadjacents], i, graph->edges[i][j]);
		      subgraph->nadjacents++;
		    }
		  else if (subgraph->nvedges[graph->edges[i][j]] == 0)
		    {
		      set_edge (subgraph->adjacents[subgraph->nadjacents], graph->edges[i][j], i);
		      subgraph->nadjacents++;
		    }
		}
	      else
		{
		  k++;
		}
	    }
	}
    }
}
