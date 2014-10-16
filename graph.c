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

graph_t * graph_alloc (size_t nv)
{
  graph_t * graph = (graph_t *) malloc (sizeof (graph_t));
  if (graph == NULL) goto fail;
  graph->max_nvertices = nv;
  graph->max_nedges = nv * nv;
  graph->nvertices = 0;
  graph->nedges = 0;
  graph->nvedges = (size_t *) malloc (sizeof (size_t) * nv);
  if (graph->nvedges == NULL) goto fail;
  graph->edges = (unsigned short **) malloc (sizeof (unsigned short *) * nv);
  if (graph->edges == NULL) goto fail;
  graph->edges[0] = (unsigned short *) malloc (sizeof (unsigned short *) * nv * nv);
  if (graph->edges[0] == NULL) goto fail;
  size_t i, j;
  for (i = 0; i < nv; i++)
    {
      graph->nvedges[i] = 0;
      if (i > 0)
	{
	  graph->edges[i] = graph->edges[0] + i * nv;
	}
      for (j = 0; j < nv; j++)
	{
	  graph->edges[i][j] = 0;
	}
    }
  graph->nbridges = 0;
  graph->bridges = (edge_t *) malloc (sizeof (edge_t) * nv * nv);
  if (graph->bridges == NULL) goto fail;
  graph->nremovable = 0;
  graph->removable = (edge_t *) malloc (sizeof (edge_t) * nv * nv);
  if (graph->removable == NULL) goto fail;
  graph->nadjacents = 0;
  graph->adjacents = (edge_t *) malloc (sizeof (edge_t) * nv * nv);
  if (graph->adjacents == NULL) goto fail;

  graph->_dfs_iter = 0;
  graph->_dfs_visited = (char *) malloc (sizeof (char) * nv);
  if (graph->_dfs_visited == NULL) goto fail;
  graph->_dfs_discovery = (unsigned short *) malloc (sizeof (unsigned short) * nv);
  if (graph->_dfs_discovery == NULL) goto fail;
  graph->_dfs_low = (unsigned short *) malloc (sizeof (unsigned short) * nv);
  if (graph->_dfs_low == NULL) goto fail;
  graph->_dfs_parent = (unsigned short *) malloc (sizeof (unsigned short) * nv);
  if (graph->_dfs_parent == NULL) goto fail;

  return graph;

 fail:
  graph_free (graph);
  return NULL;
}

void graph_free (graph_t * graph)
{
  if (graph)
    {
      if (graph->nvedges) free (graph->nvedges);
      if (graph->edges)
	{
	  if (graph->edges[0]) free (graph->edges[0]);
	  free (graph->edges);
	}
      if (graph->bridges) free (graph->bridges);
      if (graph->removable) free (graph->removable);
      if (graph->adjacents) free (graph->adjacents);
      if (graph->_dfs_visited) free (graph->_dfs_visited);
      if (graph->_dfs_discovery) free (graph->_dfs_discovery);
      if (graph->_dfs_low) free (graph->_dfs_low);
      if (graph->_dfs_parent) free (graph->_dfs_parent);
      free (graph);
    }
}

void graph_copy (const graph_t * orig, graph_t * dest)
{
  size_t i, j;

  dest->nvertices = 0;
  dest->nedges = 0;

  for (i = 0; i < orig->max_nvertices; i++)
    {
      dest->nvedges[i] = 0;
    }

  for (i = 0; i < orig->max_nvertices; i++)
    {
      for (j = 0; j < orig->nvedges[i]; j++)
	{
	  if (orig->edges[i][j] > i)
	    {
	      graph_add_edge (dest, i, orig->edges[i][j]);
	    }
	}
    }

}

void graph_add_edge (graph_t * graph, unsigned short v1, unsigned short v2)
{
  size_t i, j;

  for (i = 0; graph->edges[v1][i] < v2 && i < graph->nvedges[v1]; i++);
  for (j = graph->nvedges[v1]; j > i; j--)
    {
      graph->edges[v1][j] = graph->edges[v1][j-1];
    }
  graph->edges[v1][j] = v2;
  graph->nvedges[v1]++;

  for (i = 0; graph->edges[v2][i] < v1 && i < graph->nvedges[v2]; i++);
  for (j = graph->nvedges[v2]; j > i; j--)
    {
      graph->edges[v2][j] = graph->edges[v2][j-1];
    }
  graph->edges[v2][j] = v1;
  graph->nvedges[v2]++;

  graph->nedges++;
  if (graph->nvedges[v1] == 1) graph->nvertices++;
  if (graph->nvedges[v2] == 1) graph->nvertices++;
}

void graph_remove_edge (graph_t * graph, unsigned short v1, unsigned short v2)
{
  size_t i, j;

  for (i = 0; graph->edges[v1][i] != v2 && i < graph->nvedges[v1]; i++);
  for (j = i + 1; j < graph->nvedges[v1]; j++)
    {
      graph->edges[v1][j-1] = graph->edges[v1][j];
    }
  graph->nvedges[v1]--;

  for (i = 0; graph->edges[v2][i] != v1 && i < graph->nvedges[v2]; i++);
  for (j = i + 1; j < graph->nvedges[v2]; j++)
    {
      graph->edges[v2][j-1] = graph->edges[v2][j];
    }
  graph->nvedges[v2]--;

  graph->nedges--;
  if (graph->nvedges[v1] == 0) graph->nvertices--;
  if (graph->nvedges[v2] == 0) graph->nvertices--;
}
