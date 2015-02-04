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

#include <math.h>
#include <float.h>
#include "sample.h"
#include "random.h"

#define INIT_ADJACENTS_REMOVABLE_SWAP					\
  size_t nremovable_old, nadjacents_old;				\
  size_t nremovable_tmp, nadjacents_tmp;				\
  edge_t * removable_old = (edge_t *)					\
    malloc (sizeof (edge_t) * graph->max_nedges);			\
  if (removable_old == NULL) goto fail;					\
  edge_t * removable_tmp;						\
  edge_t * adjacents_old = (edge_t *)					\
    malloc (sizeof (edge_t) * graph->max_nedges);			\
  if (adjacents_old == NULL) goto fail;					\
  edge_t * adjacents_tmp;

#define swap_adjacents_removable()		\
  do {						\
    nadjacents_tmp = nadjacents_old;		\
    nadjacents_old = subgraph->nadjacents;	\
    subgraph->nadjacents = nadjacents_tmp;	\
    adjacents_tmp = adjacents_old;		\
    adjacents_old = subgraph->adjacents;	\
    subgraph->adjacents = adjacents_tmp;	\
    nremovable_tmp = nremovable_old;		\
    nremovable_old = subgraph->nremovable;	\
    subgraph->nremovable = nremovable_tmp;	\
    removable_tmp = removable_old;		\
    removable_old = subgraph->removable;	\
    subgraph->removable = removable_tmp;	\
  } while (0)

void wl_sample_subgraphs_EV (const graph_t * graph, graph_t * subgraph, histogramEV_t * lnhEV, double f, double flatness, unsigned long ncheck)
{
  const size_t max_x = lnhEV->max_x;
  const size_t max_y = lnhEV->max_y;

  unsigned long * visits = (unsigned long *) malloc (sizeof (unsigned long) * max_x * max_y);
  if (visits == NULL) goto fail;
  size_t i;
  for (i = 0; i < max_x * max_y; i++)
    {
      visits[i] = 0;
    }

  INIT_ADJACENTS_REMOVABLE_SWAP;
  graph_find_bridges (subgraph);
  graph_find_adjacents (graph, subgraph);

  edge_t edge;
  size_t edge_index, nedges_old, nvertices_old;
  double arg;

  unsigned long step;
  do {
    for (step = 0; step < ncheck; step++)
      {
	nedges_old = subgraph->nedges;
	nvertices_old = subgraph->nvertices;
	if (random_uniform () < 0.5)
	  { /* Attempt to add an edge. */
	    if (subgraph->nadjacents > 0)
	      {
		edge_index = random_uniform () * subgraph->nadjacents;
		set_edge (edge, subgraph->adjacents[edge_index].v1, subgraph->adjacents[edge_index].v2);
		graph_add_edge (subgraph, edge.v1, edge.v2);
		swap_adjacents_removable ();
		graph_find_bridges (subgraph);
		graph_find_adjacents (graph, subgraph);
		arg = log ((double) nadjacents_old / (double) subgraph->nremovable)
		  + histogramEV_get (lnhEV, nedges_old, nvertices_old)
		  - histogramEV_get (lnhEV, subgraph->nedges, subgraph->nvertices);
		if (arg < 0. && random_uniform () >= exp (arg))
		  { /* Reject move. */
		    /* printf ("reject add: %zu,%zu -> %zu,%zu [%g -> %g] %g\n", nedges_old, nvertices_old, subgraph->nedges, subgraph->nvertices, \ */
		    /* 	    histogramEV_get (lnhEV, nedges_old, nvertices_old), histogramEV_get (lnhEV, subgraph->nedges, subgraph->nvertices), arg); */
		    graph_remove_edge (subgraph, edge.v1, edge.v2);
		    swap_adjacents_removable ();
		  }
		else
		  {
		    /* printf ("accept add: %zu,%zu -> %zu,%zu [%g -> %g] %g\n", nedges_old, nvertices_old, subgraph->nedges, subgraph->nvertices, \ */
		    /* 	    histogramEV_get (lnhEV, nedges_old, nvertices_old), histogramEV_get (lnhEV, subgraph->nedges, subgraph->nvertices), arg); */
		  }
	      }
	  }
	else
	  { /* Attempt to remove an edge. */
	    if (subgraph->nremovable > 0 && subgraph->nedges > 1)
	      {
		edge_index = random_uniform () * subgraph->nremovable;
		set_edge (edge, subgraph->removable[edge_index].v1, subgraph->removable[edge_index].v2);
		graph_remove_edge (subgraph, edge.v1, edge.v2);
		swap_adjacents_removable ();
		graph_find_bridges (subgraph);
		graph_find_adjacents (graph, subgraph);
		arg = log ((double) nremovable_old / (double) subgraph->nadjacents)
		  + histogramEV_get (lnhEV, nedges_old, nvertices_old)
		  - histogramEV_get (lnhEV, subgraph->nedges, subgraph->nvertices);
		if (arg < 0. && random_uniform () >= exp (arg))
		  { /* Reject move. */
		    /* printf ("reject rem: %zu,%zu -> %zu,%zu [%g -> %g] %g\n", nedges_old, nvertices_old, subgraph->nedges, subgraph->nvertices, \ */
		    /* 	    histogramEV_get (lnhEV, nedges_old, nvertices_old), histogramEV_get (lnhEV, subgraph->nedges, subgraph->nvertices), arg); */
		    graph_add_edge (subgraph, edge.v1, edge.v2);
		    swap_adjacents_removable ();
		  }
		else
		  {
		    /* printf ("accept rem: %zu,%zu -> %zu,%zu [%g -> %g] %g\n", nedges_old, nvertices_old, subgraph->nedges, subgraph->nvertices, \ */
		    /* 	    histogramEV_get (lnhEV, nedges_old, nvertices_old), histogramEV_get (lnhEV, subgraph->nedges, subgraph->nvertices), arg); */
		  }
	      }
	  }
	visits[(subgraph->nedges + 1) - subgraph->nvertices + max_y * subgraph->nedges]++;
	histogramEV_inc (lnhEV, subgraph->nedges, subgraph->nvertices, f);
      }

    /* Check flatness condition. */
    unsigned long n = 0;
    unsigned long visits_sum = 0;
    unsigned long visits_max = 0;
    for (i = 0; i < max_x * max_y; i++)
      {
	if (visits[i] > 0)
	  {
	    n++;
	    visits_sum += visits[i];
	    visits_max = (visits[i] > visits_max) ? visits[i] : visits_max;
	  }
      }
    if ((double) visits_sum / (double) n > flatness * visits_max) break;
  } while (1);

  /* Tare lnhEV. */
  double lnhEV_min = DBL_MAX;
  for (i = 0; i < max_x * max_y; i++)
    {
      lnhEV_min = (lnhEV->data[i] > 0 && lnhEV->data[i] < lnhEV_min) ? lnhEV->data[i] : lnhEV_min;
    }
  for (i = 0; i < max_x * max_y; i++)
    {
      if (lnhEV->data[i] >= lnhEV_min)
	{
	  lnhEV->data[i] -= lnhEV_min;
	}
      else
	{
	  lnhEV->data[i] = -1.;
	}
    }

  graph_find_bridges (subgraph);

 fail:
  if (removable_old) free (removable_old);
  if (adjacents_old) free (adjacents_old);
  if (visits) free (visits);
}

void wl_simulate_subgraphs_EV (const graph_t * graph, graph_t * subgraph, histogramEV_t * lnhEV, double f_start, double f_target, double flatness, unsigned long ncheck)
{
  double f = f_start;
  while (f >= f_target)
    {
      printf ("WL f = %g\n", f);
      wl_sample_subgraphs_EV (graph, subgraph, lnhEV, f, flatness, ncheck);
      f *= 0.5;
    }
}

void mc_sample_subgraphs_EV (const graph_t * graph, graph_t * subgraph, const histogramEV_t * lnhEV, unsigned long nsteps)
{
  INIT_ADJACENTS_REMOVABLE_SWAP;
  graph_find_bridges (subgraph);
  graph_find_adjacents (graph, subgraph);

  edge_t edge;
  size_t edge_index, nedges_old, nvertices_old;
  double arg;

  unsigned long step;
  for (step = 0; step < nsteps; step++)
    {
      nedges_old = subgraph->nedges;
      nvertices_old = subgraph->nvertices;
      if (random_uniform () < 0.5)
	{ /* Attempt to add an edge. */
	  if (subgraph->nadjacents > 0)
	    {
	      edge_index = random_uniform () * subgraph->nadjacents;
	      set_edge (edge, subgraph->adjacents[edge_index].v1, subgraph->adjacents[edge_index].v2);
	      graph_add_edge (subgraph, edge.v1, edge.v2);
	      swap_adjacents_removable ();
	      graph_find_bridges (subgraph);
	      graph_find_adjacents (graph, subgraph);
	      arg = log ((double) nadjacents_old / (double) subgraph->nremovable)
		+ histogramEV_get (lnhEV, nedges_old, nvertices_old)
		- histogramEV_get (lnhEV, subgraph->nedges, subgraph->nvertices);
	      if (arg < 0. && random_uniform () >= exp (arg))
		{ /* Reject move. */
		  /* printf ("reject add: %zu,%zu -> %zu,%zu [%g -> %g] %g\n", nedges_old, nvertices_old, subgraph->nedges, subgraph->nvertices, \ */
		  /* 	    histogramEV_get (lnhEV, nedges_old, nvertices_old), histogramEV_get (lnhEV, subgraph->nedges, subgraph->nvertices), arg); */
		  graph_remove_edge (subgraph, edge.v1, edge.v2);
		  swap_adjacents_removable ();
		}
	      else
		{
		  /* printf ("accept add: %zu,%zu -> %zu,%zu [%g -> %g] %g\n", nedges_old, nvertices_old, subgraph->nedges, subgraph->nvertices, \ */
		  /* 	    histogramEV_get (lnhEV, nedges_old, nvertices_old), histogramEV_get (lnhEV, subgraph->nedges, subgraph->nvertices), arg); */
		}
	    }
	}
      else
	{ /* Attempt to remove an edge. */
	  if (subgraph->nremovable > 0 && subgraph->nedges > 1)
	    {
	      edge_index = random_uniform () * subgraph->nremovable;
	      set_edge (edge, subgraph->removable[edge_index].v1, subgraph->removable[edge_index].v2);
	      graph_remove_edge (subgraph, edge.v1, edge.v2);
	      swap_adjacents_removable ();
	      graph_find_bridges (subgraph);
	      graph_find_adjacents (graph, subgraph);
	      arg = log ((double) nremovable_old / (double) subgraph->nadjacents)
		+ histogramEV_get (lnhEV, nedges_old, nvertices_old)
		- histogramEV_get (lnhEV, subgraph->nedges, subgraph->nvertices);
	      if (arg < 0. && random_uniform () >= exp (arg))
		{ /* Reject move. */
		  /* printf ("reject rem: %zu,%zu -> %zu,%zu [%g -> %g] %g\n", nedges_old, nvertices_old, subgraph->nedges, subgraph->nvertices, \ */
		  /* 	    histogramEV_get (lnhEV, nedges_old, nvertices_old), histogramEV_get (lnhEV, subgraph->nedges, subgraph->nvertices), arg); */
		  graph_add_edge (subgraph, edge.v1, edge.v2);
		  swap_adjacents_removable ();
		}
	      else
		{
		  /* printf ("accept rem: %zu,%zu -> %zu,%zu [%g -> %g] %g\n", nedges_old, nvertices_old, subgraph->nedges, subgraph->nvertices, \ */
		  /* 	    histogramEV_get (lnhEV, nedges_old, nvertices_old), histogramEV_get (lnhEV, subgraph->nedges, subgraph->nvertices), arg); */
		}
	    }
	}
    }

  graph_find_bridges (subgraph);

 fail:  
  if (removable_old) free (removable_old);
  if (adjacents_old) free (adjacents_old);
}
