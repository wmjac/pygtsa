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

#ifndef __SAMPLE_H__
#define __SAMPLE_H__

#include <stdlib.h>
#include "graph.h"

typedef struct
{
  size_t max_x, max_y;
  double * data;
}
histogramEV_t;

histogramEV_t * histogramEV_alloc (const graph_t * graph);
void histogramEV_free (histogramEV_t * h);
void histogramEV_fill (histogramEV_t * h, double z);

#define histogramEV_get(h, E, V) ((h)->data[((E) + 1) - (V) + (h)->max_y * (E)])
#define histogramEV_set(h, E, V, z) do {(h)->data[((E) + 1) - (V) + (h)->max_y * (E)] = z;} while (0)
#define histogramEV_inc(h, E, V, z) do {(h)->data[((E) + 1) - (V) + (h)->max_y * (E)] += z;} while (0)

void wl_simulate_subgraphs_EV (const graph_t * graph, graph_t * subgraph, histogramEV_t * lnhEV, double f_start, double f_target, double flatness, unsigned int ncheck);
void wl_sample_subgraphs_EV (const graph_t * graph, graph_t * subgraph, histogramEV_t * lnhEV, double f, double flatness, unsigned int ncheck);

void mc_sample_subgraphs_EV (const graph_t * graph, graph_t * subgraph, const histogramEV_t * lnhEV, unsigned int nsteps);

#endif /* __SAMPLE_H__ */
