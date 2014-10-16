# Copyright (C) 2014 William M. Jacobs

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

cdef extern from "random.h":
    void random_init ()

cdef extern from "graph.h":
    ctypedef struct edge_t:
        unsigned short v1, v2

    ctypedef struct graph_t:
        size_t max_nvertices, nvertices
        size_t max_nedges, nedges
        size_t * nvedges
        unsigned short ** edges
        size_t nbridges
        edge_t * bridges
        size_t nremovable
        edge_t * removable
        size_t nadjacents
        edge_t * adjacents

    graph_t * graph_alloc (size_t nv)
    void graph_free (graph_t * graph)
    void graph_copy (const graph_t * orig, graph_t * dest)

    void graph_add_edge (graph_t * graph, unsigned short v1, unsigned short v2)
    void graph_remove_edge (graph_t * graph, unsigned short v1, unsigned short v2)
    void set_edge (edge, va, vb) # macro

    void graph_find_bridges (graph_t * graph)
    void graph_find_adjacents (const graph_t * graph, graph_t * subgraph)

cdef extern from "sample.h":
    ctypedef struct histogramEV_t:
        pass

    histogramEV_t * histogramEV_alloc (const graph_t * graph)
    void histogramEV_free (histogramEV_t * h)
    void histogramEV_fill (histogramEV_t * h, double z)

    double histogramEV_get(histogramEV_t * h, size_t E, size_t V) # macro
    void histogramEV_set(histogramEV_t * h, size_t E, size_t V, double z) # macro
    void histogramEV_inc(histogramEV_t * h, size_t E, size_t V, double z) # macro

    void wl_sample_subgraphs_EV (const graph_t * graph, graph_t * subgraph, histogramEV_t * lnhEV, \
                                 double f, double flatness, unsigned int ncheck)

    void mc_sample_subgraphs_EV (const graph_t * graph, graph_t * subgraph, const histogramEV_t * lnhEV, unsigned int nsteps)

cdef class Graph:
    cdef graph_t * _graph

