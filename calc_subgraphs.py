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

import argparse
import networkx as nx
from pygtsa.structure import Assembly
from pygtsa.cgraph import wl_simulate_EV, mc_simulate_dihedrals_adjacents_EV

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('structure', type=str, help="path to input structure file")
    parser.add_argument('--qcoord', metavar='Q', type=float, default=4., help="monomer rotation constant [4.]")
    parser.add_argument('--qdih', metavar='Q', type=float, default=3., help="dihedral angle connective constant [3.]")
    parser.add_argument('--output-prefix', metavar='PATH', type=str, default='./', help="path to output files [./]")
    clargs = parser.parse_args()

    # Initialize

    target = Assembly.read(clargs.structure)

    if clargs.output_prefix != '' and clargs.output_prefix[-1] != '/':
        clargs.output_prefix = clargs.output_prefix + '_'

    # Sanity checks

    if not nx.is_connected(target):
        print("ERROR: target assembly is disconnected")
        raise SystemExit
    if len(target.leaves()) != 0:
        print("WARNING: target assembly contains at least one leaf")
    if len(target.bridges()) != 0:
        print("WARNING: target assembly contains at least one bridge")
    if target.number_of_edges() >= 200:
        print("WARNING: target assemblies with more than 200 edges may require an "
              "unreasonable amount of computing time (with the current implementation).")

    # Run subgraphs calculation

    print("Calculating subgraph density of states...")
    print("Writing density of states to", clargs.output_prefix + 'lnhEV.dat')
    lnhEV = wl_simulate_EV(target, output=clargs.output_prefix + 'lnhEV.dat')

    # Calculate dihedrals and adjacents

    print("Calculating dihedrals and adjacents averages...")
    nsteps = 2 * target.number_of_edges()
    nsamples = 4000 * sum(sum(1 for j in range(lnhEV.h.shape[1]) if lnhEV.h[i,j] >= 0.) for i in range(lnhEV.h.shape[0]))
    print("Collecting %d samples with %d steps between samples." % (nsamples, nsteps))
    dihedrals, adjacents, visits = mc_simulate_dihedrals_adjacents_EV(target, lnhEV, clargs.qdih, nsteps, nsamples)

    print("Writing average dihedral entropy loss to", clargs.output_prefix + 'dihedrals.dat')
    with open(clargs.output_prefix + 'dihedrals.dat', 'w') as f:
        f.write("# qcoord = %g\n" % clargs.qcoord)
        f.write("# qdih = %g\n\n" % clargs.qdih)
        dihedrals.write(f)
    print("Writing average number of adjacent edges to", clargs.output_prefix + 'adjacents.dat')
    with open(clargs.output_prefix + 'adjacents.dat', 'w') as f:
        adjacents.write(f)
    print("Writing sampling histogram to", clargs.output_prefix + 'sampling_visits.dat')
    with open(clargs.output_prefix + 'sampling_visits.dat', 'w') as f:
        visits.write(f)
