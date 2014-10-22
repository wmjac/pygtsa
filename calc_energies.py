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

import argparse, random
from pygtsa.structure import Assembly
from pygtsa.histogram import EVHistogram
from pygtsa.cgraph import mc_simulate_energies_EV

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('structure', type=str, help="path to input structure file")
    parser.add_argument('lnhEV', type=str, help="path to lnhEV file")
    subparsers = parser.add_subparsers(dest='distribution')
    parser_gaussian = subparsers.add_parser('gaussian')
    parser_gaussian.add_argument('mean', metavar='MU', type=float, help="mean of Gaussian distribution")
    parser_gaussian.add_argument('stddev', metavar='SIGMA', type=float, help="standard deviation of Gaussian distribution")
    parser_gaussian.add_argument('--nsamples', type=int, default=100, help="number of independent energy samples [100]")
    parser_gaussian.add_argument('--output-prefix', metavar='PATH', type=str, default='./', help="path to output files [./]")
    clargs = parser.parse_args()

    # Initialize

    target = Assembly.read(clargs.structure)
    with open(clargs.lnhEV, 'r') as f:
        lnhEV = EVHistogram.read(f)

    if clargs.output_prefix != '' and clargs.output_prefix[-1] != '/':
        clargs.output_prefix = clargs.output_prefix + '_'

    # Sample energies

    if clargs.distribution == 'gaussian':
        energy_dicts = [{E : random.normalvariate(clargs.mean, clargs.stddev) for E in target.edges()} for i in range(clargs.nsamples)]
    print("Writing energy samples to", clargs.output_prefix + 'bonds.dat')
    with open(clargs.output_prefix + 'bonds.dat', 'w') as f:
        for E in target.edges():
            f.write("%d -- %d:" % E)
            for i in range(len(energy_dicts)):
                f.write(" %g" % energy_dicts[i][E])
            f.write("\n")

    # Calculate dihedrals and adjacents

    print("Calculating energy averages...")
    nsteps = 2 * target.number_of_edges()
    nsamples = 4000 * sum(sum(1 for j in range(lnhEV.h.shape[1]) if lnhEV.h[i,j] >= 0.) for i in range(lnhEV.h.shape[0]))
    print("Collecting %d samples with %d steps between samples." % (nsamples, nsteps))
    energies, visits = mc_simulate_energies_EV(target, lnhEV, energy_dicts, nsteps, nsamples)

    print("Writing average energies to", clargs.output_prefix + 'energies.dat')
    with open(clargs.output_prefix + 'energies.dat', 'w') as f:
        f.write("# average energies (%s distribution)\n" % clargs.distribution)
        energies.write(f)
    print("Writing sampling histogram to", clargs.output_prefix + 'sampling_visits.dat')
    with open(clargs.output_prefix + 'sampling_visits.dat', 'w') as f:
        visits.write(f)
