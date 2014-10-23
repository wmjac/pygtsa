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

import argparse, math
import numpy as np
from pygtsa.histogram import EVHistogram, read_histogram, read_meta
from pygtsa.calc_yield import load_dihedrals, find_target_indices, lnhzEV
from pygtsa.calc_fe_profile import free_energy_profile_V

def barrier_peak(F):
    Fmax = F[1:].max() - F[1]
    Vmax = F[1:].argmax() + 1
    return Vmax, Fmax

def nucleation_barrier(F):
    minima = [[]]
    for i in range(2, len(F)):
        if F[i] < F[1]:
            minima[-1].append(i)
        elif len(minima[-1]) > 0:
            minima.append([])
    if len(minima[-1]) == 0:
        minima = minima[:-1]
    if len(minima) == 0:
        return barrier_peak(F)
    elif len(minima) == 1 and minima[0][0] > 2:
        return barrier_peak(F[:minima[0][0]])
    elif len(minima) == 1 and F[-1] > F[1]:
        Fb = F[minima[0][-1]:]
        Fmax = Fb.max() - F[1:minima[0][-1]].min()
        Vmax = Fb.argmax() + minima[0][-1]
        return Vmax, Fmax
    elif len(minima) == 1:
        return None, None
    else: # len(minima) > 1:
        Fb = F[minima[0][-1]:minima[1][0] + 1]
        Fmax = Fb.max() - F[1:minima[0][-1]].min()
        Vmax = Fb.argmax() + minima[0][-1]
        return Vmax, Fmax

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('lnhEV', type=str, help="path to lnhEV file")
    parser.add_argument('dihedrals', type=str, help="path to dihedrals file")
    parser.add_argument('rho', type=float, help="dimensionless number density")
    parser.add_argument('--output', type=str, default='barrier.dat', help="path to output file [barrier.dat]")
    subparsers = parser.add_subparsers(dest='variance')
    parser_zero = subparsers.add_parser('zero-var')
    parser_const = subparsers.add_parser('const-var')
    for p in (parser_zero, parser_const):
        p.add_argument('Emin', type=float, help="minimum dimensionless energy")
        p.add_argument('Emax', type=float, help="maximum dimensionless energy")
        p.add_argument('dE', type=float, help="dimensionless energy step")
    parser_const.add_argument('distribution', type=str, help="path to energy distribution file")
    clargs = parser.parse_args()

    if clargs.variance == None:
        parser.print_usage()
        print("error: variance mode not specified")
        raise SystemExit

    with open(clargs.lnhEV, 'r') as f:
        lnhEV = read_histogram(f)
    with open(clargs.dihedrals, 'r') as f:
        dihedrals, dihedrals_meta = load_dihedrals(f)
    target_indices = find_target_indices(dihedrals)

    if clargs.variance == 'const-var':
        with open(clargs.distribution, 'r') as f:
            epsilon_expmean = read_histogram(f)
        epsilon_expmean -= epsilon_expmean[max(target_indices)]
    else:
        epsilon_expmean = np.zeros(dihedrals.shape)

    mu = math.log(clargs.rho)
    epsilon = clargs.Emin
    barrier = {}
    while epsilon < clargs.Emax + clargs.dE:
        lnhz = lnhzEV(lnhEV, dihedrals, mu, epsilon, epsilon_expmean, qcoord=dihedrals_meta['qcoord'])
        profile = free_energy_profile_V(lnhz)
        barrier[epsilon] = nucleation_barrier(profile)
        if barrier[epsilon][0] == None:
            barrier[epsilon] = (0, 0.)
        epsilon += clargs.dE

    print("Writing output to %s" % clargs.output)
    with open(clargs.output, 'w') as f:
        f.write("# epsilon rho nucleation_barrier_V nucleation_barrier_fe\n")
        for epsilon in sorted(barrier.keys()):
            f.write("%g %g %d %g\n" % (epsilon, clargs.rho, barrier[epsilon][0], barrier[epsilon][1]))
