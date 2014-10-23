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

def free_energy_profile_V(lnhz):
    F = np.zeros(lnhz.shape[0] - lnhz.shape[1] + 3)
    Vmax = 0
    for i in range(lnhz.shape[0]):
        for j in range(lnhz.shape[1]):
            if np.isfinite(lnhz[i,j]):
                E = i
                V = E - j + 1
                if V > Vmax:
                    Vmax = V
                F[V] += math.exp(lnhz[i,j])
    F[1:Vmax+1] = -np.log(F[1:Vmax+1])
    return F

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('lnhEV', type=str, help="path to lnhEV file")
    parser.add_argument('dihedrals', type=str, help="path to dihedrals file")
    parser.add_argument('rho', type=float, help="dimensionless number density")
    parser.add_argument('--output', type=str, default='fe_profile.dat', help="path to output file [fe_profile.dat]")
    subparsers = parser.add_subparsers(dest='variance')
    parser_zero = subparsers.add_parser('zero-var')
    parser_const = subparsers.add_parser('const-var')
    for p in (parser_zero, parser_const):
        p.add_argument('E', type=float, help="dimensionless energy")
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
    lnhz = lnhzEV(lnhEV, dihedrals, mu, clargs.E, epsilon_expmean, qcoord=dihedrals_meta['qcoord'])
    profile = free_energy_profile_V(lnhz)

    print("Writing output to %s" % clargs.output)
    with open(clargs.output, 'w') as f:
        f.write("# epsilon = %g\n" % clargs.E)
        f.write("# V F(V)\n")
        for i in range(1, len(profile)):
            f.write("%g %g\n" % (i, profile[i]))
