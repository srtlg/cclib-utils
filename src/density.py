"""
Calculate the particle density for a given molecule
"""
from __future__ import print_function

import argparse

from cclib.io import ccopen
from cclib.method.volume import Volume, electrondensity
from print_extent import MolecularExtent


def _parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('quantum_chemical_file')
    p.add_argument('-r', '--radius', type=float, default=1.5)
    p.add_argument('-s', '--spacing', type=float, default=0.25)
    p.add_argument('-i', '--index', type=int, default=0)
    p.add_argument('-o', '--output')
    p.add_argument('-S', '--spinidx', type=int, default=0, help='spin index 0=alpha, 1=beta')
    return p.parse_args()


def main():
    args = _parse_args()
    log_file = ccopen(args.quantum_chemical_file)
    data = log_file.parse()
    me = MolecularExtent(data, args.radius)
    me.calculate_extent(args.index)
    volume = Volume(me.lower_limit, me.upper_limit, spacing=(args.spacing, args.spacing, args.spacing))
    slice_end = data.homos[args.spinidx] + 1
    print('density of molecular orbitals [0..%d] HOMO=%d' % (slice_end, data.homos[args.spinidx]))
    density = electrondensity(data.atomcoords[args.index],
                              data.mocoeffs[args.spinidx][0:slice_end],
                              data.gbasis,
                              volume)
    print('integral', density.integrate(), density.integrate_square())
    if args.output is not None:
        density.writeascube(args.output + '.cube')


if __name__ == '__main__':
    main()
