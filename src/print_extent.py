"""
Calculate the extent of the molecule for a proper volume
"""
import argparse

import numpy as np
from cclib.io import ccopen
from cclib.parser.utils import PeriodicTable


class MolecularExtent(object):
    def __init__(self, data, radius):
        self.data = data
        self.radius = radius
        self.upper_limit = None
        self.lower_limit = None

    def calculate_extent(self, index):
        coords = self.data.atomcoords[index]
        distance = np.ones([3]) * self.radius
        self.upper_limit = np.max(coords + distance, axis=0)
        self.lower_limit = np.min(coords - distance, axis=0)

    def as_jmol(self, index, output='extent'):
        pse = PeriodicTable()
        coords = self.data.atomcoords[index]
        with open(output + '.xyz', 'w') as fout:
            fout.write('%d\n\n' % len(coords))
            for atomno, coord in zip(self.data.atomnos, coords):
                fout.write('%s %f %f %f\n' % (tuple(pse.element[atomno]) + tuple(coord)))
        with open(output + '.spt', 'w') as fout:
            fout.write('load "%s.xyz"\n' % output)
            fout.write('draw POLYGON ')
            fout.write('8 ')
            origin = self.lower_limit
            diagonal = self.upper_limit - self.lower_limit
            vertices = []
            for i in ((0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1)):
                vertex = origin + np.asarray(i) * diagonal
                vertices.append(vertex)
                fout.write('{%f %f %f} ' % tuple(vertex))
            fout.write('6 [0 1 2  3] [2 3 0  3] [0 4 5  3] [1 5 6  3] [2 6 7  3] [3 7 4  3] ')
            fout.write('mesh nofill\n')
        with open(output + '_vertices.xyz', 'w') as fout:
            fout.write('%d\n\n' % len(vertices))
            for v in vertices:
                fout.write('C %f %f %f\n' % tuple(v))


def _parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('quantum_chemical_file')
    p.add_argument('-r', '--radius', type=float, default=1.5, help='radius in Angstrom')
    p.add_argument('-i', '--index', type=int, default=0, help='index of molecule')
    return p.parse_args()


def main():
    args = _parse_args()
    log_file = ccopen(args.quantum_chemical_file)
    me = MolecularExtent(log_file.parse(), args.radius)
    me.calculate_extent(args.index)
    me.as_jmol(args.index)


if __name__ == '__main__':
    main()
