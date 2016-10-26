"""
Print the Hill Formula from a ccopen'ed log file
"""
from __future__ import print_function

import sys
from cclib.io import ccopen
from cclib.parser.utils import PeriodicTable


def main():
    pse = PeriodicTable()
    formula = {}
    log = ccopen(sys.argv[1])
    data = log.parse()
    for atom in data.atomnos:
        formula[atom] = formula.setdefault(atom, 0) + 1
    string = []
    exclude = set()
    if 6 in formula:
        string.append('C%d' % formula[6])
        exclude.add(6)
        if 1 in formula:
            string.append('H%d' % formula[1])
            exclude.add(1)
    elements = set(formula.keys()) - exclude
    for element in sorted([pse.element[atom] for atom in elements]):
        string.append(element + str(formula[pse.number[element]]))
    print(' '.join(string))


if __name__ == '__main__':
    main()
