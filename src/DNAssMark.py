

import os
import sys
import re


def verbose_usage_and_exit():
	sys.stderr.write('DNAssMark - a DNA assembly benchmarking tool.\n')
	sys.stderr.write('\n')
	sys.stderr.write('Usage:\n')
	sys.stderr.write('\t%s [mode]\n' % sys.argv[0])
	sys.stderr.write('\n')
	sys.stderr.write('\tmode:\n')
	sys.stderr.write('\t\tsetup\n')
	sys.stderr.write('\t\tbenchmark\n')
	sys.stderr.write('\n')
	exit(0)

def main():
    if len(sys.argv) != 3:
        print 'usage: ./DNAssMark.py --option file'
        print 'options:'
        print 'o1 - option 1'
        sys.exit(1)

    if sys.argv[1] == '--setup':
        f = open(sys.argv[2], 'rU')
        f.close()

    elif sys.argv[1] == '--benchmark':
        f = open(sys.argv[2], 'rU')
        f.close()

    else:
        print 'Ivalid option!'



# Standard template for the main() function
if __name__ == '__main__':
  main()
