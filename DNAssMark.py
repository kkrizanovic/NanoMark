

import os
import sys
import re


def main():
    if len(sys.argv) != 3:
        print 'usage: ./DNAssMark.py --option file'
        print 'options:'
        print 'o1 - option 1'
        sys.exit(1)

    if sys.argv[1] == '--o1':
        f = open(sys.argv[2], 'rU')
        f.close()

    else:
        print 'Ivalid option!'



# Standard template for the main() function
if __name__ == '__main__':
  main()
