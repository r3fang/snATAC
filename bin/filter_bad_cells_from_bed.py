#!/usr/bin/env python

import sys
import collections
import os


def main():
    barcode_list = [ line.strip().split()[0] for line in open(sys.argv[1]).readlines() ]
    for line in sys.stdin:
        barcode = line.split()[3].split(':')[0]
        if not barcode in barcode_list:  # if not perfectly matched
                continue    
        try:
            print line,
        except IOError:
            try:
                sys.stdout.close()
            except IOError:
                pass
            try:
                sys.stderr.close()
            except IOError:
                pass
                 
if __name__ == '__main__':
    main()
