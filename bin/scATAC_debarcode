#!/usr/bin/env python
import gzip
import bz2
import sys
import collections
import os
import operator
import os.path
import optparse

def main():
    """ main function """
    parser = optparse.OptionParser(usage='%prog [-h] [-a I1.fastq] [-b I2.fastq] [-c R1.fastq]',
                                   description='Decomplex single-cell ATAC-seq barcode allowing mismatch.')
    parser.add_option('-a',
                      dest="I1",
                      help='I1.fastq.gz'
                      )

    parser.add_option('-b',
                      dest='I2',
                      help='I2.fastq.gz'
                      )

    parser.add_option('-c',
                      dest='R1',
                      help='R1.fastq.gz'
                      )

    parser.add_option('--version',
                      dest="version",
                      default=1.0,
                      type="float",
                      )

    if len(sys.argv) < 6:
        parser.print_help()
        exit('error: too few arguments')

    args = parser.parse_args()[0]

    fi1_name = args.I1
    fi2_name = args.I2
    fr1_name = args.R1
    
    if not os.path.isfile(fi1_name): exit("error: \'%s\' not exist" % fi1_name)
    if not os.path.isfile(fi2_name): exit("error: \'%s\' not exist" % fi2_name)
    if not os.path.isfile(fr1_name): exit("error: \'%s\' not exist" % fr1_name)
    
    # check file format
    if fi1_name.endswith('.gz'):
        fi1 = gzip.open(fi1_name, 'rb')
    elif fi1_name.endswith('.bz2'):
        fi1 = bz2.BZ2File(fi1_name, 'r')
    elif fi1_name.endswith('.fastq'):
        fi1 = open(fi1_name, 'r')
    elif fi1_name.endswith('.fq'):
        fi1 = open(fi1_name, 'r')

    if fi2_name.endswith('.gz'):
        fi2 = gzip.open(fi2_name, 'rb')
    elif fi2_name.endswith('.bz2'):
        fi2 = bz2.BZ2File(fi2_name, 'r')
    elif fi2_name.endswith('.fastq'):
        fi2 = open(fi2_name, 'r')
    elif fi2_name.endswith('.fq'):
        fi2 = open(fi2_name, 'r')

    if fr1_name.endswith('.gz'):
        fr1 = gzip.open(fr1_name, 'rb')
    elif fr1_name.endswith('.bz2'):
        fr1 = bz2.BZ2File(fr1_name, 'r')
    elif fr1_name.endswith('.fastq'):
        fr1 = open(fr1_name, 'r')
    elif fr1_name.endswith('.fq'):
        fr1 = open(fr1_name, 'r')
    
    while True:
        cur_i1_name = fi1.readline().strip()[1:]
        cur_i1_read = fi1.readline().strip()
        cur_i1_plus = fi1.readline().strip()
        cur_i1_qual = fi1.readline().strip()
    
        cur_i2_name = fi2.readline().strip()[1:]
        cur_i2_read = fi2.readline().strip()
        cur_i2_plus = fi2.readline().strip()
        cur_i2_qual = fi2.readline().strip()
    
        cur_r1_name = fr1.readline().strip()[1:]
        cur_r1_read = fr1.readline().strip()
        cur_r1_plus = fr1.readline().strip()
        cur_r1_qual = fr1.readline().strip()
            
        if cur_i1_name == "" or cur_i2_name == "" or cur_r1_name == "": break        
        if not (cur_i1_name.split()[0] == cur_i2_name.split()[0] == cur_r1_name.split()[0]): sys.exit("error(main): read name not matched")        

        cur_r7 = cur_i1_read[:8]
        cur_i7 = cur_i1_read[-8:]
        cur_i5 = cur_i2_read[:8]
        cur_r5 = cur_i2_read[-8:]        
        
        cur_barcode = cur_r7 + cur_i7 + cur_i5 + cur_r5
        if cur_barcode.count('N') >= 12: continue

        try:
            print '@' + cur_barcode + ':' + cur_r1_name
            print cur_r1_read
            print '+'
            print cur_r1_qual
        except IOError:
            try:
                sys.stdout.close()
            except IOError:
                pass
            try:
                sys.stderr.close()
            except IOError:
                pass
                
    fi1.close()
    fi2.close()
    fr1.close()

if __name__ == '__main__':
    main()
