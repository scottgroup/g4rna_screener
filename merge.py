#!/usr/bin/env python

#    Identification of potential RNA G-quadruplexes by G4RNA screener.
#    Copyright (C) 2018  Jean-Michel Garant
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
import utils
import argparse
import pandas as pd

class float_range(object):
    """
    Object that defines a range of float that is authorized or returns the
    authorized range in error. Used to validate score threshold input by user.
    """
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        if self.start == None:
            return other <= self.end
        elif self.end == None:
            return self.start <= other
        else:
            return self.start <= other <= self.end
    def __repr__(self):
        return '[{0}:{1}]'.format(self.start, self.end)

def custom(sub_df):
    return list(sub_df.sequence)
    #return len(''.join(sub_df.sequence))

def merge_g4rna(df, window=60, step=10,
        cGcC=False, G4H=False, G4NN=False):
#    df = df.groupby(lambda x: group_func(df, x-1),
#            as_index=True, sort=False
#            )['start','sequence','cGcC','G4H','G4NN'].apply(
#                    custom).reset_index()
#    df = df.apply(lambda x: annotate(df, x), axis=1)
    if 'start' in df.columns:
        df['start_match'] = (
                df.start.eq(df.start.shift()+step) | 
                df.start.shift(-1).eq(df.start+step))
    if 'sequence' in df.columns:
        df['seq_match'] = (
                df.sequence.str[-50:].eq(df.sequence.str[:50].shift(-1)) | 
                df.sequence.str[-50:].shift().eq(df.sequence.str[:50]))
#    df = df.groupby(['start_match','seq_match'])
    print df
    
    #return df

def arguments():
    """
    Arguments management
    """
    # declare argument parser
    parser = argparse.ArgumentParser(formatter_class=utils.Formatter,
            prog=os.path.basename(__file__),
            description="Merge positive windows of screen.py output and "\
                    "discard windows below the threshold(s)",
            epilog="G4RNA screener  Copyright (C) 2018  Jean-Michel Garant "\
            "This program comes with ABSOLUTELY NO WARRANTY. This is free "\
            "software, and you are welcome to redistribute it under certain "\
            "conditions <http://www.gnu.org/licenses/>.")
    # TSV input from STDIN is supported by default using argument "-"
    parser.add_argument('tsv',
            type=argparse.FileType('r'),
            default=sys.stdin,
            help='TSV file (tab separated values), - for default STDIN',
            metavar='TSV')
    # Use cGcC
    parser.add_argument("--cGcC",
        type=float,
        nargs='?',
        choices=[float_range(0,None)],
        default=False,
        help="Use cGcC score threshold to determine positive windows "\
                "(default: 4.5)",
        metavar="FLOAT")
    # Use G4H
    parser.add_argument("--G4H",
        type=float,
        nargs='?',
        choices=[float_range(-4,4)],
        default=False,
        help="Use G4Hunter score threshold to determine positive windows "\
                "(default: 0.9)",
        metavar="FLOAT")
    # Use G4NN
    parser.add_argument("--G4NN",
        type=float,
        nargs='?',
        choices=[float_range(0,1)],
        default=False,
        help="Use G4NN score threshold to determine positive windows "\
                "(default: 0.5)",
        metavar="FLOAT")
    args = parser.parse_args()
    if args.cGcC == None:
        args.cGcC = 4.5
    if args.G4H == None:
        args.G4H = 0.9
    if args.G4NN == None:
        args.G4NN = 0.5
    return args

def main():
    args = arguments()
    g4rna_frame = pd.read_csv(args.tsv, sep='\t', index_col=0)
#            10,
    merge_g4rna(g4rna_frame,
            60, 10,
            args.cGcC, args.G4H, args.G4NN)#.to_csv(
#                        path_or_buf=sys.stdout, sep='\t')

if __name__ == '__main__':
    main()
