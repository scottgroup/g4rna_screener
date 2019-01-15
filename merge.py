#!/usr/bin/env python2.7

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

# temporary warning filter until packages with numpy dependancies updates
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import os
import sys
import utils
import argparse
import pandas as pd

pd.set_option('display.max_columns', None)

class float_range(object):
    """
    Object that defines a range of float that is authorized or returns the
    authorized range in error. Used to validate score threshold input by user.
    """
    # float range is a function that is needed to have any threshold of float in
    # supported range of scores
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

#def too_much_merge_g4rna(df, window=60, step=10,
#        cGcC=False, G4H=False, G4NN=False):
#    """
#    Merges consecutive windows that are scored above the provided thresholds
#    into a single hit sequence and discards windows below the thresholds.
#    Three level of verification that windows are consecutives:
#    
#    description column:
#    For window with index n, consecutive windows share the same description
#    description(n) == description(n+1)
#    
#    start column:
#    For window with index n, the consecutive windows have
#    start(n+1) - start(n) = step
#    
#    sequence column:
#    For window with index n and length l, the consecutive windows have
#    Last part of sequence(n) == first part of sequence(n+1)
#    sequence(n)[step-l:] == sequence(n+1)[:l-step]
#    """
#    verif = []
#    if 'description' in df.columns:
#        verif.append(set(df.index[
#                df.description.eq(df.description.shift(-1)) |
#                df.description.shift().eq(df.description)].tolist()))
#        #verif.append('description')
#    # start and next start are distanced by step length new column True
#    if 'start' in df.columns:
#        verif.append(set(df.index[
#                df.start.eq(df.start.shift()+step) | 
#                df.start.shift(-1).eq(df.start+step)].tolist()))
#        #verif.append('start')
#    # overlap is the length that sequential windows should share
#    overlap = window-step
#    if 'sequence' in df.columns:
#        verif.append(set(df.index[
#                df.sequence.str[-overlap:].eq(
#                    df.sequence.str[:overlap].shift(-1)) | 
#                df.sequence.str[-overlap:].shift().eq(
#                    df.sequence.str[:overlap])].tolist()))
#        #verif.append('sequence')
#    if len(verif) == 3:
#        keep = verif.pop(0).intersection(verif.pop(0),verif.pop(0))
#    elif len(verif) == 2:
#        keep = verif.pop(0).intersection(verif.pop(0))
#    elif len(verif) == 1:
#        keep = verif.pop(0)
#
#    return df

def merge_g4rna(df, window=60, step=10,
        cGcC=False, G4H=False, G4NN=False,
        score_aggregation=list):
    """
    """
    # filter windows with users thresholds
    if cGcC:
        df = df[ df.cGcC >= cGcC ].dropna()
    if G4H:
        df = df[ df.G4H >= G4H ].dropna()
    if G4NN:
        df = df[ df.G4NN >= G4NN ].dropna()
    # aggregate functions are defined in this dictionnary
    agg_fct = {
#            'description':"".join,
            'gene_symbol':'max',
            'mrnaAcc':'mode',
            'protAcc':'mode',
            'gene_stable_id':'mode',
            'transcript_stable_id':'mode',
            'full_name':'max',
            'HGNC_id':'mode',
            'identifier':'mode',
            'source':'mode',
            'genome_assembly':'mode',
            'chromosome':'mode',
            'start':'min',
            'end':'max',
            'strand':'mode',
            'range':'mode',
            'length':'max',
            'sequence':'max',
            'cGcC':score_aggregation,
            'G4H':score_aggregation,
            'G4NN':score_aggregation
            }
    # overlap is the length that sequential windows should share
    overlap = window-step
    pd.set_option('display.max_colwidth', -1)
    # the merge is dependant on the sequence columns
    if 'sequence' in df.columns:
        for ite in range(0,len(df)):#max iterations is dataframe length
            # check if any windows still needs some merging
            if any(df.sequence.str[:(overlap+step*ite)].eq(
                df.sequence.str[step:].shift(1))) is False:
                break
            # merge overlapping window by appending the [step] first
            # nucleotides to the next window
            df.loc[
                    df.sequence.str[:(overlap+step*ite)].eq(
                        df.sequence.str[step:].shift(1))
                    , 'sequence'] = df.sequence.str[:step].shift(1) + \
                            df.sequence.str[:]
        # description column is used when available which should prevent the
        # eventual problems caused by a user that submits two sequences that
        # share some subsequences
        if 'description' in df.columns:
            df_grouped = df.groupby(
                    [df.description,df.sequence.str[:overlap]],
                    sort=False,
                    as_index=False)
            # grouping the resulting windows using the first [:overlap] nts
            return df_grouped.agg(
                    {k:agg_fct[k] for k in df.columns.drop(['description'])}
                    ).reindex(columns=df.columns)
        else:
            df_grouped = df.groupby(
                    df.sequence.str[-overlap:],
                    sort=False,
                    as_index=False)
            return df_grouped.agg(
                    {k:agg_fct[k] for k in df.columns},
                    ).reindex(df.columns, axis=1)
    else:
        sys.stderr.write("UsageError: 'sequence' column must be provided\n")
        sys.exit()

def arguments():
    """
    Arguments management
    """
    # declare argument parser
    parser = argparse.ArgumentParser(formatter_class=utils.Formatter,
            prog=os.path.basename(__file__),
            description="[WORK IN PROGRESS] [DO NOT DISTRIBUTE] Merge positive windows of screen.py output and "\
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
    # windows length
    parser.add_argument("-w", "--window",
        type=int,
        default=60,
        help="Windows length of input",
        metavar="INT")
    # step length
    parser.add_argument("-s","--step",
        type=int,
        default=10,
        help="Step lenght of input",
        metavar="INT")
    # aggregation function for scores
    parser.add_argument("-a", "--aggregation",
        #nargs="+",
        choices=["max","min","median","mean","std","sem",list],
        default=list,
        help="Aggregation function to pool scores of merged windows "\
                "(default: list)",
        metavar="STR")
    parser.add_argument("-V", "--version",
            action="version",
            version="%(prog)s 0.3")
    # useful for debug, not meant for users
    parser.add_argument("-e", "--error",
            action="store_true",
            default=False,
            help="Raise errors and exceptions")
    # needed to to have default help display
    if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()
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
    try:
        merge_g4rna(g4rna_frame,
                args.window, args.step,
                args.cGcC, args.G4H, args.G4NN,
                args.aggregation).to_csv(
                            path_or_buf=sys.stdout, sep='\t')
    except:
        # raise python error calls if -e --error is used
        if args.error:
            raise
        # custom error message
        else:
            sys.stderr.write(parser.prog+': error: '\
            'An option is missing, incorrect or not authorized\n')

if __name__ == '__main__':
    main()
