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

from g4base import *
import os
import sys
import argparse

def apply_network(ann,
        fasta,
        columns,
        wdw_len,
        wdw_step,
        bedgraph=None,
        verbose=False):
    """
    Wrapping function
    Uses a provided ANN in pickled format (.pkl) and retrieves a dataframe of
    sequences to obtain their G4NN score. Saves the complete values in a .csv
    file.
    
    Wrapping functions don't return values but combine functions to achieve
    something.
    """
    if "all" in columns:
        columns = ['gene_symbol','mrnaAcc','protAcc','gene_stable_id',
                'transcript_stable_id','full_name','HGNC_id','identifier',
                'source','genome_assembly','chromosome','start','end','strand',
                'length','sequence','cGcC','G4H','G4NN']
    #else:
    #    columns = regex.split(",", columns.strip("[]"))
    columns_to_drop = []
    for essential in ['length', 'sequence', 'g4']:
        if essential not in columns:
            columns.append(essential)
            columns_to_drop.append(essential)
    if str(fasta)[0] == '>':
        RNome_df = gen_G4RNA_df(fasta_str_fetcher(fasta, verbose=verbose),
                columns, 1, int(wdw_len), int(wdw_step), verbose=verbose)
    elif str(fasta)[-3:] == '.fa'\
    or str(fasta)[-4:] in ['.fas', '.txt']\
    or str(fasta)[-6:] == '.fasta'\
    or fasta == "/dev/stdin":
        RNome_df = gen_G4RNA_df(fasta_fetcher(fasta, 0, 0, verbose=verbose),
                    columns, 1, int(wdw_len), int(wdw_step), verbose=verbose)
    else:
        screen_usage(52, 'fasta input not specified or not supported')
    if 'G4NN' in columns:
    #    network_file = open(ann,'r')
        ann = pickle.load(ann)
    #    network_file.close()
        RNome_trans_df = trimer_transfo(RNome_df, 'sequence', verbose=verbose)
#        RNome_trans_df = kmer_transfo(RNome_df, 3, 'length', 'sequence', 'g4',
#                int(wdw_len), jellyfish=False, overlapped=True,
#                verbose=verbose)
        RNome_df = submit_seq(ann, RNome_trans_df.drop('G4NN',axis=1),
                [c for c in columns if c != 'G4NN'], "G4NN",
                verbose=verbose)
    if bedgraph:
        sys.stdout.write('browser position %s:%d-%d\n'%(
            RNome_df['chromosome'].iloc[0],
            RNome_df[
                RNome_df.chromosome == RNome_df[
                    'chromosome'].iloc[0]].start.min(),
            RNome_df[
                RNome_df.chromosome == RNome_df[
                    'chromosome'].iloc[0]].end.max()))
        sys.stdout.write('track type=bedGraph name=%s visibility=full \
color=200,100,0\n'%RNome_df.drop(columns_to_drop, axis=1).columns[-1])
    return RNome_df.drop(columns_to_drop, axis=1)

def screen_usage(error_value=False, error_message=False):
    """
    Provide the user with instructions to use screen.py.
    """
    print "Usage: PATH/TO/screen.py [OPTIONS...]"
    print "Use -? or --help to show this message"
    print "Use -V or --version to show program version\n"
    print "Apply options:"
    print "  -a, --ann       \tSupply a pickled ANN (.pkl format)"
    print "  -f, --fasta     \tSupply a fasta file (.fa .fas format)"
    print "  -w, --window    \tWindow length"
    print "  -s, --step      \tStep length between windows"
    print "  -b, --bedgraph  \tDisplay output as bedGraph, user must "\
            "provide columns"
    print "  -c, --columns   \tColumns to display: gene_symbol,sequence,..."
    print "                  \tTo browse available columns use: -c list\n"
    if "-c" and "list" in sys.argv or "--columns" and "list" in sys.argv:
        print "Available columns:"
        print "  description\t\tDescription as available in fasta (Default)"
        print "  all        \t\tAll of the following except description\n"
        print "  gene_symbol\t\tGene symbol"
        print "  mrnaAcc    \t\tRefSeq mRNA accession number"
        print "  protAcc    \t\tRefseq protein accession number"
        print "  gene_stable_id\tEnsembl gene stable ID"
        print "  transcript_stable_id\tEnsembl transcript stable ID"
        print "  full_name  \t\tGene full name (From HGNC)"
        print "  HGNC_id    \t\tHGNC numeric ID"
        print "  identifier \t\tIdentifier"
        print "  source     \t\tSource of the data"
        print "  genome_assembly\tGenome build version"
        print "  chromosome \t\tChromosome"
        print "  start      \t\tStart position"
        print "  end        \t\tEnd position"
        print "  strand     \t\tCoding strand"
        print "  range      \t\tInitial chromosomic range"
        print "  length     \t\tLength of sequence analyzed"
        print "  sequence   \t\tSequence analyzed"
        print "  cGcC       \t\tcGcC score"
        print "  G4H        \t\tG4Hunter score"
        print "  G4NN       \t\tG4NN score of similitude"
        print "             \t\t(must be specified to use ANN)\n"
    print "Other options:"
    print "  -v, --verbose   \tVerbose output with timed operations"
    print "  -e, --error     \tRaise errors and exceptions\n"
    if error_value and error_message:
        sys.stderr.write("UsageError: "+error_message+"\n\n")
        sys.exit(error_value)
    else:
        sys.exit(0)

def main():
    """
    Handles arguments.
    """
    global start_time
    start_time = time.time()
    #Default values here in option_dict
    option_dict = {"--columns":"description,sequence,start,cGcC,G4H,G4NN",
            "--ann":os.path.dirname(__file__)+"/G4RNA_2016-11-07.pkl",
            "--window":60,
            "--step":10,
            "--fasta":"STDIN"}
    for no, arg in enumerate(sys.argv):
        if arg[0] == "-":
            if arg in ["-?","--help"]:
                screen_usage()
            elif arg in ["-V","--version"]:
                print "Version: G4RNA screener 0.3"
                sys.exit(0)
            elif arg in ["-b","--bedgraph",
                    "-v","--verbose",
                    "-e","--error"]:
                option_dict[arg] = True
            elif arg in ["-a","--ann",
                    "-f","--fasta",
                    "-c","--columns",
                    "-w","--window",
                    "-s","--step"]:
                try:
                    option_dict[arg] = sys.argv[no+1]
                except:
                    if "-e" in option_dict.keys() or "--error" \
                    in option_dict.keys():
                        raise
                    else:
                        screen_usage(51, 'No value provided for option "%s"'%arg)
            else:
                screen_usage(51, 'Argument "%s" not recognized'%arg)
    if ("-c" in option_dict.keys() and option_dict["-c"] == "list") \
    or (option_dict["--columns"] == "list"):
        screen_usage()
    if len(sys.argv) == 1 and sys.stdin.isatty():
        screen_usage(51, "no arguments detected")
    if ("-b" in option_dict.keys() or "--bedgraph" in option_dict.keys()):
        if "-c" in  option_dict.keys():
            column_str = "-c"
        elif "--columns" in option_dict.keys():
            column_str = "--columns"
        if len(option_dict.get(column_str).split(',')) != 4 \
        or (['chromosome','start','end'] <= option_dict.get(
            column_str).split(',')) == False \
        or (set(['G4NN','cGcC','G4H']).isdisjoint(option_dict.get(
            column_str).split(','))):
            screen_usage(51, 'bedGraph format requires 4 columns: '\
                    'chromosome,start,end,[SCORE]\n'\
                    '               where [SCORE] is either cGcC, G4H or G4NN')
    if "-f" in option_dict.keys() and option_dict['-f'] == "STDIN":
        option_dict['-f'] = "/dev/stdin"
    elif option_dict['--fasta'] == "STDIN":
        option_dict['--fasta'] = "/dev/stdin"
    try:
        apply_network(option_dict.get("-a") or option_dict.get("--ann"),
                option_dict.get("-f") or option_dict.get("--fasta"),
                option_dict.get("-c") or option_dict.get("--columns"),
                option_dict.get("-w") or option_dict.get("--window"),
                option_dict.get("-s") or option_dict.get("--step"),
                option_dict.get("-b") or option_dict.get("--bedgraph"),
                verbose=option_dict.get("-v") or option_dict.get("--verbose")
                ).to_csv(
                        path_or_buf=sys.stdout, sep='\t',
                        index=(option_dict.get("-b")==None and
                            option_dict.get("--bedgraph")==None),
                        header=(option_dict.get("-b")==None and
                            option_dict.get("--bedgraph")==None))
    except:
        if "-e" in option_dict.keys() or "--error" in option_dict.keys():
            raise
        else:
            screen_usage(50, 'An option is missing, incorrect or not authorized')

class Formatter(argparse.HelpFormatter):
    # use defined argument order to display usage
    def _format_usage(self, usage, actions, groups, prefix):
        if prefix is None:
            prefix = 'usage: '

        # if usage is specified, use that
        if usage is not None:
            usage = usage % dict(prog=self._prog)

        # if no optionals or positionals are available, usage is just prog
        elif usage is None and not actions:
            usage = '%(prog)s' % dict(prog=self._prog)
        elif usage is None:
            prog = '%(prog)s' % dict(prog=self._prog)
            # build full usage string
            action_usage = self._format_actions_usage(actions, groups) # NEW
            usage = ' '.join([s for s in [prog, action_usage] if s])
            # omit the long line wrapping code
        # prefix with 'usage:'
        return '%s%s\n\n' % (prefix, usage)

def arguments():
    """
    Handles arguments
    """
    parser = argparse.ArgumentParser(formatter_class=Formatter,
            prog="screen.py",
            description="Identification of potential RNA G-quadruplexes",
            epilog="G4RNA screener  Copyright (C) 2018  Jean-Michel Garant "\
            "This program comes with ABSOLUTELY NO WARRANTY. This is free "\
            "software, and you are welcome to redistribute it under certain "\
            "conditions <http://www.gnu.org/licenses/>.")
    parser.add_argument('FASTA',
            type=str,
            help='FASTA file (.fa .fas)')
    parser.add_argument("-a", "--ann",
            type=argparse.FileType('r'),
            default=os.path.dirname(__file__)+"/G4RNA_2016-11-07.pkl",
            help="Supply a picled ANN (default: G4RNA_2016-11-07.pkl)")
    parser.add_argument("-w", "--window",
        type=int,
        default=60,
        help="Window length (default: 60)",
        metavar="INT")
    parser.add_argument("-s", "--step",
            type=int,
            default=10,
            help="Step length between windows (default: 10)",
            metavar="INT")
    parser.add_argument("-b", "--bedgraph",
            action="store_true",
            default=False,
            help="Display output as BedGraph, user must provides columns")
            ## TODO use choices of three scores as bedgraph options which will
            ## select columns for the user, must include verifications
    parser.add_argument("-c", "--columns",
            nargs="+",
            choices=["list",
                "all",
                "description",
                "gene_symbol",
                "mrnaAcc",
                "protAcc",
                "gene_stable_id",
                "transcript_stable_id",
                "full_name",
                "HGNC_id",
                "identifier",
                "source",
                "genome_assembly",
                "chromosome",
                "start",
                "end",
                "strand",
                "range",
                "length",
                "sequence",
                "cGcC",
                "G4H",
                "G4NN",
                ],
            default=["description","sequence","start","cGcC","G4H","G4NN"],
            help="Columns to display (default: description). "\
                    "To browse available columns use: -c list",
            metavar="")
    parser.add_argument("-v", "--verbose",
            action="store_true",
            default=False,
            help="Verbose output with timed operations")
    parser.add_argument("-e", "--error",
            action="store_true",
            default=False,
            help="Raise errors and exceptions")

    return parser

def to_replace_main():
    """
    Functions calls
    """
    parser = arguments()
    args = parser.parse_args()
    if args.columns == ["list"]:
        splitted_help = parser.format_help().split(
        ". To browse available columns use:\n\
                        -c list (default: description)")
        print("\n\t".join([splitted_help[0],
                "Available columns:",
                "description\t\tDescription as available in fasta (Default)",
                "all        \t\tAll of the following except description\n",
                "gene_symbol\t\tGene symbol",
                "mrnaAcc    \t\tRefSeq mRNA accession number",
                "protAcc    \t\tRefseq protein accession number",
                "gene_stable_id\tEnsembl gene stable ID",
                "transcript_stable_id\tEnsembl transcript stable ID",
                "full_name  \t\tGene full name (From HGNC)",
                "HGNC_id    \t\tHGNC numeric ID",
                "identifier \t\tIdentifier",
                "source     \t\tSource of the data",
                "genome_assembly\tGenome build version",
                "chromosome \t\tChromosome",
                "start      \t\tStart position",
                "end        \t\tEnd position",
                "strand     \t\tCoding strand",
                "range      \t\tInitial chromosomic range",
                "length     \t\tLength of sequence analyzed",
                "sequence   \t\tSequence analyzed",
                "cGcC       \t\tcGcC score",
                "G4H        \t\tG4Hunter score",
                "G4NN       \t\tG4NN score of similitude",
                "           \t\t(must be specified to use ANN)",
                splitted_help[1]]))
    if args.bedgraph and len(args.columns) != 4 and args.columns[-1] not in [
            'cGcC', 'G4H', 'G4NN']:
        parser.print_usage()
        sys.stderr.write(parser.prog+': error: '\
                'BedGraph format requires 4 columns: '\
                'chromosome start end [SCORE]\n'\
                'where [SCORE] is either cGcC, G4H or G4NN\n')
        sys.exit()
    print args
    try:
        apply_network(args.ann,
                args.FASTA,
                args.columns,
                args.window,
                args.step,
                args.bedgraph,
                args.verbose
                ).to_csv(
                        path_or_buf=sys.stdout, sep='\t',
                        index=(args.bedgraph==False),
                        header=(args.bedgraph==False))
    except:
        if args.error:
            raise
        else:
            screen_usage(50, 'An option is missing, incorrect or not authorized')


if __name__ == '__main__':
#    main()
    to_replace_main()
