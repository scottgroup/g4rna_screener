#!/usr/bin/env python

from g4base import *

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
                'source','genome_build','chromosome','start','end','strand',
                'length','sequence','cGcC','G4H','G4NN']
    else:
        columns = regex.split(",", columns.strip("[]"))
    columns_to_drop = []
    for essential in ['length', 'sequence', 'g4']:
        if essential not in columns:
            columns.append(essential)
            columns_to_drop.append(essential)
    if type(fasta) == 'str' and fasta[0] == '>':
        RNome_df = gen_G4RNA_df(fasta_str_fetcher(fasta, verbose=verbose),
                columns, 1, int(wdw_len), int(wdw_step), verbose=verbose)
    else:
        RNome_df = gen_G4RNA_df(fasta_fetcher(fasta, 0, 0, verbose=verbose),
                columns, 1, int(wdw_len), int(wdw_step), verbose=verbose)
    if 'G4NN' in columns:
        network_file = open(ann,'r')
        ann = pickle.load(network_file)
        network_file.close()
        RNome_trans_df = kmer_transfo(RNome_df, 3, 'length', 'sequence', 'g4',
                int(wdw_len), jellyfish=False, overlapped=True,
                verbose=verbose)
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
    return RNome_df.drop(columns_to_drop, axis=1).to_csv(
            path_or_buf=sys.stdout, sep='\t',
            index=(bedgraph==None), header=(bedgraph==None))

def screen_usage(error_message=False):
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
    print "  -b, --bedgraph  \tDisplay output as bedGraph, user must \
provide columns"
    print "  -c, --columns   \tColumns to display: gene,sequence,..."
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
        print "  genome_build\t\tGenome build version"
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
    if error_message:
        print "UsageError:", error_message
        sys.exit(500)
    else:
        sys.exit(0)

def main():
    """
    Handles arguments.
    """
    global start_time
    start_time = time.time()
    #Default values here in option_dict
    option_dict = {"--columns":"description,sequence,cGcC,G4H,G4NN",
            "--ann":"G4RNA_2016-11-07.pkl",
            "--window":60,
            "--step":10,
            "--fasta":"STDIN"}
    for no, arg in enumerate(sys.argv):
        if arg[0] == "-":
            if arg in ["-?","--help"]:
                screen_usage()
            elif arg in ["-V","--version"]:
                print "Version: G4RNA screener 0.1"
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
                        screen_usage('No value provided for option "%s"'%arg)
            else:
                screen_usage('Argument "%s" not recognized'%arg)
    if ("-c" in option_dict.keys() and option_dict["-c"] == "list") \
    or (option_dict["--columns"] == "list"):
        screen_usage()
    if len(sys.argv) == 1 and sys.stdin.isatty():
        screen_usage("no arguments detected")
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
            screen_usage('bedGraph format requires 4 columns: \
chrome,start,end,[SCORE]\n\
        where [SCORE] is either cGcC, G4H or G4NN')
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
                verbose=option_dict.get("-v") or option_dict.get("--verbose"))
    except:
        if "-e" in option_dict.keys() or "--error" in option_dict.keys():
            raise
        else:
            screen_usage('An option is missing, incorrect or not authorized')


if __name__ == '__main__':
    main()
