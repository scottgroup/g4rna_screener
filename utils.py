#!/usr/bin/env python

import time
import sys
import regex
import pickle
import pandas as pd
from collections import Counter, OrderedDict


def hms_string(time_elapsed):
    """
    Format time in hours, minutes, seconds.
    
    Must be included in a string to format a time value.
    """
    h = int(time_elapsed / (60 * 60))
    m = int((time_elapsed % (60 * 60)) / 60)
    s = time_elapsed % 60.
    return "{}:{:>02}:{:>05.2f}".format(h, m, s)

def verbosify(verbose, message, time_it=False):
    """
    Take care of the verbosity for the user.
    
****Time_it requires a global variable -> start_time****
    
    Supports both Boolean value of verbose and numerical level of verbose.
    Either print or flush message.
    """
    if (verbose == True or verbose > 0) and time_it == True:
        print message+"\t"*(4-len(message)/8)+\
        "{}".format(hms_string(time.time() - start_time))
    elif (verbose == True or verbose > 0) and time_it == False:
        print message
    elif verbose == None:
        pass
    else:
        if "\n" in message:
            print message
        else:
            sys.stdout.write(message+"..."+" "*(77-len(message))+"\r")
            sys.stdout.flush()

def connect_psql(host, user, passwd, db, number_of_cursor):# schema,
    """
    Connects to a PostGreSQL database and generates cursors.
    
    Returns a list.
    [0] connection object
    [1] first cursor
    [2] second cursor
    ...
    [n] last cursor
    """
    import psycopg2
    mydb = psycopg2.connect(None, db, user, passwd, host)
    cursors = [mydb]
    for cursor_nb in range(number_of_cursor):
        cursors.append(mydb.cursor())
#        cursors[cursor_nb+1].execute("SET search_path TO %s"%schema)
    return cursors

def retrieve_RefSeq(mrna_accession=None, prot_accession=None):
    """
    Retrieve cross-reference information from UCSC database using a
    RefSeq accession.
    
    Returns a list.
    [0] gene symbol
    [1] gene complete name
    [2] mRNA accession
    [3] protein accession
    """
    import mysql.connector
    mydb = mysql.connector.connect(
            host='genome-mysql.cse.ucsc.edu', user='genome', db='hgFixed')
    cursor = mydb.cursor(buffered=True)
    cursor.execute(
            'SELECT mrnaAcc, protAcc, name, product FROM refLink WHERE \
mrnaAcc = "%s" OR protAcc = "%s"'%(mrna_accession, prot_accession))
    if mrna_accession == None and prot_accession == None:
        return [None,None,None,None]
    else:
        return list(cursor.fetchone())

def retrieve_xref_Ensembl(stable_id=None,mrnaAcc=None,
        gene_acronym=None):
    """
    Retrieve cross-reference information from Ensembl.
    
    Returns a list.
    [0] gene stable id
    [1] transcript stable id
    [2] accession
    [3] gene id
    [4] transcript id
    [5] gene acronym (HGNC)
    [6] gene description
    """
    import mysql.connector
    mydb = mysql.connector.connect(
            host='ensembldb.ensembl.org', user='anonymous',
            db='homo_sapiens_core_86_38')
    cursor = mydb.cursor(buffered=True)
    base_sql = 'SELECT gene.stable_id, transcript.stable_id, xref.dbprimary_acc\
, gene.gene_id, transcript.transcript_id, CASE WHEN \
gene_attrib.attrib_type_id = 4 THEN gene_attrib.value ELSE NULL END AS \
gene_symbol, gene.description FROM xref JOIN object_xref ON xref.xref_id = \
object_xref.xref_id JOIN transcript ON transcript.transcript_id = \
object_xref.ensembl_id JOIN gene ON transcript.gene_id = gene.gene_id JOIN \
gene_attrib ON gene_attrib.gene_id = transcript.gene_id WHERE '
    try:
        cursor.execute(base_sql+'(dbprimary_acc LIKE "N_\_%%" OR dbprimary_acc \
LIKE "X_\_%%") AND (gene.stable_id = "%s" OR transcript.stable_id = "%s") \
ORDER BY attrib_type_id'%(stable_id,stable_id))
        return list(cursor.fetchone())
    except:
        try:
            cursor.execute(base_sql+'dbprimary_acc = "%s" ORDER BY attrib_type_id'%
                    mrnaAcc)
            return list(cursor.fetchone())
        except:
            cursor.execute(base_sql+'(dbprimary_acc LIKE "N_\_%%" OR \
dbprimary_acc LIKE "X_\_%%") AND value = "%s"'%gene_acronym)
            return list(cursor.fetchone())

def format_description(fas_description, verbose=None):
    '''
    Takes a fasta description line and try to retrieve informations out of it
    if it has a known format.
    
    Returns a dictionnary of available informations
    '''
    infos = {}
    try:
        infos = regex.match("(?<description>(?<genome_build>\D\D\d+)?\
(?:_(?<source>[^_]*))?_?(?:(?P<mrnaAcc>[N|X][M|R]_\d+)|(?P<protAcc>[N|X]P_\d+))\
(?: range=(?<range>(?<chromosome>chr.*):(?<start>(\d*))-(?<end>(\d*))))?\
(?: 5'pad=(?<pad5>\d*))?(?: 3'pad=(?<pad3>\d*))?(?: strand=(?<strand>.))?\
(?: repeatMasking=(?<repeatMasking>.*))?)",
                fas_description).groupdict()
    except:
        verbosify(verbose,"RefSeq not recognised for %s"%fas_description)
    try:
        infos = regex.match("(?<description>(?:(?<identifier>[^ ]*) )?(?<stable_id>ENS[T|G]\d*)?(?:\.\d)?(?: (?<info>[^ ]*))?(?: chromosome:(?<genome_build>[^:]*):(?<range>(?<chromosome>[^:]*):(?<start>\d*):(?<end>\d*)):(?<strand>.)(?:.*)))", 
                fas_description).groupdict()
        infos['source'] = "Ensembl"
    except:
        verbosify(verbose,"Ensembl not recognised for %s"%fas_description)
    if 'description' not in infos.keys() or infos['description'] == '':
        infos['description'] = fas_description
    if 'chromosome' in infos.keys() and infos['range'][:3] != 'chr':
        infos['chromosome'] = 'chr' + infos['chromosome']
    if 'range' in infos.keys() and infos['range'][:3] != 'chr':
        infos['range'] = 'chr' + infos['range']
    if 'strand' in infos.keys() and infos['strand'] == '1':
        infos['strand'] = '+'
    return infos

def fasta_fetcher(
        fasta_file,
        number_to_fetch,
        seq_size,
        verbose=False):
    """
    Fetch for sequences from a fasta file and returns a defined number of 
    random sequences or random window from random sequences if seq_size is
    not 0.
    
    number_to_fetch = 0 takes everything
    seq_size = 0 takes full length sequences
    
    Return a dictionnary.
    {Description:sequence}
    """
    from Bio import SeqIO
    fas_dic = OrderedDict()
    for seq in SeqIO.parse(fasta_file, 'fasta'):
        if len(seq.seq) > seq_size and seq_size != 0:
            r_int = np.random.randint(0, len(seq.seq)-seq_size)
            fas_dic[seq.description] = str(seq.seq)[r_int:r_int+seq_size]
        else:
            fas_dic[seq.description] = str(seq.seq)
    dic = {}
    verbosify(verbose, "File fetched")
    if number_to_fetch == 0:
        return fas_dic
    else:
        randomize = np.random.permutation(len(fas_dic))
        for i in range(number_to_fetch):
            dic[fas_dic.keys()[randomize[i]]] = fas_dic[
                    fas_dic.keys()[randomize[i]]].strip('N').strip('n')
        return dic

def fasta_str_fetcher(fasta_string, verbose=False):
    """
    Fetch for sequences in a fasta file presented as a string.
    
    Return a dictionnary.
    {Description:sequence}
    """
    fas_dic = OrderedDict()
    for instance in regex.split(r'\\r\\n>|\\n>|>', fasta_string)[1:]:
        [description, seq] = regex.split(r'\\r\\n|\\n', instance, maxsplit=1)
        fas_dic[description] = regex.sub(r'\\r\\n|\\n','', seq)
    return fas_dic

def kmer_transfo(
        df_,
        depth,
        sort_column,
        sequence_column,
        target_column,
        window,
        jellyfish=False,
        overlapped=True,
        verbose=False):
    """
    Define sequences by their kmers proportions and returns a bigger
    dataframe containing it.
    
    jellyfish = True uses jellyfish command to get the kmers
    overlapped = True allows kmers to overlapped each other (almost
    always the case)
    
    Return pandas dataframe.
    """
    df = df_.copy()
    nts = ['A','U','C','G']
    di_nts = []
    if depth == 1:
        di_nts = nts
    else:
        for nt1 in nts:
            for nt2 in nts:
                if depth == 2:
                    di_nts.append(nt1+nt2)
                elif depth == 3:
                    for nt3 in nts:
                        di_nts.append(nt1+nt2+nt3)
                elif depth == 4:
                    for nt3 in nts:
                        for nt4 in nts:
                            di_nts.append(nt1+nt2+nt3+nt4)
                else:
                    print "This kmer length isn't available"
                    break
    for each in di_nts:
        df[each] = .0
    lst_rows = []
    pos_dframe = pd.DataFrame(columns=df.columns)
    for row in df.index:
        lst_rows.append(row)
        if jellyfish is True:
            di_nt_cnts = {}
            for line_out in subprocess.check_output(u'echo ">0\n%s" | \
sed "s/U/T/g" | jellyfish count -m 3 -s 100 -o /dev/stdout /dev/stdin | \
jellyfish dump -ct /dev/stdin | sed "s/T/U/g"'%(df.ix[row,sequence_column].
    upper().replace('T','U')), shell=True).split('\n')[:-1]:
                di_nt_cnts[line_out.split('\t')[0]] = int(line_out.split('\t')[1])
        else:
            di_nt_lst = regex.findall(
                    '.{%d}'%depth,df.ix[row,sequence_column].upper().replace('T','U'),
                    overlapped=True)
            di_nt_cnts = Counter(di_nt_lst)
        if len([di_nt_cnts[x] for x in di_nt_cnts]) > 4**depth:
            print "Crap! There's a non AUCG nt in the sequence at row %d"%row
            break
        total_di_nt = sum(di_nt_cnts.values())
        di_nt_freqs = [(str(di_nt), float(di_nt_cnts[di_nt])/total_di_nt)
                for di_nt in di_nt_cnts if "N" not in di_nt]
        for di_ntd, freq in di_nt_freqs:
            df.ix[row,di_ntd] = freq
    verbosify(verbose, "Kmer transformed")
    return df