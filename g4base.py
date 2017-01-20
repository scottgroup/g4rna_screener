#!/usr/bin/env python

from utils import *
from pybrain.datasets import ClassificationDataSet


def cgcc_scorer(sequence):
    """
    Returns the cGcC scores of a sequence.
    
****NOTE: The cGcC score does not consider flanking G/C like the G4Hunter****
    
    Returns a float.
    """
    G_score = 0
    C_score = 0
    for i in xrange(1,len(sequence)+1):
        G_score+=len(regex.findall(ur'[gG]{%d}'%i,sequence,overlapped=True))*10*i
        C_score+=len(regex.findall(ur'[cC]{%d}'%i,sequence,overlapped=True))*10*i
    if C_score == 0:
        C_score = 1
    return float(G_score)/float(C_score)

def g4hunter(sequence, output_map=False):
    """
    Apply a G4Hunter score calculation to the sequence.
    
    Default output is a mean value for the sequence, it should be used to
    analyze a sequence as whole.
    Map output is avaible providing G4H values for each nt, it allows the G/Cs
    at start or end of a window to have a value which takes into account the
    possible flanking G/Cs since they would have an impact on the G4H score at
    those positions.
    
    Returns a float between -4 and 4. [DEFAULT]
    Returns a list of floats between -4 and 4. [output_map=True]
    """
    g4h_list = [gr.group(1) for gr in regex.finditer(u'(?i)((?P<nt>[a-z])(?P=nt)*)',
        sequence, overlapped=False)]
    g4h_map = [(float(value) if -4 <= value <=4
        else (4.0 if value > 4 else -4.0))
        for match in [[len(x)]*len(x) if x[0]=='G' or x[0]=='g'
            else ([-len(x)]*len(x) if x[0]=='C' or x[0]=='c' else [0]*len(x))
            for x in g4h_list] for value in match]
    if output_map:
        return g4h_map
    else:
        return sum(g4h_map)/float(len(g4h_map))

def gen_G4RNA_df(
        seq_dict,
        columns,
        first_id,
        window_fragment=False,
        window_step=False,
        verbose=False):
    """
    Generate a dataframe like the one used to train the ANN from G4RNA to
    eventually concatenate them usually to G4RNA data. The dataframe can be
    filled with any controls (RNomic or scrambled) provided in a dictionnary.
    
    When window_fragment and window_step are defined, the dataframe is
    populated with substring from dictionnary instead of whole sequences.
    
    Return pandas dataframe.
    """
    data = []
    for key in seq_dict.keys():
        infos = format_description(key)
        content = {}
        for ke in infos.keys():
            content[ke] = infos[ke]
        if "gene_symbol" in columns or "mrnaAcc" in columns\
        or "protAcc" in columns or "gene_stable_id" in columns\
        or "transcript_stable_id" in columns or "full_name" in columns\
        or "HGNC_id" in columns:
            try:
                [content["gene_stable_id"],
                        content["transcript_stable_id"],
                        content["mrnaAcc"],
                        content["gene_id"],
                        content["transcript_id"],
                        content["gene_symbol"],
                        content["gene_description"]
                        ] = retrieve_xref_Ensembl(infos.get('stable_id'),
                                infos.get('mrnaAcc'), retrieve_RefSeq(
                                    infos.get('mrnaAcc'),infos.get('protAcc')
                                    )[2])
                [content["full_name"],content["HGNC_id"]] = list(
                        regex.search(
                        '(?<full_name>.*) \[Source:HGNC Symbol;Acc:HGNC:\
(?<HGNC_    id>\d+)\]',content["gene_description"]).group('full_name','HGNC_id'))
            except:
                if 'stable_id' in infos.keys() and infos.get('stable_id')[3] == 'T':
                    content['transcript_stable_id'] = infos.get('stable_id')
                elif 'stable_id' in infos.keys() and infos.get('stable_id')[3] == 'G':
                    content['gene_stable_id'] = infos.get('stable_id')
                else:
                    pass
            try:
                [content['mrnaAcc'],
                        content['protAcc'],
                        content['gene_symbol'],
                        content['product']
                        ] = retrieve_RefSeq(infos["mrnaAcc"],infos["protAcc"])
            except:
                try:
                    [content['mrnaAcc'],
                            content['protAcc'],
                            content['gene_symbol'],
                            content['product']
                            ] = retrieve_RefSeq(content["mrnaAcc"])
                except:
                    pass
        for ke in infos.keys():
            if infos[ke] != '' and infos[ke] is not None:
                content[ke] = infos[ke]
        if window_fragment == False or window_fragment >= len(seq_dict[key]):
            seq = [seq_dict[key]]
        else:
            seq = [seq_dict[key][w:w+window_fragment] for w in range(
                0,len(seq_dict[key])-window_fragment+window_step,window_step)]
            g4h_map = g4hunter(seq_dict[key], output_map=True)
        try:
            strt = int(infos.get('start'))
        except:
            strt = None
        for no, s in enumerate(seq):
            row = []
            content['id'] = first_id+len(data)
            if strt != None:
                content['start'] = int(strt+no*window_step)
                content['end'] = int(strt+no*window_step+len(s)-1)
            content['length'] = len(s)
            content['sequence'] = s
            content['cGcC'] = cgcc_scorer(s)
            if window_fragment == False or window_fragment >= len(seq_dict[key]):
                content['G4H'] = g4hunter(s)
            else:
                content['G4H'] = sum(g4h_map[
                    no*window_step:no*window_step+window_fragment])/float(len(
                        g4h_map[no*window_step:no*window_step+window_fragment]))
            content['g4'] = 'N/A'
            for column in columns:
                try:
                    row.append(content[column])
                except:
                    row.append(None)
            data.append(row)
    verbosify(verbose, "DataFrame built")
    return pd.DataFrame(data=data, index=range(
        first_id, first_id+len(data)), columns=columns)

def submit_seq(
        ann,
        df_,
        except_columns,
        score_name,
        verbose=False):
    """
    Submit sequences from a dataframe to an ANN to be scored from 0 to 1.
    
    Return pandas dataframe with a new column "score_name".
    """
    df=df_.copy()
    alldata_tst = ClassificationDataSet(len(df.columns)-len(except_columns),
            1, nb_classes=2)
    trans_set_tst = df.drop(except_columns, axis=1)
    for n in df.index:
        if str(df.loc[n].g4) in ['yes', 'True']:
            target = 1.
        elif str(df.loc[n].g4) in ['no', 'False', 'N/A']:
            target = 0.
        alldata_tst.addSample(trans_set_tst.loc[n][:], target)
    alldata_tst._convertToOneOfMany()
    test_results = ann.activateOnDataset(alldata_tst)
    final_df = df_.copy()
    final_df[score_name]= pd.Series(test_results[:,1], index=final_df.index)
    verbosify(verbose, 'Sequence submitted')
    return final_df.drop(
            [c for c in final_df.columns if c not in except_columns][:-1],
            axis=1)
