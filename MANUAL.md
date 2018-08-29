<Use a Markdown document viewer to display this file as an HTML file in your>
<internet browser>

**G4RNA SCREENER MANUAL**
=========================

1. [Name](#name)
2. [Synopsis](#synopsis)
3. [Screen.py](#screen.py)
    1. [Usage](#screen-usage)
    2. [Description](#screen-description)
    3. [Arguments](#screen-arguments)
    4. [Fasta Example](#screen-fasta-example) 
4. [Merge.py](#merge.py)
    1. [Usage](#merge-usage)
    2. [Description](#merge-description)
    3. [Arguments](#merge-arguments)
    4. [Example](#merge-example)

## **NAME**

G4RNA Screener - The nucleic acid screener for RNA G-quadruplexes

## **SYNOPSIS**

```
./screen.py [-h|--help] [-V|--version]

./screen.py <path>[.fa|.fas|.fasta]
            [-a|--ann <path>.pkl]
            [-c|--columns <list>]
            [-w|--window <int>]
            [-s|--step <int>]
            [-b|--bedgraph]
            [-e|--error]
            [-v|--verbose]
```

```
./merge.py [-h|--help] [-V|--version]

./merge.py <path>[.tsv|.csv|.txt]
           [--cGcC [<float>]]
           [--G4H [<float>]]
           [--G4NN [<float>]]
           [-w|--window <int>]
           [-s|--step <int>]
           [-a|--aggregate <str>]
```

> Combination using pipeline
> ```./screen.py <path>.fa | ./merge.py -```
> example:
> ```./screen.py example.fasta | ./merge.py - --G4NN > example.tsv```

**SCREEN.PY**
=============

## **SCREN USAGE**

    screen.py [-h] FASTA [-a ANN] [-w INT] [-s INT] [-c  [...]] [-b] [-v] [-V] [-e]

## **SCREEN DESCRIPTION**

Score nucleic acid using an artificial neural network classifier (G4NN) that was
trained on sequences found in the [G4RNA database](http://scottgroup.med.
usherbrooke.ca/G4RNA/). It also provides the previously described: 
[G4Hunter score](https://www.ncbi.nlm.nih.gov/pubmed/26792894) and 
[cG/cC score](https://www.ncbi.nlm.nih.gov/pubmed/24121682) if specified.

## **SCREEN ARGUMENTS**

**_FASTA Positional Argument_** = - (default STDIN)

> Path to the fasta file to analyze. Will support string value as long as it
respects the fasta format. Use "-" to feed standard input.

**_-h, --help_**

> Show usage and exit.

**_-a, --ann_** = G4RNA_2016-11-07.pkl

> Path to the pickled ANN (.pkl) which will provide the program a particular
pattern to evaluate each sequences or windows of sequences.

**_-w, --window_** = 60

> Length of the sliding window that is used to analyze long sequences.

**_-s, --step_** = 10

> Step size that moves the window along the long sequences.

**_-c, --columns_** = description

> List of columns to be displayed in the output, use "all" to display them all.
The information must be available in the fasta description of the sequence in
order for the program to retrieve it. It currently supports RefSeq (refGene)
description format and Ensembl (Chromosomic sequences) format. The list of 
columns must be supplied using single commas without spaces:
gene,description,chromosome,strand,...
>
> + description:    _Sequence description as provided by the fasta file
                    annotation (Default display)_
>
> + gene_symbol:    _Gene symbol (Automatically retrieved from UCSC database
                    if the provided description is supported)_
>
> + mrnaAcc:        _mRNA Accession number (Automatically retrieved from UCSC
                    database if the provided description is supported)_
>
> + protAcc:        _Protein Accession number (Automatically retrieved from UCSC
                    database if the provided description is supported)_
>
> + gene_stable_id: _Gene stable ID (Automaticlly retrieved from Ensembl
                    database if the provided description is supported)_
>
> + transcript_stable_id: _Transcript stable ID (Automaticlly retrieved
                          from Ensembl database if the provided description
                          is supported)_
>
> + full_name:      _Gene full name (Automatically retrieved from Ensembl
                    database if the provided description is supported)_
>
> + identifier:     _Identifier provided by user with Ensembl range_
>
> + source:         _Source of the sequence_
>
> + genome_assembly: _Genome assembly_
>
> + start:          _Start position on genomic positive strand_
>
> + end:            _End position on genomic positive strand_
>
> + strand:         _Coding strand_
>
> + range:          _Chromosomic range as provided_
>
> + length:         _Length of sequence analyzed_
>
> + sequence:       _Sequence analyzed_
>
> + cGcC:           _cGcC score_
>
> + G4H:            _G4Hunter score_
>
> + G4NN:           _Score obtained through G4NN (Must be specified when
                    enumerating columns in list)_

**_-v, --verbose_**

> Rudimental verbose

**_-V, --version_**

> Show version information and exit.

**_-e, --error_**

> Raise and display error and warnings

###### Screen Fasta Examples

Any fasta file can support columns: description,length,sequence,cGcC,G4H,G4NN
```html
>Description1
SEQUENCE1
>Description2 *Make sure that the descriptions are unique
SEQUENCE2
```

To use the other columns, the user must provide a fasta file with one of the
following description formats:

**RefSeq**

The identification of a RefSeq description is provided by a regular expression
that requires an NM_######, NR_######, XM_######, XR_######. It also support
protein accession number NP_###### or XP_######. The regular expression will
capture the information provided in the description if it is presented in the
correct order such as the examples below:
```html
>genome-build_source_Accession range=chromosome:start-end 5'pad= 3'pad= strand= repeatMasking=
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>hg38_refGene_NM_001276352 range=chr1:67092176-67134971 5'pad=0 3'pad=0 strand=- repeatMasking=none
CGAGUAACCCGCGGAGCCAGAAGAGGGAGGAAAGGAGAUGAGguuuguuu
>NP_852042
CGAGTAACCCGCGGAGCCAGAAGAGGGAGGAAAGGAGATGAGgtttgttt
```

**Ensembl**

The identification of an Ensembl description is provided by a regular
expression that requires an optional identifier as first argument, an Ensembl
stable ID and a large range like argument starting by literal "chromosome" and
ending by strand (-1 or 1) such as the example below:
```html
>identifier stable_id alphabet:chromosome chromosome:genome-build:chromosome:start:end:strand
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>1 ENST00000380152 dna:chromosome chromosome:GRCh38:11:106073501:106100710:1
tgacggctattatcgcatgatcGGCGGCGGCGGGGcggcggcgatatcgt
>ENST00000380152 chromosome:GRCh38:11:106073501:106100710:1
ACGTTTACGGCTATTCGTATGCTTGTGCTATTGCATGTACTG
```

**MERGE.PY**
============

## **MERGE USAGE**

    merge.py [-h] TSV [--cGcC [FLOAT]] [--G4H [FLOAT]] [--G4NN [FLOAT]] [-w INT] [-s INT] [-a STR]

## **MERGE DESCRIPTION**

Transforms the tabular output of `screen.py` in order to merge overlapping
regions that are scored above threshold(s). The omission of using --cGcC, --G4H
and --G4NN argument wil merge all overlapping windows back into the original
sequences. The -w,--window and -s,--step arguments must be the same used for
the `screen.py`. 

## **MERGE ARGUMENTS**

**_TSV Positional Argument_** = - (default STDIN)

> Path to the tabular separated values file to analyze. Use "-" to feed
standard input.

**_-h, --help_**

> Show usage and exit.

**_--cGcC_** = 4.5

> Use consecutive G over consecutive C score to determine positive windows
above threshold and merge them.

**_--G4H_** = 0.9

> Use G4Hunter score to determine positive windows above threshold and merge
them.

**_--G4NN_** = 0.5

> Use G4 neural network score to determine positive windows above threshold and
merge them.

**_-w, --window_** = 60

> Length of the sliding window that was used to analyze long sequences.

**_-a, --aggregation_** = list (only available with python 2.7.15rc1 +)

> Scores of merged windows will be aggregated using an aggregation function:
(choose from 'max', 'min', 'median', 'mean', 'std', 'sem')

**_-s, --step_** = 10

> Step size that defines the overlap between the windows sequences.

**_-v, --verbose_**

> Show version information and exit.

**_-e, --error_**

> Raise and display error and warnings


###### Merge Examples

`screen.py` output
```
./screen.py spinach.fas -w 60 -s 10
	description	sequence	start	cGcC	G4H	G4NN
1	Spinach aptamer (false negative example)	GGACGCGACCGAAAUGGUGAAGGACGGGUCCAGUGCGAAACACGCACUGUUGAGUAGAGU	1	2.1875	0.3	0.209915966961
2	Spinach aptamer (false negative example)	GAAAUGGUGAAGGACGGGUCCAGUGCGAAACACGCACUGUUGAGUAGAGUGUGAGCUCCG	11	2.2	0.283333333333	0.146798576587
3	Spinach aptamer (false negative example)	AGGACGGGUCCAGUGCGAAACACGCACUGUUGAGUAGAGUGUGAGCUCCGUAACUGGUCG	21	1.88235294118	0.233333333333	0.154917456609
4	Spinach aptamer (false negative example)	CAGUGCGAAACACGCACUGUUGAGUAGAGUGUGAGCUCCGUAACUGGUCGCGUC	31	1.26666666667	0.0555555555556	0.118000916901
```

`screen.py | merge.py` output
```
./screen.py spinach.fas -w 60 -s 10 | ./merge.py - -w 60 -s 10
1	Spinach aptamer (false negative example)	GGACGCGACCGAAAUGGUGAAGGACGGGUCCAGUGCGAAACACGCACUGUUGAGUAGAGUGUGAGCUCCGUAACUGGUCGCGUC	1	[2.1875, 2.2, 1.88235294118, 1.2666666666700002]	[0.3, 0.283333333333, 0.233333333333, 0.0555555555556]	[0.209915966961, 0.14679857658700002, 0.154917456609, 0.118000916901]
```
