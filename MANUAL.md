<Use a Markdown document viewer to display this file as an HTML file in your>
<internet browser:>
<Markdown Viewer 1.12 (Firefox extension) *recommended>
<Markdown Reader 1.0.12 (Chrome extension)>

**G4RNA SCREENER MANUAL**
=========================

## **NAME**

G4RNA Screener - The nucleic acid screener for RNA G-quadruplexes

## **SYNOPSIS**

    ./screen.py [-?|--help] [-V|--version]

    ./screen.py [-a|--ann <path>.pkl] [-c|--columns <list>] [-w|--window <int>]
                [-s|--step <int>] [-b|--bedgraph]
                [-f|--fasta <path>.fas] [-e|--error] [-v|--verbose]

## **DESCRIPTION**

Score nucleic acid using an artificial neural network classifier (G4NN) that was
trained on sequences found in the [G4RNA database](http://scottgroup.med.
usherbrooke.ca/G4RNA/). It also provides the previously described: 
[G4Hunter score](https://www.ncbi.nlm.nih.gov/pubmed/26792894) and 
[cG/cC score](https://www.ncbi.nlm.nih.gov/pubmed/24121682) if specified.

## **OPTIONS**

**_-?, --help_**

> Show usage and exit.

**_-V, --version_**

> Show version information and exit.

**_-a, --ann_** = G4RNA_2016-11-07.pkl

> Path to the pickled ANN (.pkl) which will provide the program a particular
pattern to evaluate each sequences or windows of sequences.

**_-f, --fasta_** = STDIN

> Path to the fasta file to analyze. Will support string value as long as it
respects the fasta format. Use "STDIN" to feed standard input to -f argument.

**_-w, --window_**

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

###### Fasta example

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
