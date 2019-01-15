<Use a Markdown document viewer to display this file as an HTML file in your>
<internet browser:>

**README**
==========

## **USER INFO**

#### Instigator

name:               _Jean-Michel Garant_

email:              _jean-michel.garant@usherbrooke.ca_

gitlab username:    _J-Michel_


## **REPO INFO**

gitlab url:         *gitlabscottgroup.med.usherbrooke.ca/J-Michel/g4rna_screener*

Version: G4RNA screener 0.3

**Please consider cloning/downloading a stable branch**


## **DEPENDENCIES**

Here are listed dependencies in format:
```bash
library_name (recommended version)
```

##### From environment

```bash
python2.7 (2.7.15rc1)
```

##### From Python

Available from pip

```bash
biopython (1.68)
numpy (1.11.0)
pandas (0.18.1)
PyBrain (0.3)
regex (2016.9.22)
scipy (0.18.1)
```

> ```mysql-connector-python (2.1.4)``` is also highly recommended for
> automated retrieval of information on UCSC and Ensembl databases using
> -c, --columns arguments. It is not available through pip but here are the
> steps to follow to install it:
> ```bash
> sudo -i
> cd PATH/TO/PYTHON/dist-packages/    or    cd PATH/TO/PYTHON/site-packages/
> wget https://dev.mysql.com/get/Downloads/Connector-Python/mysql-connector-python-2.1.4.tar.gz --no-check-certificate
> tar -xzf mysql-connector-python-2.1.4.tar.gz
> cd mysql-connector-python-2.1.4
> python setup.py install
> ```

**Consider adding G4RNA screener to your environment path**

Example for a bash terminal
```bash
cd PATH/TO/g4rna_screener
echo "" >> ~/.bashrc
echo "# add G4RNA screener to PATH" >> ~/.bashrc
echo "export PATH=\"\$PATH:$(pwd)\"" >> ~/.bashrc
source ~/.bashrc
```

## **GLOSSARY**

* ANN: _Artificial Neural Network_

* AUC: _Area Under [ROC](https://en.wikipedia.org/wiki/Receiver_operating_characteristic)
 Curve_

* cGcC: _consecutive guanine over consecutive cytosine, usually expressed as a
[score](http://www.ncbi.nlm.nih.gov/pubmed/24121682)_

* csv: _Comma separated values. A tabular text file readable by spreadsheet
softwares such as Microsoft Excel and LibreOffice Calc. Files used in this
project are actually .tsv ([tab separated values](https://en.wikipedia.org/wiki/Comma-separated_values#Standardization)) 
for better visualization in the terminal._

* cv: _Cross-validation_

* db: _Database_

* DnB: _[Dot'n'Bracket](http://ultrastudio.org/en/Dot-Bracket_Notation) notation
for secondary structure_

* [Ensembl](http://www.ensembl.org/index.html): Joint project between [European
Bioinformatics Institute](http://www.ebi.ac.uk/) (EBI) and the [Wellcome Trust Sanger Institute](http://www.sanger.ac.uk/) (WTSI)

* fasta: _[File format for biological sequences](http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml)_

* G4: _G-quadruplex_

* [G4RNA](http://scottgroup.med.usherbrooke.ca/G4RNA/): _G-quadruplex RNA
database_

* [NCBI](http://www.ncbi.nlm.nih.gov/): _National Center for Biotechnology Information_

* mfe: _Minimum free energy (kcal/mol)_

* mp2: _[Mammouth parallel 2](https://wiki.calculquebec.ca/w/Mp2)_

* [MySQL](http://www.mysql.com/): _Open Source SQL database management system 
from Oracle Corporation_

* nt: _Nucleotide_

* PG4: _Potential G-quadruplex_

* PRAC: _Pavillon de Recherche Appliquée sur le Cancer de [l'Université de
Sherbrooke](http://www.usherbrooke.ca/) (Pavilion of Applied Research on Cancer)_

* [PSQL](https://www.postgresql.org): _Open source object-relational database system_

* [RefSeq](http://www.ncbi.nlm.nih.gov/refseq/): _The Reference Sequence
database of NCBI_

* regex: _[Regular expression](http://dobromirivanov.net/wp-content/uploads/2013/06/Oreilly.Introducing.Regular.Expressions.Jul.2012.pdf)_

* ROC: _[Receiver Operating Characteristic](https://en.wikipedia.org/wiki/Receiver_operating_characteristic)_

* RNA: _RiboNucleic acid_

* UCSC: _University of California in Santa Cruz_
