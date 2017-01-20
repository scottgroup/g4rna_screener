<Use a Markdown document viewer to display this file as an HTML file in your>
<internet browser:>
<Markdown Viewer 1.12 (Firefox extension) *recommended>
<Markdown Reader 1.0.12 (Chrome extension)>

**README**
==========

## **USER INFO**

#### Instigator

name:               _Jean-Michel Garant_

email:              _jean-michel.garant@usherbrooke.ca_

gitlab username:    _J-Michel_


## **REPO INFO**

gitlab url:         *gitlabscottgroup.med.usherbrooke.ca/J-Michel/g4rna-screener*

Version: G4RNA screener Alpha


## **DEPENDENCIES**

Here are listed dependencies in format:
```bash
library_name (recommended version)
```

##### From environment

```bash
python2.7 (2.7.12)
```

##### From Python

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
> cd PATH/TO/PYTHON/dist-packages/    or    cd PATH/TO/PYTHON/site-packages/
> sudo -i
> wget https://dev.mysql.com/get/Downloads/Connector-Python/mysql-connector-python-2.1.4.tar.gz
> tar -xzf mysql-connector-python-2.1.4.tar.gz
> cd mysql-connector-python-2.1.4
> python setup.py install
> ```
