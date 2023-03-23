# Mishpokhe is a self-supervised algorithm to discover functional spatial clusters

It predicts the spatial functional clusters in ordered proteomes based on the unordered set of functionally related proteins.

## Usage

Mishpokhe is a command line tool

`python3 py_mishpokhe2.py -q {QUERYDB} -t {TARGETDB} -r {RES} -tf {TARGETFA} -qf {QUERYFA}`

The help can be called with 

`python3 py_mishpokhe2.py -h`

# Options

  * `-h`, `--help`            show this help message and exit
  * `-q QUERYDB`, `--querydb QUERYDB`
                        Provide path to query MMseqs2 database (protein,
                        unordered)
  * `-t TARGETDB`, `--targetdb TARGETDB`
                        Provide path to target MMseqs2 database (protein,
                        ordered)
  * `-r RES`, `--res RES`     Specify the name to be given to the results of the
                        search
  * `-tf TARGETFA`, `--targetfa TARGETFA`
                        Provide path to target sequence (fasta (linearized!))
			That is the initial fasta containing proteomes (ordered) to use as a target
  * `-qf QUERYFA`, `--queryfa QUERYFA`
                        Provide path to query sequence (fasta)
			That is the initial fasta containing proteins (unordered) to use as query
  * `-i ITER`, `--iter ITER`  Give the number of iterations, default is 1
			Recommended to set not more than 3
  * `-s SINGLETON`, `--singleton SINGLETON`
                        Set to 1 if you have singleton queries, default is 0
		        Only set in quite special cases


## Tips:

MMseqs2 can be installed via conda. 

`conda install -c bioconda mmseqs2`

MMseqs2 database can be made with:

`mmseqs createdb yourseq.fasta yourseq_db`

*For more options on the installation and usage of MMseqs2, please, visit the GitHub  [MMseqs2](https://github.com/soedinglab/MMseqs2/ "MMseqs2")


Fasta can be linearized with:

`awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' your_seqs.fasta > your_seqs.faa_lin`

It is important to keep mishpokhe together with other accessory scripts in one directory. 




