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
                        <br /> &emsp;&emsp;&emsp;Provide path to query MMseqs2 database (protein,
                        unordered)
  * `-t TARGETDB`, `--targetdb TARGETDB`
                        <br /> &emsp;&emsp;&emsp;Provide path to target MMseqs2 database (protein,
                        ordered, provided in Prodigal .faa format)
  * `-r RES`, `--res RES`
    			<br /> &emsp;&emsp;&emsp;Specify the name to be given to the results of the
                        search
  * `-tf TARGETFA`, `--targetfa TARGETFA`
                        <br /> &emsp;&emsp;&emsp;Provide path to target sequence (fasta (linearized!))
			<br /> &emsp;&emsp;&emsp;That is the initial fasta containing proteomes (ordered) to use as a target, provided in Prodigal .faa format
  * `-qf QUERYFA`, `--queryfa QUERYFA`
                        <br /> &emsp;&emsp;&emsp;Provide path to query sequence (fasta)
			<br /> &emsp;&emsp;&emsp;That is the initial fasta containing proteins (unordered) to use as query
  * `-i ITER`, `--iter ITER`  Give the number of iterations, default is 1
			<br /> &emsp;&emsp;&emsp;Recommended to set not more than 3
  * `-s SINGLETON`, `--singleton SINGLETON`
                        <br /> &emsp;&emsp;&emsp;Set to 1 if you have singleton queries, default is 0
		        <br /> &emsp;&emsp;&emsp;Only set in quite special cases
  *  `-u EVALFILTERUSE`, `--evalfilteruse EVALFILTERUSE`
                        <br /> &emsp;&emsp;&emsp;Set to 1 if you want to use mishpokhe clusters e-value filter, default is 1
  * `e EVAL`, `--eval EVAL`,
                        <br /> &emsp;&emsp;&emsp;Specify the e-value threshold, default is 1
  * `-f FRAC_OCC_MIN`, `--frac_occ_min FRAC_OCC_MIN`,
                        Specify the threshold for the fraction of the cluster
                        matches in which each sequence cluster occurs, default
                        is >= 0
  * `-if MIN_FRAC_INSIDE`, `--min_frac_inside MIN_FRAC_INSIDE`,
                        Specify the threshold for the fraction of the cluster
                        matches in which each sequence cluster occurs to all
                        the matches, default is >= 0
  * `-c SEARCH_COV`, `--search_cov SEARCH_COV`,
                        Specify the coverage threshold for the mmseqs search
  * `-af ARC_FILTER`, `--arc_filter ARC_FILTER`,
                        Set to 0 if you do NOT want to use architecture
                        clustering and filtering, default is 1
  * `-b BIAS`, `--bias BIAS`,
    			Specify bias for enrichment scores calculation,
                        default is 0

## How to read the output

There are numerous files that Mishpokhe outputs (including MMseqs2-related), but the the most important ones are:

* {res}..._stat, files starting with the name given by user in the argument -r, --res and ending with "stat". They store unfiltered clusters obtained in each iteration.
* {res}..._stat_filtered, files starting with the name given by user in the argument -r, --res and ending with "stat_filtered". They store clusters filtered by e-value from each iteration. These files are only made in case if -u, --evalfilteruse option is not switched off.
* {res}..._stat_filtered_clu_filter, files starting with the name given by user in the argument -r, --res and ending with "stat_filtered_clu_filter". They store clusters filtered by architecture clustering to the initial queries from each iteration.
  


## Tips:

MMseqs2 can be installed via conda. 

`conda install -c bioconda mmseqs2`

MMseqs2 database can be made with:

`mmseqs createdb yourseq.fasta yourseq_db`

*For more options on the installation and usage of MMseqs2, please, visit the GitHub  [MMseqs2](https://github.com/soedinglab/MMseqs2/ "MMseqs2")


Fasta can be linearized with:

`awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' your_seqs.fasta > your_seqs.faa_lin`

It is important to keep mishpokhe together with other accessory scripts in one directory. 

If you run mishpokhe again on the same dataset, it is important to delete the "tmp" folder as well as query database starting from the first iteration (simply do `rm {QUERYDB}1*`, it removes also the later databases) and all results files (do `rm {RES}*`).


# Credits

https://github.com/vlcc/basplice

