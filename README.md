# FSDA: Fragment Size Distribution Analysis for non-invasive prenatal CNV prediction
##

## Requirements
* Python 2.6 or newer
* NumPy & SciPy
* samtools
* bcftools

## Pre-process
Use the preprocess.sh script in the preprocess directory to create the required files for prediction. Based on the size of the sample and reference, this can take hours.

Usage:

```
./preprocess.sh [sample_dir] [ref_dir] [prep_dir] [ref_genome] [dbsnp_file]
```

* **sample_dir** is the directory containing the .bam and the .bai files for the sample cfDNA data
* **ref_dir** is the directory containing the .bam and the .bai files for the reference cfDNA data
* **prep_dir** is the directory where the script puts the created files inside
* **ref_genome** is path to the reference genome fasta file (e.g. hg19.fa) 
* **dbsnp_file** is path to a list of commons snps for the reference genome, such as snp142Common.txt, obtainable from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp142Common.txt.gz

## Prediction
The script predict_region.sh in the directory prediction can be used to predict the copy number of a region. 

Usage:

```
./predict_region.sh [sample_dir] [ref_dir] [prep_dir] [chr] [begin] [end] [ks_threshold] [neighbours_threshold] [fraction] 
```

* **sample_dir** is the the directory containing the .bam and the .bai files for the sample cfDNA data
* **ref_dir** is the the directory containing the .bam and the .bai files for the reference cfDNA data
* **prep_dir** is the the directory containing the files created in the pre-process
* **chr** is the chromosome of the target region
* **begin** is the beginning position of the target region
* **end** is the end position of the target region
* **ks_threshold** is the threshold for control similarity
* **neighbours_threshold** is the threshold for the required minimum number of neighbours
* **fetal_fraction** is the fraction of the cfDNA in the sample that is fetal origin

Example:

```
./predict_region.sh ~/sample_dir ~/ref_dir ~/prep_dir chr3 10000000 11000000 0.002 5 0.13
```

