# SELMA_data_analysis

## Figure 1, Supplementary Figure 1, Supplementary Figure 2:
### To generate the scatter plots in Figure 1 d-k, Supplementary Figure 1 a-c:
#### 1. estimate naive k-mer bias 
script: ATACseqbias_fromBED.py, example command line: 
```sh
python ATACseqbias_fromBED.py -p hg38.chrom.sizes -t ATACseq_reads.bed -f 5 -s hg38.2bit -o ATACseq_naive10merBias.txt
```
\# hg38.chrom.sizes and hg38.2bit are the chromosome sizes and genome sequence of human hg38 genome (or mouse mm10 genome) downloaded from UCSC genome browser or NCBI.<br>
ATACseq_reads.bed is the aligned reads in bed format (input file, support both SE or PE data). <br>
ATACseq_naive10merBias.txt is the naive k-mer bias matrix (output files, 1st column for k-mer, 2nd column for the naive k-mer bias). 
#### 2. estimate SELMA k-mer bias from naive k-mer bias
script: Seqbias_compare_8mer_encoding_pred_obs.py, example command line: 
```sh
python Seqbias_compare_8mer_encoding_pred_obs.py -b ATACseq_naive10merBias.txt -o ATACseq_SELMA10merBias.txt
```
where ATACseq_naive10merBias.txt is the naive k-mer bias (input file of this step, output of step1) and ATACseq_SELMA10merBias.txt is the SELMA k-mer bias (output file, with header, 1st column for k-mer, 2nd column for naive k-mer bias, 3rd column for SELMA k-mer bias). 
#### 3. generate scatter plots
read in the naive or SELMA k-mer bias matrix into R and generate scatter plots with any R functions by comparing the bias from two conditions (two datasets). The Pearson correlation coefficients were calculated with the "cor" function in R. The correlation coefficients were also plotted as bar plots in Supplementary Figure 1d. 

### To generate the genome-wide signal correlation between observed and expected cuts in Figure 1 l-o, Supplementary Figure 2:
#### 1. 



## Requirements
### python version 2.7, R >= 3.0
### Required modules for scripts: 
ATACseqbias_fromBED.py: twobitreader <br>
Seqbias_compare_8mer_encoding_pred_obs.py: numpy <br>
