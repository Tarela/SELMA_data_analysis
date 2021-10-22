# SELMA_data_analysis

## Section 1: Bias estimation (Figure 1-3, Supplementary Figure 1-4):
### To generate the scatter plots in Figure 1 d-k, Supplementary Figure 1 a-c:
#### 1. estimate naive k-mer bias 
script: ATACseqbias_fromBED.py, example command line: 
```sh
python ATACseqbias_fromBED.py -p hg38.chrom.sizes -t ATACseq_reads.bed -f 5 -s hg38.2bit -o ATACseq_naive10merBias.txt
```
\# hg38.chrom.sizes and hg38.2bit are the chromosome sizes and genome sequence of human hg38 genome (or mouse mm10 genome) downloaded from UCSC genome browser or NCBI.<br>
\# ATACseq_reads.bed is the aligned reads in bed format (input file, support both SE or PE data). <br>
\# ATACseq_naive10merBias.txt is the naive k-mer bias matrix (output files, 1st column for k-mer, 2nd column for the naive k-mer bias). 
#### 2. estimate SELMA k-mer bias from naive k-mer bias
script: Seqbias_compare_8mer_encoding_pred_obs.py, example command line: 
```sh
python Seqbias_compare_8mer_encoding_pred_obs.py -b ATACseq_naive10merBias.txt -o ATACseq_SELMA10merBias.txt
```
\# ATACseq_naive10merBias.txt is the naive k-mer bias (input file of this step, output of step1)
\# ATACseq_SELMA10merBias.txt is the SELMA k-mer bias (output file, with header, 1st column for k-mer, 2nd column for naive k-mer bias, 3rd column for SELMA k-mer bias). 
#### 3. generate scatter plots
read in the naive or SELMA k-mer bias matrix into R and generate scatter plots with any R functions by comparing the bias from two conditions (two datasets). The Pearson correlation coefficients were calculated with the "cor" function in R. The correlation coefficients were also plotted as bar plots in Supplementary Figure 1d. 

### To generate the genome-wide signal correlation between observed and expected cuts in Figure 1 l-o, Supplementary Figure 2:
\# this step was also implemented in the SELMA bulk mode. Users can run SELMA bulk mode directly to get the observed and bisExpected cuts on peak regions.
#### 1. pileup genome-wide cleavage signal
note that the aligned reads bed file were first splitted and transformed to +/- cleavage sites. <br>
example for paired end data: (note that paired end reads bed file contain 3 columns for chromosome, fragment start, and fragment end). 
```sh
awk '{OFS="\t";print $1,$2,$2+1,".",".","+"}' ATACseq_reads.bed > ATACseq_plus.bed
awk '{OFS="\t";print $1,$3-1,$3,".",".","-"}' ATACseq_reads.bed > ATACseq_minus.bed
```
use macs2 pileup function to generate genome-wide +/- cleavage signal in bedGraph/bigWig format 
```sh
macs2 pileup -i ATACseq_plus.bed -o ATACseq_obsCuts_plus.bdg --extsize 1 -f BED
macs2 pileup -i ATACseq_minus.bed -o ATACseq_obsCuts_minus.bdg --extsize 1 -f BED
```
then transform the ATACseq_obsCuts_plus.bdg and ATACseq_obsCuts_minus.bdg to bigWig format (ATACseq_obsCuts_plus.bw, ATACseq_obsCuts_minus.bw) with UCSC tools (bedGraphToBigWig)

#### 2. call peaks and extend to summit +/- 200bp as target regions
use macs2 callpeak function to call ATAC/DNase-seq peaks and extend +/- 200bp from the peak summits
```sh
macs2 callpeak -q 0.01 --keep-dup 1 -f BEDPE -g hs -t ATACseq_reads.bed -n ATACseq
awk '{OFS="\t";if ($2>=200) print $1,$2-200,$2+200,$4,$5;}' ATACseq_summits.bed  > ATACseq_summit200.bed
```
#### 3. scan and calculate observed and biasExpected cleavages on peak regions
script: scan_cuts_bias_region.py, example command line: 
```sh
python /sfs/qumulo/qproject/CPHG/ZANG/sh8tv/Script/ATAC/scan_cuts_bias_region.py  --Cspan 25  -t flank  -i ATACseq_summit200.bed -o ATACseq_obsExpCuts.txt -b ATACseq_SELMA10merBias.txt -p ATACseq_obsCuts_plus.bw -n ATACseq_obsCuts_minus.bw
Rscript pred_obs_cmp.r ATACseq_obsExpCuts.txt ATACseq_obsExpCutsCor.txt
```
\# ATACseq_summit200.bed, ATACseq_SELMA10merBias.txt, ATACseq_obsCuts_plus.bw, and ATACseq_obsCuts_minus.bw are the intermediate results generated from the above steps
\# the parameter "Cspan" represent the length of flanking background region considered for each bp
\# the parameter "flank" represent the 5' only mode of bias consideration. In the following section it will be replaced by "-t fxr" for SELMA suggested methods. 
\# ATACseq_obsExpCuts.txt is the output files for observed and expected cleavages on peaks. 
\# the Rscript pred_obs_cmd.r will return a tiny file (ATACseq_obsExpCutsCor.txt) for the correlation coefficient between observed and expected cleavages, the correlation coefficient was used as the height of the bar in the barplots (Figure 1 l-o, Supplementary Figure 2)


## Section 2: bulk footprint analysis (Figure 4, Supplementary Figure 5-6)

## Section 3: single cell ATAC-seq analysis (Figure 5, Supplementary Figure 7)




## Requirements
### python version 2.7, R >= 3.0
### Required modules for scripts: 
ATACseqbias_fromBED.py: twobitreader <br>
Seqbias_compare_8mer_encoding_pred_obs.py: numpy <br>
scan_cuts_bias_region.py: bx, twobitreader, numpy <br>
