# SELMA_data_analysis

## Figure 1, Supplementary Figure 1, Supplementary Figure 2:
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
\# this step was also implemented in the SELMA bulk mode. 
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
then transform the ATACseq_obsCuts_plus.bdg and ATACseq_obsCuts_minus.bdg to bigWig format with UCSC tools (bedGraphToBigWig)

#### 2. call peaks and extend to summit +/- 200bp as target regions
use macs2 callpeak function to call ATAC/DNase-seq peaks and extend +/- 200bp from the peak summits
```sh
macs2 callpeak -q 0.01 --keep-dup 1 -f BEDPE -g hs -t ATACseq_reads.bed -n ATACseq
awk '{OFS="\t";if ($2>=200) print $1,$2-200,$2+200,$4,$5;}' ATACseq_summits.bed  > ATACseq_summit200.bed
```
#### 3. estimate biasExpected cuts on peak regions

python /sfs/qumulo/qproject/CPHG/ZANG/sh8tv/Script/ATAC/scan_cuts_bias_region.py  --Cspan 25  -t ${T}  -i /nv/vol190/zanglab/sh8tv/Project/scATAC/Data/ChIP_DNase_ATAC_sametissue/difftype_bias_8mer/use_top50k_summits/${C}_${D}_summit200.bed -o ${C}_${D}_${K}mer_${T}_obsExpCuts.txt -b /sfs/qumulo/qproject/CPHG/ZANG/sh8tv/ATAC/BiasMat/upper_letters/${B} -p /nv/vol190/zanglab/sh8tv/Project/scATAC/Data/ChIP_DNase_ATAC_sametissue/combine_bed/${C}_${D}_plus.bw  -n /nv/vol190/zanglab/sh8tv/Project/scATAC/Data/ChIP_DNase_ATAC_sametissue/combine_bed/${C}_${D}_minus.bw



## Requirements
### python version 2.7, R >= 3.0
### Required modules for scripts: 
ATACseqbias_fromBED.py: twobitreader <br>
Seqbias_compare_8mer_encoding_pred_obs.py: numpy <br>
