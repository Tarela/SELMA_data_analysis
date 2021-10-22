# SELMA_data_analysis

## Section 1: Bias estimation (Figure 1-3, Supplementary Figure 1-4):
### To generate the scatter plots for k-mers' bias score (Figure 1 d-k, Supplementary Figure 1 a-c, Supplementary Figure 3 k-n):
#### 1. estimate naive k-mer bias 
script: ATACseqbias_fromBED_upperletter.py, example command line: 
```sh
python ATACseqbias_fromBED_upperletter.py -p hg38.chrom.sizes -t ATACseq_reads.bed -f 5 -s hg38.2bit -o ATACseq_naive10merBias.txt
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
\# ATACseq_summit200.bed, ATACseq_SELMA10merBias.txt, ATACseq_obsCuts_plus.bw, and ATACseq_obsCuts_minus.bw are the intermediate results generated from the above steps.<br>
\# the parameter "Cspan" represent the length of flanking background region considered for each bp. <br>
\# the parameter "flank" represent the 5' only mode of bias consideration. For example to use 10-mer model, the --flank parameter was set to 5 (k/2 of the k-mer) <br>
\# ATACseq_obsExpCuts.txt is the output files for observed and expected cleavages on peaks. <br>
\# the Rscript pred_obs_cmd.r will return a tiny file (ATACseq_obsExpCutsCor.txt) for the correlation coefficient between observed and expected cleavages, the correlation coefficient was used as the height of the bar in the barplots (Figure 1 l-o) or the color in the heatmap (Supplementary Figure 2).<br>

### To calculate the correlation coefficient with different bias estimation method (Figure 2d-i, Supplementary Figure 3 c-j)
The correlation coefficient for all these plots were calculated in exactly the same way as for figure 1. The only difference is that the bias expected cleavages in Figure 2 were calculated with different methods (indicated by the x-axis of the barplots)
#### 1. 5' only method
The "5' only" method of bias estimation was estimated with the method described in the above section but the k-mer length was fixed at 10-mer ("--flank 5").
#### 2. SELMA method
In the The "SELMA" method of bias for each 10-mer was estimated with the same way as "5'only" method. The bias expected cleavages were estimated similarly using the script: scan_cuts_bias_region.py with parameter "-t fxr" (short for forwardXreverse), example command line: 
```sh
python /sfs/qumulo/qproject/CPHG/ZANG/sh8tv/Script/ATAC/scan_cuts_bias_region.py  --Cspan 25  -t fxr  -i ATACseq_summit200.bed -o ATACseq_obsExpCuts.txt -b ATACseq_SELMA10merBias.txt -p ATACseq_obsCuts_plus.bw -n ATACseq_obsCuts_minus.bw
```
#### 3. published method
\# note that for the published method we didn't add the SELMA bias improvement (i.e., the script Seqbias_compare_8mer_encoding_pred_obs.py was no longer used). 
\# for the DNase-seq analsyis (Supplementary Figure 3 c-d) we used the same scripts and parameters for ATACseq data and just switch the data to DNase-seq data. 
a. For "Martins2017" method the bias for each k-mer was estimated using the script ATACseqbias_fromBED_upperletter.py with the similar parameters as for the Figure 1 but set the "--biastype" parameter to "sob" (short for seqOutBias). Then the observed and bias expected cleavages were estimated by the example command line:
```sh
python /sfs/qumulo/qproject/CPHG/ZANG/sh8tv/Script/ATAC/scan_cuts_bias_region.py  --Cspan 25  -t sob  -i ATACseq_summit200.bed -o ATACseq_obsExpCuts.txt -b ATACseq_sob10merBias.txt -p ATACseq_obsCuts_plus.bw -n ATACseq_obsCuts_minus.bw
```
where the ATACseq_sob10merBias.txt was the bias score estimated by ATACseqbias_fromBED_upperletter.py with the parameter "--biastype sob".<br>

b. For "Baek2017" method the bias for each k-mer was estimated using the script ATACseqbias_fromBED_upperletter.py with the similar parameters as for the Figure 1 but set the "--biastype" parameter to "bagfoot" (short for bagfootr). Then the observed and bias expected cleavages were estimated by the example command line:
```sh
python /sfs/qumulo/qproject/CPHG/ZANG/sh8tv/Script/ATAC/scan_cuts_bias_region.py  --Cspan 25  -t bag55  -i ATACseq_summit200.bed -o ATACseq_obsExpCuts.txt -b ATACseq_bagfoot10merBias.txt -p ATACseq_obsCuts_plus.bw -n ATACseq_obsCuts_minus.bw
```
where the ATACseq_bagfoot10merBias.txt was the bias score estimated by ATACseqbias_fromBED_upperletter.py with the parameter "--biastype bagfoot".<br>

c. For "Calviello2019" method the observed and bias expected cleavages were estimated by the example command line:
```sh
python /sfs/qumulo/qproject/CPHG/ZANG/sh8tv/Script/ATAC/scan_cuts_bias_region.py  --Cspan 25  -t repFoot  -i ATACseq_summit200.bed -o ATACseq_obsExpCuts.txt -b ATACseq_repfoot10merBias.txt -p ATACseq_obsCuts_plus.bw -n ATACseq_obsCuts_minus.bw
```
where the ATACseq_bagfoot10merBias.txt was the bias score downloaded from the original study.<br>

### To compare the bp-resolution cleavage pattern on motif centers (Supplementary Figure 3 a-b)
we first extract the observed cleavages from the ATAC-seq data (using strand-splitted bigWig files calculated as described above: plus.bw, minus.bw) using the script: get_cleavage_pattern_motif_fixlen.py, example command line: 
```sh
python /nv/vol190/zanglab/sh8tv/Script/ATAC/get_cleavage_pattern_motif_fixlen.py -i CTCF_motif.bed -o CTCF_motif_ATAC_cleavage.bed -f /PATH/cleavage_bw/ -n useATACseqBW --ext 100
```
\# CTCF_motif is the genome-wide CTCF motif sites (see the method section of SELMA paper for motif sites scanning).
\# /PATH/cleavage_bw/ is the path for the cleavage bigwig files (calculated as described above), there should be two bigWig files, one for plus strand bigWig and the other for minus strand. And the bigWig files should have the name as useATACseqBW_plus.bw and useATACseqBW_minus.bw. The prefix of the file name is specificed by the "-n" parameter. 
\# In the output file (i.e., CTCF_motif_ATAC_cleavage.bed) there are 400 extra columns added to the end of the original 6 motif columns, represent the bp-resolution cleavages (1-200 for plus and 201-400 for minus). Then the file is read into R and the center +/- 50bp of the signal for differet motifs were aligned and take average for the aggregate plot. The plus and minus strand signal are separated for the red and blue lines (labeled as 5' and 3' respectively in the Supplementary Figure). To plot the 9bp shifted version, the plus strand signal is kept the same and the minus strand signal is selected 9bp left to the plus strand signal. 



## Requirements
### python version 2.7, R >= 3.0
### Required modules for scripts: 
ATACseqbias_fromBED_upperletter.py: twobitreader <br>
Seqbias_compare_8mer_encoding_pred_obs.py: numpy <br>
scan_cuts_bias_region.py: bx, twobitreader, numpy <br>
get_cleavage_pattern_motif_fixlen.py: bx <br>



