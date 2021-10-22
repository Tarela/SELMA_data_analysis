#!/usr/bin/env python
#Time-stamp:<Tarela>
"""
Description: 

"""

# ------------------------------------
# Python Modual
# ------------------------------------

import os,sys,re
from optparse import OptionParser
import logging
import string
import math,time
try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()
import twobitreader
import numpy
#from numpy import linalg as la

#import scipy.stats.distributions

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def rev(seq):
    revseq = ""
    for i in seq[::-1]:
        if i == 'A':
            r = 'T'
        elif i == 'T':
            r = 'A'
        elif i == 'C':
            r = 'G'
        elif i == 'G':
            r = 'C'
        else:
            r=i#print i
        revseq += r
    return revseq

#def readBG(bgmatrix):
#    BGraw = {}
#    BGenc = {}
#    #nBG = {}
#    BGenc_list = []
#    BGraw_list = []
#    inf = open(bgmatrix)
#    for line in inf:
#        if line.startswith("seqtype"):
#            continue
#        ll = line.split()
#        name = ll[0]
#        if ll[1] == "NA":
#            BGraw[name] = 0
#            BGraw_list.append(0)
#        else:
#            BGraw[name] = float(ll[1])#pow(numpy.e,float(ll[1]))
#            BGraw_list.append(float(ll[2]))
#
#        BGenc[name] = float(ll[2])#pow(numpy.e,float(ll[2]))
#        BGenc_list.append(float(ll[2]))
#    inf.close()
#    seqlen = len(name)
#
#    mean_BGraw = numpy.mean(BGraw_list)
#    sd_BGraw = numpy.std(BGraw_list)
#
#    mean_BGenc = numpy.mean(BGenc_list)
#    sd_BGenc = numpy.std(BGenc_list)
#
#    BG_raw = {}
#    BG_enc = {}
#    for kmer in BGenc.keys():
#        BG_raw_normalize[kmer] = pow(numpy.e, (BGraw[kmer]-mean_BGraw)/sd_BGraw)
#        BG_enc_normalize[kmer] = pow(numpy.e, (BGenc[kmer]-mean_BGenc)/sd_BGenc)
#
#    return BG_raw_normalize,BG_enc_normalize,seqlen

def readBG(bgmatrix):
    BGraw = {}
    BGenc = {}
    #nBG = {}
    inf = open(bgmatrix)
    total_raw = 0
    total_enc = 0
    for line in inf:
        if line.startswith("seqtype"):
            continue
        ll = line.split()
        name = ll[0]
        if ll[1] == "NA":
            BGraw[name] = 0
        else:
            BGraw[name] = pow(numpy.e,float(ll[1]))
            total_raw += pow(numpy.e,float(ll[1]))
        BGenc[name] = pow(numpy.e,float(ll[2]))
        total_enc += pow(numpy.e,float(ll[2]))
    inf.close()

    out_raw = {}
    out_enc = {}

    for kmer in BGraw.keys():
        out_raw[kmer] = BGraw[kmer] / total_raw
        out_enc[kmer] = BGenc[kmer] / total_enc

    seqlen = len(name)
    return out_raw,out_enc,seqlen

def readBG_repfoot(bgmatrix):
    BGraw = {}
    #nBG = {}
    inf = open(bgmatrix)
    total_raw = 0
    for line in inf:
        if line.startswith("seqtype"):
            continue
        ll = line.split()
        name = ll[0]
        BGraw[name] = ll
        if ll[1] == "NA":
            BGraw[name] = 0
        else:
            BGraw[name] = pow(numpy.e,float(ll[1]))
            total_raw += pow(numpy.e,float(ll[1]))
        BGenc[name] = pow(numpy.e,float(ll[2]))
        total_enc += pow(numpy.e,float(ll[2]))
    inf.close()

    out_raw = {}
    out_enc = {}

    for kmer in BGraw.keys():
        out_raw[kmer] = BGraw[kmer] / total_raw
        out_enc[kmer] = BGenc[kmer] / total_enc

    seqlen = len(name)
    return out_raw,seqlen

def SOBfetchSEQ(fullseq,seqlen):
    if seqlen == 11:
        outseq = fullseq[:3]+fullseq[4]+fullseq[8:10]+fullseq[11:15]+fullseq[16]
    elif seqlen == 10:
        outseq = fullseq[:3]+fullseq[4]+fullseq[8:10]+fullseq[11:14]+fullseq[16]
    elif seqlen == 8:
        outseq = fullseq[0]+fullseq[2]+fullseq[8:10]+fullseq[11:14]+fullseq[16]
    elif seqlen == 6:
        outseq = fullseq[0]+fullseq[2]+fullseq[8:10]+fullseq[12:14]
    elif seqlen == 4:
        outseq = fullseq[2]+fullseq[8]+fullseq[12:14]
    else:
        outseq = "NA"
    return outseq

def bias_scan_flank(peak,outfile,w_plus,w_minus,biasMat,Cspan,Gen):

#    t = time.time()
    genome = twobitreader.TwoBitFile(Gen)
    
    BGraw,BGenc,Nmer = readBG(biasMat)
    flank = int(Nmer)/2
    w_plus_H=BigWigFile(open(w_plus, 'rb'))
    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    
    inf = open(peak)
#    w_plus_H=BigWigFile(open(w_plus, 'rb'))
#    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    outf = open(outfile,'w')
    
    for line in inf:
        ll = line.split()
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])

        if start-Cspan < 0:
            continue

        plus_obj = w_plus_H.summarize(chrm,start-Cspan,end+Cspan,(end-start+2*Cspan))
        minus_obj = w_minus_H.summarize(chrm,start-Cspan,end+Cspan,(end-start+2*Cspan))
        if not plus_obj  :
            plus_vector = numpy.array([0]*(end-start+2*Cspan))
        else:
            plus_vector = plus_obj.sum_data 
        if not minus_obj  :
            minus_vector = numpy.array([0]*(end-start+2*Cspan))
        else:
            minus_vector = minus_obj.sum_data 
       
        plus_single_bias_raw_vector = []
        plus_single_bias_enc_vector = []
        #plus_cb_bias_vector = []
        minus_single_bias_raw_vector = []
        minus_single_bias_enc_vector = []
        #minus_cb_bias_vector = []
    
        for pos in range(start-Cspan,end+Cspan):
            plus_seq = genome[chrm][(pos-flank):(pos+flank)].upper()
#            plus_reverse_seq = rev(genome[chrm][(pos+offset-flank):(pos+offset+flank)].upper())
            minus_seq = rev(genome[chrm][(pos-flank+1):(pos+flank+1)].upper())
#            minus_reverse_seq = genome[chrm][(pos-offset+1-flank):(pos-offset+1+flank)].upper()

            ### calculate bias value for each bases
            
            if len(plus_seq) == Nmer and not "N" in plus_seq:
                plus_bias_raw = BGraw[plus_seq]
                plus_bias_enc = BGenc[plus_seq]
            else:
                plus_bias_raw = 0
                plus_bias_enc = 0
                
            #if len(plus_reverse_seq) == Nmer and not "N" in plus_reverse_seq:
            #    plus_reverse_bias = pow(numpy.e,BG[plus_reverse_seq])
            #else:
            #    plus_reverse_bias = 0
            
            if len(minus_seq) == Nmer and not "N" in minus_seq:
                minus_bias_raw = BGraw[minus_seq]
                minus_bias_enc = BGenc[minus_seq]
            else:
                minus_bias_raw = 0
                minus_bias_enc = 0

            #if len(minus_reverse_seq) == Nmer and not "N" in minus_reverse_seq:
            #    minus_reverse_bias = pow(numpy.e,BG[minus_reverse_seq])
            #else:
            #    minus_reverse_bias = 0

            #if plus_bias != "NA" and plus_reverse_bias != "NA":
            #    plus_cb_bias = (plus_bias + plus_reverse_bias ) / 2
            #else:
            #    plus_cb_bias = -10
            
            #if minus_bias != "NA" and minus_reverse_bias != "NA":
            #    minus_cb_bias = (minus_bias + minus_reverse_bias ) / 2
            #else:
            #    minus_cb_bias = -10

            plus_single_bias_raw_vector.append(plus_bias_raw)
            plus_single_bias_enc_vector.append(plus_bias_enc)
            #plus_cb_bias_vector.append(plus_cb_bias)
            minus_single_bias_raw_vector.append(minus_bias_raw)
            minus_single_bias_enc_vector.append(minus_bias_enc)
            #minus_cb_bias_vector.append(minus_cb_bias)

        #if "NA" in plus_single_bias_vector or "NA" in plus_cb_bias_vector or "NA" in minus_single_bias_vector or "NA" in minus_cb_bias_vector:
        #    print ll
        #    continue        
        #### construct bias vector and transform to linear
        Plus_Single_rawBias = numpy.array(plus_single_bias_raw_vector)
        Plus_Single_encBias = numpy.array(plus_single_bias_enc_vector)
        #Plus_Combine_Bias_log = numpy.array(plus_cb_bias_vector)
        Minus_Single_rawBias = numpy.array(minus_single_bias_raw_vector)
        Minus_Single_encBias = numpy.array(minus_single_bias_enc_vector)
        #Minus_Combine_Bias_log = numpy.array(minus_cb_bias_vector)
##### now here
        #Plus_Single_Bias_linear = pow(numpy.e,Plus_Single_Bias_log)
        #Plus_Combine_Bias_linear = pow(numpy.e,Plus_Combine_Bias_log)
        #Minus_Single_Bias_linear = pow(numpy.e,Minus_Single_Bias_log)
        #Minus_Combine_Bias_linear = pow(numpy.e,Minus_Combine_Bias_log)

        #roundN = 4

        #### assign bias to bp and proportion
        for outpos in range(Cspan,(end-start+Cspan)):
                    
            this_plus_single_raw = Plus_Single_rawBias[outpos]
            this_minus_single_raw = Minus_Single_rawBias[outpos]
            this_plus_sum_raw = sum(Plus_Single_rawBias[(outpos-Cspan):(outpos+Cspan)])
            this_minus_sum_raw = sum(Minus_Single_rawBias[(outpos-Cspan):(outpos+Cspan)])

            this_plus_single_enc = Plus_Single_encBias[outpos]
            this_minus_single_enc = Minus_Single_encBias[outpos]
            this_plus_sum_enc = sum(Plus_Single_encBias[(outpos-Cspan):(outpos+Cspan)])
            this_minus_sum_enc = sum(Minus_Single_encBias[(outpos-Cspan):(outpos+Cspan)])

            #this_plus_single_prop = round(Plus_Single_Bias_linear[outpos]/sum(Plus_Single_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)
            #this_plus_cb_prop = round(Plus_Combine_Bias_linear[outpos]/sum(Plus_Combine_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)
            #this_minus_single_prop = round(Minus_Single_Bias_linear[outpos]/sum(Minus_Single_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)
            #this_minus_cb_prop = round(Minus_Combine_Bias_linear[outpos]/sum(Minus_Combine_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)

            this_plus = plus_vector[outpos]
            this_minus = minus_vector[outpos]
            this_plus_cuts_sum = sum(plus_vector[(outpos-Cspan):(outpos+Cspan)])
            this_minus_cuts_sum = sum(minus_vector[(outpos-Cspan):(outpos+Cspan)])

            out_chrm = chrm
            out_start = start + outpos - Cspan
            out_end = out_start+1
            
            if this_plus > 0:
                expcut_plus_raw = this_plus_cuts_sum * (this_plus_single_raw/this_plus_sum_raw)
                expcut_plus_enc = this_plus_cuts_sum * (this_plus_single_enc/this_plus_sum_enc)
                outf.write("\t".join(map(str,[out_chrm+"_"+str(out_start)+"_"+str(out_end)+"_+",this_plus,format(expcut_plus_raw,".3e"),format(expcut_plus_enc,".3e") ]))+"\n")

            if this_minus > 0:
                expcut_minus_raw = this_minus_cuts_sum * (this_minus_single_raw/this_minus_sum_raw)
                expcut_minus_enc = this_minus_cuts_sum * (this_minus_single_enc/this_minus_sum_enc)
                outf.write("\t".join(map(str,[out_chrm+"_"+str(out_start)+"_"+str(out_end)+"_-",this_minus,format(expcut_minus_raw,".3e"),format(expcut_minus_enc,".3e") ]))+"\n")

            #outf_rawPlus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_single]))+"\n")
            #outf_rawMinus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_single]))+"\n")
            #outf_cbPlus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_cb]))+"\n")
            #outf_cbMinus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_cb]))+"\n")
            #
            #outf_rawPlus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_single_prop]))+"\n")
            #outf_rawMinus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_single_prop]))+"\n")
            #outf_cbPlus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_cb_prop]))+"\n")
            #outf_cbMinus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_cb_prop]))+"\n")
    outf.close()
    #outf_rawPlus.close()
    #outf_rawMinus.close()
    #outf_cbPlus.close()
    #outf_cbMinus.close()
    #outf_rawPlus_prop.close()
    #outf_rawMinus_prop.close()
    #outf_cbPlus_prop.close()
    #outf_cbMinus_prop.close()

    inf.close()


def bias_scan_fxr(peak,outfile,w_plus,w_minus,biasMat,Cspan,Gen,offset):

#    t = time.time()
    genome = twobitreader.TwoBitFile(Gen)
    
    BGraw,BGenc,Nmer = readBG(biasMat)
    flank = int(Nmer)/2
    
    inf = open(peak)
    w_plus_H=BigWigFile(open(w_plus, 'rb'))
    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    outf = open(outfile,'w')
    
    for line in inf:
        ll = line.split()
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])

        plus_obj = w_plus_H.summarize(chrm,start-Cspan,end+Cspan,(end-start+2*Cspan))
        minus_obj = w_minus_H.summarize(chrm,start-Cspan,end+Cspan,(end-start+2*Cspan))
        if not plus_obj  :
            plus_vector = numpy.array([0]*(end-start+2*Cspan))
        else:
            plus_vector = plus_obj.sum_data 
        if not minus_obj  :
            minus_vector = numpy.array([0]*(end-start+2*Cspan))
        else:
            minus_vector = minus_obj.sum_data 
    
        plus_single_bias_fxr_vector = []
        #plus_cb_bias_vector = []
        minus_single_bias_fxr_vector = []
        #minus_cb_bias_vector = []

        for pos in range(start-Cspan,end+Cspan):
            plus_seq = genome[chrm][(pos-flank):(pos+flank)].upper()
            plus_reverse_seq = rev(genome[chrm][(pos+offset-flank):(pos+offset+flank)].upper())
            minus_seq = rev(genome[chrm][(pos-flank+1):(pos+flank+1)].upper())
            minus_reverse_seq = genome[chrm][(pos-offset+1-flank):(pos-offset+1+flank)].upper()

            ### calculate bias value for each bases
            
            if len(plus_seq) == Nmer and not "N" in plus_seq:
                plus_bias = BGenc[plus_seq]
            else:
                plus_bias = 0
                
            if len(plus_reverse_seq) == Nmer and not "N" in plus_reverse_seq:
                plus_reverse_bias = BGenc[plus_reverse_seq]
            else:
                plus_reverse_bias = 0
            
            if len(minus_seq) == Nmer and not "N" in minus_seq:
                minus_bias = BGenc[minus_seq]
            else:
                minus_bias = 0

            if len(minus_reverse_seq) == Nmer and not "N" in minus_reverse_seq:
                minus_reverse_bias = BGenc[minus_reverse_seq]
            else:
                minus_reverse_bias = 0

            #if plus_bias != "NA" and plus_reverse_bias != "NA":
            plus_cb_bias = numpy.sqrt(plus_bias * plus_reverse_bias ) 
            #else:
            #    plus_cb_bias = -10
            
            #if minus_bias != "NA" and minus_reverse_bias != "NA":
            minus_cb_bias = numpy.sqrt(minus_bias * minus_reverse_bias)
            #else:
            #    minus_cb_bias = -10

            plus_single_bias_fxr_vector.append(plus_cb_bias)
            #plus_cb_bias_vector.append(plus_cb_bias)
            minus_single_bias_fxr_vector.append(minus_cb_bias)
            #minus_cb_bias_vector.append(minus_cb_bias)

        #if "NA" in plus_single_bias_vector or "NA" in plus_cb_bias_vector or "NA" in minus_single_bias_vector or "NA" in minus_cb_bias_vector:
        #    print ll
        #    continue        
        #### construct bias vector and transform to linear
        Plus_Single_fxrBias = numpy.array(plus_single_bias_fxr_vector)
        #Plus_Combine_Bias_log = numpy.array(plus_cb_bias_vector)
        Minus_Single_fxrBias = numpy.array(minus_single_bias_fxr_vector)
        #Minus_Combine_Bias_log = numpy.array(minus_cb_bias_vector)
##### now here
        #Plus_Single_Bias_linear = pow(numpy.e,Plus_Single_Bias_log)
        #Plus_Combine_Bias_linear = pow(numpy.e,Plus_Combine_Bias_log)
        #Minus_Single_Bias_linear = pow(numpy.e,Minus_Single_Bias_log)
        #Minus_Combine_Bias_linear = pow(numpy.e,Minus_Combine_Bias_log)

        #roundN = 4

        #### assign bias to bp and proportion
        for outpos in range(Cspan,(end-start+Cspan)):
            this_plus_single = Plus_Single_fxrBias[outpos]
            this_minus_single = Minus_Single_fxrBias[outpos]
            this_plus_sum = sum(Plus_Single_fxrBias[(outpos-Cspan):(outpos+Cspan)])
            this_minus_sum = sum(Minus_Single_fxrBias[(outpos-Cspan):(outpos+Cspan)])

            #this_plus_single_prop = round(Plus_Single_Bias_linear[outpos]/sum(Plus_Single_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)
            #this_plus_cb_prop = round(Plus_Combine_Bias_linear[outpos]/sum(Plus_Combine_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)
            #this_minus_single_prop = round(Minus_Single_Bias_linear[outpos]/sum(Minus_Single_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)
            #this_minus_cb_prop = round(Minus_Combine_Bias_linear[outpos]/sum(Minus_Combine_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)

            this_plus = plus_vector[outpos]
            this_minus = minus_vector[outpos]
            this_plus_cuts_sum = sum(plus_vector[(outpos-Cspan):(outpos+Cspan)])
            this_minus_cuts_sum = sum(minus_vector[(outpos-Cspan):(outpos+Cspan)])

            out_chrm = chrm
            out_start = start + outpos - Cspan
            out_end = out_start+1
            
            if this_plus > 0:
                expcut_plus = this_plus_cuts_sum * (this_plus_single/this_plus_sum)
                outf.write("\t".join(map(str,[out_chrm+"_"+str(out_start)+"_"+str(out_end)+"_+",this_plus,format(expcut_plus,".3e") ]))+"\n")

            if this_minus > 0:
                expcut_minus = this_minus_cuts_sum * (this_minus_single/this_minus_sum)
                outf.write("\t".join(map(str,[out_chrm+"_"+str(out_start)+"_"+str(out_end)+"_-",this_minus,format(expcut_minus,".3e") ]))+"\n")

            
            #outf.write("\t".join(map(str,[out_chrm+":"+str(out_start)+"-"+str(out_end),this_plus_single,this_plus_sum,this_minus_single,this_minus_sum]))+"\n")

            #outf_rawPlus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_single]))+"\n")
            #outf_rawMinus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_single]))+"\n")
            #outf_cbPlus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_cb]))+"\n")
            #outf_cbMinus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_cb]))+"\n")
            #
            #outf_rawPlus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_single_prop]))+"\n")
            #outf_rawMinus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_single_prop]))+"\n")
            #outf_cbPlus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_cb_prop]))+"\n")
            #outf_cbMinus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_cb_prop]))+"\n")
    outf.close()
    #outf_rawPlus.close()
    #outf_rawMinus.close()
    #outf_cbPlus.close()
    #outf_cbMinus.close()
    #outf_rawPlus_prop.close()
    #outf_rawMinus_prop.close()
    #outf_cbPlus_prop.close()
    #outf_cbMinus_prop.close()

    inf.close()




def bias_scan_sob(peak,outfile,w_plus,w_minus,biasMat,Cspan,Gen):

#    t = time.time()
    genome = twobitreader.TwoBitFile(Gen)
    
    BGraw,BGenc,Nmer = readBG(biasMat)
    
    inf = open(peak)
    w_plus_H=BigWigFile(open(w_plus, 'rb'))
    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    outf = open(outfile,'w')
    
    for line in inf:
        ll = line.split()
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])

        plus_obj = w_plus_H.summarize(chrm,start-Cspan,end+Cspan,(end-start+2*Cspan))
        minus_obj = w_minus_H.summarize(chrm,start-Cspan,end+Cspan,(end-start+2*Cspan))
        if not plus_obj  :
            plus_vector = numpy.array([0]*(end-start+2*Cspan))
        else:
            plus_vector = plus_obj.sum_data 
        if not minus_obj  :
            minus_vector = numpy.array([0]*(end-start+2*Cspan))
        else:
            minus_vector = minus_obj.sum_data 
    
        plus_single_bias_raw_vector = []
        plus_single_bias_enc_vector = []
        #plus_cb_bias_vector = []
        minus_single_bias_raw_vector = []
        minus_single_bias_enc_vector = []
        #minus_cb_bias_vector = []

        for pos in range(start-Cspan,end+Cspan):
            plus_seq_tmp = genome[chrm][(pos-6):(pos+11)].upper()
#            plus_reverse_seq = rev(genome[chrm][(pos+offset-flank):(pos+offset+flank)].upper())
            minus_seq_tmp = rev(genome[chrm][(pos-11+1):(pos+6+1)].upper())
#            minus_reverse_seq = genome[chrm][(pos-offset+1-flank):(pos-offset+1+flank)].upper()
            if len(plus_seq_tmp) != 17 or len(minus_seq_tmp) != 17:
                continue

            plus_seq = SOBfetchSEQ(plus_seq_tmp,Nmer)
            minus_seq = SOBfetchSEQ(minus_seq_tmp,Nmer)
            ### calculate bias value for each bases
            
            if len(plus_seq) == Nmer and not "N" in plus_seq:
                plus_bias_raw = BGraw[plus_seq]
                plus_bias_enc = BGenc[plus_seq]
            else:
                plus_bias_raw = 0
                plus_bias_enc = 0
                
            #if len(plus_reverse_seq) == Nmer and not "N" in plus_reverse_seq:
            #    plus_reverse_bias = pow(numpy.e,BG[plus_reverse_seq])
            #else:
            #    plus_reverse_bias = 0
            
            if len(minus_seq) == Nmer and not "N" in minus_seq:
                minus_bias_raw = BGraw[minus_seq]
                minus_bias_enc = BGenc[minus_seq]
            else:
                minus_bias_raw = 0
                minus_bias_enc = 0

            #if len(minus_reverse_seq) == Nmer and not "N" in minus_reverse_seq:
            #    minus_reverse_bias = pow(numpy.e,BG[minus_reverse_seq])
            #else:
            #    minus_reverse_bias = 0

            #if plus_bias != "NA" and plus_reverse_bias != "NA":
            #    plus_cb_bias = (plus_bias + plus_reverse_bias ) / 2
            #else:
            #    plus_cb_bias = -10
            
            #if minus_bias != "NA" and minus_reverse_bias != "NA":
            #    minus_cb_bias = (minus_bias + minus_reverse_bias ) / 2
            #else:
            #    minus_cb_bias = -10

            plus_single_bias_raw_vector.append(plus_bias_raw)
            plus_single_bias_enc_vector.append(plus_bias_enc)
            #plus_cb_bias_vector.append(plus_cb_bias)
            minus_single_bias_raw_vector.append(minus_bias_raw)
            minus_single_bias_enc_vector.append(minus_bias_enc)
            #minus_cb_bias_vector.append(minus_cb_bias)

        #if "NA" in plus_single_bias_vector or "NA" in plus_cb_bias_vector or "NA" in minus_single_bias_vector or "NA" in minus_cb_bias_vector:
        #    print ll
        #    continue        
        #### construct bias vector and transform to linear
        Plus_Single_rawBias = numpy.array(plus_single_bias_raw_vector)
        Plus_Single_encBias = numpy.array(plus_single_bias_enc_vector)
        #Plus_Combine_Bias_log = numpy.array(plus_cb_bias_vector)
        Minus_Single_rawBias = numpy.array(minus_single_bias_raw_vector)
        Minus_Single_encBias = numpy.array(minus_single_bias_enc_vector)
        #Minus_Combine_Bias_log = numpy.array(minus_cb_bias_vector)
##### now here
        #Plus_Single_Bias_linear = pow(numpy.e,Plus_Single_Bias_log)
        #Plus_Combine_Bias_linear = pow(numpy.e,Plus_Combine_Bias_log)
        #Minus_Single_Bias_linear = pow(numpy.e,Minus_Single_Bias_log)
        #Minus_Combine_Bias_linear = pow(numpy.e,Minus_Combine_Bias_log)

        #roundN = 4

        #### assign bias to bp and proportion
        for outpos in range(Cspan,(end-start+Cspan)):
            this_plus_single_raw = Plus_Single_rawBias[outpos]
            this_minus_single_raw = Minus_Single_rawBias[outpos]
            this_plus_sum_raw = sum(Plus_Single_rawBias[(outpos-Cspan):(outpos+Cspan)])
            this_minus_sum_raw = sum(Minus_Single_rawBias[(outpos-Cspan):(outpos+Cspan)])

            this_plus_single_enc = Plus_Single_encBias[outpos]
            this_minus_single_enc = Minus_Single_encBias[outpos]
            this_plus_sum_enc = sum(Plus_Single_encBias[(outpos-Cspan):(outpos+Cspan)])
            this_minus_sum_enc = sum(Minus_Single_encBias[(outpos-Cspan):(outpos+Cspan)])

            #this_plus_single_prop = round(Plus_Single_Bias_linear[outpos]/sum(Plus_Single_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)
            #this_plus_cb_prop = round(Plus_Combine_Bias_linear[outpos]/sum(Plus_Combine_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)
            #this_minus_single_prop = round(Minus_Single_Bias_linear[outpos]/sum(Minus_Single_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)
            #this_minus_cb_prop = round(Minus_Combine_Bias_linear[outpos]/sum(Minus_Combine_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)
            this_plus = plus_vector[outpos]
            this_minus = minus_vector[outpos]
            this_plus_cuts_sum = sum(plus_vector[(outpos-Cspan):(outpos+Cspan)])
            this_minus_cuts_sum = sum(minus_vector[(outpos-Cspan):(outpos+Cspan)])

            out_chrm = chrm
            out_start = start + outpos - Cspan
            out_end = out_start+1
            
            if this_plus > 0:
                expcut_plus_raw = this_plus_cuts_sum * (this_plus_single_raw/this_plus_sum_raw)
                expcut_plus_enc = this_plus_cuts_sum * (this_plus_single_enc/this_plus_sum_enc)
                outf.write("\t".join(map(str,[out_chrm+"_"+str(out_start)+"_"+str(out_end)+"_+",this_plus,format(expcut_plus_raw,".3e"),format(expcut_plus_enc,".3e") ]))+"\n")

            if this_minus > 0:
                expcut_minus_raw = this_minus_cuts_sum * (this_minus_single_raw/this_minus_sum_raw)
                expcut_minus_enc = this_minus_cuts_sum * (this_minus_single_enc/this_minus_sum_enc)
                outf.write("\t".join(map(str,[out_chrm+"_"+str(out_start)+"_"+str(out_end)+"_-",this_minus,format(expcut_minus_raw,".3e"),format(expcut_minus_enc,".3e") ]))+"\n")

            #outf_rawPlus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_single]))+"\n")
            #outf_rawMinus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_single]))+"\n")
            #outf_cbPlus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_cb]))+"\n")
            #outf_cbMinus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_cb]))+"\n")
            #
            #outf_rawPlus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_single_prop]))+"\n")
            #outf_rawMinus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_single_prop]))+"\n")
            #outf_cbPlus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_cb_prop]))+"\n")
            #outf_cbMinus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_cb_prop]))+"\n")
    outf.close()
    #outf_rawPlus.close()
    #outf_rawMinus.close()
    #outf_cbPlus.close()
    #outf_cbMinus.close()
    #outf_rawPlus_prop.close()
    #outf_rawMinus_prop.close()
    #outf_cbPlus_prop.close()
    #outf_cbMinus_prop.close()

    inf.close()
 


def bias_scan_bagfoot(peak,outfile,w_plus,w_minus,biasMat,Cspan,Gen,shiftbp_plus,shiftbp_minus):

#    t = time.time()
    genome = twobitreader.TwoBitFile(Gen)
    
    BGraw,BGenc,Nmer = readBG(biasMat)
    flank = int(Nmer)/2
    w_plus_H=BigWigFile(open(w_plus, 'rb'))
    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    
    inf = open(peak)
#    w_plus_H=BigWigFile(open(w_plus, 'rb'))
#    w_minus_H=BigWigFile(open(w_minus, 'rb'))
    outf = open(outfile,'w')
    
    for line in inf:
        ll = line.split()
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])

        if start-Cspan < 0:
            continue

        plus_obj = w_plus_H.summarize(chrm,start-Cspan,end+Cspan,(end-start+2*Cspan))
        minus_obj = w_minus_H.summarize(chrm,start-Cspan,end+Cspan,(end-start+2*Cspan))
        if not plus_obj  :
            plus_vector = numpy.array([0]*(end-start+2*Cspan))
        else:
            plus_vector = plus_obj.sum_data 
        if not minus_obj  :
            minus_vector = numpy.array([0]*(end-start+2*Cspan))
        else:
            minus_vector = minus_obj.sum_data 
       
        plus_single_bias_raw_vector = []
        plus_single_bias_enc_vector = []
        #plus_cb_bias_vector = []
        minus_single_bias_raw_vector = []
        minus_single_bias_enc_vector = []
        #minus_cb_bias_vector = []
    
        for pos in range(start-Cspan,end+Cspan):
            plus_seq = genome[chrm][(pos-flank+shiftbp_plus):(pos+flank+shiftbp_plus)].upper()
#            plus_reverse_seq = rev(genome[chrm][(pos+offset-flank):(pos+offset+flank)].upper())
            minus_seq = genome[chrm][(pos-flank+1-shiftbp_minus):(pos+flank+1-shiftbp_minus)].upper()
#            minus_reverse_seq = genome[chrm][(pos-offset+1-flank):(pos-offset+1+flank)].upper()
            #if len(plus_seq) != 2*flank or len(minus_seq) != 2*flank:
            #    continue

            ### calculate bias value for each bases
            
            if len(plus_seq) == Nmer and not "N" in plus_seq:
                plus_bias_raw = BGraw[plus_seq]
                plus_bias_enc = BGenc[plus_seq]
            else:
                plus_bias_raw = 0
                plus_bias_enc = 0
                
            #if len(plus_reverse_seq) == Nmer and not "N" in plus_reverse_seq:
            #    plus_reverse_bias = pow(numpy.e,BG[plus_reverse_seq])
            #else:
            #    plus_reverse_bias = 0
            
            if len(minus_seq) == Nmer and not "N" in minus_seq:
                minus_bias_raw = BGraw[minus_seq]
                minus_bias_enc = BGenc[minus_seq]
            else:
                minus_bias_raw = 0
                minus_bias_enc = 0

            #if len(minus_reverse_seq) == Nmer and not "N" in minus_reverse_seq:
            #    minus_reverse_bias = pow(numpy.e,BG[minus_reverse_seq])
            #else:
            #    minus_reverse_bias = 0

            #if plus_bias != "NA" and plus_reverse_bias != "NA":
            #    plus_cb_bias = (plus_bias + plus_reverse_bias ) / 2
            #else:
            #    plus_cb_bias = -10
            
            #if minus_bias != "NA" and minus_reverse_bias != "NA":
            #    minus_cb_bias = (minus_bias + minus_reverse_bias ) / 2
            #else:
            #    minus_cb_bias = -10

            plus_single_bias_raw_vector.append(plus_bias_raw)
            plus_single_bias_enc_vector.append(plus_bias_enc)
            #plus_cb_bias_vector.append(plus_cb_bias)
            minus_single_bias_raw_vector.append(minus_bias_raw)
            minus_single_bias_enc_vector.append(minus_bias_enc)
            #minus_cb_bias_vector.append(minus_cb_bias)

        #if "NA" in plus_single_bias_vector or "NA" in plus_cb_bias_vector or "NA" in minus_single_bias_vector or "NA" in minus_cb_bias_vector:
        #    print ll
        #    continue        
        #### construct bias vector and transform to linear
        Plus_Single_rawBias = numpy.array(plus_single_bias_raw_vector)
        Plus_Single_encBias = numpy.array(plus_single_bias_enc_vector)
        #Plus_Combine_Bias_log = numpy.array(plus_cb_bias_vector)
        Minus_Single_rawBias = numpy.array(minus_single_bias_raw_vector)
        Minus_Single_encBias = numpy.array(minus_single_bias_enc_vector)
        #Minus_Combine_Bias_log = numpy.array(minus_cb_bias_vector)
##### now here
        #Plus_Single_Bias_linear = pow(numpy.e,Plus_Single_Bias_log)
        #Plus_Combine_Bias_linear = pow(numpy.e,Plus_Combine_Bias_log)
        #Minus_Single_Bias_linear = pow(numpy.e,Minus_Single_Bias_log)
        #Minus_Combine_Bias_linear = pow(numpy.e,Minus_Combine_Bias_log)

        #roundN = 4

        #### assign bias to bp and proportion
        for outpos in range(Cspan,(end-start+Cspan)):
                    
            this_plus_single_raw = Plus_Single_rawBias[outpos]
            this_minus_single_raw = Minus_Single_rawBias[outpos]
            this_plus_sum_raw = sum(Plus_Single_rawBias[(outpos-Cspan):(outpos+Cspan)])
            this_minus_sum_raw = sum(Minus_Single_rawBias[(outpos-Cspan):(outpos+Cspan)])

            this_plus_single_enc = Plus_Single_encBias[outpos]
            this_minus_single_enc = Minus_Single_encBias[outpos]
            this_plus_sum_enc = sum(Plus_Single_encBias[(outpos-Cspan):(outpos+Cspan)])
            this_minus_sum_enc = sum(Minus_Single_encBias[(outpos-Cspan):(outpos+Cspan)])

            #this_plus_single_prop = round(Plus_Single_Bias_linear[outpos]/sum(Plus_Single_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)
            #this_plus_cb_prop = round(Plus_Combine_Bias_linear[outpos]/sum(Plus_Combine_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)
            #this_minus_single_prop = round(Minus_Single_Bias_linear[outpos]/sum(Minus_Single_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)
            #this_minus_cb_prop = round(Minus_Combine_Bias_linear[outpos]/sum(Minus_Combine_Bias_linear[(outpos-Cspan):(outpos+Cspan)]),roundN)

            this_plus = plus_vector[outpos]
            this_minus = minus_vector[outpos]
            this_plus_cuts_sum = sum(plus_vector[(outpos-Cspan):(outpos+Cspan)])
            this_minus_cuts_sum = sum(minus_vector[(outpos-Cspan):(outpos+Cspan)])

            out_chrm = chrm
            out_start = start + outpos - Cspan
            out_end = out_start+1
            
            if this_plus > 0:
                expcut_plus_raw = this_plus_cuts_sum * (this_plus_single_raw/this_plus_sum_raw)
                expcut_plus_enc = this_plus_cuts_sum * (this_plus_single_enc/this_plus_sum_enc)
                outf.write("\t".join(map(str,[out_chrm+"_"+str(out_start)+"_"+str(out_end)+"_+",this_plus,format(expcut_plus_raw,".3e"),format(expcut_plus_enc,".3e") ]))+"\n")

            if this_minus > 0:
                expcut_minus_raw = this_minus_cuts_sum * (this_minus_single_raw/this_minus_sum_raw)
                expcut_minus_enc = this_minus_cuts_sum * (this_minus_single_enc/this_minus_sum_enc)
                outf.write("\t".join(map(str,[out_chrm+"_"+str(out_start)+"_"+str(out_end)+"_-",this_minus,format(expcut_minus_raw,".3e"),format(expcut_minus_enc,".3e") ]))+"\n")

            #outf_rawPlus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_single]))+"\n")
            #outf_rawMinus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_single]))+"\n")
            #outf_cbPlus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_cb]))+"\n")
            #outf_cbMinus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_cb]))+"\n")
            #
            #outf_rawPlus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_single_prop]))+"\n")
            #outf_rawMinus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_single_prop]))+"\n")
            #outf_cbPlus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_cb_prop]))+"\n")
            #outf_cbMinus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_cb_prop]))+"\n")
    outf.close()
    #outf_rawPlus.close()
    #outf_rawMinus.close()
    #outf_cbPlus.close()
    #outf_cbMinus.close()
    #outf_rawPlus_prop.close()
    #outf_rawMinus_prop.close()
    #outf_cbPlus_prop.close()
    #outf_cbMinus_prop.close()

    inf.close()
# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--interval",dest="interval",type="str",
                         help="5 column summit200bp file ,sorted")
    optparser.add_option("-o","--outname",dest="outname",type="str",
                         help="name of output file")
    optparser.add_option("-b","--biasMat",dest="bgmatrix",type="str",default = "/scratch/sh8tv/Project/scATAC/Data/Summary_Data/bias_matrix/summary36bp/singleMat/NakedYeast_ATAC_Enc8mer.txt",
                         help="bias matrix in giver N-mer, default is 8-mer estimated from Yeast naked DNA (simplex encoding)")
    optparser.add_option("-t","--biastype",dest="biastype",type="choice",choices=["sob","flank","fxr","bag44","bag45","bag55","repFoot"],default = "flank",
                         help="bias type, choices from flank, sob and fxr")
    optparser.add_option("-p","--plusbw",dest="plusbw",type="str",
                         help="plus cuts in bw format")
    optparser.add_option("-n","--minusbw",dest="minusbw",type="str",
                         help="minus cuts in bw format")

#========minor options=============
    optparser.add_option("--Cspan",dest="Cspan",type="int",default = 25,
                         help="region for get total signal in single bp, default = 25 means +-25bp(total 50bp) signal as total for each bp")
    optparser.add_option("--genome",dest="genome",type="str",default = "/sfs/qumulo/qproject/CPHG/ZANG/sh8tv/Data/Genome/hg38/hg38.2bit",
                         help="genome sequence in 2bit format")
    optparser.add_option("--offset",dest="offset",type="int",default = 9,
                         help="offset related, distance of pair of +/- related cut,default = 9")


    (options,args) = optparser.parse_args()

    if not options.interval:
        optparser.print_help()
        sys.exit(1)

    interval = options.interval
    outname = options.outname
    bgmatrix = options.bgmatrix
    gen = options.genome
    Cspan = options.Cspan
    offset = options.offset
    biastype = options.biastype
    plusbw = options.plusbw
    minusbw = options.minusbw

    if biastype == "flank":
        bias_scan_flank(interval,outname,plusbw,minusbw,bgmatrix,Cspan,gen)
    if biastype == "repFoot":
        bias_scan_flank(interval,outname,plusbw,minusbw,bgmatrix,Cspan,gen)
    if biastype == "sob":
        bias_scan_sob(interval,outname,plusbw,minusbw,bgmatrix,Cspan,gen)
    if biastype == "fxr":
        bias_scan_fxr(interval,outname,plusbw,minusbw,bgmatrix,Cspan,gen,offset)
    if biastype == "bag44":
        bias_scan_bagfoot(interval,outname,plusbw,minusbw,bgmatrix,Cspan,gen,4,4)
    if biastype == "bag45":
        bias_scan_bagfoot(interval,outname,plusbw,minusbw,bgmatrix,Cspan,gen,4,5)
    if biastype == "bag55":
        bias_scan_bagfoot(interval,outname,plusbw,minusbw,bgmatrix,Cspan,gen,5,5)

    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

