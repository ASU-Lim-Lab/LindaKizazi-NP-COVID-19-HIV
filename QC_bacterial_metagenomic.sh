## For QC of bacterial metagenomic samples 
# Adapted from ViromeWorkflow2020_v4.py
# Input sample names file: can use file created from excel or text editor; if excel used, save in tab delimited text (.txt) file format  
# First line of the input file must be the path to the directory of the run; path to where the fastq files and new directories will be made/populated

import sys
import csv
import re
import string
import os.path
import fileinput
import math

try: 
	input_name=str(sys.argv[1])+"_name_of_script.sh"
except (ValueError,IndexError):
	print("Please enter the file with full sample names, including the .gz if files zipped")
	sys.exit(1)
output=open(input_name,'w')
sampleNames=open("SampleNames_"+sys.argv[1], 'w')
file=sys.argv[1]
# read1=open("read1.txt",'w')
# read2=open("read2.txt", 'w')
names = []

with open(file, 'r') as f: 
	reads=[[],[]] ## creates two lists which have reads1 full name (including .gz) and same for reads2 fastq sample files
	lines=f.readlines()
	i=0
	for each in lines:
		if each.startswith("Path:"):
			path=str(each)
			path=path.strip().split(":")
			path=path[1]
		elif re.match("^[a-zA-Z]+.*", each):
			reads[i].append(each.strip())
			i ^= 1
			Each = each.strip().split("_S")
			# print each 
			# print Each[0]
			SamplesID= "_S"+Each[1]
			# print SamplesID
			Names=Each[0]
			if Names not in names: 
				names.append(Names)
				sampleNames.write(Names+'\n')
c = str(path) 
if c[-1]!= "/":
	print("Please add '/' to the end of path")
	sys.exit(1)

# print names
# print reads[1]
reading1=[]
reading2=[]
for read1 in reads[0]: 
	# read1.write(r1+'\n')
	read1=str(read1)
	read1=read1.strip().split('.gz')
	read1=read1[0]
	reading1.append(read1)
# print reading1
for read2 in reads[1]: 
	# read2.write(r2+'\n')
	read2=str(read2)
	read2=read2.strip().split('.gz')
	read2=read2[0]
	reading2.append(read2)
	# print r2[0]
	# r2=str(r2[0])
	# print r2
	# read1.write(r2+'\n') # creates the reads1 files without the .gz

output.write("#!/bin/bash\n")
#output.write("mkdir "+path+"hostRemoved_2022;\n")
output.write("mkdir "+path+"qc_fasta_2022;\n")
output.write("mkdir "+path+"qc_fastq_2022;\n")
output.write("mkdir "+path+"log_files_2022;\n")
output.write("mkdir "+path+"fastqc_raw;\n")
output.write("mkdir "+path+"fastqc_raw/logFiles;\n")
output.write("mkdir "+path+"fastqc_clean;\n")
output.write("mkdir "+path+"fastqc_clean/logFiles;\n")


try:
	for run,r1,r2 in zip (names, reading1, reading2):
		output.write("fastqc "+path+r1+".gz --outdir="+path+"fastqc_raw/ -t _ 1> "+path+"fastqc_raw/logFiles/"+r1+"_fastqc_raw.log.txt 2>&1;\n")
		output.write("fastqc "+path+r2+".gz --outdir="+path+"fastqc_raw/ -t  1> "+path+"fastqc_raw/logFiles/"+r2+"_fastqc_raw.log.txt 2>&1;\n")
except (ValueError, NameError): 
		print("No path given in sample file")
		sys.exit(1)

try:
	for run,r1,r2 in zip (names, reading1, reading2):
			output.write("gunzip "+path+r1+".gz;\n")
			output.write("gunzip "+path+r2+".gz;\n")
			output.write("mkdir "+path+run+";\n")
			output.write("cutadapt -a file:/path/to/contaminants/fasta.fa -A file:/path/to/contaminants/fasta.fa -o "+path+run+"/"+run+"_cutadapt_adaptTrim_R1_2022.fastq -p "+path+run+"/"+run+"_cutadapt_adaptTrim_R2_2022.fastq "+path+r1+" "+path+r2+" -O 7 -j 88 1> "+path+run+"/"+run+"_cutadapt_adaptTrim_2022.log.txt 2>&1;\n")
			output.write("bbduk.sh in="+path+run+"/"+run+"_cutadapt_adaptTrim_R1_2022.fastq in2="+path+run+"/"+run+"_cutadapt_adaptTrim_R2_2022.fastq out="+path+run+"/"+run+"-cutadapt_QC_R1_2022.fastq out2="+path+run+"/"+run+"-cutadapt_QC_R2_2022.fastq qtrim=rl trimq=30 minlength=75 minavgquality=20 removeifeitherbad=f tpe=t overwrite=t 1> "+path+run+"/"+run+"-cutadapt_QC_2022.log.txt 2>&1;\n")
			output.write("bbduk.sh in="+path+run+"/"+run+"-cutadapt_QC_R1_2022.fastq in2="+path+run+"/"+run+"-cutadapt_QC_R2_2022.fastq "+"ref=/path/to/phix174_ill.ref.fa.gz out="+path+run+"/"+run+"-R1-cutadapt_phixRemoved_2022.fastq out2="+path+run+"/"+run+"-R2-cutadapt_phixRemoved_2022.fastq k=31 hdist=1 overwrite=t 1> "+path+run+"/"+run+"_cutadapt_phixRemoved_2022.log.txt 2>&1;\n")
			output.write("bbmap.sh minid=.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 -Xmx64g path=/path/to/Human in="+path+run+"/"+run+"-R1-cutadapt_phixRemoved_2022.fastq in2="+path+run+"/"+run+"-R2-cutadapt_phixRemoved_2022.fastq outu="+path+run+"/"+run+"_cutadapt_hostRemoved_2022.fastq outm="+path+run+"/"+run+"_cutadapt_hostMatched_2022.fastq 1>"+path+run+"/"+run+"_cutadapt_hostRemoval.log_2022.txt 2>&1;\n")			
			output.write("dedupe.sh in="+path+run+"/"+run+"_cutadapt_hostRemoved_2022.fastq out="+path+run+"/"+run+"_cutadapt_firstDeduplication_2022.fastq outd="+path+run+"/"+run+"_cutadapt_firstDuplication_2022.fastq csf=dedupe.cluster.stats overwrite=t minidentity=99 1> "+path+run+"/"+run+"_cutadapt_firstDeduplication.log_2022.txt 2>&1;\n")
			output.write("bbmerge.sh in="+path+run+"/"+run+"_cutadapt_firstDeduplication_2022.fastq out="+path+run+"/"+run+"_cutadapt_firstDeduplicationMerged_2022.fastq outu="+path+run+"/"+run+"_cutadapt_firstDeduplicationUnMerged_2022.fastq 1>"+path+run+"/"+run+"_cutadapt_firstDeduplicationMerged.log_2022.txt 2>&1;\n")
			output.write("cat "+path+run+"/"+run+"_cutadapt_firstDeduplicationMerged_2022.fastq "+path+run+"/"+run+"_cutadapt_firstDeduplicationUnMerged_2022.fastq > "+path+run+"/"+run+"_cutadapt_firstDeduplicationMerged_UnMerged_2022.fastq;\n")
			output.write("dedupe.sh in="+path+run+"/"+run+"_cutadapt_firstDeduplicationMerged_UnMerged_2022.fastq out="+path+run+"/"+run+"_cutadapt_secondDeduplication_2022.fastq outd="+path+run+"/"+run+"_cutadapt_secondDuplication_2022.fastq csf=dedupe.cluster.stats overwrite=t minidentity=100 ac=f 1> "+path+run+"/"+run+"_cutadapt_secondDeduplication.log_2022.txt 2>&1;\n")
			output.write("bbduk.sh in="+path+run+"/"+run+"_cutadapt_secondDeduplication_2022.fastq out="+path+run+"/"+run+"_cutadapt_secondDeduplication_filtered_2022.fastq minlength=75 overwrite=t 1>"+path+run+"/"+run+"_cutadapt_secondDeduplication_filtered.log_2022.txt 2>&1; \n")
			output.write("sed -n '1~4s/^@/>/p;2~4p' "+path+run+"/"+run+"_cutadapt_secondDeduplication_filtered_2022.fastq > "+path+run+"/"+run+"_cutadapt_secondDeduplication_filtered_2022.fasta;\n")
			output.write("fastqc "+path+run+"/"+run+"_cutadapt_secondDeduplication_filtered_2022.fastq --outdir="+path+"fastqc_clean/ 1> "+path+"fastqc_clean/logFiles/"+run+"_fastqc_clean.log.txt 2>&1;\n")
#			output.write("mv "+path+run+"/"+run+"_cutadapt_hostRemoved_2022.fastq "+path+"hostRemoved_2022;\n")
			output.write("mv "+path+run+"/"+run+"_cutadapt_secondDeduplication_filtered_2022.fastq "+path+"qc_fastq_2022;\n")
			output.write("mv "+path+run+"/"+run+"_cutadapt_secondDeduplication_filtered_2022.fasta "+path+"qc_fasta_2022;\n")
			output.write("mv "+path+run+"/"+run+"_cutadapt_adaptTrim_2022.log.txt "+path+"log_files_2022;\n")
			output.write("mv "+path+run+"/"+run+"-cutadapt_QC_2022.log.txt "+path+"log_files_2022;\n")
			output.write("mv "+path+run+"/"+run+"_cutadapt_phixRemoved_2022.log.txt "+path+"log_files_2022;\n")
			output.write("mv "+path+run+"/"+run+"_cutadapt_hostRemoval.log_2022.txt "+path+"log_files_2022;\n")
			output.write("mv "+path+run+"/"+run+"_cutadapt_firstDeduplication.log_2022.txt "+path+"log_files_2022;\n")
			output.write("mv "+path+run+"/"+run+"_cutadapt_firstDeduplicationMerged.log_2022.txt "+path+"log_files_2022;\n")
			output.write("mv "+path+run+"/"+run+"_cutadapt_secondDeduplication.log_2022.txt "+path+"log_files_2022;\n")
			output.write("mv "+path+run+"/"+run+"_cutadapt_secondDeduplication_filtered.log_2022.txt "+path+"log_files_2022;\n")
			#output.write("blastx -db /scratch/ekaelin/RefSeqPlusNeighborSeq_2020ViralNR/RefSeqPlusNeighbor2020 -query "+path+run+"/"+run+"_cutadapt_secondDeduplication_filtered_2022.fasta -evalue 1e-3 -num_threads 56 -out "+path+run+"/"+run+"_cutadapt_secondDeduplication_filtered.blastx_2022.out;\n")

except (ValueError, NameError): 
		print("No path given in sample file")
		sys.exit(1)
