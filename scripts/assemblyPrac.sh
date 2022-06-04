1. Genome Assembly
Practical

In this practical, we will perform the assembly of Escherichia coli str. K-12 substr. W3110. Data is obtained from the study published in 2020 by Bleichert P et al., in Appl Environ Microbiol, 2020 Dec 17;87(1).

Getting the data

Escherichia coli was sequenced using the MiSeq platform (2 x 250bp). The reads were deposited in the NCBI SRA database under the accession SRR12793243.

You need to have SRAToolkit program installed in order to download the fastq files which are associated with the accession number. Once you have installed SRAToolkit, you can run the following command on the command line.

## fastq-dump -I SRR12793243  --split-files  (just for your information PLEASE don't run that command !!)

This command will download two files: SRR12793243_1.fastq.gz and SRR12793243_2.fastq.gz. For the sake of easiness, we have already downloaded the data.

Check the files by typing (or better copy&paste) the following command: 
ls /home/Epicass_workshop/Day3/Data/illumina/

Question

How many reads are in the files?
 

Data Quality Check

The downloaded fastq files are raw data, therefore we have to do some quality checks and quality trimming of the data ourselves!

Now we create indir and XX variables assigning the PATH and the accession number to it. By doing this, we can provide the input files in a short and easy way. 

indir="/home/Epicass_workshop/Day3/Data/illumina/"

XX="SRR12793243"

Now run the FastQC program to see the quality of raw data. You may store the output in /home/Epicass_workshop/Day3/myResults directory. We must create the outdir variable and then create the directories:

outdir="/home/Epicass_workshop/Day3/myResults"

mkdir ${outdir}

mkdir ${outdir}/fastqc

mkdir ${outdir}/fastp



fastqc -t 4 ${indir}/${XX}_*.fastq.gz -o ${outdir}/fastqc

Please check the results of the above command by downloading the *.html file to your local computer and looking at it with your browser. 

Now we will run fast to trim and filter the bad quality reads with the following command.

fastp -i ${indir}/${XX}_1.fastq.gz -I ${indir}/${XX}_2.fastq.gz -o ${outdir}/fastp/${XX}_1trim.fastq.gz -O ${outdir}/fastp/${XX}_2trim.fastq.gz -j ${outdir}/fastp/${XX}_fastp.json  -h  ${outdir}/fastp/ ${XX}_fastp.html --thread 4 --trim_poly_g -l 250 -q 20 -n 0

The quality-filtered files are stored as SRR12793243_1trim.fastq.gz and SRR12793243_2trim.fastq.gz in the fastp folder under results directory.

Run the FastQC again on the quality trimmed data and you will see quality metrics of the trimmed data.

fastqc -t 4 ${outdir}/fastp/${XX}_*trim.fastq.gz -o ${outdir}/fastqc

You can check the output by downloading the *.html files in your local computer and open them with your browser.

Now that we have performed QC and filtering of the data, we can now start the genome assembly.




de novo assembly with Illumina short paired-end reads
We will be using the Spades de novo assembler to assemble our bacterium. Run the following command to start the assembly: 

spades.py --careful -m 16 -t 8 -k 21,33,55,77,99,127 -1 ${outdir}/fastp/${XX}_1trim.fastq.gz -2 ${outdir}/fastp/${XX}_2trim.fastq.gz -o ${outdir}/asm_spades

This will take approximately 10 minutes to complete the assembly.  

The result of the assembly is in the path /Epicass_workshop/day3/results/asm_spades containing scaffolds.fasta and  assembly_graph_with_scaffolds.gfa files.


Quality of the Assembly

QUAST is software that evaluates the quality of genome assemblies by computing various metrics.

Run Quast on your assembly by typing:

quast.py ${outdir}/asm_spades/scaffolds.fasta -o ${outdir}/quast_spades_report

and take a look at the text report

cat ${outdir}/quast_spades_report/report.txt

You should see a summary of stats about your assembly.

**Quast Report**

Assembly                    scaffolds
# contigs (>= 0 bp)         116      
# contigs (>= 1000 bp)      62       
# contigs (>= 5000 bp)      53       
# contigs (>= 10000 bp)     49       
# contigs (>= 25000 bp)     45       
# contigs (>= 50000 bp)     30       
Total length (>= 0 bp)      4577414  
Total length (>= 1000 bp)   4558637  
Total length (>= 5000 bp)   4539510  
Total length (>= 10000 bp)  4505373  
Total length (>= 25000 bp)  4435770  
Total length (>= 50000 bp)  3904856  
# contigs                   72       
Largest contig              364193   
Total length                4564990  
GC (%)                      50.77    
N50                         130538   
N75                         70023    
L50                         11       
L75                         23       
# N's per 100 kbp           6.57


Additionally, you can download the file quast_spades_report/report.html and view it in your own web browser.

 

Assembly Completeness

Although quast outputs a range of metrics to assess how contiguous our assembly is, having a long N50 does not guarantee a good assembly, it could be riddled with misassemblies!

We will run busco to try to find marker genes in our assembly. Marker genes are conserved across a range of species and finding intact conserved genes in our assembly would be a good indication of the assembly quality.

First, we need to download and unpack the bacterial datasets used by the busco:

busco --list-datasets

We will use enterobacterales_odb10 dataset to identify how many marker genes are present in our assembly. We can run the busco as following :

busco -i ${outdir}/asm_spades/scaffolds.fasta -l enterobacterales_odb10 -o ${outdir}/busco_spades -m genome

### In case Busco is not working in the VM, please see the pre-calculated results in this directory:

ls /home/Epicass_workshop/Day3/Results/Busco/


Question

How many marker genes have busco found?


de novo assembly using pacbio Hifi long reads
For de novo assembly using long Hifi reads, we will use two different long-read assemblers Flye and Hifi-asm to see which one of them produces better genome assembly.

We already have QC processed data available in /home/Epicass_workshop/Day3/Data/pacbio directory which you should use directly to run a genome assembly. 

Hifi data will take a longer time to finish the assembly. So, for this tutorial, we have reduced the data to 40,000 reads only to finish the assembly process quickly.

However, if you want to process the raw reads later yourself, you can download the data and clean it like this:

## ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR109/019/SRR10971019/SRR10971019_subreads.fastq.gz with wget command and repeat the QC steps as shown below. 

## XX="/change/path/to/read/file/SRR10971019_subreads"

## fastqc -t 4 ${XX}.fastq -o fastqc_result

## Fastp command is slight different for hifi reads:

## fastp -i ${XX}.fastq -W 10 -5 -3 -M 10 -l 10000 -G -A -w 4 -f 10 -j ${XX}_fastp.json -h ${XX}_fastp.html -o ${XX}_clean.fastq

(just for your information PLEASE don't run these commands !!)



Run Flye assembler

You can check the existing processed data here:

ls /home/Epicass_workshop/Day3/Data/pacbio/

indir="/home/Epicass_workshop/Day3/Data/pacbio/"

Now run the Flye assembler:

flye --pacbio-hifi ${indir}/reducedPB_clean.fastq --out-dir ${outdir}/asm_flye  --genome-size 4.5m --threads 8

The assembly files are stored under asm_flye directory with the following main files

assembly.fasta
assembly_graph.gfa


INFO: Assembly statistics:

Total length: 4642499
Fragments: 1
Fragments N50: 4642499
Largest frg: 4642499
Scaffolds: 0
Mean coverage: 31


Run assembly with hifiasm 0.16.1

mkdir -p ${outdir}/asm_hifi_reducedPB

hifiasm -o ${outdir}/asm_hifi_reducedPB/asm_hifi_reducedPB -t8 ${indir}/reducedPB_clean.fastq 2> ${outdir}/asm_hifi_reducedPB/hifiasm.log

perl -ane 'if ($F[0]eq"S") { print ">".$F[1]."\n".$F[2]."\n"}' ${outdir}/asm_hifi_reducedPB/asm_hifi_reducedPB.bp.p_ctg.gfa > ${outdir}/asm_hifi_reducedPB/asm_hifi_reducedPB.p_ctg.fa

Check the assembly output files in /home/Epicass_workshop/Day3/myResults/asm_hifi_reducedPB directory. They should look like:

asm_hifi_reducedPB.p_ctg.fa
asm_hifi_reducedPB.bp.p_ctg.noseq.gfa

plus many other files.


Quality of the Assembly with long reads
We will run Quast again to see assembly statistics. But this time we will provide those assemblies to Quast that we generated with other assemblers, to visualize the quality of all assemblies simultaneously.

Run Quast on Illumina, Flye and hifi-asm assemblies with the following command

res_path="/home/Epicass_workshop/Day3/myResults"

quast.py ${res_path}/asm_spades/scaffolds.fasta  ${res_path}/asm_flye/assembly.fasta ${res_path}/asm_hifi_reducedPB/asm_hifi_reducedPB.p_ctg.fa -o ${res_path}/quast_all_assemblies_report


Now take a look at the text report

cat ${res_path}/quast_all_assemblies_report/report.txt

Download the file ${res_path}/quast_all_assemblies _report/report.html and open it in your own web browser to compare the quality of different assemblies.
 

Assembly completeness of long reads
Run busco on Flye:

busco -i ${res_path}/asm_flye/assembly.fasta -l enterobacterales_odb10 -o ${res_path}/busco_flye -m genome

Run busco on hifi-asm assembly:

busco -i ${res_path}/asm_hifi_reducedPB/asm_hifi_reducedPB.p_ctg.fa -l enterobacterales_odb10 -o ${res_path}/busco_hifi -m genome


Question

How many marker genes have busco found in long-read assemblies?  
 

Visualizing assembly graphs 
Bandage is a program for visualizing de novo assembly graphs. It displays connections that are not visible in the contigs file.


We will now visualize our bacterial assembly graphs which were generated by Spades, Flye, and Hifi-asm during the genome assembly process.


Bandage image ${res_path}/asm_spades/assembly_graph_with_scaffolds.gfa ${res_path}/spades_graph.png

Bandage image ${res_path}/asm_flye/assembly_graph.gfa ${res_path}/flye_graph.png --bases 100000

Bandage image ${res_path}/asm_hifi_reducedPB/asm_hifi_reducedPB.bp.p_ctg.noseq.gfa ${res_path}/hifi_graph.png --bases 100000



This is not working on the VM:

Similarly, you can check the information of each assembly by running the following commands:

Bandage info ${res_path}/asm_spades/assembly_graph_with_scaffolds.gfa

Bandage info ${res_path}/asm_flye/assembly_graph.gfa

Bandage info ${res_path}/asm_hifi_reducedPB/asm_hifi_reducedPB.bp.p_ctg.noseq.gfa



flye assembly graph produced with Bandage Hifi assembly graph - Bandage


