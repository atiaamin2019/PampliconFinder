
<img width="1443" alt="Screenshot 2025-05-18 at 7 56 27 AM" src="https://github.com/user-attachments/assets/6c142b6d-4b50-48af-bcab-207825d9cca4" />



# What is AmpliconFinder

AmpliconFinder is a Python-based bioinformatics pipeline for the automated detection of extrachromosomal circular DNA from long read DNA sequencing data. 

Extrachromosomal circular DNA (ecDNA) amplicons have been investigated across diverse organisms. For example, in fungi, different plants, and parasites. In cancer, ecDNAs are recognized as key drivers of oncogene amplification, tumor evolution, and drug resistance, and potential diagnostic biomarkers. Using the long read DNA sequencing such as Oxford Nanopore Technology, these amplicons can be identified within a single read. We built AmpliconFinder to use long reads to identify and characterize the amplicons. This is a Python based bioinformatics pipeline, where several bioinformatics tools are used such as [minimap2](https://github.com/lh3/minimap2), [samtools](https://github.com/samtools), [bedtoos](https://github.com/arq5x/bedtools2), [BLASTn](https://github.com/enormandeau/ncbi_blast_tutorial), [IGV](https://github.com/igvteam/igv), [metaflye](https://github.com/mikolmogorov/Flye), [FeatureCount](https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html), [CD-hit](https://github.com/weizhongli/cdhit), [Seqkit](https://github.com/shenwei356/seqkit).

**Python Packages Used**
• [BioPython](https://github.com/biopython/biopython) ( Key modules used:
SeqIO, SeqRecord, Seq)
• [Seaborn](https://github.com/mwaskom/seaborn) (Used for creating
statistical plots such as histograms,
violin plots)
• Pandas (Used for data analysis)
• Matplotlib (pyplot was used for
visualizing genome-wide coverage)
• Numpy (For data analysis)
• [Pysam](https://github.com/pysam-developers/pysam) (For coverage calculation)
• [Bioinfokit](https://github.com/reneshbedre/bioinfokit) (Used for drawing
Manhattan plots)

AmpliconFinder pipeline can be run in two ways. One way is a targeted approach where the target gene within the amplicon is known. Two in-silico primer like DNA fragments are required to run this approach.
And the second approach is for identifying the amplicons throughout the whole genome. This approach requires clustering of the ONT reads and de novo assembly of the amplicons.  

### Targeted Approach for Amplicon Identification:
![Screenshot 2025-05-13 at 4 06 52 PM](https://github.com/user-attachments/assets/f14991b4-919a-4a96-9f51-3a054edb4902)

**Step1**: The following script takes as input the raw ONT data as ONT_reads.fastq.gz files, identifies the reads carrying amplicon junctions, and trim the reads for the probe start site and end site. It generates an output fasta file with the trimmed reads that starts with one probe and end with another probe. 

```
python identifying_junction_reads.py
```


**Step2**: Run SeqKit to get the trimmed read length from the output fasta file of the previous step. 

```
seqkit fx2tab --length --name --header-line  trimmed_amplicon_reads.fasta > trimmed_read_length.txt
```

The resulting trimmed_read_length.txt file will have two columns, first one carrying the read ID and the second one has the trimmed read length. A histogram of the read length will show the junction distribution. 

**Step3**: Run the following script which takes as input the raw ONT_reads.fastq.gz file and the trimmed_read_length.txt file from Seqkit and gives as output the full amplicon reads in a separate fasta file. This output fasta file can be used for amplicon assembly pipeline. 

```
python extracting_amplicon_reads.py
```

**Step4**: Alignment of the extracted fasta file to the reference genome using minimap2.

```
minimap2 -a reference.fasta extracted_reads.fasta > alignment_output.sam
```

**Step5**: Alignment output needs to be modified and filtered using Samtools

```
samtools view -bS alignment_output.sam | samtools view -h -F 0x900 - | samtools sort -o alignment_output_sorted.bam - && samtools index alignment_output_sorted.bam 
```

### Whole genome based approach for amplicon identification
![Screenshot 2025-05-13 at 4 00 11 PM](https://github.com/user-attachments/assets/1f0ab64c-35e6-4ba7-9bb3-c4bd914d511e)

**Step1**: From the alignment output BAM files calculate the normalized mean coverage using the following script.

```
python calculate_coverage_bedtools.py 
```

**Step2**: Use the following script to process the output.csv file from the previous step to process the bedtools coverage output.

```
python bedtool_cov_data_modification.py
```

**Step3**: Use the following script to process the normalized mean coverage and generate a histogram showing the coverage increase.

```
python ploidy_calculation.py
```
