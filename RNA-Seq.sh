cd /DATA4T/Brassica_napus_RNA_Seq/0genome
hisat2-build youcai.genome.fasta youcai.genome.fasta -p 10
cd ../
genome_index="/DATA4T/Brassica_napus_RNA_Seq/0genome/youcai.genome.fasta"
for fq in $(ls 1fastq/*_R1.fastq.gz); do prefix=$(basename 1fastq/$fq |sed 's/_R1.fastq.gz//g'| uniq ); if [ ! -e 2hisat2/${prefix}.sam ]; then echo $fq; hisat2 -x $genome_index -1 1fastq/${prefix}_R1.fastq.gz -2 1fastq/${prefix}_R2.fastq.gz -S 2hisat2/${prefix}.sam  --dta-cufflinks  --no-mixed --no-discordant --new-summary --summary-file 2hisat2/${prefix}.txt -p 12;fi;done

### featurecounts
for m in $(ls 2hisat2/*.sam); do echo $m; prefix=$(basename ${m%.sam}); if [ ! -e 3featurecount/${prefix} ]; then featureCounts -a 0genome/youcai.genome.gff -F GTF -p -g Parent -o 3featurecount/$prefix $m -T 10; fi; done

# stringtie
mkdir 3sorted_bam
for m in $(ls 2hisat2/*.sam); do echo $m; samtools view -u $m | samtools sort -o 3sorted_bam/$(basename ${m%.sam}).sorted.bam -@ 10; done
for m in $(ls 3sorted_bam/*.bam); do echo $m; stringtie $m -G 0genome/youcai.genome.gff -e -B -A 4stringtie/$(basename ${m%.sorted.bam}).genes.tsv -o 4stringtie/$(basename ${m%.sorted.bam}).gff -p 10; done

