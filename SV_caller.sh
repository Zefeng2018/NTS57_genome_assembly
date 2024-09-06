cd /public1/home/GSAUWuzf/MyResearch/1Brassica_napus/0pangenome/minimap2NTS57  # 
ssh xnode20
vi minimap2.sh
for sp in Darmor ganganF73 no2127 quintaA shengli3 tapidor3 zheyou73 westar ZS11
do
   minimap2 -ax  asm5 --eqx -t 40 ../0genome/youcai.chr.fasta ../0genome/${sp}.chr.fasta | samtools view -b > ${sp}2NTS57.bam 
done
################## run syri ###############################
conda activate syri
cd /public1/home/GSAUWuzf/MyResearch/1Brassica_napus/0pangenome/minimap2NTS57 
for sp in Darmor ganganF73 no2127 quintaA shengli3 tapidor3 zheyou73 westar zs11
do 
    syri -c ${sp}2NTS57.bam -r ../0genome/youcai.chr.fasta  -q ../0genome/${sp}*.chr.fasta -k -F B --prefix ${sp}2NTS57 --nc 40  # -F B bam格式输入
done

(syri) [GSAUWuzf@xnode17 minimap2NTS57]$ nohup sh minimap2.sh &
[1] 107400

# SV to vcf
## INV
#! /bin/sh
for i in Darmor ganganF73 no2127 quintaA shengli3 tapidor3 zheyou73 westar zs11
do
cat ${i}2NTS57invOut.txt | grep '#' > ${i}2NTS57_INV.xls
        rm ${i}2NTS57_syri_INV.vcf
        cat ${i}2NTS57_INV.xls | while read n
        do
                chr=`echo $n | awk '{print $2}'`
                pos=`echo $n | awk '{print $3}'`
                ref_pos=`echo $n | awk '{print $2":"$3"-"$4}'`
                alt_pos=`echo $n | awk '{print $6":"$7"-"$8}'`
                ref_seq=`samtools faidx ../0genome/youcai.chr.fasta $ref_pos | grep -v '>' | sed ':a;N;s/\n//g;ta'`
                alt_seq=`samtools faidx ../0genome/${i}*chr.fasta $alt_pos | grep -v '>' | sed ':a;N;s/\n//g;ta'`
                echo -e "${chr}\t${pos}\t.\t${ref_seq}\t${alt_seq}\t30\tPASS\t.\tGT\t1/1" >> convert_vcf/${i}_NTS57_syri_INV.vcf
        done
done
nohup sh syri2vcf.sh &


### DEL INS
for i in Darmor ganganF73 no2127 quintaA shengli3 tapidor3 zheyou73 westar zs11
do
cat ${i}2NTS57*syri.vcf|egrep 'INS|DEL'|awk 'length($4)-length($5)>50' > convert_vcf/${i}_NTS57_syri_DEL.vcf
cat ${i}2NTS57*syri.vcf|egrep 'INS|DEL'|awk 'length($5)-length($4)>50' > convert_vcf/${i}_NTS57_syri_INS.vcf
done


## TRANS
for i in Darmor ganganF73 no2127 quintaA shengli3 tapidor3 zheyou73 westar zs11
do
./output_translocation_from_vcf.py ${i}2NTS57*syri.vcf> ${i}_NTS57_TRANS.xls
        rm ${i}_NTS57_syri_TRANS.vcf
        cat ${i}_NTS57_TRANS.xls | while read n
        do
                chr=`echo $n | awk '{print $2}'`
                pos=`echo $n | awk '{print $3}'`
                ref_pos=`echo $n | awk '{print $2":"$3"-"$4}'`
                alt_pos=`echo $n | awk '{print $6":"$7"-"$8}'`
                ref_seq=`samtools faidx ../0genome/youcai.chr.fasta $ref_pos | grep -v '>' | sed ':a;N;s/\n//g;ta'`
                alt_seq=`samtools faidx ../0genome/${i}*chr.fasta $alt_pos | grep -v '>' | sed ':a;N;s/\n//g;ta'`
                echo -e "${chr}\t${pos}\t.\t${ref_seq}\t${alt_seq}\t30\tPASS\t.\tGT\t1/1" >> convert_vcf/${i}_NTS57_syri_TRANS.vcf
        done
        rm ${i}_NTS57_TRANS.xls
done

## DUP
for i in Darmor ganganF73 no2127 quintaA shengli3 tapidor3 zheyou73 westar zs11
do
cat ${i}2NTS57dupOut.txt | grep '#' > ${i}2NTS57_dup.xls
        rm ${i}2NTS57_syri_dup.vcf
        cat ${i}2NTS57_dup.xls | while read n
        do
                chr=`echo $n | awk '{print $2}'`
                pos=`echo $n | awk '{print $3}'`
                ref_pos=`echo $n | awk '{print $2":"$3"-"$4}'`
                alt_pos=`echo $n | awk '{print $6":"$7"-"$8}'`
                ref_seq=`samtools faidx ../0genome/youcai.chr.fasta $ref_pos | grep -v '>' | sed ':a;N;s/\n//g;ta'`
                alt_seq=`samtools faidx ../0genome/${i}*chr.fasta $alt_pos | grep -v '>' | sed ':a;N;s/\n//g;ta'`
                echo -e "${chr}\t${pos}\t.\t${ref_seq}\t${alt_seq}\t30\tPASS\t.\tGT\t1/1" >> convert_vcf/${i}_NTS57_syri_dup.vcf
        done
done
# circos 做图
bedtools coverage -a circos/genome.windows -b DEL.bed | cut -f 1-4 > circos/DEL_num.txt
bedtools coverage -a circos/genome.windows -b INS.bed | cut -f 1-4 > circos/INS_num.txt
bedtools coverage -a circos/genome.windows -b INV.bed | cut -f 1-4 > circos/INV_num.txt
bedtools coverage -a circos/genome.windows -b TRANS.bed | cut -f 1-4 > circos/TRANS_num.txt
bedtools coverage -a circos/genome.windows -b DUP.bed | cut -f 1-4 > circos/DUP_num.txt

# 统计SV大小，基于上述xls文件统计
for i in Darmor ganganF73 no2127 quintaA shengli3 tapidor3 zheyou73 westar zs11; do cat ${i}2NTS57*syri.vcf|egrep 'INS|DEL'|awk 'length($4)-length($5)>50' | awk 'OFS="\t"{print $1,$2,$2+length($4)-length($5)}'> ${i}_NTS57_DEL.xls;done
for i in Darmor ganganF73 no2127 quintaA shengli3 tapidor3 zheyou73 westar zs11; do cat ${i}2NTS57*syri.vcf|egrep 'INS|DEL'|awk 'length($5)-length($4)>50' | awk 'OFS="\t"{print $1,$2,$2+length($5)-length($4)}'> ${i}_NTS57_INS.xls;done
