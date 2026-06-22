#!/bin/bash
vcf=$1
pop=$2
nchr=$3

[ $# -lt 3 ] && echo $0 [vcf] [pop] [nchr] && exit 0

for i in `awk '{print $2}' $pop | sort -u`; do
	# cut vcf by pop
	awk -v p=$i '$2==p {print $1}' $pop >$i.ind.txt
    
    # sub vcf
    bcftools view --threads 8 -W -S "${i}.ind.txt" "$vcf" -Oz -o "${i}.sub.vcf.gz"
    
	# impute
	java -Xmx256g -jar /opt/bio/populations/beagle.jar \
	gt=$i.sub.vcf.gz out=$i.impute
    bcftools index --threads 8 $i.impute.vcf.gz
    # recovery
    bcftools annotate --threads 8 \
    -a $i.impute.vcf.gz \
    -c FMT/GT -W \
    ${i}.sub.vcf.gz \
    -Oz -o $i.vcf.gz
    # clean
    rm ${i}.sub.vcf.gz 
done
