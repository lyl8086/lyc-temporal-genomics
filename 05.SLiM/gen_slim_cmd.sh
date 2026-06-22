:>run_slim.cmd
for i in chr{01..24};do
    echo slim -d \"vcf=\'VCF/$i.vcf\'\" -d \"geno=\'CHR/$i.fa\'\" -d \"chr=\'$i\'\" SLiM.final.txt >>run_slim.cmd
done
