ref="Larimichthys_crocea_chr.fa"
angsd="/opt/bio/populations/angsd/angsd"
chr=$1
minIND=120

$angsd -out $chr -b bam.list -rf $chr.good.txt \
-nThreads 8 -minInd $minIND -C 50 \
-uniqueOnly 1 -minMapQ 20 -only_proper_pairs 0 \
--ignore-RG 0 -ref $ref -anc $ref \
-doBcf 1 -docounts 1 -dogeno 1 \
-geno_minDepth 3 -geno_maxDepth 40 \
-gl 1 -dopost 1 -domajorminor 1 -domaf 1 

