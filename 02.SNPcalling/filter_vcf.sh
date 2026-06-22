chr=$1

perl VCF_filter_PL.pl \
--i $chr.bcf --o $chr.vcf.gz \
--P popmap \
--c 12 --C 0.9 \
--m 4 --avgDP 4:40 \
--l 0.2 --g 0.05 \
--H 0.85 --f \
--T 16 --s
