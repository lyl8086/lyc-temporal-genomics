
vcf=$1
rep=20
thr=80
cb=$2

out=${vcf%%.*}
plink --vcf $vcf --real-ref-alleles --geno 0.1 --mac 2 --hwe 0.01 --recode 'vcf-iid' --out $out --allow-extra-chr --autosome-num 77
mkdir -p $out

:>$out.GONE2.cmd
for i in `seq 1 $rep`; do
    echo " gone2 $3 -r $cb -t $thr -o $out/${out}_$i $out.vcf" >>$out.GONE2.cmd
done

echo "
library(ggpubr)
library(data.table)

x<-fread(\"${out}_1_GONE2_Ne\", header=T)
names(x) <- c(\"Gen\",\"Ne\")
p<-ggplot(data=x[1:min(100, nrow(x)),], aes(x=Gen, y=Ne))+
geom_line(alpha=0.2)+theme_pubr()+
labs(title = \"${out}, cMMb=${cb}\")+
theme(plot.title=element_text(color=\"red\",hjust=0.5))
a<-fread(\"${out}_1_GONE2_Ne\", header=T)
names(a)<-c(\"Gen\",\"rep1\")

for (i in seq(2,${rep})) {
	x<-fread(paste0(\"${out}_\", i, \"_GONE2_Ne\"), header=T)
	names(x)<-c(\"Gen\",\"Ne\")
    a[x, on=\"Gen\", paste0(\"rep\",i) := Ne]
	p<-p+geom_line(data=x[1:min(100, nrow(x)),],aes(x=Gen, y=Ne),alpha=0.2)
}

mm <- function(x) { exp(mean(log(x[x>0]), na.rm=TRUE))}
a[,Ne := apply(a[,2:(${rep}+1), with=FALSE], 1, mm)]
fwrite(a,\"${out}.xls\",sep=\"\t\",quote=FALSE)
p<-p+geom_line(data=a[1:min(100, nrow(a)),],aes(x=Gen, y=Ne),color=\"red\",alpha=0.5, size=1)+
labs(x='Generation',y='Ne')
ggsave(\"${out}.GONE.pdf\",p)
" >$out/${out}.plot.R

