#!/bin/bash


scriptdir=/nas/HARPONE/QTL_Analysis/mQTL/scripts
plink=/nas/software

wdir=$1
gwas_file=$2
species=$3
pgen_file=$4
bed_file=$5

cd $wdir/GWAS/$species
cp $scriptdir/calculate_PCs.sh .

awk '{print $1,$1}' Pheno.$species.txt | awk 'NR>1' > $species.txt
cat <(echo IID) <(awk '{print $1}' Pheno.$species.txt | awk 'NR>1') > IIDs_total.$species.txt

bash calculate_PCs.sh $species.txt
cp $species/PCs-$species.txt .
rm -r $species

### UNIR PCs a matriz de covariables. Este comando de awk permite sort saltandose el header. Arreglar este join (y coger IIDs de los que no tengan NA en el pheno)
join -1 1 -2 1 \
<(fgrep -wf IIDs_total.$species.txt $wdir/Covs_noPCs_$gwas_file |  awk 'NR<2{print $0;next}{print $0| "sort -r"}') \
<(fgrep -wf IIDs_total.$species.txt PCs-$species.txt |  awk 'NR<2{print $2,$3,$4,$5,$6;next}{print $2,$3,$4,$5,$6| "sort -r"}') \
| sed 's/\ /\t/g' > Covs.$species.txt

for x in {1..22}; do
   $plink/plink2 \
   --pfile ${pgen_file}_chr${x} \
   --glm no-firth hide-covar cols=+a1freq --covar-variance-standardize Age_LP \
   --pheno Pheno.$species.txt \
   --covar Covs.$species.txt \
   --mac 20 \
   --out chr$x
 done

if [[ -f $species.glm.linear ]]; then rm $species.glm.linear; fi

 CHR_TO_PROCESS="$(seq 1 22)"
 THREADS=6
 /nas/software/parallel -j ${THREADS} grep -w ADD chr{}.$species.glm.linear >> $species.glm.linear ::: ${CHR_TO_PROCESS}

 # Poner header y ordenar
 cat <(head -1 chr22.$species.glm.linear  | sed 's/#//g') <(sort -nk1 -nk2 $species.glm.linear) > tmp; mv tmp $species.glm.linear
 #Poner alelo 2
 paste <(cut -f 1-6 $species.glm.linear) <(awk '{if ($6 == $4) print $5; else if ($6 == $5) print $4; else if ($6 == "A1") print "A2"}' $species.glm.linear) <(cut -f 7- $species.glm.linear) > tmp; mv tmp $species.glm.linear
 
 # Clumping. Default parameters for clump-p1, clump-p2, clump-r2.
 
 for x in {1..22}; do 
   $plink/plink \
   --bfile ${bed_file}_chr${x} \
   --clump chr$x.$species.glm.linear \
   --clump-snp-field ID \
   --clump-p1 0.000001 \
   --clump-p2 0.00001 \
   --clump-r2 0.001 \
   --clump-kb 250 \
   --out Clumped_chr$x.$species.glm.linear
 done
 
 # Join all clumps
 
 head -1 Clumped_chr22.$species.glm.linear.clumped > Clumped.$species.glm.linear
 for x in {1..22}; do awk 'NR>1' Clumped_chr$x.$species.glm.linear.clumped >> Clumped.$species.glm.linear; done
 sed '/^[[:space:]]*$/d' Clumped.$species.glm.linear > tmp; mv tmp Clumped.$species.glm.linear # Eliminar lineas vacias
 
 ## Remove intermediate files
 for x in {1..22}; do rm chr$x.$species.glm.linear; done
 for x in {1..22}; do rm Clumped_chr$x.$species.glm.linear.clumped Clumped_chr$x.$species.glm.linear.log Clumped_chr$x.$species.glm.linear.nosex; done
 
