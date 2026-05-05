#!/bin/bash
for i in Common_Snipe Far_Eastern Great_Knot Grey-tailed Red-necked Ruddy_Turnstone Spoon-billed Terek Bar-tailed
do
        for j in {1..10}
        do
                shuf -n 7 $i > $i.$j
        PopLDdecay -InVCF ../03cleanvcf/use/auto.nomissing.recode.vcf.recode.vcf -SubPop $i.$j -OutStat $i.$j
        gunzip $i.$j.stat.gz
        done
done
