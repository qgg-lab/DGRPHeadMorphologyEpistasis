bash /var/dgrp2/www/data/freeze2March11/pipeline.sh --plink /var/dgrp2/www/data/freeze2March11/plink --fastlmm /var/dgrp2/www/data/freeze2March11/fastlmmc --ldplot /var/dgrp2/www/data/freeze2March11/LDheatmap.R --annot /var/dgrp2/www/data/freeze2March11/freeze2.biallelic.all.annot.txt --bfile /var/dgrp2/www/data/freeze2March11/freeze2.common --nbfile /var/dgrp2/www/data/freeze2March11/freeze2.common.nonsnp.recoded --rel /var/dgrp2/www/data/freeze2March11/freeze2.common.rel.mat --qqplot /var/dgrp2/www/data/freeze2March11/qqplot.R --adjust /var/dgrp2/www/data/freeze2March11/adjustData.RData --maf 4 --cut 0.00001 --mpheno jingFace.mpheno --fpheno jingFace.fpheno --out jingFace --final jingFaceFinal > jingFace.out 2>&1 &

bash /var/dgrp2/www/data/freeze2March11/pipeline.sh --plink /var/dgrp2/www/data/freeze2March11/plink --fastlmm /var/dgrp2/www/data/freeze2March11/fastlmmc --ldplot /var/dgrp2/www/data/freeze2March11/LDheatmap.R --annot /var/dgrp2/www/data/freeze2March11/freeze2.biallelic.all.annot.txt --bfile /var/dgrp2/www/data/freeze2March11/freeze2.common --nbfile /var/dgrp2/www/data/freeze2March11/freeze2.common.nonsnp.recoded --rel /var/dgrp2/www/data/freeze2March11/freeze2.common.rel.mat --qqplot /var/dgrp2/www/data/freeze2March11/qqplot.R --adjust /var/dgrp2/www/data/freeze2March11/adjustData.RData --maf 4 --cut 0.00001 --mpheno invFace.mpheno --fpheno invFace.fpheno --out invFace --final invFaceFinal > invFace.out 2>&1 &

bash /var/dgrp2/www/data/freeze2March11/pipeline.sh --plink /var/dgrp2/www/data/freeze2March11/plink --fastlmm /var/dgrp2/www/data/freeze2March11/fastlmmc --ldplot /var/dgrp2/www/data/freeze2March11/LDheatmap.R --annot /var/dgrp2/www/data/freeze2March11/freeze2.biallelic.all.annot.txt --bfile /var/dgrp2/www/data/freeze2March11/freeze2.common --nbfile /var/dgrp2/www/data/freeze2March11/freeze2.common.nonsnp.recoded --rel /var/dgrp2/www/data/freeze2March11/freeze2.common.rel.mat --qqplot /var/dgrp2/www/data/freeze2March11/qqplot.R --adjust /var/dgrp2/www/data/freeze2March11/adjustData.RData --maf 4 --cut 0.00001 --mpheno jingHead.mpheno --fpheno jingHead.fpheno --out jingHead --final jingHeadFinal > jingHead.out 2>&1 &

bash /var/dgrp2/www/data/freeze2March11/pipeline.sh --plink /var/dgrp2/www/data/freeze2March11/plink --fastlmm /var/dgrp2/www/data/freeze2March11/fastlmmc --ldplot /var/dgrp2/www/data/freeze2March11/LDheatmap.R --annot /var/dgrp2/www/data/freeze2March11/freeze2.biallelic.all.annot.txt --bfile /var/dgrp2/www/data/freeze2March11/freeze2.common --nbfile /var/dgrp2/www/data/freeze2March11/freeze2.common.nonsnp.recoded --rel /var/dgrp2/www/data/freeze2March11/freeze2.common.rel.mat --qqplot /var/dgrp2/www/data/freeze2March11/qqplot.R --adjust /var/dgrp2/www/data/freeze2March11/adjustData.RData --maf 4 --cut 0.00001 --mpheno invHead.mpheno --fpheno invHead.fpheno --out invHead --final invHeadFinal > invHead.out 2>&1 &

# make raleigh plots
# ============================================================

Rscript manhattan.plot.R jingFaceFinal/output jingFaceWidth.png "CSB-jing Face Width" &
Rscript manhattan.plot.R invFaceFinal/output invFaceWidth.png "CSB-inv Face Width" &
Rscript manhattan.plot.R jingHeadFinal/output jingHeadWidth.png "CSB-jing Head Width" &
Rscript manhattan.plot.R invHeadFinal/output invHeadWidth.png "CSB-inv Head Width" &