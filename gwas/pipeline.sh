#! /bin/bash

args=`getopt -o "f:,m:,s:,o:,b:,R:,p:,a:,m:,r:,n:,l:,z:,q:,j:,u:,c:" -l "fpheno:,mpheno:,spheno:,out:,bfile:,Rscript:,plink:,annot:,fastlmm:,rel:,nbfile:,ldplot:,final:,qqplot:,adjust:,maf:,cut:" -- "$@"`

echo "command arguments given: $args"

eval set -- "$args" 

plink="plink"
Rscript="Rscript"
fastlmm="fastlmm"
outDir="./plinkOut/output/"
resDir="./plinkOut/"
finalDir="./finalDir"
maf=0.05
cut=0.00001

# parse arguments
while true;
do
  case $1 in
    -f|--fpheno)
       femalePhenoFile=$2
       shift 2;;

    -u|--maf)
       maf=$2
       shift 2;;

    -c|--cut)
       cut=$2
       shift 2;;

    -m|--mpheno)
       malePhenoFile=$2
       shift 2;;

    -s|--spheno)
       singlePhenoFile=$2
       shift 2;;

    -o|--out)
       outDir="$2/output/"
       resDir="$2/"
       shift 2;;

    -j|--adjust)
       adjustData=$2
       shift 2;;

    -z|--final)
       finalDir="$2/"
       shift 2;;

    -b|--bfile)
       bFile=$2
       shift 2;;

    -R|--Rscript)
       Rscript=$2
       shift 2;;

    -p|--plink)
       plink=$2
       shift 2;;

    -m|--fastlmm)
       fastlmm=$2
       shift 2;;

    -a|--annot)
       annot=$2
       shift 2;;

    -l|--ldplot)
       ldplot=$2
       shift 2;;

    -r|--rel)
       rel=$2
       shift 2;;

    -b|--nbfile)
       nbFile=$2
       shift 2;;
    -q|--qqplot)
       qqplot=$2
       shift 2;;

    --)
       shift
       break;;
  esac
done

# merge female, male, phenotype and 
# assume that the files have been parsed (John Jack?) to become line mean csv
# and contain only complete phenotypes without header

# check if necessary programs exist
RscriptCheck=`command -v $Rscript`
if [ $RscriptCheck ]
then
  echo "using $RscriptCheck"
else
  echo "cannot find Rscript"
  exit 1 # 1 = command not found
fi

plinkCheck=`command -v $plink`
if [ $plinkCheck ]
then
  echo "using $plinkCheck"
else
  echo "cannot find plink"
  exit 1
fi

fastlmmCheck=`command -v $fastlmm`
if [ $fastlmmCheck ]
then
  echo "using $fastlmmCheck"
else
  echo "cannot find fastlmm"
  exit 1
fi

# check if all files exist

if [ $malePhenoFile ]
then
  if [[ ! -e $malePhenoFile ]]
  then
    echo "$malePhenoFile cannot be found"
    exit 2 # 2 = file not found
  fi
fi

if [[ ! -e $ldplot ]]
  then
    echo "$ldplot cannot be found"
    exit 2
fi

if [[ ! -e $qqplot ]]
  then
    echo "$qqplot cannot be found"
    exit 2
fi

if [ $femalePhenoFile ]
then
  if [[ ! -e $femalePhenoFile ]]
  then
    echo "$femalePhenoFile cannot be found"
    exit 2
  fi
fi

if [ $singlePhenoFile ]
then
  if [[ ! -e $singlePhenoFile ]]
  then
    echo "$singlePhenoFile cannot be found"
    exit 2
  fi
fi

if [ $annot ]
then
  if [[ ! -e $annot ]]
  then
    echo "$annot cannot be found"
    exit 2
  fi
else
  echo "annotation file not given"
  exit 2
fi

if [ $adjustData ]
then
  if [[ ! -e $adjustData ]]
  then
    echo "$adjustData not found"
    exit 2
  fi
else
  echo "adjustment data not given"
  exit 2
fi

if [[ ! -e $bFile.bim || ! -e $bFile.bed || ! -e $bFile.fam ]]
then
  echo "plink formatted genotypes $bFile.bed/bim/fam not found"
  exit 2
fi

if [[ ! -e $nbFile.bim || ! -e $nbFile.bed || ! -e $nbFile.fam ]]
then
  echo "plink formatted genotypes $nbFile.bed/bim/fam for the indel recoded not found"
  exit 2
fi

if [ $rel ]
then
  if [[ ! -e $rel ]]
  then
    echo "$rel cannot be found"
    exit 2
  fi
else
  echo "relationship file not given"
  exit 2
fi

# make sure that only a single file or both sexes are given, but not both

if [ $singlePhenoFile ]
then
  if [[ $femalePhenoFile || $malePhenoFile ]]
  then
    echo "provide (both female and male) or (a single) phenotype files"
    exit 3 # 3 non-conforming format
  fi
else
  if [[ ! $femalePhenoFile || ! $malePhenoFile ]]
  then
    echo "provide (both female and male) or (a single) phenotype files"
    exit 3 # 3 non-conforming format
  fi
fi

# check if directory exist
# This program assumes that the resDir exists
# this is part of the preprocessing process by JJ
# when running the pipeline offline, need to make the directory if it does not exist already
if [[ ! -e $resDir ]]
then
  echo "$resDir does not exist, making one"
  mkdir $resDir
fi

if [ -e $outDir ]
then
  echo "$outDir exists"
  exit 3
else
  mkdir $outDir
  if [ $? == 1 ]
  then
    echo "cannot make $outDir"
    exit 3
  fi
fi

# if there are two phenotypes, make the average and difference
# and perform association

if [[ $femalePhenoFile && $malePhenoFile ]]
then
  date
  echo "both female and male phenotypes are given"
  echo "preparing phenotypes and genotypes"
  sed 's/,/\t/' $femalePhenoFile | awk '{print "line_"$0}' | sort -k1,1 | perl -pe 's/\r\n$/\n/g' > $outDir"tmp.female.pheno"
  sed 's/,/\t/' $malePhenoFile | awk '{print "line_"$0}' | sort -k1,1 | perl -pe 's/\r\n$/\n/g' > $outDir"tmp.male.pheno"
  echo $outDir | Rscript -e 'out.dir <- scan(pipe("cat /dev/stdin"), what = "", quiet = T); closeAllConnections(); female.pheno <- read.table(paste(out.dir, "/tmp.female.pheno", sep = ""), header = F, as.is = T); female.pheno[, 2] <- suppressWarnings(as.numeric(female.pheno[, 2])); male.pheno <- read.table(paste(out.dir, "/tmp.male.pheno", sep = ""), header = F, as.is = T); male.pheno[, 2] <- suppressWarnings(as.numeric(male.pheno[, 2])); merge.data <- merge(female.pheno, male.pheno, by.x = 1, by.y = 1, all = T); merge.data <- cbind(merge.data[, 1], merge.data[, 1], (merge.data[, 2] + merge.data[, 3])/2, merge.data[, 2] - merge.data[, 3], merge.data[, 2:3]); write.table(merge.data, file = paste(out.dir, "/tmp.pheno", sep = ""), sep = " ", row.names = F, col.names = c("FID", "IID", "avg", "diff", "female", "male"), quote = F);'
  # extract IDs
  tail -n+2 $outDir"/tmp.pheno" | grep -v NA | cut -d " " -f 1,2 > $outDir"tmp.nonmiss.id"
  # correct phenotype
  # phenotype before correction
  cp $outDir"/tmp.pheno" $outDir"/tmp.before.pheno"
  echo "line_id avg_raw diff_raw female_raw male_raw avg_adjusted diff_adjusted female_adjusted male_adjusted" > $outDir"/raw.adjusted.pheno.txt"
  echo $outDir"/tmp.pheno" $adjustData | $Rscript -e 'pheno.file.name = scan(pipe("cat /dev/stdin"), what = "", quiet = T); closeAllConnections(); load(pheno.file.name[2]); pheno <- read.table(pheno.file.name[1], header = T, as.is = T); pheno[, 5] <- adjustPheno(pheno[, c(1, 5)], "Female"); pheno[, 6] <- adjustPheno(pheno[, c(1, 6)], "Male"); pheno[, 3] <- adjustPheno(pheno[, c(1, 3)], "Avg"); pheno[, 4] <- adjustPheno(pheno[, c(1, 4)], "Diff"); write.table(pheno, file = pheno.file.name[1], col.names = T, row.names = F, sep = " ", quote = F);' > $outDir"/pheno.adjust.txt"
  paste -d " " $outDir"/tmp.before.pheno" $outDir"/tmp.pheno" | tail -n+2 | cut -d " " -f 1,3-6,9-12 >> $outDir"/raw.adjusted.pheno.txt"
  
  # extract relevant columns from the relationship matrix, unnecessary
  # echo $rel $outDir"tmp.nonmiss.id" $outDir | $Rscript -e ' files <- scan(con <- pipe("cat /dev/stdin"), what = "", quiet = TRUE); close(con); rel <- read.table(files[1], header = FALSE, as.is = TRUE, sep = "\t"); rel.ids <- unlist(strsplit(unlist(rel[1, -1]), split = " ")); rel.ids <- rel.ids[seq(1, length(rel.ids), 2)]; ids <- read.table(files[2], header = FALSE, as.is = TRUE)[, 1]; rel.select <- which(rel.ids %in% ids); rel <- rel[c(1, 1+rel.select), c(1, 1+rel.select)]; write.table(rel, file = paste(files[3], "tmp.rel.mat", sep = ""), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE);'

  # filter on genotypes
  # if $maf < 1, use $maf as a cutoff
  # if $maf > 1, use $maf as a cutoff of minimum number of lines
  # always enforce call_rate > 0.8

  date
  echo "calculating allele frequency"
  $plink --bfile $bFile --silent --noweb --keep $outDir"tmp.nonmiss.id" --freq --out $outDir"tmp.freq"

  # check if $maf is less than 1
  nind=`wc -l $outDir/tmp.nonmiss.id | cut -d " " -f 1 | awk '{print int($1*0.8)}'`
  uu=`echo $maf | awk '$1 < 1' | wc -l`
  # if $maf < 1
  
  if [ $uu -eq 1 ]
  then
    tail -n+2 $outDir"tmp.freq.frq" | awk '$5 >= '"$maf"' && ($6 + 1)/2 >= '"$nind"' && ($5*$6 + 1)/2 >= 4 && ($6 + 1)/2 >= 30 {print $2}' | sort > $outDir"tmp.test.variant"
  else # when $maf > 1
    tail -n+2 $outDir"tmp.freq.frq" | awk '($5*$6 + 1)/2 >= '"$maf"' && ($6 + 1)/2 >= '"$nind"' && ($5*$6 + 1)/2 >= 4 && ($6 + 1)/2 >= 30 {print $2}' | sort > $outDir"tmp.test.variant"
  fi
  
  tail -n+2 $outDir/tmp.freq.frq | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"int($6/2*$5 + 0.5)"\t"int($6/2*(1-$5) + 0.5)}' | sort -k1,1 | join $outDir"tmp.test.variant" - -t $'\t' > $outDir/tmp.freq.sorted

  
  # association
  date
  echo "performing association analyses"
  $plink --bfile $bFile --silent --noweb --extract $outDir"tmp.test.variant" --assoc --pheno $outDir"tmp.pheno" --missing-phenotype NA --all-pheno --out $outDir"tmp.assoc"

  # fastlmm
  date
  echo "performing mixed model GWAS"
  $fastlmm -bfile $nbFile -extract $outDir"tmp.test.variant" -pheno $outDir"tmp.pheno" -sim $rel -missingPhenotype NA -mpheno 1 -out $outDir/tmp.lmm.avg.assoc > $outDir/tmp.lmm.avg.assoc.log 2>&1 
  tail -n+2 $outDir/tmp.lmm.avg.assoc | awk '{print $1"\t"$5}' | sort -k1,1 > $outDir/tmp.lmm.avg.assoc.sorted
  rm $outDir/tmp.lmm.avg.assoc

  $fastlmm -bfile $nbFile -extract $outDir"tmp.test.variant" -pheno $outDir"tmp.pheno" -sim $rel -missingPhenotype NA -mpheno 2 -out $outDir/tmp.lmm.diff.assoc > $outDir/tmp.lmm.diff.assoc.log 2>&1 
  tail -n+2 $outDir/tmp.lmm.diff.assoc | awk '{print $1"\t"$5}' | sort -k1,1 > $outDir/tmp.lmm.diff.assoc.sorted
  rm $outDir/tmp.lmm.diff.assoc

  $fastlmm -bfile $nbFile -extract $outDir"tmp.test.variant" -pheno $outDir"tmp.pheno" -sim $rel -missingPhenotype NA -mpheno 3 -out $outDir/tmp.lmm.female.assoc > $outDir/tmp.lmm.female.assoc.log 2>&1 
  tail -n+2 $outDir/tmp.lmm.female.assoc | awk '{print $1"\t"$5}' | sort -k1,1 > $outDir/tmp.lmm.female.assoc.sorted
  rm $outDir/tmp.lmm.female.assoc

  $fastlmm -bfile $nbFile -extract $outDir"tmp.test.variant" -pheno $outDir"tmp.pheno" -sim $rel -missingPhenotype NA -mpheno 4 -out $outDir/tmp.lmm.male.assoc > $outDir/tmp.lmm.male.assoc.log 2>&1 
  tail -n+2 $outDir/tmp.lmm.male.assoc | awk '{print $1"\t"$5}' | sort -k1,1 > $outDir/tmp.lmm.male.assoc.sorted
  rm $outDir/tmp.lmm.male.assoc

  # assemble results
  date
  echo "assembling results"
  paste $outDir/tmp.lmm.female.assoc.sorted $outDir/tmp.lmm.male.assoc.sorted $outDir/tmp.lmm.avg.assoc.sorted $outDir/tmp.lmm.diff.assoc.sorted | cut -f 1,2,4,6,8 | sort -k1,1 > $outDir"tmp.lmm.all.assoc"
  paste $outDir"tmp.assoc.female.qassoc" $outDir"tmp.assoc.male.qassoc" $outDir"tmp.assoc.avg.qassoc" $outDir"tmp.assoc.diff.qassoc" | tail -n+2 | awk '{print $2"\t"$1"\t"$3"\t"$5"\t"$6"\t"$9"\t"$14"\t"$15"\t"$18"\t"$23"\t"$24"\t"$27"\t"$32"\t"$33"\t"$36}' | sort -k1,1 > $outDir"tmp.gwas.all.assoc"

  echo "ID MinorAllele MajorAllele MAF MinorAlleleCount MajorAlleleCount FemalePval FemaleMixedPval MalePval MaleMixedPval AvgPval AvgMixedPval DiffPval DiffMixedPval" > $outDir"gwas.all.assoc"
  join $outDir/tmp.freq.sorted $outDir"tmp.gwas.all.assoc" -t $'\t' | join - $outDir"tmp.lmm.all.assoc" -t $'\t' | sort -k7,7 -k8,8n | awk '{print $1" "$2" "$3" "$4" "$5" "$6" "$11" "$21" "$14" "$22" "$17" "$23" "$20" "$24}' >> $outDir"gwas.all.assoc"
  
  # annotate top snps
  date
  echo "annotating top snps"
  echo "ID MinorAllele MajorAllele RefAllele MAF MinorAlleleCount MajorAlleleCount FemaleEff FemalePval FemaleMixedPval MaleEff MalePval MaleMixedPval AvgEff AvgPval AvgMixedPval DiffEff DiffPval DiffMixedPval GeneAnnotation RegulationAnnotation" > $outDir/gwas.top.annot
  cat $outDir"tmp.gwas.all.assoc" | join $outDir/tmp.freq.sorted - -t $'\t' | join - $outDir"tmp.lmm.all.assoc" -t $'\t' | join -e "-" -a 1 -o '0,1.2,1.3,2.2,1.4,1.5,1.6,1.9,1.11,1.21,1.12,1.14,1.22,1.15,1.17,1.23,1.18,1.20,1.24,2.3,2.4' - $annot | awk '$9 < '$cut' || $10 < '$cut' || $12 < '$cut' || $13 < '$cut' || $15 < '$cut' || $16 < '$cut' || $18 < '$cut' || $19 < '$cut'' | perl -wne 'use List::Util qw(min); chomp $_; @line = split / /, $_; $minp = min(@line[(9,12,15,18)]); print $minp, " ", $_, "\n";' | sort -t " " -k1,1g | head -n 1000 | cut -d " " -f 2- | awk -F " " '{print $1" "$2" "$3" "$4" "$5" "$6" "$7" "0-$8" "$9" "$10" "0-$11" "$12" "$13" "0-$14" "$15" "$16" "0-$17" "$18" "$19" "$20" "$21}' >> $outDir/gwas.top.annot
  # make LD heatmap
  # need to be careful when there is no snp among the top
  date
  echo "making LD heatmap"
  tail -n+2 $outDir/gwas.top.annot | cut -d " " -f 1 > $outDir/tmp.ld.top.snp.id
  if [ `wc -l $outDir/tmp.ld.top.snp.id | cut -d " " -f 1` -ge 2 ]
  then
    # extract snp genotypes
    $plink --bfile $bFile --silent --noweb --extract $outDir/tmp.ld.top.snp.id --recode12 --out $outDir/tmp.top.snp.geno
    $Rscript $ldplot $outDir/tmp.top.snp.geno.ped $outDir/tmp.top.snp.geno.map $outDir/LDheat.eps $outDir/tmp.pheno $outDir
  fi
  # clean up
  date
  echo "analyses finished, cleaning up"
  rm $outDir/tmp.*

fi

if [ $singlePhenoFile ]
then
  date
  echo "only a single phenotype is given"
  echo "preparing phenotypes and genotypes"
  
  sed 's/,/\t/' $singlePhenoFile | awk '{print "line_"$0}' | sort -k1,1 | perl -pe 's/\r\n$/\n/g' > $outDir"tmp.pheno"
  echo $outDir | Rscript -e 'out.dir <- scan(pipe("cat /dev/stdin"), what = "", quiet = T); closeAllConnections(); single.pheno <- read.table(paste(out.dir, "/tmp.pheno", sep = ""), header = F, as.is = T); single.pheno[, 2] <- suppressWarnings(as.numeric(single.pheno[, 2])); write.table(cbind(single.pheno[, 1], single.pheno), file = paste(out.dir, "/tmp.pheno", sep = ""), sep = " ", row.names = F, col.names = c("FID", "IID", "single"), quote = F);'
  # extract IDs
  tail -n+2 $outDir"/tmp.pheno" | grep -v NA | cut -d " " -f 1,2 > $outDir"tmp.nonmiss.id"
  
  # extract relevant columns from the relationship matrix, unnecessary
  # echo $rel $outDir"tmp.nonmiss.id" $outDir | $Rscript -e ' files <- scan(con <- pipe("cat /dev/stdin"), what = "", quiet = TRUE); close(con); rel <- read.table(files[1], header = FALSE, as.is = TRUE, sep = "\t"); rel.ids <- unlist(strsplit(unlist(rel[1, -1]), split = " ")); rel.ids <- rel.ids[seq(1, length(rel.ids), 2)]; ids <- read.table(files[2], header = FALSE, as.is = TRUE)[, 1]; rel.select <- which(rel.ids %in% ids); rel <- rel[c(1, 1+rel.select), c(1, 1+rel.select)]; write.table(rel, file = paste(files[3], "tmp.rel.mat", sep = ""), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE);'
  
  cp $outDir"/tmp.pheno" $outDir"/tmp.before.pheno"
  echo $outDir"/tmp.pheno" $adjustData | $Rscript -e 'pheno.file.name = scan(pipe("cat /dev/stdin"), what = "", quiet = T); closeAllConnections(); load(pheno.file.name[2]); pheno <- read.table(pheno.file.name[1], header = T, as.is = T); pheno[, 3] <- adjustPheno(pheno[, c(1, 3)], "Single"); write.table(pheno, file = pheno.file.name[1], col.names = T, row.names = F, sep = " ", quote = F);' > $outDir"/pheno.adjust.txt"
  echo "line_id single_raw single_adjusted" > $outDir"/raw.adjusted.pheno.txt"
  paste -d " " $outDir"/tmp.before.pheno" $outDir"/tmp.pheno" | tail -n+2 | cut -d " " -f 1,3,6 >> $outDir"/raw.adjusted.pheno.txt"

  # filter on genotypes
  date
  echo "calculating allele frequency"
  $plink --bfile $bFile --silent --noweb --keep $outDir"tmp.nonmiss.id" --freq --out $outDir"tmp.freq"


  nind=`wc -l $outDir/tmp.nonmiss.id | cut -d " " -f 1 | awk '{print int($1*0.8)}'`
  uu=`echo $maf | awk '$1 < 1' | wc -l`
  # if $maf < 1
  
  if [ $uu -eq 1 ]
  then
    tail -n+2 $outDir"tmp.freq.frq" | awk '$5 >= '"$maf"' && ($6 + 1)/2 >= '"$nind"' && ($5*$6 + 1)/2 >= 4 && ($6 + 1)/2 >= 30 {print $2}' | sort > $outDir"tmp.test.variant"
  else # when $maf > 1
    tail -n+2 $outDir"tmp.freq.frq" | awk '($5*$6 + 1)/2 >= '"$maf"' && ($6 + 1)/2 >= '"$nind"' && ($5*$6 + 1)/2 >= 4 && ($6 + 1)/2 >= 30 {print $2}' | sort > $outDir"tmp.test.variant"
  fi
  
  tail -n+2 $outDir/tmp.freq.frq | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"int($6/2*$5 + 0.5)"\t"int($6/2*(1-$5) + 0.5)}' | sort -k1,1 | join $outDir"tmp.test.variant" - -t $'\t' > $outDir/tmp.freq.sorted

  # association
  date
  echo "performing association analyses"
  $plink --bfile $bFile --silent --noweb --extract $outDir"tmp.test.variant" --assoc --pheno $outDir"tmp.pheno" --missing-phenotype NA --all-pheno --out $outDir"tmp.assoc"

  # fastlmm
  date
  echo "performing mixed model GWAS"
  $fastlmm -bfile $nbFile -extract $outDir"tmp.test.variant" -pheno $outDir"tmp.pheno" -sim $rel -missingPhenotype NA -mpheno 1 -out $outDir/tmp.lmm.single.assoc > $outDir/tmp.lmm.single.assoc.log 2>&1 
  tail -n+2 $outDir/tmp.lmm.single.assoc | awk '{print $1"\t"$5}' | sort -k1,1 > $outDir/tmp.lmm.single.assoc.sorted
  rm $outDir/tmp.lmm.single.assoc

  # assemble results
  date
  echo "assembling results"
  cat $outDir"tmp.assoc.single.qassoc" | tail -n+2 | awk '{print $2"\t"$1"\t"$3"\t"$5"\t"$6"\t"$9}' | sort -k1,1 > $outDir"tmp.gwas.all.assoc"

  echo "ID MinorAllele MajorAllele MAF MinorAlleleCount MajorAlleleCount SingleEff SinglePval SingleMixedPval" > $outDir"gwas.all.assoc"
  join $outDir/tmp.freq.sorted $outDir"tmp.gwas.all.assoc" -t $'\t' | join - $outDir"tmp.lmm.single.assoc.sorted" -t $'\t' | sort -k7,7 -k8,8n | awk '{print $1" "$2" "$3" "$4" "$5" "$6" "$9" "$11" "$12}' >> $outDir"gwas.all.assoc"
  
  # annotate top snps
  date
  echo "annotating top snps"
  echo "ID MinorAllele MajorAllele RefAllele MAF MinorAlleleCount MajorAlleleCount SingleEff SinglePval SingleMixedPval GeneAnnotation RegulationAnnotation" > $outDir/gwas.top.annot
  cat $outDir"tmp.gwas.all.assoc" | join $outDir/tmp.freq.sorted - -t $'\t' | join - $outDir"tmp.lmm.single.assoc.sorted" -t $'\t' | join -e "-" -a 1 -o '0,1.2,1.3,2.2,1.4,1.5,1.6,1.9,1.11,1.12,2.3,2.4' - $annot | awk '$9 < '$cut' || $10 < '$cut'' | sort -t " " -k10,10g | head -n 1000 | awk -F " " '{print $1" "$2" "$3" "$4" "$5" "$6" "$7" "0-$8" "$9" "$10" "$11" "$12}' >> $outDir/gwas.top.annot

  # make LD heatmap
  # need to be careful when there is no snp among the top
  date
  echo "making LD heatmap"
  tail -n+2 $outDir/gwas.top.annot | cut -d " " -f 1 > $outDir/tmp.ld.top.snp.id
  if [ `wc -l $outDir/tmp.ld.top.snp.id | cut -d " " -f 1` -ge 2 ]
  then
    # extract snp genotypes
    $plink --bfile $bFile --silent --noweb --extract $outDir/tmp.ld.top.snp.id --recode12 --out $outDir/tmp.top.snp.geno
    $Rscript $ldplot $outDir/tmp.top.snp.geno.ped $outDir/tmp.top.snp.geno.map $outDir/LDheat.eps $outDir/tmp.pheno $outDir
  fi

  # clean up
  date
  echo "analyses finished, cleaning up"
  rm $outDir/tmp.*
  
fi

# make qqplots
echo "making qqplots"
$Rscript $qqplot $outDir/gwas.all.assoc $outDir/
zip -j $outDir/qqplots.zip $outDir/*.png
rm $outDir/*.png

# move to final_dir
echo "moving files to final_dir"
mkdir $finalDir
mv $resDir/* $finalDir
rmdir $resDir

