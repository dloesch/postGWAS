#set up directories
mkdir -p inputs
mkdir -p BASE
mkdir -p FINAL

#summary stats are grouped by chromosome
#script it set up in case there are multiple loci per chromosome
files=$(ls data/ | sort -V) 

#specify type of annotations
#S. Nigra
type=Nigra
prefix=SN
#Brain 
#type=Brain
#prefix=BRAIN

for f in $files; do

	LOCUS=./data/${f} #file is in format chr4.txt

	chr=$(basename $f ".txt" | sed 's/chr//g') #basename will extract chr4

	#subset ref data, speeds up utility
	REF=ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	newREF=temp.chr$chr.vcf.gz

	#position column is 4 if not lifted, 2 if lifted
	start=$(awk 'NR>1 {print $2}' $LOCUS | sort -n | head -1)
	end=$(awk 'NR>1 {print $2}' $LOCUS | sort -n | tail -1)
	bcftools view -r $chr:${start}-${end} -M2 -o $newREF -Oz $REF
	
	#run paintor utility to create LD matrix and processed sumstats file
	MAP=integrated_call_samples_v3.20130502.ALL.panel
	POP=AMR

	dir=./PAINTOR_V3.0/PAINTOR_Utilities
	e=alt
	a=ref
	pos=pos
	Z=Score.Stat #algebraically equivalent to the Zscore obtained from beta/betaSE
	c=chr
	out=$(basename $f ".txt")

	python $dir/CalcLD_1KG_VCF.py -l $LOCUS -r $newREF --m $MAP -i $pos -e $e -a $a -p $POP -z $Z -o $out

	#annotations
		
	filepath=/mnt/hdd/Annotation_Paths
	LOCUS=$out.processed

	#get only desired annotations
	cat $filepath | grep $type > $prefix.txt
	annot=$prefix.txt		
	echo $filepath
	python $dir/AnnotateLocus.py -i $annot -l $LOCUS -o $out.annotations -c $c -p $pos
	

	#organize files
	mkdir -p inputs
	mv $out.processed inputs/$out
	mv $out.ld inputs
	mv $out.annotations inputs
	

done

ls inputs | grep -v "ld" | grep -v annot | sort -V > locus.files

input=locus.files
PAINTOR -input $input -Zhead $Z -LDname ld -in ./inputs -out ./BASE -enumerate 2 -Gname Enrich.BASE -Lname BF.BASE

#set up outputs for annnotations
cat $filepath | grep $type | cut -f6 -d"/" > ${prefix}_data.txt
mkdir -p $prefix

#set up while loop for number of annotations
N=$(cat ${prefix}_data.txt | wc -l)
i=1

while [ $i -le $N ]; do 

	anno=$(head -$i ${prefix}_data.txt| tail -1)
	R=Results.$i
	G=Enrich.$i

	mkdir -p SN
	PAINTOR -input $input -Zhead $Z -LDname ld -in ./inputs  -out ./SN -enumerate 3 -Gname $G -RESname $R -annotations $anno -Lname BF.$i
	
	#if this annotation has significant enrichment, just stop loop and use it
	Rscript eval_annot_run.R $i 57
	
	status=$(cat runfile.txt)
	
	if [[ $status == "STOP" ]]; then break ; fi
	i=$((i+1))
done

#if no annotation reached the Bonferonni significance level, evaluate all and keep top 3
N=$((N+1))
if [[ $i == $N ]]; then
	
	Rscript eval_annot_all.R chr$chr.annotations

	#final model
	input=locus.files
	temp=$(cat ${prefix}_top_results.txt | grep -v '#' | tr '\n' ',')
	anno=$(echo "${temp%?}")

	PAINTOR -input $input -Zhead $Z -LDname ld -in ./inputs  -out ./FINAL -enumerate 3  -Gname Enrich.FINAL -RESname Results.FINAL -annotations $anno -Lname BF.FINAL

	#filter by 95% credible set
	Rscript filter_sumstats.R $chr
else
	cp ${prefix}/chr$chr.Results.$i ./FINAL
	Rscript filter_sumstats.R $chr

fi

