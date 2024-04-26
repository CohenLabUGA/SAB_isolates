#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --mem=10G
#SBATCH --job-name=koNames
#SBATCH --output=%x.out
#SBATCH --error=%x.err

#while read line; do hit=$(curl -# http://rest.kegg.jp/list/pathway | grep -w "${line/ko/map}"); echo -e "${line}\t${hit}" |cut -f1,3 >> all.ko.txt; done < all.ko.txt

source all.ko.sh
#source all.brite.txt
#source all.path.txt
#source all.modules.txt

for s in `echo $ko`; do
	echo ${s}
	hit=$(curl https://rest.kegg.jp/find/ko/${s})
	echo ${hit} >> ../kegg_names/ko_def.csv
	echo ${hit}
done


#for m in `echo $module`; do
# hit=$(curl https://rest.kegg.jp/find/module/${m})
#    echo ${hit} >> ../kegg_names/module_def.csv
#	echo ${hit}

#done
#for p in `echo $path`; do
# hit=$(curl https://rest.kegg.jp/find/pathway/${p})
#    echo ${hit} >> ../kegg_names/pathway_def.csv
#        echo ${hit}

#done

#for b in `echo $brite`; do
# hit=$(curl https://rest.kegg.jp/find/brite/${b})
  #  echo ${hit} >> ../kegg_names/brite_def.csv
 #       echo ${hit}

#done
