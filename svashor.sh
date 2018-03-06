#!/bin/bash
# Generate reports with Structural Variants of Alpha Satellite Higher-Order Repeats

set -o pipefail
set -o errexit

if [[ $# -ne 1 ]]
then
    echo "usage: $0 path_to_the_input_folder" 1>&2
    exit 1
fi

if [ ! -d $1 ]
then
    echo "$1 is not a directory!" 1>&2
    exit 1
fi

dir="$1"

scrdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

hsdir="$dir/AS_HOR_structure"

asmFrags="$scrdir"/hg38-frags-sorted.bed
prcnAS="$scrdir"/human-GRC-hg38-M1SFsv2.2-sorted.bed

mkdir "$hsdir"

# Count columns in file asmFrags
asmNF="$( awk '(!asmNF)&&($1!~/^#/){asmNF=NF;print asmNF}' "$asmFrags" )"

# http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
# bedtools intersect screens for overlaps between two sets:
# asmFrags and multiple files in directory
bedtools intersect -sortout -wb -filenames \
-a "$asmFrags" \
-b "$dir"/*.bed |
awk -v asmNF="$asmNF" '!(once++){
        n=split($0,a);
        for(i=1;i<=n;i++){
            if(a[i]~/\.bed/) fd=i+1;
        }
        if(!fd) fd=asmNF+1
    }
    {
    for(i=fd;i<=NF-1;i++)
        printf("%s\t", $i) >> "'$hsdir'/"$1"_"$4".bed";
        printf("%s\n", $NF) >> "'$hsdir'/"$1"_"$4".bed"
    }'

cd "$hsdir"

"$scrdir"/bedHORContigMap.sh "$hsdir"
"$scrdir"/bedHORContigStats.sh "$hsdir"

for horStatFile in *-stat.txt ; do
    if [[ -s $horStatFile ]]; then
        fnwosuff=$(basename "$horStatFile" -stat.txt)
        fragName=${fnwosuff%%-*}
        horAbbr=${fnwosuff#*-}
        unEscHorAbbr=$( echo "$horAbbr" | sed -e 's%_%\/%g' )
        horStatA=()
        while IFS= read -r line; do
            horStatA+=( "$line" )
        # Parameter cn define a copy number threshold for variants
        # Parameter propTH define number of сonsecutive monomers in proper HORs
        done < <(awk -F"\t" -v cn=10 -v propTH=2 -f "$scrdir"/horstat2list.awk "$horStatFile")

        varN=1
        for i in "${!horStatA[@]}"; do

            awk -F"\t" -v OFS="\t" '
            f{
                end=$2;
                print chr,start,end-1,pNR;
                f=0
            }
            {
                if($3 ~ /^'"${horStatA[$i]}"'$/){
                    f=1;
                    chr=$1;
                    start=$2;
                    pNR=NR
                }
            }' "$fnwosuff"-map.txt > tmp-svashor-"$fnwosuff"-var"${varN}".bed

            bedtools intersect -sortout -wb \
            -a tmp-svashor-"$fnwosuff"-var"${varN}".bed \
            -b "$prcnAS" > tmp-svashor-intersectBed-"$fnwosuff"-var"${varN}".bed

            # Test if fragName present in prcnAS
            if [[ -s tmp-svashor-intersectBed-"$fnwosuff"-var"${varN}".bed ]]; then

                cut -f5- tmp-svashor-intersectBed-"$fnwosuff"-var"${varN}".bed > tmp-svashor-intersectBed-"$fnwosuff"-var"${varN}"-cutf5.bed

                bedtools map -c 4 -o collapse -delim "_" \
                -a tmp-svashor-"$fnwosuff"-var"${varN}".bed \
                -b tmp-svashor-intersectBed-"$fnwosuff"-var"${varN}"-cutf5.bed > tmp-svashor-mapBed-"$fnwosuff"-var"${varN}".bed

                awk '{a[$5]++}END{for(i in a)printf("%s\t%d\n",i,a[i])}' tmp-svashor-mapBed-"$fnwosuff"-var"${varN}".bed > tmp-svashor-allvars-"$fnwosuff"-var"${varN}".txt

                sort -t"$(printf "\t")" -k2nr tmp-svashor-allvars-"$fnwosuff"-var"${varN}".txt > tmp-svashor-allvars-"$fnwosuff"-var"${varN}"-sorted.txt

                sed -e 's/_/\t/g' tmp-svashor-allvars-"$fnwosuff"-var"${varN}"-sorted.txt > tmp-svashor-allvars-"$fnwosuff"-var"${varN}"-sorted-fmt1.txt

                unEscHorStatAi=$( echo "${horStatA[$i]}" | sed -e 's%\\%%g' )

                # number after parameter -n defines how many PERCON variants
                # should be printed out, -n1 prints most high copy variant
                head -n1 tmp-svashor-allvars-"$fnwosuff"-var"${varN}"-sorted-fmt1.txt |
                awk -v s="$unEscHorStatAi" '
                    BEGIN{
                        split(s,sA)
                    }
                    {
                        for(i=1;i<NF;i++){
                            split($i,mA,",");
                            printf("%s\t%s\t%s%s\t%s%s\t%s\t%s\t%s\t%s\n",
                            "'$fragName'","'$unEscHorAbbr'","var","'${varN}'","var",NR,i,mA[2],mA[1],sA[i])
                        }
                    }' - >> "$fnwosuff".tsv

            else echo "$fragName" was not found in "$prcnAS"
            fi

            varN=$(expr $varN + 1)

        done

        if [[ -s "$hsdir/$fnwosuff".tsv ]]; then
            # Add column titles to tsv file
            sed -i '1s/^/Frag\tHORabbr\tHVar\tPVar\tMon\tType\tClass\tNum\n/' "$hsdir/$fnwosuff".tsv

            "$scrdir"/svashor.R "$hsdir/$fnwosuff".tsv "$hsdir/$horStatFile" "$hsdir/$fnwosuff".html
        fi

    fi
done

# Concatenate HTML outputs by HOR name
find . -name "*.html" -print | cut -d- -f2- | sort | uniq | while read -r line; do find . -name "*$line*" -exec cat {} > "../$line" \; ; done

cp "$scrdir"/style.css .
cp "$scrdir"/style.css ./..
find "$hsdir" -name "tmp-svashor-*" -delete
