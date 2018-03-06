#!/bin/bash
# Part of svashor which generate reports with Structural Variants of
# Alpha Satellite Higher-Order Repeats
# Make HOR map for all BED files in directory
# It is intended to use for long contigs after the HMMER HOR analysis

while [ "$#" -gt 0 ]; do
  case "$1" in
    -h) horabbr="$2"; shift 2;;

    --hor=*) horabbr="${1#*=}"; shift 1;;
    --hor) echo "$1 requires an argument" >&2; exit 1;;

    -*) echo "unknown option: $1" >&2; exit 1;;
    *) dir="$1"; shift 1;;
  esac
done

if [ ! -d "$dir" ]
then
    echo "$dir is not a directory!" 1>&2
    exit 1
fi

export horabbr="$horabbr"

for bedFile in "$dir"/*.bed; do
    # Select HORs containing more than 'montr' number of monomers
    tmpHN=()
    if [ -z "$horabbr" ]; then
        while IFS= read -r line; do
            tmpHN+=( "$line" )
        done < <(
            awk -v montr=1 '{
                n = index($4, ".");
                s = substr($4, 1, n-1);
                a[s]++
            }
            END{
                for(i in a) if(a[i]>montr) print i
            }' "$bedFile"
        );
    else
        tmpHN="$horabbr"
    fi
    # Print maps for each selected HOR
    for i in "${!tmpHN[@]}"; do
        escHN=$( echo "${tmpHN[$i]}" | sed -e "s%\/%\\\/%g" )
        sanHN=$( echo "${tmpHN[$i]}" | sed -e "s/[^A-Za-z0-9._-]/_/g" )
        # Remove text before dot in HOR name and
        # print names of the monomers starting from the first
        awk -v maxHORLenght=50 -v wos=0 '
            /^chr|^AC|^GJ/{
                if(!once++) printf("%s\t%s\t", $1, $2);
                sub(/^'"$escHN"'\./, "", $4);
                if(wos != 0){
                    sub(/\/.+?$/, "", $4);
                    sub(/-.+?$/, "", $4)
                }
                if(($4 ~ /^1$/)||(n==maxHORLenght)){
                    printf("\n");
                    printf("%s\t%s\t", $1, $2);
                    n=0
                }
                printf("%s ", $4);
                n++
            }
            END{
                printf("\n%s\t%s\t\n", $1, $3)
            }' "$bedFile" > "$dir"/$(basename "$bedFile" .bed)-$sanHN-map.txt
    done
done
