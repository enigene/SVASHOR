NR==1{ topCn=$2 }
{
    hvar=$1;
    propHOR=0;
    propPass=0;
    n=split(hvar,a," ");
    for(i=1;i<=n;i++){
        if(a[i]!~/[A-Za-z]/){
            propHOR++
        } else {
            propHOR=0
        }
        if(propHOR>=propTH){
            propPass=1;
            break
        }
    }
    if((propPass)&&(topCn<cn)&&($2<cn)){
        if(split(hvar,a," ")>1){
            gsub(/\//,"\\/",hvar);
            print hvar
        }
    }
    if((propPass)&&($2>=cn)){
        gsub(/\//,"\\/",hvar);
        print hvar
    }
}
