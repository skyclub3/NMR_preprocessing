#!/bin/csh

#replaced by rawdata_parsing.pl for handling chain ; (not yet) 'or' string in file (ambiguity hydrogens ; old version treat them as seperate one) 
#082621 excerpt from /data/jinhyuk/package/STAProtamer_refine/data
#082621 consider ambigous restraints ; but the ambigous restraints has the same number in the first column ; already done ; so no change

set input    = $1 # rawdata (X-PLOR type) downloaded from BMRB
set libawk   = /data2/yycho/Protein/yycho/opt/src_packages/lib.awk

if ( $#argv == 1 ) then
else
    echo $0 "[BMRBnoerawdata]"
    exit
endif

#awk -i ~/bin/lib.awk 'BEGIN{RS=" assi";FS="\n"}{j=0;for(i=1;i<=NF;i++){if(length(trim($i))&&trim($i)!="or"){j++;print NR,j,remquot(trimparen(trim($i)))}}}' /home/josh617/NMR_Refinement/JSP_NMR_refinement_01/1a24/mknmrdata/rawdata

awk -i /data2/yycho/Protein/yycho/opt/src_packages/lib.awk 'BEGIN{RS=" assi";FS="\n"}\
{\
    j=0;\
    for(i=1;i<=NF;i++){\
	if(length(trim($i))&&trim($i)!="or"){\
	    j++;\
	    #print NR,j,remquot(trimparen(trim($i)));\
	    if(j==3){\
		dist=trim($i);\
		nc=split(dist,c," ");\
		if(c[3]<0){ \
		  dist=sprintf("%s %s %s",c[1]-c[2],0,c[2]+c[3]) \
		}else{ \
		} \
		#print dist;\
	    }else{\
		#print NR,j,remquot(trimparen(trim($i)));\
		str[j]=(NR remquot(trimparen(trim($i))));\
	    }\
	}\
    }\
    #tj=(j>3?j-1:j);\
    tj=j;\
    for(i=1;i<tj;i=i+2){\
	k=(i>=3?i+1:i);\
	na=split(str[k],a," ");nb=split(str[k+1],b," ");\
	print a[1],a[3],a[6],"null","null",a[9],b[3],b[6],"null","null",b[9],dist;\
    }\
}' $input
