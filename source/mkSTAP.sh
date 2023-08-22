#!/bin/csh

#072420 two versions of STAP, STAP1 (old one, 201? ; /home/jinhyuk/package/cspline_interpolation/redund_rm_set/), STAP2 (new one, 2020 ; /data/jinhyuk/package/STAPpdb ; only 'no' SI version is available) 

set tarpdb   = $1
set dbtype   = $2 # 30 40 50 70 90 95 100 no no2
set forc     = $3
set dunbrack = 0  # 1: do not consider PRE residue because dunbrack set has no PRE or if u don't use PRE ; 0: default
set dele     = 1

#122320 include an option for analysis calculation energy by a residue
set analyres = 0  # 1: activate ; 0: no (hidden option ; activated manually)

#temporal files
set tmppdb   = tmp.pdb.$$
set tmpstr   = tmp.str.$$
set resresn  = resresn.tmp.$$
set ttmp     = t.tmp.$$

if ( $#argv == 3 ) then
else
    echo $0 "[PDB] [DBtype] [Force_constant] ,where [PDB]: PDBfile ; [DBtype]: 30 40 50 70 90 95 100 no no2, no: old STAP1, no2: new STAP2 ('no' SI ; all PDBs are included) ; [Force_constant]: STAP force constant ; if u want to control this from charmm input, use @IN1,2 (analysis residue number or STAP force constant)"
    exit
endif

awk 'function trim(str){sub(/^[ \t]+/,"",str);sub(/[ \t]+$/,"",str);return str}BEGIN{FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4"}{if(length($10)){printf "%s%4s%s\n",substr($0,1,22),trim($9 $10),(" " substr($0,28))}else{print $0}}' $tarpdb > $tmppdb

#set fres     = `grep -i " CA " $tarpdb | head -n 1 | awk '{print $6}'`
#set lres     = `grep -i " CA " $tarpdb | tail -n 1 | awk '{print $6}'`
set fres     = `awk 'BEGIN {FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4"}{if($4==" CA "){print $4,$9}}' $tmppdb | head -n 1 | awk '{print $2}'`
set lres     = `awk 'BEGIN {FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4"}{if($4==" CA "){print $4,$9}}' $tmppdb | tail -n 1 | awk '{print $2}'`

#set dbp      = /home/jinhyuk/package/cspline_interpolation/redund_rm_set/${dbtype}   # phipsi phichi1 psichi1 chi1chi2; extension version
#072420 add no2 for STAP2
if ( $dbtype == "no2" ) then

   #set dbp      = /home/jinhyuk/package/STAPpdb/no                                      # phipsi phichi1 psichi1 chi1chi2; extension version
   #122221 debug for PRE (dunbrack)
   if ($dunbrack) then
       set dbp      = `pwd`/source/no                                   # phipsi phichi1 psichi1 chi1chi2; extension version
   else
       set dbp      = `pwd`/source/no.pre                               # phipsi phichi1 psichi1 chi1chi2; extension version
   endif

endif

set potlist  = ( phipsi phichi1 psichi1 chi1chi2 ) # original extension version

######################################################################
#phipsi
set poten    = $potlist[1]
set npos     = `awk '{if(NR==1){print $1;exit}}' ${dbp}/${poten}/ARG.inp`
set xmin     = `awk '{if(NR==1){print $2;exit}}' ${dbp}/${poten}/ARG.inp`
set xmax     = `awk '{if(NR==1){print $3;exit}}' ${dbp}/${poten}/ARG.inp`
set xint     = `awk '{if(NR==1){print $4;exit}}' ${dbp}/${poten}/ARG.inp`
set ymin     = `awk '{if(NR==1){print $5;exit}}' ${dbp}/${poten}/ARG.inp`
set ymax     = `awk '{if(NR==1){print $6;exit}}' ${dbp}/${poten}/ARG.inp`
set yint     = `awk '{if(NR==1){print $7;exit}}' ${dbp}/${poten}/ARG.inp`

if ($dunbrack) then
    awk -v fr=$fres -v lr=$lres 'BEGIN {FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4"}{if($4==" CA "&&$9!=fr&&$9!=lr){print $6,$9}else{}}' $tmppdb | sed 's/HSD/HIS/g' | sed 's/HSP/HIS/g' | sed 's/HSE/HIS/g' | awk '{if(NF==0){}else{print $0}}' > $resresn
else
    awk -v fr=$fres -v lr=$lres 'BEGIN {FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4"}{if($4==" CA "&&$9!=fr&&$9!=lr){print $6,$9}else{}}' $tmppdb | sed 's/HSD/HIS/g' | sed 's/HSP/HIS/g' | sed 's/HSE/HIS/g' | awk '{if(length(pl)&&$1=="PRO"){printf "PRE %s\n",substr(pl,4)}else{printf "%s\n",pl};pl=$0}END{printf "%s\n",pl}' | awk '{if(NF==0){}else{print $0}}' > $resresn
endif

#122320 include an option for analysis calculation energy by a residue
if ( $analyres ) then
   awk -v forc=$forc -v ppp=${dbp}/${poten} -v npos=$npos -v xmin=$xmin -v xmax=$xmax -v xint=$xint -v ymin=$ymin -v ymax=$ymax -v yint=$yint 'BEGIN{print ""}{print "if @IN1 .eq. ",$2," then\nmod\nrots\nassign sele resid ",$2-1,".and. type C end -\n       sele resid ",$2,".and. type N end -\n       sele resid ",$2,".and. type CA end -\n       sele resid ",$2,".and. type C end -\n       sele resid ",$2,".and. type N end -\n       sele resid ",$2,".and. type CA end -\n       sele resid ",$2,".and. type C end -\n       sele resid ",$2+1,".and. type N end -\nforc",forc,"npos",npos,"-\nxmin",xmin,"xmax",xmax,"xint",xint,"-\nymin",ymin,"ymax",ymax,"yint",yint,"-\npos -"; pot=sprintf("%s/%s.inp",ppp,$1); sss=""; while((getline < pot)>0){sss=sprintf("%s %s",sss,$0)}; na=split(sss,a); for(i=8;i<=na;i++){printf "%s%s",a[i],(((i-7)%6==0)?" -\n":"  ");close(pot)};printf "\n\nend\nendif\n\n"}END{print "\nreturn"}' $resresn > ${tarpdb:r}_${poten}.str
else

   awk -v forc=$forc -v ppp=${dbp}/${poten} -v npos=$npos -v xmin=$xmin -v xmax=$xmax -v xint=$xint -v ymin=$ymin -v ymax=$ymax -v yint=$yint 'BEGIN{print "mod\nrots\n"}{print "assign sele resid ",$2-1,".and. type C end -\n       sele resid ",$2,".and. type N end -\n       sele resid ",$2,".and. type CA end -\n       sele resid ",$2,".and. type C end -\n       sele resid ",$2,".and. type N end -\n       sele resid ",$2,".and. type CA end -\n       sele resid ",$2,".and. type C end -\n       sele resid ",$2+1,".and. type N end -\nforc",forc,"npos",npos,"-\nxmin",xmin,"xmax",xmax,"xint",xint,"-\nymin",ymin,"ymax",ymax,"yint",yint,"-\npos -"; pot=sprintf("%s/%s.inp",ppp,$1); sss=""; while((getline < pot)>0){sss=sprintf("%s %s",sss,$0)}; na=split(sss,a); for(i=8;i<=na;i++){printf "%s%s",a[i],(((i-7)%6==0)?" -\n":"  ");close(pot)};printf "\n\n"}END{print "\nend\nreturn"}' $resresn > ${tarpdb:r}_${poten}.str

endif

######################################################################
#phixi1
set poten     = $potlist[2]

awk -v fr=$fres 'BEGIN {FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4"}{if($4==" CA "&&$9!=fr){print $6,$9}else{}}' $tmppdb | sed 's/HSD/HIS/g' | sed 's/HSP/HIS/g' | sed 's/HSE/HIS/g' | awk '{if($1=="ALA"||$1=="GLY"){}else{print $0}}' > $resresn

#122320 include an option for analysis calculation energy by a residue
if ( $analyres ) then
   awk -v forc=$forc -v ppp=${dbp}/${poten} -v npos=$npos -v xmin=$xmin -v xmax=$xmax -v xint=$xint -v ymin=$ymin -v ymax=$ymax -v yint=$yint 'BEGIN{print ""}{print "if @IN1 .eq. ",$2," then\nmod\nrots\nassign sele resid ",$2-1,".and. type C end -\n       sele resid ",$2,".and. type N end -\n       sele resid ",$2,".and. type CA end -\n       sele resid ",$2,".and. type C end -\n       sele resid ",$2,".and. type N end -\n       sele resid ",$2,".and. type CA end -\n       sele resid ",$2,".and. type CB end -\n       sele resid ",$2,".and. type",(($1=="CYS")?"SG":(($1=="ILE")?"CG1":(($1=="SER")?"OG":(($1=="THR")?"OG1":(($1=="VAL")?"CG1":"CG"))))),"end -\nforc",forc,"npos",npos,"-\nxmin",xmin,"xmax",xmax,"xint",xint,"-\nymin",ymin,"ymax",ymax,"yint",yint,"-\npos -"; pot=sprintf("%s/%s.inp",ppp,$1); sss=""; while((getline < pot)>0){sss=sprintf("%s %s",sss,$0)}; na=split(sss,a); for(i=8;i<=na;i++){printf "%s%s",a[i],(((i-7)%6==0)?" -\n":"  ");close(pot)};printf "\n\nend\nendif\n\n"}END{print "\nreturn"}' $resresn > ${tarpdb:r}_${poten}.str
else

   awk -v forc=$forc -v ppp=${dbp}/${poten} -v npos=$npos -v xmin=$xmin -v xmax=$xmax -v xint=$xint -v ymin=$ymin -v ymax=$ymax -v yint=$yint 'BEGIN{print "mod\nrots\n\n"}{print "assign sele resid ",$2-1,".and. type C end -\n       sele resid ",$2,".and. type N end -\n       sele resid ",$2,".and. type CA end -\n       sele resid ",$2,".and. type C end -\n       sele resid ",$2,".and. type N end -\n       sele resid ",$2,".and. type CA end -\n       sele resid ",$2,".and. type CB end -\n       sele resid ",$2,".and. type",(($1=="CYS")?"SG":(($1=="ILE")?"CG1":(($1=="SER")?"OG":(($1=="THR")?"OG1":(($1=="VAL")?"CG1":"CG"))))),"end -\nforc",forc,"npos",npos,"-\nxmin",xmin,"xmax",xmax,"xint",xint,"-\nymin",ymin,"ymax",ymax,"yint",yint,"-\npos -"; pot=sprintf("%s/%s.inp",ppp,$1); sss=""; while((getline < pot)>0){sss=sprintf("%s %s",sss,$0)}; na=split(sss,a); for(i=8;i<=na;i++){printf "%s%s",a[i],(((i-7)%6==0)?" -\n":"  ");close(pot)};printf "\n\n"}END{print "\nend\nreturn"}' $resresn > ${tarpdb:r}_${poten}.str

endif

######################################################################
#psixi1
set poten     = $potlist[3]

awk -v lr=$lres 'BEGIN {FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4"}{if($4==" CA "&&$9!=lr){print $6,$9}else{}}' $tmppdb | sed 's/HSD/HIS/g' | sed 's/HSP/HIS/g' | sed 's/HSE/HIS/g' | awk '{if($1=="ALA"||$1=="GLY"){}else{print $0}}' > $resresn

#122320 include an option for analysis calculation energy by a residue
if ( $analyres ) then
   awk -v forc=$forc -v ppp=${dbp}/${poten} -v npos=$npos -v xmin=$xmin -v xmax=$xmax -v xint=$xint -v ymin=$ymin -v ymax=$ymax -v yint=$yint 'BEGIN{print ""}{print "if @IN1 .eq. ",$2," then\nmod\nrots\nassign sele resid ",$2,".and. type N end -\n       sele resid ",$2,".and. type CA end -\n       sele resid ",$2,".and. type C end -\n       sele resid ",$2+1,".and. type N end -\n       sele resid ",$2,".and. type N end -\n       sele resid ",$2,".and. type CA end -\n       sele resid ",$2,".and. type CB end -\n       sele resid ",$2,".and. type",(($1=="CYS")?"SG":(($1=="ILE")?"CG1":(($1=="SER")?"OG":(($1=="THR")?"OG1":(($1=="VAL")?"CG1":"CG"))))),"end -\nforc",forc,"npos",npos,"-\nxmin",xmin,"xmax",xmax,"xint",xint,"-\nymin",ymin,"ymax",ymax,"yint",yint,"-\npos -"; pot=sprintf("%s/%s.inp",ppp,$1); sss=""; while((getline < pot)>0){sss=sprintf("%s %s",sss,$0)}; na=split(sss,a); for(i=8;i<=na;i++){printf "%s%s",a[i],(((i-7)%6==0)?" -\n":"  ");close(pot)};printf "\n\nend\nendif\n\n"}END{print "\nreturn"}' $resresn > ${tarpdb:r}_${poten}.str
else

   awk -v forc=$forc -v ppp=${dbp}/${poten} -v npos=$npos -v xmin=$xmin -v xmax=$xmax -v xint=$xint -v ymin=$ymin -v ymax=$ymax -v yint=$yint 'BEGIN{print "mod\nrots\n\n"}{print "assign sele resid ",$2,".and. type N end -\n       sele resid ",$2,".and. type CA end -\n       sele resid ",$2,".and. type C end -\n       sele resid ",$2+1,".and. type N end -\n       sele resid ",$2,".and. type N end -\n       sele resid ",$2,".and. type CA end -\n       sele resid ",$2,".and. type CB end -\n       sele resid ",$2,".and. type",(($1=="CYS")?"SG":(($1=="ILE")?"CG1":(($1=="SER")?"OG":(($1=="THR")?"OG1":(($1=="VAL")?"CG1":"CG"))))),"end -\nforc",forc,"npos",npos,"-\nxmin",xmin,"xmax",xmax,"xint",xint,"-\nymin",ymin,"ymax",ymax,"yint",yint,"-\npos -"; pot=sprintf("%s/%s.inp",ppp,$1); sss=""; while((getline < pot)>0){sss=sprintf("%s %s",sss,$0)}; na=split(sss,a); for(i=8;i<=na;i++){printf "%s%s",a[i],(((i-7)%6==0)?" -\n":"  ");close(pot)};printf "\n\n"}END{print "\nend\nreturn"}' $resresn > ${tarpdb:r}_${poten}.str

endif

######################################################################
#xi1xi2
set poten     = $potlist[4]

awk 'BEGIN {FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4"}{if($4==" CA "){print $6,$9}else{}}' $tmppdb | sed 's/HSD/HIS/g' | sed 's/HSP/HIS/g' | sed 's/HSE/HIS/g' | awk '{if($1=="ALA"||$1=="GLY"||$1=="CYS"||$1=="SER"||$1=="THR"||$1=="VAL"){}else{print $0}}' > $resresn

#122320 include an option for analysis calculation energy by a residue
if ( $analyres ) then
   awk -v forc=$forc -v ppp=${dbp}/${poten} -v npos=$npos -v xmin=$xmin -v xmax=$xmax -v xint=$xint -v ymin=$ymin -v ymax=$ymax -v yint=$yint 'BEGIN{print ""}{print "if @IN1 .eq. ",$2," then\nmod\nrots\nassign sele resid ",$2,".and. type N end -\n       sele resid ",$2,".and. type CA end -\n       sele resid ",$2,".and. type CB end -\n       sele resid ",$2,".and. type",(($1=="CYS")?"SG":(($1=="ILE")?"CG1":(($1=="SER")?"OG":(($1=="THR")?"OG1":(($1=="VAL")?"CG1":"CG"))))),"end -\n       sele resid ",$2,".and. type CA end -\n       sele resid ",$2,".and. type CB end -\n       sele resid ",$2,".and. type",(($1=="CYS")?"SG":(($1=="ILE")?"CG1":(($1=="SER")?"OG":(($1=="THR")?"OG1":(($1=="VAL")?"CG1":"CG"))))),"end -\n       sele resid ",$2,".and. type",(($1=="ARG")?"CD":(($1=="ASN")?"OD1":(($1=="ASP")?"OD1":(($1=="GLN")?"CD":(($1=="GLU")?"CD":(($1=="HIS")?"ND1":(($1=="ILE")?"CD":(($1=="LEU")?"CD1":(($1=="LYS")?"CD":(($1=="MET")?"SD":(($1=="PHE")?"CD1":(($1=="PRO")?"CD":(($1=="TRP")?"CD1":"CD1"))))))))))))),"end -\nforc",forc,"npos",npos,"-\nxmin",xmin,"xmax",xmax,"xint",xint,"-\nymin",ymin,"ymax",ymax,"yint",yint,"-\npos -"; pot=sprintf("%s/%s.inp",ppp,$1); sss=""; while((getline < pot)>0){sss=sprintf("%s %s",sss,$0)}; na=split(sss,a); for(i=8;i<=na;i++){printf "%s%s",a[i],(((i-7)%6==0)?" -\n":"  ");close(pot)};printf "\n\nend\nendif\n\n"}END{print "\nreturn"}' $resresn > ${tarpdb:r}_${poten}.str
else

   awk -v forc=$forc -v ppp=${dbp}/${poten} -v npos=$npos -v xmin=$xmin -v xmax=$xmax -v xint=$xint -v ymin=$ymin -v ymax=$ymax -v yint=$yint 'BEGIN{print "mod\nrots\n\n"}{print "assign sele resid ",$2,".and. type N end -\n       sele resid ",$2,".and. type CA end -\n       sele resid ",$2,".and. type CB end -\n       sele resid ",$2,".and. type",(($1=="CYS")?"SG":(($1=="ILE")?"CG1":(($1=="SER")?"OG":(($1=="THR")?"OG1":(($1=="VAL")?"CG1":"CG"))))),"end -\n       sele resid ",$2,".and. type CA end -\n       sele resid ",$2,".and. type CB end -\n       sele resid ",$2,".and. type",(($1=="CYS")?"SG":(($1=="ILE")?"CG1":(($1=="SER")?"OG":(($1=="THR")?"OG1":(($1=="VAL")?"CG1":"CG"))))),"end -\n       sele resid ",$2,".and. type",(($1=="ARG")?"CD":(($1=="ASN")?"OD1":(($1=="ASP")?"OD1":(($1=="GLN")?"CD":(($1=="GLU")?"CD":(($1=="HIS")?"ND1":(($1=="ILE")?"CD":(($1=="LEU")?"CD1":(($1=="LYS")?"CD":(($1=="MET")?"SD":(($1=="PHE")?"CD1":(($1=="PRO")?"CD":(($1=="TRP")?"CD1":"CD1"))))))))))))),"end -\nforc",forc,"npos",npos,"-\nxmin",xmin,"xmax",xmax,"xint",xint,"-\nymin",ymin,"ymax",ymax,"yint",yint,"-\npos -"; pot=sprintf("%s/%s.inp",ppp,$1); sss=""; while((getline < pot)>0){sss=sprintf("%s %s",sss,$0)}; na=split(sss,a); for(i=8;i<=na;i++){printf "%s%s",a[i],(((i-7)%6==0)?" -\n":"  ");close(pot)};printf "\n\n"}END{print "\nend\nreturn"}' $resresn > ${tarpdb:r}_${poten}.str

endif

if ( $dele ) rm -rf $tmpstr $resresn $ttmp >& /dev/null
