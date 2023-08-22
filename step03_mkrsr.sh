#!/bin/csh


# Setting the 'prog' variable to the name of the script without the path or extension.
set prog = `echo $0:r`

# Appending ".proc" to the 'prog' variable to create the 'pproc' variable.
set pproc = ${prog}".proc"

# Setting various parameters.
set dele = 1
set noersr = rsr/noersr
set dihrsr = rsr/dihrsr
set rawnoe = rawdata_parsing
set rawdih = rawdata_dihe
set rawdtm = dihe.$$
set backbonewidth = 1.0
set sidechainwidth = 0.0
set dcut = 0.5

# Setting to handle ambiguous restraints.
set ambiguous = 1   # 1: ambiguous ; 0: original (all individual)
set numambiguous = numambiguous.tmp
set noetype = noetype.dat
set typeambiguous = ambiguous.type

# Setting various paths and variables.
set workpath = `pwd`
set READ = ${workpath}/pReads
set SOFTPATH = ${workpath}/source 
set mol = $1
set tmpf = ${prog}.$$
set tmpdih = tmp.dihe
set model = `echo $mol | cut -d. -f 1`
set modell = `echo $model | awk '{print tolower($1)}'`
set modelu = `echo $modell | awk '{print toupper($1)}'`
set RESP = ${workpath}/pResults/${modelu}
set chain = `cat ${READ}/$mol | grep ^ATOM | awk '{print $5}' | sort | uniq | head -n 2 | tail -n 1`
set native = ${RESP}/$mol:r_n1.pdb
set targetb = `echo $mol | tr '[a-z]' '[A-Z]'`

# Checking if the ssbond.str file exists to handle ambiguous restraints.
if (-z $RESP/ssbond.str) then
    set ssbond = 1
else
    set ssbond = 0
endif

# Converting the 'chain' variable to lowercase.
set chain = `echo $chain | awk '{print tolower($0)}'`

# Creating a directory if it doesn't exist.
mkdir -p ${RESP}/rsr >& /dev/null

# Parsing raw data and converting them compatible with CHARMM topology.
${SOFTPATH}/rawdata_parsing.sh ${RESP}/rawdata > ${RESP}/$rawnoe
# Modifying atom names: HG -> HG1 in SER, HD1# -> HD# in ILE, HG -> HG1 in CYS, CD1 -> CD in ILE.
# These operations seem to be modifying atom names based on specific amino acid types.
# The 'num' variable stores residue numbers from the grep commands.

# Removing NOEs outside residue number range and specific atom names.
# The 'fres' and 'lres' variables store the first and last residue numbers.
# This section seems to process and filter NOE data.

# HG -> HG1 in SER
foreach num ( `grep -i "HG1 SER" $native | awk '{print $6}'` )
    awk -v num=$num '{if($3==num && $6=="HG"){print $1,$2,$3,$4,$5," HG1 ",$7,$8,$9,$10,$11,$12,$13,$14}else{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}}' ${RESP}/$rawnoe | awk -v num=$num '{if($8==num && $11=="HG"){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10," HG1 ",$12,$13,$14}else{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}}' > ${RESP}/$tmpf
    yes | cp -rf ${RESP}/$tmpf ${RESP}/$rawnoe >& /dev/null
end
# HD1# -> HD# in ILE
foreach num ( `grep -i "CD  ILE" $native | awk '{print $6}'` )
    awk -v num=${num} '{if($3==num && $6=="HD1#"){print $1,$2,$3,$4,$5," HD# ",$7,$8,$9,$10,$11,$12,$13,$14}else{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}}' ${RESP}/$rawnoe | awk -v num=${num} '{if($8==num && $11=="HD1#"){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10," HD# ",$12,$13,$14}else{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}}' > ${RESP}/$tmpf
    yes | cp -rf ${RESP}/$tmpf ${RESP}/$rawnoe >& /dev/null
end
# HG -> HG1 in CYS
foreach num (`grep -i "HG1 CYS" $native | awk '{print $6}'`)
    awk -v num=${num} '{if($3==num && $6=="HG"){print $1,$2,$3,$4,$5," HG1 ",$7,$8,$9,$10,$11,$12,$13,$14}else{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}}' ${RESP}/$rawnoe |  awk -v num=${num} '{if($8==num && $11=="HG"){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10," HG1 ",$12,$13,$14}else{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}}' > ${RESP}/$tmpf
    yes | cp -rf ${RESP}/$tmpf ${RESP}/$rawnoe >& /dev/null
end
# CD1 -> CD in ILE
foreach num (`grep -i "CD  ILE" $native | awk '{print $6}'`)
    awk -v num=${num} '{if($3==num && $6=="CD1"){print $1,$2,$3,$4,$5," CD ",$7,$8,$9,$10,$11,$12,$13,$14}else{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}}' ${RESP}/$rawnoe | awk -v num=${num} '{if($8==num && $11=="CD1"){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10," CD ",$12,$13,$14}else{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}}' > ${RESP}/$tmpf
    yes | cp -rf ${RESP}/$tmpf ${RESP}/$rawnoe >& /dev/null
end
# remove the NOEs outside residue number
set fres = `grep -i " CA " $native | head -n 1 | awk '{print $6}'`
set lres = `grep -i " CA " $native | tail -n 1 | awk '{print $6}'`
awk -v fn=${fres} -v ln=${lres} '{if($3>=fn && $3<=ln && $8>=fn && $8<=ln){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}else{}}' ${RESP}/$rawnoe | awk -v fn=${fres} -v ln=${lres} '{if(($3==fn && $6=="HN") || ($8==fn && $11=="HN")){}else{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}}'  > ${RESP}/$tmpf
yes | cp -rf ${RESP}/$tmpf ${RESP}/$rawnoe >& /dev/null

echo "Rawdata parsing done"


# add handling ambiguous restraints
if ( $ambiguous ) then

   # Count occurrences of NOE numbers and save to a file
   awk '{print $1}' ${RESP}/$rawnoe | uniq -c > ${RESP}/$numambiguous   # numbers noe_number

   # Calculate types of ambiguous restraints
   awk -i ${SOFTPATH}/lib.awk -v numam=${RESP}/$numambiguous 'BEGIN{while((getline < numam)>0){num[$2]=$1};resid1="";type1="";resid2="";type2="";j=0}{number=num[$1];if(number==1){print sprintf("%16d %25s %16d %25s %3sx%-3s",$3,$6,$8,$11,"1","1")}else{resid1=sprintf("%s%c%s",resid1,(resid1==""?" ":"_"),$3);type1=sprintf("%s%c%s",type1,(type1==""?" ":"_"),$6);resid2=sprintf("%s%c%s",resid2,(resid2==""?" ":"_"),$8);type2=sprintf("%s%c%s",type2,(type2==""?" ":"_"),$11);residtype1[sprintf("%s%s",$3,$6)]=0;residtype2[sprintf("%s%s",$8,$11)]=0;j++;if(j==number){print sprintf("%16s %25s %16s %25s %3dx%-3d",resid1,type1,resid2,type2,alen(residtype1),alen(residtype2));resid1="";type1="";resid2="";type2="";j=0;delete residtype1;delete residtype2}else{}}}' ${RESP}/$rawnoe > ${RESP}/$noetype

   awk '{print $NF}' ${RESP}/$noetype | sort | uniq -c > ${RESP}/$typeambiguous

   # Check for multiple x multiple assignment ; if so, stop process before CHARMM support this ; hard coding is required ; in this test, skip this target
   set mulmul = `awk '{na=split($2,a,"x");if(a[1]>=2&&a[2]>=2){existmulti=1}else{existmulti=0}}END{if(existmulti){print 1}else{print 0}}' ${RESP}/$typeambiguous`
   if ( $mulmul ) then
      echo "NOE restraints includes multiple x multiple type so it requires CHARMM new assignment by hard code"
      echo "$modelu" >> multimulti_rsrlist.txt
      exit
   else
      #all
      awk -i ${SOFTPATH}/lib.awk -v numam=${RESP}/$numambiguous -v chain=$chain 'BEGIN{print "*NOE \n*\n\nif @?rexp eq 0   set rexp = -0.166666666666667 \nif @?kmin eq 0   set kmin = 1 \nif @?kmax eq 0   set kmax = 1 \nif @?fmax eq 0   set fmax = 1 \nif @?rswi eq 0   set rswi = 3 \nif @?sexp eq 0   set sexp = 1 \nNOE";n=0;while((getline < numam)>0){num[$2]=$1};resid1="";type1="";resid2="";type2="";j=0}{number=num[$1];if(number==1){if(tolower($2)==chain&&tolower($7)==chain){print "assi sele resid ",$3," .and. type ",$6," end sele resid ",$8," .and. type  ",$11," end -\n"," rmin ",$12-$13," rmax ",$12+$14,"-\n"," rexp @rexp fmax @fmax rswi @rswi sexp @sexp kmin @kmin kmax @kmax SUMR ";n=n+1}}else{resid1=sprintf("%s%c%s",resid1,(resid1==""?" ":"_"),$3);type1=sprintf("%s%c%s",type1,(type1==""?" ":"_"),$6);resid2=sprintf("%s%c%s",resid2,(resid2==""?" ":"_"),$8);type2=sprintf("%s%c%s",type2,(type2==""?" ":"_"),$11);residtype1[sprintf("%s_%s",$3,$6)]=0;residtype2[sprintf("%s_%s",$8,$11)]=0;j++;if(j==number){tnum=0;for(i in residtype1){tnum++;na=split(i,a,"_");fres[tnum]=a[1];ftype[tnum]=a[2]};fnum=tnum;fstr="";for(i=1;i<=fnum;i++){fstr=sprintf("%s ( resid %s .and. type %s )",(length(fstr)==0?"":sprintf("%s .or. ",fstr)),fres[i],ftype[i])};tnum=0;for(i in residtype2){tnum++;na=split(i,a,"_");sres[tnum]=a[1];stype[tnum]=a[2]};snum=tnum;sstr="";for(i=1;i<=snum;i++){sstr=sprintf("%s ( resid %s .and. type %s )",(length(sstr)==0?"":sprintf("%s .or. ",sstr)),sres[i],stype[i])};if(tolower($2)==chain&&tolower($7)==chain){print "assi sele ",fstr," end -\n sele ",sstr," end - \n"," rmin ",$12-$13," rmax ",$12+$14,"-\n"," rexp @rexp fmax @fmax rswi @rswi sexp @sexp kmin @kmin kmax @kmax SUMR !",fnum,"x",snum;n=n+1;resid1="";type1="";resid2="";type2="";j=0;delete residtype1;delete residtype2}else{}}}}END{print "END\nreturn\n""! Number of NOEs:"n}' ${RESP}/$rawnoe > ${RESP}/${noersr}_ambi


      # Generate NOE restraints for different restraint classes
      # Intra-residue
      awk -i ${SOFTPATH}/lib.awk -v numam=${RESP}/$numambiguous -v chain=$chain 'BEGIN{print "*NOE \n*\n\nif @?rexp eq 0   set rexp = -0.166666666666667 \nif @?kmin eq 0   set kmin = 1 \nif @?kmax eq 0   set kmax = 1 \nif @?fmax eq 0   set fmax = 1 \nif @?rswi eq 0   set rswi = 3 \nif @?sexp eq 0   set sexp = 1 \nNOE";n=0;while((getline < numam)>0){num[$2]=$1};resid1="";type1="";resid2="";type2="";j=0}{if(abs($3-$8)==0){number=num[$1];if(number==1){if(tolower($2)==chain&&tolower($7)==chain){print "assi sele resid ",$3," .and. type ",$6," end sele resid ",$8," .and. type  ",$11," end -\n"," rmin ",$12-$13," rmax ",$12+$14,"-\n"," rexp @rexp fmax @fmax rswi @rswi sexp @sexp kmin @kmin kmax @kmax SUMR ";n=n+1}}else{resid1=sprintf("%s%c%s",resid1,(resid1==""?" ":"_"),$3);type1=sprintf("%s%c%s",type1,(type1==""?" ":"_"),$6);resid2=sprintf("%s%c%s",resid2,(resid2==""?" ":"_"),$8);type2=sprintf("%s%c%s",type2,(type2==""?" ":"_"),$11);residtype1[sprintf("%s_%s",$3,$6)]=0;residtype2[sprintf("%s_%s",$8,$11)]=0;j++;if(j==number){tnum=0;for(i in residtype1){tnum++;na=split(i,a,"_");fres[tnum]=a[1];ftype[tnum]=a[2]};fnum=tnum;fstr="";for(i=1;i<=fnum;i++){fstr=sprintf("%s ( resid %s .and. type %s )",(length(fstr)==0?"":sprintf("%s .or. ",fstr)),fres[i],ftype[i])};tnum=0;for(i in residtype2){tnum++;na=split(i,a,"_");sres[tnum]=a[1];stype[tnum]=a[2]};snum=tnum;sstr="";for(i=1;i<=snum;i++){sstr=sprintf("%s ( resid %s .and. type %s )",(length(sstr)==0?"":sprintf("%s .or. ",sstr)),sres[i],stype[i])};if(tolower($2)==chain&&tolower($7)==chain){print "assi sele ",fstr," end -\n sele ",sstr," end - \n"," rmin ",$12-$13," rmax ",$12+$14,"-\n"," rexp @rexp fmax @fmax rswi @rswi sexp @sexp kmin @kmin kmax @kmax SUMR !",fnum,"x",snum;n=n+1;resid1="";type1="";resid2="";type2="";j=0;delete residtype1;delete residtype2}else{}}}}else{}}END{print "END\nreturn\n""! Number of NOEs:"n}' ${RESP}/$rawnoe > ${RESP}/${noersr}_intra_ambi
      
      # Sequential (|i-j|=1)
      awk -i ${SOFTPATH}/lib.awk -v numam=${RESP}/$numambiguous -v chain=$chain 'BEGIN{print "*NOE \n*\n\nif @?rexp eq 0   set rexp = -0.166666666666667 \nif @?kmin eq 0   set kmin = 1 \nif @?kmax eq 0   set kmax = 1 \nif @?fmax eq 0   set fmax = 1 \nif @?rswi eq 0   set rswi = 3 \nif @?sexp eq 0   set sexp = 1 \nNOE";n=0;while((getline < numam)>0){num[$2]=$1};resid1="";type1="";resid2="";type2="";j=0}{if(abs($3-$8)==1){number=num[$1];if(number==1){if(tolower($2)==chain&&tolower($7)==chain){print "assi sele resid ",$3," .and. type ",$6," end sele resid ",$8," .and. type  ",$11," end -\n"," rmin ",$12-$13," rmax ",$12+$14,"-\n"," rexp @rexp fmax @fmax rswi @rswi sexp @sexp kmin @kmin kmax @kmax SUMR ";n=n+1}}else{resid1=sprintf("%s%c%s",resid1,(resid1==""?" ":"_"),$3);type1=sprintf("%s%c%s",type1,(type1==""?" ":"_"),$6);resid2=sprintf("%s%c%s",resid2,(resid2==""?" ":"_"),$8);type2=sprintf("%s%c%s",type2,(type2==""?" ":"_"),$11);residtype1[sprintf("%s_%s",$3,$6)]=0;residtype2[sprintf("%s_%s",$8,$11)]=0;j++;if(j==number){tnum=0;for(i in residtype1){tnum++;na=split(i,a,"_");fres[tnum]=a[1];ftype[tnum]=a[2]};fnum=tnum;fstr="";for(i=1;i<=fnum;i++){fstr=sprintf("%s ( resid %s .and. type %s )",(length(fstr)==0?"":sprintf("%s .or. ",fstr)),fres[i],ftype[i])};tnum=0;for(i in residtype2){tnum++;na=split(i,a,"_");sres[tnum]=a[1];stype[tnum]=a[2]};snum=tnum;sstr="";for(i=1;i<=snum;i++){sstr=sprintf("%s ( resid %s .and. type %s )",(length(sstr)==0?"":sprintf("%s .or. ",sstr)),sres[i],stype[i])};if(tolower($2)==chain&&tolower($7)==chain){print "assi sele ",fstr," end -\n sele ",sstr," end - \n"," rmin ",$12-$13," rmax ",$12+$14,"-\n"," rexp @rexp fmax @fmax rswi @rswi sexp @sexp kmin @kmin kmax @kmax SUMR !",fnum,"x",snum;n=n+1;resid1="";type1="";resid2="";type2="";j=0;delete residtype1;delete residtype2}else{}}}}else{}}END{print "END\nreturn\n""! Number of NOEs:"n}' ${RESP}/$rawnoe > ${RESP}/${noersr}_seq_ambi

      # Medium range (2<=|i-j|<5)
      awk -i ${SOFTPATH}/lib.awk -v numam=${RESP}/$numambiguous -v chain=$chain 'BEGIN{print "*NOE \n*\n\nif @?rexp eq 0   set rexp = -0.166666666666667 \nif @?kmin eq 0   set kmin = 1 \nif @?kmax eq 0   set kmax = 1 \nif @?fmax eq 0   set fmax = 1 \nif @?rswi eq 0   set rswi = 3 \nif @?sexp eq 0   set sexp = 1 \nNOE";n=0;while((getline < numam)>0){num[$2]=$1};resid1="";type1="";resid2="";type2="";j=0}{if(abs($3-$8)>=2&&abs($3-$8)<=4){number=num[$1];if(number==1){if(tolower($2)==chain&&tolower($7)==chain){print "assi sele resid ",$3," .and. type ",$6," end sele resid ",$8," .and. type  ",$11," end -\n"," rmin ",$12-$13," rmax ",$12+$14,"-\n"," rexp @rexp fmax @fmax rswi @rswi sexp @sexp kmin @kmin kmax @kmax SUMR ";n=n+1}}else{resid1=sprintf("%s%c%s",resid1,(resid1==""?" ":"_"),$3);type1=sprintf("%s%c%s",type1,(type1==""?" ":"_"),$6);resid2=sprintf("%s%c%s",resid2,(resid2==""?" ":"_"),$8);type2=sprintf("%s%c%s",type2,(type2==""?" ":"_"),$11);residtype1[sprintf("%s_%s",$3,$6)]=0;residtype2[sprintf("%s_%s",$8,$11)]=0;j++;if(j==number){tnum=0;for(i in residtype1){tnum++;na=split(i,a,"_");fres[tnum]=a[1];ftype[tnum]=a[2]};fnum=tnum;fstr="";for(i=1;i<=fnum;i++){fstr=sprintf("%s ( resid %s .and. type %s )",(length(fstr)==0?"":sprintf("%s .or. ",fstr)),fres[i],ftype[i])};tnum=0;for(i in residtype2){tnum++;na=split(i,a,"_");sres[tnum]=a[1];stype[tnum]=a[2]};snum=tnum;sstr="";for(i=1;i<=snum;i++){sstr=sprintf("%s ( resid %s .and. type %s )",(length(sstr)==0?"":sprintf("%s .or. ",sstr)),sres[i],stype[i])};if(tolower($2)==chain&&tolower($7)==chain){print "assi sele ",fstr," end -\n sele ",sstr," end - \n"," rmin ",$12-$13," rmax ",$12+$14,"-\n"," rexp @rexp fmax @fmax rswi @rswi sexp @sexp kmin @kmin kmax @kmax SUMR !",fnum,"x",snum;n=n+1;resid1="";type1="";resid2="";type2="";j=0;delete residtype1;delete residtype2}else{}}}}else{}}END{print "END\nreturn\n""! Number of NOEs:"n}' ${RESP}/$rawnoe > ${RESP}/${noersr}_med_ambi

      # Long range (|i-j|>=5)
      awk -i ${SOFTPATH}/lib.awk -v numam=${RESP}/$numambiguous -v chain=$chain 'BEGIN{print "*NOE \n*\n\nif @?rexp eq 0   set rexp = -0.166666666666667 \nif @?kmin eq 0   set kmin = 1 \nif @?kmax eq 0   set kmax = 1 \nif @?fmax eq 0   set fmax = 1 \nif @?rswi eq 0   set rswi = 3 \nif @?sexp eq 0   set sexp = 1 \nNOE";n=0;while((getline < numam)>0){num[$2]=$1};resid1="";type1="";resid2="";type2="";j=0}{if(abs($3-$8)>=5){number=num[$1];if(number==1){if(tolower($2)==chain&&tolower($7)==chain){print "assi sele resid ",$3," .and. type ",$6," end sele resid ",$8," .and. type  ",$11," end -\n"," rmin ",$12-$13," rmax ",$12+$14,"-\n"," rexp @rexp fmax @fmax rswi @rswi sexp @sexp kmin @kmin kmax @kmax SUMR ";n=n+1}}else{resid1=sprintf("%s%c%s",resid1,(resid1==""?" ":"_"),$3);type1=sprintf("%s%c%s",type1,(type1==""?" ":"_"),$6);resid2=sprintf("%s%c%s",resid2,(resid2==""?" ":"_"),$8);type2=sprintf("%s%c%s",type2,(type2==""?" ":"_"),$11);residtype1[sprintf("%s_%s",$3,$6)]=0;residtype2[sprintf("%s_%s",$8,$11)]=0;j++;if(j==number){tnum=0;for(i in residtype1){tnum++;na=split(i,a,"_");fres[tnum]=a[1];ftype[tnum]=a[2]};fnum=tnum;fstr="";for(i=1;i<=fnum;i++){fstr=sprintf("%s ( resid %s .and. type %s )",(length(fstr)==0?"":sprintf("%s .or. ",fstr)),fres[i],ftype[i])};tnum=0;for(i in residtype2){tnum++;na=split(i,a,"_");sres[tnum]=a[1];stype[tnum]=a[2]};snum=tnum;sstr="";for(i=1;i<=snum;i++){sstr=sprintf("%s ( resid %s .and. type %s )",(length(sstr)==0?"":sprintf("%s .or. ",sstr)),sres[i],stype[i])};if(tolower($2)==chain&&tolower($7)==chain){print "assi sele ",fstr," end -\n sele ",sstr," end - \n"," rmin ",$12-$13," rmax ",$12+$14,"-\n"," rexp @rexp fmax @fmax rswi @rswi sexp @sexp kmin @kmin kmax @kmax SUMR !",fnum,"x",snum;n=n+1;resid1="";type1="";resid2="";type2="";j=0;delete residtype1;delete residtype2}else{}}}}else{}}END{print "END\nreturn\n""! Number of NOEs:"n}' ${RESP}/$rawnoe > ${RESP}/${noersr}_long_ambi
    endif

   # Extract and count unique NOE numbers
   awk -v chain=$chain 'BEGIN{print "*NOE \n*\n\nif @?rexp eq 0   set rexp = -0.166666666666667 \nif @?kmin eq 0   set kmin = 1 \nif @?kmax eq 0   set kmax = 1 \nif @?fmax eq 0   set fmax = 1 \nif @?rswi eq 0   set rswi = 3 \nif @?sexp eq 0   set sexp = 1 \nNOE";n=0}{if(tolower($2)==chain&&tolower($7)==chain){print "assi sele resid ",$3," .and. type ",$6," end sele resid ",$8," .and. type  ",$11," end - \n"," rmin ",$12-$13," rmax ",$12+$14,"- \n "," rexp @rexp fmax @fmax rswi @rswi sexp @sexp kmin @kmin kmax @kmax SUMR !",n+1;n=n+1}}END{print "END\nreturn\n""! Number of NOEs:"n}' ${RESP}/$rawnoe > ${RESP}/$noersr

   #intra
   awk -v chain=$chain 'BEGIN{print "*noe \n*\n\nif @?rexp eq 0   set rexp = -0.166666666666667 \nif @?kmin eq 0   set kmin = 1 \nif @?kmax eq 0   set kmax = 1 \nif @?fmax eq 0   set fmax = 1 \nif @?rswi eq 0   set rswi = 3 \nif @?sexp eq 0   set sexp = 1 \nNOE";n=0}{if(tolower($2)==chain&&tolower($7)==chain){if($3-$8==0){print "assi sele resid ",$3," .and. type ",$6," end sele resid ",$8," .and. type  ",$11," end - \n"," rmin ",$12-$13," rmax ",$12+$14,"- \n "," rexp @rexp fmax @fmax rswi @rswi sexp @sexp kmin @kmin kmax @kmax SUMR ";n=n+1}else{}}}END{print "END\nreturn\n""! Number of NOEs:"n}' ${RESP}/$rawnoe > ${RESP}/${noersr}_intra

   #sequential
   awk -v chain=$chain 'BEGIN{print "*noe \n*\n\nif @?rexp eq 0   set rexp = -0.166666666666667 \nif @?kmin eq 0   set kmin = 1 \nif @?kmax eq 0   set kmax = 1 \nif @?fmax eq 0   set fmax = 1 \nif @?rswi eq 0   set rswi = 3 \nif @?sexp eq 0   set sexp = 1 \nNOE";n=0}{if(tolower($2)==chain&&tolower($7)==chain){if(($3-$8)==1||($3-$8)==-1){print "assi sele resid ",$3," .and. type ",$6," end sele resid ",$8," .and. type  ",$11," end - \n"," rmin ",$12-$13," rmax ",$12+$14,"- \n "," rexp @rexp fmax @fmax rswi @rswi sexp @sexp kmin @kmin kmax @kmax SUMR ";n=n+1}else{}}}END{print "END\nreturn\n""! Number of NOEs:"n}' ${RESP}/$rawnoe > ${RESP}/${noersr}_seq

   #medium
   awk -v chain=$chain 'BEGIN{print "*noe \n*\n\nif @?rexp eq 0   set rexp = -0.166666666666667 \nif @?kmin eq 0   set kmin = 1 \nif @?kmax eq 0   set kmax = 1 \nif @?fmax eq 0   set fmax = 1 \nif @?rswi eq 0   set rswi = 3 \nif @?sexp eq 0   set sexp = 1 \nNOE";n=0}{if(tolower($2)==chain&&tolower($7)==chain){if((($3-$8)>=2&&($3-$8)<5)||(($3-$8)>-5&&($3-$8)<=-2)){print "assi sele resid ",$3," .and. type ",$6," end sele resid ",$8," .and. type  ",$11," end - \n"," rmin ",$12-$13," rmax ",$12+$14,"- \n "," rexp @rexp fmax @fmax rswi @rswi sexp @sexp kmin @kmin kmax @kmax SUMR ";n=n+1}else{}}}END{print "END\nreturn\n""! Number of NOEs:"n}' ${RESP}/$rawnoe > ${RESP}/${noersr}_med

   #long
   awk -v chain=$chain 'BEGIN{print "*noe \n*\n\nif @?rexp eq 0   set rexp = -0.166666666666667 \nif @?kmin eq 0   set kmin = 1 \nif @?kmax eq 0   set kmax = 1 \nif @?fmax eq 0   set fmax = 1 \nif @?rswi eq 0   set rswi = 3 \nif @?sexp eq 0   set sexp = 1 \nNOE";n=0}{if(tolower($2)==chain&&tolower($7)==chain){if(($3-$8)>=5||($3-$8)<=-5){print "assi sele resid ",$3," .and. type ",$6," end sele resid ",$8," .and. type  ",$11," end - \n"," rmin ",$12-$13," rmax ",$12+$14,"- \n "," rexp @rexp fmax @fmax rswi @rswi sexp @sexp kmin @kmin kmax @kmax SUMR ";n=n+1}else{}}}END{print "END\nreturn\n""! Number of NOEs:"n}' ${RESP}/$rawnoe > ${RESP}/${noersr}_long

endif


# Check if the input file exists and is not empty
if ( -s ${RESP}/$rawdih ) then

    awk 'BEGIN{RS=" ASSIGN"}{if(NR!=1){printf(" ASSIGN %s resid %s and name %s resid %s and name %s resid %s and name %s resid %s and name %s %s %s %s %s\n",substr($3,1,length($3)-1),$6,$9,$16,$19,$26,$29,$36,$39,$41,$42,$43,$44,$45)}}' ${RESP}/$rawdih > ${RESP}/$rawdtm

    # Change CD2 in ILE to CD compatible with CHARMM topology
    foreach num ( `grep -i "CD  ILE" $native | awk '{print $6}'` )

    # Process data based on conditions and copy to tmpf
    awk -v num=${num} '{if($4==num && $7=="CD1"){print $1,$2,$3,$4,$5,$6," CD ",$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27}else{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}}' ${RESP}/$rawdtm |  awk -v num=${num} '{if($9==num && $12=="CD1"){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11," CD ",$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}else{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}}' |  awk -v num=${num} '{if($14==num && $17=="CD1"){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16," CD ",$18,$19,$20,$21,$22,$23,$24,$25,$26}else{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}}' |  awk -v num=${num} '{if($19==num && $22=="CD1"){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21," CD ",$23,$24,$25,$25,$26}else{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}}' > ${RESP}/$tmpf
        yes | cp -rf ${RESP}/$tmpf ${RESP}/$rawdtm >& /dev/null
    end

    # Process data based on chain, residue numbers, and types
    # Remove the NOEs outside residue number
set chain = `echo $chain | awk '{print toupper($0)}'`

                 set fres = `grep -i " CA " $native | head -n 1 | awk '{print $6}'`
                 set lres = `grep -i " CA " $native | tail -n 1 | awk '{print $6}'`
     echo $chain
    rm -rf   ${RESP}/$tmpf
    awk -v fn=${fres} -v ln=${lres} -v c=$chain '{if($2 == c && $4>=fn && $4<=ln && $9>=fn && $9<=ln && $14>=fn && $14<=ln && $19>=fn && $19<=ln ){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}else{}}' ${RESP}/$rawdtm >> ${RESP}/$tmpf

    awk 'BEGIN{pi=3.1415926535897932384626433832795;print "*dihe \n*\n\nMOD\nmper\n\n";n=0}{print "assign sele segid","pro0",".and. resid ",$4," .and. type ",$7," end -\n       sele segid","pro0",".and. resid ",$9," .and. type  ",$12," end -\n       sele segid","pro0",".and. resid ",$14," .and. type  ",$17, " end -\n       sele segid","pro0",".and. resid ",$19," .and. type  ",$22, " end -\n","      forc @kdihe nmin 1 fbsp -\n       min ",$24*pi/180.0," -\n       sigma ",$25*pi/180.0,"\n";n=n+1}END{print "\nEND\nreturn\n""! Number of DIHEs:"n}' ${RESP}/$tmpf > ${RESP}/$dihrsr
endif


# Make STAP restraints
${SOFTPATH}/mkSTAP.sh $native no2 1.0 #
yes | mv ${RESP}/${modelu}_n1_phichi1.str  ${RESP}/rsr >& /dev/null
yes | mv ${RESP}/${modelu}_n1_phipsi.str   ${RESP}/rsr >& /dev/null
yes | mv ${RESP}/${modelu}_n1_psichi1.str  ${RESP}/rsr >& /dev/null
yes | mv ${RESP}/${modelu}_n1_chi1chi2.str ${RESP}/rsr >& /dev/null

if ( $dele ) rm -rf ${RESP}/$rawnoe ${RESP}/$numambiguous ${workpath}/tmp.pdb.$$ >& /dev/null

if ( ( -e ${RESP}/rsr/${modelu}_n1_phichi1.str || -e ${RESP}/$dihrsr || -e ${RESP}/$noersr ) ) then
    echo "${prog}:passed 1" > ${RESP}/$pproc
else
    echo "${prog}:failed 0" > ${RESP}/$pproc
endif

echo "$modelu step03_mkrsr.sh done"

