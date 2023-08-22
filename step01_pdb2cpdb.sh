#!/bin/csh

# This script processes PDB files to modify atom names and coordinates based on specific criteria.
# It also deals with residue substitutions and SS bond patches.




# Set the working path where the script operates.
set workpath = `pwd` 
set softpath = ${workpath}/source
set READ = ${workpath}/pReads




# Set script-related variables.
set prog = `echo $0:r`            # Extract the script name.
set pproc = ${prog}".proc"         # Create a process file.
set dele = 1                       # A flag indicating whether to delete temporary files.
set mol = $1                       # The PDB file to process.
set hetero = 0                     # Flag to consider hetero atoms (0: no, 1: yes).




# Check if a PDB filename is provided.
if ( $mol == "" ) then
    echo $prog" [PDBfilename]"
    exit
endif



# Extract the chain information from the PDB file.
set chain  = `cat ${READ}/$mol |grep ^ATOM | awk '{print $5}' | uniq | head -n 1 `



# Extract the model information from the PDB filename.
set model  = `echo $mol | cut -d. -f 1`
set modell  = `echo $model  | awk '{print tolower($1)}'`
set modelu  = `echo $modell | awk '{print toupper($1)}'`



# Extract the number of models in the PDB file.
set tmol = `awk 'BEGIN{n=0}{if($1=="MODEL"){n=n+1}else{}}END{if(n==0){print 1}else{print n}}' ${READ}/$mol`
set tmppdb  = tmppdb.$$
set tmppdb2 = tmppdb2.$$
set tmptxt  = tmptxt.$$
set RESP = ${workpath}/pResults/${modelu}




# Create necessary directories.
mkdir -p ${workpath}/pResults/tmp
mkdir -p ${workpath}/pResults/${modelu}
grep -i "^MODRES" ${READ}/$mol > ${RESP}/modres.txt
set nmod = `cat ${RESP}/modres.txt | wc -l | awk '{print $1}'`
#mkdir native >& /dev/null




# Display model count.
echo "$mol number of model: $tmol"



# Loop through each model.
foreach num (`seq -s" " 1 ${tmol}`)



    # Check if a flag is set to consider hetero atoms and perform necessary transformations.
    if ( $tmol == 1 ) then

        # Transform atom names for various residues.
        # Note: This section contains a series of 'sed' commands for replacing specific atom names.
        if ( $hetero ) then
           $softpath/toolset/perl/convpdb.pl -segnames -chain $chain < ${READ}/$mol | sed -e "s/HB3 ARG/HB1 ARG/g"  | sed -e "s/HG  CYS/HG1 CYS/g" | sed -e "s/HB3 CYS/HB1 CYS/g" | sed -e "s/HB3 ASN/HB1 ASN/g" | sed -e "s/HB3 ASP/HB1 ASP/g" | sed -e "s/HB3 GLN/HB1 GLN/g"  | sed -e "s/HB3 GLU/HB1 GLU/g" | sed -e "s/HB3 HSD/HB1 HSD/g"  | sed -e "s/HB3 HSP/HB1 HSP/g"  | sed -e "s/HB3 HSE/HB1 HSE/g"  | sed -e "s/HB3 LEU/HB1 LEU/g" | sed -e "s/HB3 LYS/HB1 LYS/g" | sed -e "s/HB3 MET/HB1 MET/g"  | sed -e "s/HB3 PHE/HB1 PHE/g" | sed -e "s/HB3 PRO/HB1 PRO/g"  | sed -e "s/HB3 SER/HB1 SER/g"  | sed -e "s/HB3 TRP/HB1 TRP/g" | sed -e "s/HB3 TYR/HB1 TYR/g"  | sed -e "s/HG3 ARG/HG1 ARG/g" | sed -e "s/HG3 GLN/HG1 GLN/g" | sed -e "s/HG3 GLU/HG1 GLU/g" | sed -e "s/HG3 LYS/HG1 LYS/g"  | sed -e "s/HG3 MET/HG1 MET/g" | sed -e "s/HG3 PRO/HG1 PRO/g"  | sed -e "s/HD3 ARG/HD1 ARG/g" | sed -e "s/HD3 LYS/HD1 LYS/g" | sed -e "s/HD3 PRO/HD1 PRO/g" | sed -e "s/HE3 LYS/HE1 LYS/g" | sed -e "s/HG13 ILE/HG11 ILE/g" | sed -e "s/HD11 ILE/ HD1 ILE/g" | sed -e "s/HD12 ILE/ HD2 ILE/g" | sed -e "s/HD13 ILE/ HD3 ILE/g"  | sed -e "s/HA3 GLY/HA1 GLY/g" | sed -e "s/HG  SER/HG1 SER/g" | sed -e "s/PROA/PRO0/g" | sed -e "s/P001/PRO0/g" | sed -e "s/P002/PRO0/g" | sed -e "s/P003/PRO0/g" | awk 'BEGIN {FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4"}{printf "%-6s%5s%1s%4s%1s%3s%1s%1s%4s%1s%3s%8s%8s%8s%6s%6s%6s%4s\n",$1,$2,$3,($4==" H  "?" HN ":$4),$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' | sed -e "s/HD22 ASN/HD23 ASN/g" | sed -e "s/HD21 ASN/HD22 ASN/g" | sed -e "s/HD23 ASN/HD21 ASN/g" > ${RESP}/$tmppdb
        else
            $softpath/toolset/perl/convpdb.pl -nohetero -segnames -chain $chain < ${READ}/$mol | sed -e "s/HB3 ARG/HB1 ARG/g"  | sed -e "s/HG  CYS/HG1 CYS/g" | sed -e "s/HB3 CYS/HB1 CYS/g" | sed -e "s/HB3 ASN/HB1 ASN/g" | sed -e "s/HB3 ASP/HB1 ASP/g" | sed -e "s/HB3 GLN/HB1 GLN/g"  | sed -e "s/HB3 GLU/HB1 GLU/g" | sed -e "s/HB3 HSD/HB1 HSD/g"  | sed -e "s/HB3 HSP/HB1 HSP/g"  | sed -e "s/HB3 HSE/HB1 HSE/g"  | sed -e "s/HB3 LEU/HB1 LEU/g" | sed -e "s/HB3 LYS/HB1 LYS/g" | sed -e "s/HB3 MET/HB1 MET/g"  | sed -e "s/HB3 PHE/HB1 PHE/g" | sed -e "s/HB3 PRO/HB1 PRO/g"  | sed -e "s/HB3 SER/HB1 SER/g"  | sed -e "s/HB3 TRP/HB1 TRP/g" | sed -e "s/HB3 TYR/HB1 TYR/g"  | sed -e "s/HG3 ARG/HG1 ARG/g" | sed -e "s/HG3 GLN/HG1 GLN/g" | sed -e "s/HG3 GLU/HG1 GLU/g" | sed -e "s/HG3 LYS/HG1 LYS/g"  | sed -e "s/HG3 MET/HG1 MET/g" | sed -e "s/HG3 PRO/HG1 PRO/g"  | sed -e "s/HD3 ARG/HD1 ARG/g" | sed -e "s/HD3 LYS/HD1 LYS/g" | sed -e "s/HD3 PRO/HD1 PRO/g" | sed -e "s/HE3 LYS/HE1 LYS/g" | sed -e "s/HG13 ILE/HG11 ILE/g" | sed -e "s/HD11 ILE/ HD1 ILE/g" | sed -e "s/HD12 ILE/ HD2 ILE/g" | sed -e "s/HD13 ILE/ HD3 ILE/g"  | sed -e "s/HA3 GLY/HA1 GLY/g" | sed -e "s/HG  SER/HG1 SER/g" | sed -e "s/PROA/PRO0/g" | sed -e "s/P001/PRO0/g" | sed -e "s/P002/PRO0/g" | sed -e "s/P003/PRO0/g" | awk 'BEGIN {FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4"}{printf "%-6s%5s%1s%4s%1s%3s%1s%1s%4s%1s%3s%8s%8s%8s%6s%6s%6s%4s\n",$1,$2,$3,($4==" H  "?" HN ":$4),$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' | sed -e "s/HD22 ASN/HD23 ASN/g" | sed -e "s/HD21 ASN/HD22 ASN/g" | sed -e "s/HD23 ASN/HD21 ASN/g" > ${RESP}/$tmppdb
        endif
    else
        if ( $hetero ) then
           $softpath/toolset/perl/convpdb.pl -segnames -model ${num} -chain $chain < ${READ}/$mol | sed -e "s/HB3 ARG/HB1 ARG/g"  | sed -e "s/HG  CYS/HG1 CYS/g" | sed -e "s/HB3 CYS/HB1 CYS/g" | sed -e "s/HB3 ASN/HB1 ASN/g" | sed -e "s/HB3 ASP/HB1 ASP/g" | sed -e "s/HB3 GLN/HB1 GLN/g"  | sed -e "s/HB3 GLU/HB1 GLU/g" | sed -e "s/HB3 HSD/HB1 HSD/g"  | sed -e "s/HB3 HSP/HB1 HSP/g"  | sed -e "s/HB3 HSE/HB1 HSE/g"  | sed -e "s/HB3 LEU/HB1 LEU/g" | sed -e "s/HB3 LYS/HB1 LYS/g" | sed -e "s/HB3 MET/HB1 MET/g"  | sed -e "s/HB3 PHE/HB1 PHE/g" | sed -e "s/HB3 PRO/HB1 PRO/g"  | sed -e "s/HB3 SER/HB1 SER/g"  | sed -e "s/HB3 TRP/HB1 TRP/g" | sed -e "s/HB3 TYR/HB1 TYR/g"  | sed -e "s/HG3 ARG/HG1 ARG/g" | sed -e "s/HG3 GLN/HG1 GLN/g" | sed -e "s/HG3 GLU/HG1 GLU/g" | sed -e "s/HG3 LYS/HG1 LYS/g"  | sed -e "s/HG3 MET/HG1 MET/g" | sed -e "s/HG3 PRO/HG1 PRO/g"  | sed -e "s/HD3 ARG/HD1 ARG/g" | sed -e "s/HD3 LYS/HD1 LYS/g" | sed -e "s/HD3 PRO/HD1 PRO/g" | sed -e "s/HE3 LYS/HE1 LYS/g" | sed -e "s/HG13 ILE/HG11 ILE/g" | sed -e "s/HD11 ILE/ HD1 ILE/g" | sed -e "s/HD12 ILE/ HD2 ILE/g" | sed -e "s/HD13 ILE/ HD3 ILE/g"  | sed -e "s/HA3 GLY/HA1 GLY/g" | sed -e "s/HG  SER/HG1 SER/g" | sed -e "s/PROA/PRO0/g" | sed -e "s/P001/PRO0/g" | sed -e "s/P002/PRO0/g" | sed -e "s/P003/PRO0/g" | awk 'BEGIN {FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4"}{printf "%-6s%5s%1s%4s%1s%3s%1s%1s%4s%1s%3s%8s%8s%8s%6s%6s%6s%4s\n",$1,$2,$3,($4==" H  "?" HN ":$4),$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' | sed -e "s/HD22 ASN/HD23 ASN/g" | sed -e "s/HD21 ASN/HD22 ASN/g" | sed -e "s/HD23 ASN/HD21 ASN/g" > ${RESP}/$tmppdb
        else
           $softpath/toolset/perl/convpdb.pl -nohetero -segnames -model ${num} -chain $chain < ${READ}/$mol | sed -e "s/HB3 ARG/HB1 ARG/g"  | sed -e "s/HG  CYS/HG1 CYS/g" | sed -e "s/HB3 CYS/HB1 CYS/g" | sed -e "s/HB3 ASN/HB1 ASN/g" | sed -e "s/HB3 ASP/HB1 ASP/g" | sed -e "s/HB3 GLN/HB1 GLN/g"  | sed -e "s/HB3 GLU/HB1 GLU/g" | sed -e "s/HB3 HSD/HB1 HSD/g"  | sed -e "s/HB3 HSP/HB1 HSP/g"  | sed -e "s/HB3 HSE/HB1 HSE/g"  | sed -e "s/HB3 LEU/HB1 LEU/g" | sed -e "s/HB3 LYS/HB1 LYS/g" | sed -e "s/HB3 MET/HB1 MET/g"  | sed -e "s/HB3 PHE/HB1 PHE/g" | sed -e "s/HB3 PRO/HB1 PRO/g"  | sed -e "s/HB3 SER/HB1 SER/g"  | sed -e "s/HB3 TRP/HB1 TRP/g" | sed -e "s/HB3 TYR/HB1 TYR/g"  | sed -e "s/HG3 ARG/HG1 ARG/g" | sed -e "s/HG3 GLN/HG1 GLN/g" | sed -e "s/HG3 GLU/HG1 GLU/g" | sed -e "s/HG3 LYS/HG1 LYS/g"  | sed -e "s/HG3 MET/HG1 MET/g" | sed -e "s/HG3 PRO/HG1 PRO/g"  | sed -e "s/HD3 ARG/HD1 ARG/g" | sed -e "s/HD3 LYS/HD1 LYS/g" | sed -e "s/HD3 PRO/HD1 PRO/g" | sed -e "s/HE3 LYS/HE1 LYS/g" | sed -e "s/HG13 ILE/HG11 ILE/g" | sed -e "s/HD11 ILE/ HD1 ILE/g" | sed -e "s/HD12 ILE/ HD2 ILE/g" | sed -e "s/HD13 ILE/ HD3 ILE/g"  | sed -e "s/HA3 GLY/HA1 GLY/g" | sed -e "s/HG  SER/HG1 SER/g" | sed -e "s/PROA/PRO0/g" | sed -e "s/P001/PRO0/g" | sed -e "s/P002/PRO0/g" | sed -e "s/P003/PRO0/g" | awk 'BEGIN {FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4"}{printf "%-6s%5s%1s%4s%1s%3s%1s%1s%4s%1s%3s%8s%8s%8s%6s%6s%6s%4s\n",$1,$2,$3,($4==" H  "?" HN ":$4),$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' | sed -e "s/HD22 ASN/HD23 ASN/g" | sed -e "s/HD21 ASN/HD22 ASN/g" | sed -e "s/HD23 ASN/HD21 ASN/g" > ${RESP}/$tmppdb
        endif
    endif

    #The lastest MMTSB discriminate three types of HIS
    #if HSD with HD1, HD2, and HE1 must be HSD; if HSD with HE1, HE2, and HD2 must be HSE; if HSD with HE1, HE2, HD1, and HD2 must be HSP
    awk '{if($4=="HSD" && $3=="HD1"){print $6}else{}}' ${RESP}/$tmppdb > ${RESP}/$tmptxt

    foreach tmp (`awk '{if($4=="HSD" && $3=="HE2"){print $6}else{}}' ${RESP}/$tmppdb`)
        set exist = `grep "${tmp}" ${RESP}/$tmptxt | wc -l`
        if ( $exist ) then
            #echo "Residue:"${tmp}"is HSP"
            awk -v val=${tmp} 'BEGIN {FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4"}{if($9==val){printf "%-6s%5s%1s%4s%1s%3s%1s%1s%4s%1s%3s%8s%8s%8s%6s%6s%6s%4s\n",$1,$2,$3,$4,$5,"HSP",$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}else{print }}' ${RESP}/$tmppdb > ${RESP}/$tmppdb2
        else
            #echo "Residue:"${tmp}"is HSE"
            awk -v val=${tmp} 'BEGIN {FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4"}{if($9==val){printf "%-6s%5s%1s%4s%1s%3s%1s%1s%4s%1s%3s%8s%8s%8s%6s%6s%6s%4s\n",$1,$2,$3,$4,$5,"HSE",$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}else{print }}' ${RESP}/$tmppdb > ${RESP}/$tmppdb2
        endif
        yes | mv ${RESP}/$tmppdb2 ${RESP}/$tmppdb >& /dev/null
    end




# Enginneering residue substitution according to MODRES remark
    if ( $nmod ) then
        foreach mmm ( `seq -s" " 1 $nmod` )
            set vvvv = `awk -v ll=$mmm '{if(NR==ll){print $5"_"$6}}' ${RESP}/modres.txt`
            set renum = `echo $vvvv | cut -d_ -f 1`
            set renam = `echo $vvvv | cut -d_ -f 2`
            awk -v rnu=${renum} -v rna=$renam 'BEGIN {FIELDWIDTHS=" 6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4"}{if($9==rnu){printf "%-6s%5s%1s%4s%1s%3s%1s%1s%4s%1s%3s%8s%8s%8s%6s%6s%6s%4s\n",$1,$2,$3,$4,$5,rna,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}else{print }}' ${RESP}/$tmppdb > ${RESP}/$tmppdb2
            yes | mv ${RESP}/$tmppdb2 ${RESP}/$tmppdb >& /dev/null
        end
    endif



    # Copy modified PDB file to a new filename
    yes | cp ${RESP}/$tmppdb ${RESP}/$mol:r_n${num}.pdb >& /dev/null



    # Check if specific conditions are met and update the process status
    if ( `grep -i "SSBOND" ${RESP}/$mol:r_n${num}.pdb | wc -l` ) then
        #echo "SSBOND exists"
    else
        #echo "No SSBOND exists"
    endif
end



# Create ssbond.str file containing specific patches and instructions
set nss = `grep -i "SSBOND" ${RESP}/$mol:r_n1.pdb | wc -l`


cat <<EOF > ${RESP}/ssbond.str
EOF

if ( $nss ) then
    foreach lnum ( `seq -s" " 1 ${nss}` )
        set fnum = `grep -i "SSBOND" ${RESP}/$mol:r_n1.pdb | head -n ${lnum} | tail -n 1 | awk '{print $5}' `
        set snum = `grep -i "SSBOND" ${RESP}/$mol:r_n1.pdb | head -n ${lnum} | tail -n 1 | awk '{print $8}' `
        echo "patch disu PRO0 ${fnum} PRO0 ${snum} setup warn" >> ${RESP}/ssbond.str
        echo "autogenerate angels dihedrals" >> ${RESP}/ssbond.str
    end
    echo "return" >> ${RESP}/ssbond.str
endif



# Perform further modifications on the PDB files
foreach num (`seq -s" " 1 ${tmol}`)
sed -i "s/P01A/PRO0/g" ${RESP}/$mol:r_n${num}.pdb
end



# Check if the modified PDB files exist and update process status
if ( -s ${RESP}/$mol:r_n1.pdb ) then
    echo "${prog}:passed 1" > ${RESP}/$pproc
else
    echo "${prog}:failed 0" > ${RESP}/$pproc
endif



# Check if the modified PDB files exist and update process status
echo "step01_npdb2cpdb.sh $mol  finished"
if ( $dele ) rm -rf ${RESP}/$tmppdb ${RESP}/$tmppdb2 ${RESP}/$tmptxt >& /dev/null
