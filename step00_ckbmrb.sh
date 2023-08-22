#!/bin/csh

# This script fetches various types of NMR-related data from the BMRB database for a given model (PDBID).
# It checks the availability of general distance restraints, hydrogen bond data, dihedral angles, and dipolar coupling information.

# Check if the correct number of arguments is provided
if ( $#argv != 1 ) then
    echo "Usage: $0 [PDBID]"
    exit
endif

set model = $1
set prog  = $0:r
set tmpf  = ${prog}.tmp.$$
set dele  = 1

# BMRB data download site URL
set bmrbhost = https://restraintsgrid.bmrb.io/NRG/

# Check availability of general distance restraints data
wget --spider --no-check-certificate "${bmrbhost}MRGridServlet?block_text_type=3-converted-DOCR&db_username=wattos1&file_detail=3-converted-DOCR&format=ambi&pdb_id=${model}&program=XPLOR%2FCNS&request_type=archive&subtype=general+distance&type=distance" >& $tmpf
set val = `cat $tmpf | grep -i "Length:" | awk '{str=sprintf("%s",substr($4,2,length($4)-2));if(str=="application/zip"){print 1}else{print 0}}'`
if ( $val > 0 ) then
    set gdist = 1
else
    set gdist = 0
endif

# Check availability of hydrogen bond data
wget --spider --no-check-certificate "${bmrbhost}MRGridServlet?block_text_type=3-converted-DOCR&db_username=wattos1&file_detail=3-converted-DOCR&format=ambi&pdb_id=${model}&request_type=archive&subtype=hydrogen+bond&type=distance" >& $tmpf
set val = `cat $tmpf | grep -i "Length:" | awk '{str=sprintf("%s",substr($4,2,length($4)-2));if(str=="application/zip"){print 1}else{print 0}}'`
if ( $val > 0 ) then
    set hdist = 1
else
    set hdist = 0
endif

# Check availability of dihedral angles data
wget --spider --no-check-certificate "${bmrbhost}MRGridServlet?block_text_type=3-converted-DOCR&db_username=wattos1&file_detail=3-converted-DOCR&format=n%2Fa&pdb_id=${model}&request_type=archive&subtype=n%2Fa&type=dihedral+angle" >& $tmpf
set val = `cat $tmpf | grep -i "Length:" | awk '{str=sprintf("%s",substr($4,2,length($4)-2));if(str=="application/zip"){print 1}else{print 0}}'`
if ( $val > 0 ) then
    set dihe = 1
else
    set dihe = 0
endif

# Check availability of dipolar coupling (RDC) data
wget --spider --no-check-certificate "${bmrbhost}MRGridServlet?block_text_type=3-converted-DOCR&db_username=wattos1&file_detail=3-converted-DOCR&format=n%2Fa&pdb_id=${model}&request_type=archive&subtype=n%2Fa&type=dipolar+coupling" >& $tmpf
set val = `cat $tmpf | grep -i "Length:" | awk '{str=sprintf("%s",substr($4,2,length($4)-2));if(str=="application/zip"){print 1}else{print 0}}'`
if ( $val > 0 ) then
    set rdc = 1
else
    set rdc = 0
endif

# Clean up the temporary file
if ( $dele ) rm -rf $tmpf >& /dev/null

# Append model identifier and data availability flags to the output file
echo "$model $gdist $hdist $dihe $rdc" >> NMR_ID.sort.ckbmrb.txt

