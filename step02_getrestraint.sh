#!/bin/csh

# Set script parameters
set model = $1
set noeg = $2
set noeh = $3
set dihe = $4
set rdc = $5   # Not used
set prog = `echo $0:r`
set pproc = ${prog}".proc"
set tmpf = /tmp/${prog}$$
set dele = 1

# BMRB data download site
set bmrbhost = https://restraintsgrid.bmrb.io/NRG/
set READ = `pwd`/pReads 
set RESP = `pwd`/pResults

# Check if correct number of arguments is provided
if ( $#argv == 5 ) then
else
    echo $0 "[PDBID] [gennoe] [HBnoe] [Dihe] [RDC], where [gennoe], [HBnoe], [Dihe], and [RDC] are 1 (download) or 0 (no download)."
    exit
endif

# Convert model to lowercase and uppercase
set modell = `echo $model  | awk '{print tolower($1)}'`
set modelu = `echo $modell | awk '{print toupper($1)}'`

# Set working directory and create a subdirectory for the model
set dir = `pwd`
set twowd = `echo $modell | cut -c2-3`
mkdir -p ${RESP}/${modelu}

# Download and process NOE restraints
if ( $noeg || $noeh ) then

    # Download general distance restraints
    wget -T 7 --no-check-certificate "${bmrbhost}MRGridServlet?block_text_type=3-converted-DOCR&db_username=wattos1&file_detail=3-converted-DOCR&format=ambi&pdb_id=${model}&program=XPLOR%2FCNS&request_type=archive&subtype=general+distance&type=distance" -O ${RESP}/${modelu}/noeg.zip >& /dev/null

    # Download hydrogen bond type restraints
    wget -T 7 --no-check-certificate "${bmrbhost}MRGridServlet?block_text_type=3-converted-DOCR&db_username=wattos1&file_detail=3-converted-DOCR&format=ambi&pdb_id=${model}&request_type=archive&subtype=hydrogen+bond&type=distance" -O ${RESP}/${modelu}/noeh.zip >& /dev/null

    # Unzip downloaded files and concatenate data
    foreach tmp (`ls ${RESP}/${modelu}/*.zip*`)
        yes | unzip ${tmp} -d ${RESP}/${modelu} >& /dev/null
    end
    cat ${RESP}/${modelu}/block_*.tbl > ${RESP}/${modelu}/rawdata

    # Copy data for AQUA
    yes | cp -rf ${RESP}/${modelu}/rawdata ${RESP}/${modelu}/rawdata_foraqua >& /dev/null

    # Clean up extracted files
    rm -rf ${RESP}/${modelu}/*.tbl ${RESP}/${modelu}/*.zip* >& /dev/null

endif

# Download and process dihedral angle restraints
if ( $dihe ) then

    # Download dihedral angles restraints
    wget -T 7 --no-check-certificate "${bmrbhost}MRGridServlet?block_text_type=3-converted-DOCR&db_username=wattos1&file_detail=3-converted-DOCR&format=n%2Fa&pdb_id=${model}&request_type=archive&subtype=n%2Fa&type=dihedral+angle" -O ${RESP}/${modelu}/dihe.zip >& /dev/null

    # Unzip downloaded files and concatenate data
    foreach tmp (`ls ${RESP}/${modelu}/*.zip*`)
        yes | unzip ${tmp} -d ${RESP}/${modelu} >& /dev/null
    end
    cat ${RESP}/${modelu}/block_*.tbl > ${RESP}/${modelu}/rawdata_dihe

    # Clean up extracted files
    rm -rf ${RESP}/${modelu}/*.tbl ${RESP}/${modelu}/*.zip* >& /dev/null

endif

# Check if rawdata files exist and update process status
cd ${RESP}/${modelu}
if ( ( -e rawdata || -e rawdata_foraqua || -e rawdata_dihe ) && ( ! -z rawdata || ! -z rawdata_foraqua || ! -z rawdata_dihe ) ) then
    echo "${prog}:passed 1" > $pproc
else
    echo "${prog}:failed 0" > $pproc
endif

