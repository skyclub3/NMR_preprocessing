#!/bin/csh

# This line specifies that the script should be executed using the C Shell interpreter.

set mol = $1
set model = `echo $mol | cut -d. -f 1`
set modell = `echo $model | awk '{print tolower($1)}'`
set modelu = `echo $modell | awk '{print toupper($1)}'`
# The above lines set various variables based on the input argument ($1) provided when running the script.
# $1 is expected to be a file name, and the script extracts different versions of the file name.

set workpath = `pwd` 
set softpath = ${workpath}/source
set READ = ${workpath}/pReads
set RESP = ${workpath}/pResults/${modelu}
# These lines set the paths to various directories used in the script.

echo $model $modell $modelu
# Print out the values of the extracted variables for debugging purposes.

######################################################################
#2. NOE check
######################################################################

# The following lines perform checks related to nuclear Overhauser effect (NOE) data.

#general
set noenum = `awk 'END{na=split($NF,a,":");print a[2]}' ${RESP}/rsr/noersr_ambi`

if ( $noenum == 0 ) then
    set noerun = 0
else
    set noerun = 1
endif
# The above lines extract and process some value using awk and check whether noenum is zero.
# Based on the condition, it sets the value of noerun to either 0 or 1.

if ( $noerun == 0 ) then
    set noestate = "noNOE"
else
    #aqua
    mkdir -p ${RESP}/aqua >& /dev/null 
    yes |cp ${softpath}/aqua/{noe.emacs,aqua.sh,hist.sh,run_emacs.bash}  ${RESP}/aqua  >& /dev/null
    yes |cp ${RESP}/rawdata_foraqua ${RESP}/aqua/rawdata >& /dev/null

    ${RESP}/aqua/aqua.sh ${modelu}.pdb > ${RESP}/aqua/aqua.out
    set noestate = `${RESP}/aqua/hist.sh ${modelu}.pdb`

    cd ../
endif
# The above lines execute certain operations related to the "aqua" process if noerun is not zero.
# It involves copying and moving files, running the "aqua.sh" script, and executing "hist.sh".

echo $modelu $noerun $noestate > ${RESP}/noeinfo.dat

echo "step04 aqua check $modelu finished"
# Print out information about the completion of the aqua check process.

