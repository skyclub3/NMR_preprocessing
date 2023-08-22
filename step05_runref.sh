#!/bin/bash



CHARMM=$PWD/source/charmm_modeller/bin/charmm
workpath=$PWD
softpath=${workpath}/source
mol=$1
model=`echo $mol | cut -d. -f 1`
modell=`echo $model  | awk '{print tolower($1)}'`
modelu=`echo $modell | awk '{print toupper($1)}'`
READ=${workpath}/pReads
RESP=${workpath}/pResults/${modelu}
NUM=`awk 'BEGIN{n=0}{if($1=="MODEL"){n=n+1}else{}}END{if(n==0){print 1}else{print n}}' ${READ}/$mol`
name=`basename $mol`
fileName="${name%.*}"


for num in $(seq 1 $NUM)
do

input=${fileName}_n${num}_refine.inp
noeinput=${fileName}_n${num}_noe_calc.inp
noenati=${fileName}_n${num}_noe_nati.inp
diheinput=${fileName}_n${num}_diheviol.inp
diheinputrefine=${fileName}_n${num}_diheviol_refine.inp

tmpsub=sub.tmp.$$


seed=51494
TYPE=noe
noeampl=1
kdihe=100
distancewidth=0.0
backbonewidth=0.0
sidechainwidth=0.0
analysis=1 # 0: analysis only ; 1: refine str and analysis 
#ensfirst=`ls native/*.pdb | cut -d/ -f2 | cut -d. -f1 | cut -c2- | sort -n | head -n1`


echo "$fileName n: $num"

if [ -s $RESP/ssbond.str ] ; then 
ssbond=1
else
ssbond=0
fi


echo $ssbond



#if ( $# == 8 ) then
#else
#    echo $prog" [number] [ noe | dist | dihe | nati ] [distancewidth:A] [qsubname:B] [analysis] [ssbond]"
#    echo ", where A: any-name in noe, dihe and nati ; B: any-name ; [analysis]: 0:only_analy ; 1:refine_analy ; [ssbond] 1: run ssbond ; 0: run no ssbond"
#    exit
#endif


#==========================================================================
cat <<EOF > $RESP/$input
* Refine
*

set timestep = 0.001
set hnstep   = 3200
set anstep   = 4000
set cnstep   = 8000

! Read topology and parameter files
open read card unit 10  name "${softpath}/top_all36_prot.rtf"
read  rtf card unit 10

open read card unit 20 name "${softpath}/par_all36m_prot.prm"
read para card unit 20

! Read PROA
open read card unit 10 name "${RESP}/${fileName}_n${NUM}.pdb"
read sequence pdb unit 10
generate pro0 setup warn 

open read card unit 10 name "${RESP}/${fileName}_n${NUM}.pdb"
read coor pdb  unit 10 resid

set ssbond = $ssbond
if @ssbond .eq. 1 then
    stream "${RESP}/ssbond.str"
    set ssstr = ss
else
    set ssstr = 
endif

ic para
ic build
hbuild


EOF
#==========================================================================



if [ $type=="noe" -o $type=="NOE" ]; then 
cat<<EOF >> $RESP/$input
set kdihe = $kdihe
set scale = $noeampl
NOE
  RESET
END
set fmax  = 0.01
set kmax = 20
prnl 0
EOF
fi




#==========================================================================
if [ $type=="noe" -a -s $RESP/rawdata_dihe ]; then
cat <<EOF >> $RESP/$input
stream "${RESP}/rsr/dihrsr"
stream "${RESP}/rsr/noersr_intra_ambi"
stream "${RESP}/rsr/noersr_seq_ambi" 
stream "${RESP}/rsr/noersr_med_ambi"
stream "${RESP}/rsr/noersr_long_ambi"
stream "${RESP}/rsr/${modelu}_n1_phipsi.str"
stream "${RESP}/rsr/${modelu}_n1_phichi1.str" 
stream "${RESP}/rsr/${modelu}_n1_psichi1.str"
stream "${RESP}/rsr/${modelu}_n1_chi1chi2.str" 
prnl 5
NOE
  SCALE @scale
END
EOF

elif [ $type=="noe" -a ! -s $RESP/rawdata_dihe ]; then
cat <<EOF >> $RESP/$input
stream "${RESP}/rsr/noersr_intra_ambi" 
stream "${RESP}/rsr/noersr_seq_ambi" 
stream "${RESP}/rsr/noersr_med_ambi"
stream "${RESP}/rsr/noersr_long_ambi" 
stream "${RESP}/rsr/${modelu}_n1_phipsi.str" 
stream "${RESP}/rsr/${modelu}_n1_phichi1.str" 
stream "${RESP}/rsr/${modelu}_n1_psichi1.str" 
stream "${RESP}/rsr/${modelu}_n1_chi1chi2.str" 
prnl 5
NOE
  SCALE @scale
END
EOF
fi



#==========================================================================

cat<<EOF >> $RESP/$input
!solvation (GBSW)
!
! GBSW setup
!

prnlev 0
stream "${softpath}/radii_prot_na.str"
stream "${softpath}/radius_gbsw.str"
prnlev 5

!!double check if some heavy atoms have a zero radius
scalar wmain statistics select .not. type H* end
define check select ( .not. type H* ) .and. ( property wmain .eq. 0.0 ) show end
if ?nsel ne 0 stop

gbsw sgamma 0.005 nang 50

nbond atom switch cdie vdw vswitch ctonnb 16.0 ctofnb 16.0 cutnb 20.0

ener

set enerbeforemini = ?ener

energy

!short minimization
mini sd   nstep 100
mini abnr tolgrd 0.0001 nstep 500 nprint 1000

shake bonh para tol 1e-8

!!heating
! 100 -> 1000 K (hnstep = 1600 steps)
dynamics leap verlet start timestep @timestep nstep @hnstep -
     nprint 1000 iprfrq 1000 echeck 2e99 -
     iunread -1 iunwri -1 iuncrd -1 kunit -1 iseed @seed -
     nsavc 0 nsavv 0 iunvel -1 isvfrq 2000 -
     firstt 100 finalt 1000 twindl -10.0 twindh 10.0 -
     ichecw 1 teminc 5.0 ihtfrq 20

!!annealing
!! 500 K (anstep = 2000 steps)
dynamics leap verlet start timestep @timestep nstep @anstep -
     nprint 1000 iprfrq 1000 echeck 1e99 -
     iunread -1 iunwri -1 iuncrd -1 kunit -1 iseed @seed -
     nsavc 0 nsavv 0 iunvel -1 isvfrq 1000 -
     firstt 1000 finalt 1000 twindh 25 twindl -25 ieqfrq 200

!!cooling
!! 500 -> 25 K (cnstep = 4000 steps)
dynamics leap verlet start timestep @timestep nstep @cnstep -
     nprint 1000 iprfrq 1000 echeck 2e99 -
     iunread -1 iunwri -1 iuncrd -1 kunit -1 iseed @seed -
     nsavc 0 nsavv 0 iunvel -1 isvfrq 2000 -
     firstt 1000 finalt 25 twindl -10.0 twindh 10.0 -
     ichecw 1 teminc -5.0 ihtfrq 40

!short minimization
mini sd   nstep 100
mini abnr tole 0.001 tolgrd 0.001 nstep 3000 nprint 1000


energy

EOF

#==========================================================================


if [ $type=="noe" -o $type=="NOE" ] ; then
cat<<EOF >> $RESP/$input
open writ card unit 10 name "${RESP}/${modelu}_n${num}_refine_noe.pdb"
writ coor pdb  unit 10
* ?ener
*
stop
EOF
fi


#==========================================================================
#ANALYSIS




cat<<EOF > $RESP/$noeinput
* noe_refinement
*

! Read topology and parameter files
open read card unit 10  name "${softpath}/top_all36_prot.rtf"
read  rtf card unit 10

open read card unit 20 name "${softpath}/par_all36m_prot.prm"
read para card unit 20

if @ssbond .eq. 1 then
    set ssstr = ss
else
    set ssstr = 
endif


! Read PROA
open read card unit 10 name "${RESP}/${modelu}_n${num}_refine_noe.pdb"
read sequence pdb unit 10
generate pro0 setup warn

open read card unit 10 name "${RESP}/${modelu}_n${num}_refine_noe.pdb"
read coor pdb  unit 10 resid

if @ssbond .eq. 1 then
    stream "${RESP}/ssbond.str"
    ic para
    ic build
    hbuild
else
    ic para
    ic build
    hbuild
endif

open write card unit 13 name "${RESP}/${modelu}_n${num}_refine_noe_viol_rmsd.dat"

!stream ana_noe_bin.str 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!all
set noetype = all
NOE
  RESET
END
set kdihe = $kdihe
set scale = $noeampl
set fmax  = 0.01
set kmax = 20
prnlev 0
stream "${RESP}/rsr/noersr_intra_ambi"
stream "${RESP}/rsr/noersr_seq_ambi"
stream "${RESP}/rsr/noersr_med_ambi"
stream "${RESP}/rsr/noersr_long_ambi"

prnlev 6
NOE
  SCALE @scale
END

!stream ana_noe_bin.str 

set dcut    = 0.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
*  TYPE   dcut  Nviol  Maxnoe  RMSD       Num of Vio. Perc   NOE ener.
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU @violper  ?NOE
*
endif

set dcut    = 0.01
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU @violper  ?NOE
*
endif

set dcut    = 0.5
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 1.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 2.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!intra
set noetype = intra
NOE
  RESET
END
set scale = $noeampl
set fmax  = 0.01
set kmax = 20
prnlev 0
stream "${RESP}/rsr/noersr_intra_ambi"
prnlev 5
NOE
  SCALE @scale
END

!stream ana_noe_bin.str 

set dcut    = 0.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum 
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 0.01
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU @violper  ?NOE
*
endif

set dcut    = 0.5
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 1.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 2.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!seq
set noetype = seq
NOE
  RESET
END
set scale = $noeampl
set fmax  = 0.01
prnlev 0
stream "${RESP}/rsr/noersr_seq_ambi"
prnlev 5
NOE
  SCALE @scale
END

!stream ana_noe_bin.str 

set dcut    = 0.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 0.01
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU @violper  ?NOE
*
endif

set dcut    = 0.5
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 1.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 2.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!med
set noetype = med
NOE
  RESET
END
set scale = $noeampl
set fmax  = 0.01
set kmax = 20
prnlev 0
stream "${RESP}/rsr/noersr_med_ambi"
prnlev 5
NOE
  SCALE @scale
END

!stream ana_noe_bin.str 

set dcut    = 0.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 0.01
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU @violper  ?NOE
*
endif

set dcut    = 0.5
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 1.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 2.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!long
set noetype = long
NOE
  RESET
END
set scale = $noeampl
set fmax  = 0.01
set kmax = 20
prnlev 0
stream "${RESP}/rsr/noersr_long_ambi"
prnlev 5
NOE
  SCALE @scale
END

!stream ana_noe_bin.str 

set dcut    = 0.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 0.01
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU @violper  ?NOE
*
endif

set dcut    = 0.5
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 1.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 2.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

stop
EOF




cat<<EOF > $RESP/$noenati
* noe_native
*

! Read topology and parameter files
open read card unit 10  name "${softpath}/top_all36_prot.rtf"
read  rtf card unit 10

open read card unit 20 name "${softpath}/par_all36m_prot.prm"
read para card unit 20

if @ssbond .eq. 1 then
    set ssstr = ss
else
    set ssstr =
endif


! Read PROA
open read card unit 10 name "${RESP}/${fileName}_n${NUM}.pdb"
read sequence pdb unit 10
generate pro0 setup warn

open read card unit 10 name "${RESP}/${fileName}_n${NUM}.pdb"
read coor pdb  unit 10 resid

if @ssbond .eq. 1 then
    stream "${RESP}/ssbond.str"
    ic para
    ic build
    hbuild
else
    ic para
    ic build
    hbuild
endif

open write card unit 13 name "${RESP}/${modelu}_n${num}_nati_noe_viol_rmsd.dat"

!stream ana_noe_bin.str

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!all
set noetype = all
NOE
  RESET
END
set kdihe = $kdihe
set scale = $noeampl
set fmax  = 0.01
set kmax = 20
prnlev 0
stream "${RESP}/rsr/noersr_intra_ambi"
stream "${RESP}/rsr/noersr_seq_ambi"
stream "${RESP}/rsr/noersr_med_ambi"
stream "${RESP}/rsr/noersr_long_ambi"

prnlev 6
NOE
  SCALE @scale
END

!stream ana_noe_bin.str

set dcut    = 0.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
*  TYPE   dcut  Nviol  Maxnoe  RMSD       Num of Vio. Perc   NOE ener.
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU @violper  ?NOE
*
endif

set dcut    = 0.01
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU @violper  ?NOE
*
endif

set dcut    = 0.5
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 1.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 2.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!intra
set noetype = intra
NOE
  RESET
END
set scale = $noeampl
set fmax  = 0.01
set kmax = 20
prnlev 0
stream "${RESP}/rsr/noersr_intra_ambi"
prnlev 5
NOE
  SCALE @scale
END

!stream ana_noe_bin.str

set dcut    = 0.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 0.01
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU @violper  ?NOE
*
endif

set dcut    = 0.5
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 1.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 2.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!seq
set noetype = seq
NOE
  RESET
END
set scale = $noeampl
set fmax  = 0.01
prnlev 0
stream "${RESP}/rsr/noersr_seq_ambi"
prnlev 5
NOE
  SCALE @scale
END

!stream ana_noe_bin.str

set dcut    = 0.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 0.01
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU @violper  ?NOE
*
endif

set dcut    = 0.5
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 1.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 2.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!med
set noetype = med
NOE
  RESET
END
set scale = $noeampl
set fmax  = 0.01
set kmax = 20
prnlev 0
stream "${RESP}/rsr/noersr_med_ambi"
prnlev 5
NOE
  SCALE @scale
END

!stream ana_noe_bin.str

set dcut    = 0.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 0.01
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU @violper  ?NOE
*
endif

set dcut    = 0.5
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 1.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 2.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!long
set noetype = long
NOE
  RESET
END
set scale = $noeampl
set fmax  = 0.01
set kmax = 20
prnlev 0
stream "${RESP}/rsr/noersr_long_ambi"
prnlev 5
NOE
  SCALE @scale
END

!stream ana_noe_bin.str

set dcut    = 0.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 0.01
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU @violper  ?NOE
*
endif

set dcut    = 0.5
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 1.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

set dcut    = 2.0
NOE
  PRINT ANAL CUT @dcut
END

if ?nnum .eq. 0 then
write title unit 13
* @noetype @dcut ?nnum
*
else
!set nviolintra2 = ?VIOLNU
!calc violintra = sqrt ( ?VIOL * ?viol / ?nnum )
! RMSD of NOE
calc violintra   = ?RMSN
set nviolintra   = ?VIOLNU
set nintra       = ?NOENUM
energy
set enernoeintra = ?NOE
calc violper     = ?violnu / ?nnum * 100
write title unit 13
* @noetype @dcut ?nnum ?mxvn  ?RMSN      ?VIOLNU  @violper    ?NOE
*
endif

stop
EOF

cat<<EOF > $RESP/$diheinput
* noe_dihe
*

! Read topology and parameter files
open read card unit 10  name "${softpath}/top_all36_prot.rtf"
read  rtf card unit 10

open read card unit 20 name "${softpath}/par_all36m_prot.prm"
read para card unit 20

if @type .eq. noe then
    set pdbfile = "${RESP}/${fileName}_n${num}.pdb"
    set outfile = "${RESP}/noe_viol_dihe_rmsd_noe_${num}.dat"
endif


! Read PROA
open read card unit 10 name @pdbfile
read sequence pdb unit 10
generate pro0 setup warn

open read card unit 10 name @pdbfile
read coor pdb  unit 10 resid

if @ssbond .eq. 1 then
    stream "${RESP}/ssbond.str"
    ic para
    ic build
    hbuild
else
    ic para
    ic build
    hbuild
endif

open write card unit 13 name @outfile

set kdihe = @kdihe
set scale = 50
set fmax  = 0.01
prnlev 0
stream "${RESP}/rsr/dihrsr" @kdihe
prnlev 6
NOE
  SCALE @scale
END

energy
stop
EOF

cat<<EOF > $RESP/$diheinputrefine
* noe_dihe_refine
*

! Read topology and parameter files
open read card unit 10  name "${softpath}/top_all36_prot.rtf"
read  rtf card unit 10

open read card unit 20 name "${softpath}/par_all36m_prot.prm"
read para card unit 20

if @type .eq. noe then
    set pdbfile = "${RESP}/${fileName}_n${num}_refine_noe.pdb"
    set outfile = "${RESP}/noe_viol_dihe_rmsd_noe_${num}_refine.dat"
endif


! Read PROA
open read card unit 10 name @pdbfile
read sequence pdb unit 10
generate pro0 setup warn

open read card unit 10 name @pdbfile
read coor pdb  unit 10 resid

if @ssbond .eq. 1 then
    stream "${RESP}/ssbond.str"
    ic para
    ic build
    hbuild
else
    ic para
    ic build
    hbuild
endif

open write card unit 13 name @outfile

set kdihe = @kdihe
set scale = 50
set fmax  = 0.01
prnlev 0
stream "${RESP}/rsr/dihrsr" @kdihe
prnlev 6
NOE
  SCALE @scale
END

energy
stop
EOF

done
exit
