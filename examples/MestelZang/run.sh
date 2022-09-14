#!/bin/bash
#PBS -S /bin/sh
#PBS -N MestelUnstable_1.0
#PBS -o ./log/MestelUnstable_1.0.log
#PBS -e ./log/MestelUnstable_1.0.err
#PBS -l nodes=1:ppn=128,mem=1000g,walltime=10:00:00

module purge
module load julia/1.7.2

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=128
export JULIA_CPU_THREADS=128

JULIA=julia
PROJECT=/home/roule/Disc_2D/JuliaCallAResponse
CODE=./runbash.jl

PREFIX=/home/roule/Disc_2D/JuliaCallAResponse/examples/MestelZang
cd ${PREFIX}

# Basis parameters
basisname="CluttonBrock"
G=1.0
rb=5.0
kKA=7
lharmonic=2
nmax=100

# Potential parameters
modelname="Mestel"
R0=20.0
V0=1.0

# DF parameters
dfname="Zang"
q0=6
Rin=1.0
Rout=11.5
Rmax=20.0
xi=1.0
mu=5
nu=6

# Integration parameters
rmin=0.1
rmax=20.0
Ku=100
Kv=100
kw=200
n1max=50

# Frequency computations
omgmin=0.5
omgmax=1.5
etamin=-0.1
etamax=0.5
nomg=100
neta=100

# Other parameters
verbose=1
overwrite=false


${JULIA} --project=${PROJECT} --threads $JULIA_NUM_THREADS ${CODE} # --G $G --basisname $basisname --rb $rb --kKA $kKA --lharmonic $lharmonic --nmax $nmax --modelname $modelname --R0 $R0 --V0 $V0 --DFname $dfname --q0 $q0 --Rin $Rin --Rout $Rout --Rmax $Rmax --xi $xi --mu $mu --nu $nu --rmin $rmin --rmax $rmax --Ku $Ku --Kv $Kv --Kw $Kw --n1max $n1max --omgmin $omgmin --omgmax $omgmax --etamin $etamin --etamax $etamax --nomg $nomg --neta $neta --verbose $verbose --overwrite $overwrite
