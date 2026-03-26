#!/bin/bash

val=( '1e-6' '5e-6' '1e-5' '5e-5' )
AnalysisName="MutationRateSensitivityAnalysis"

for Mu in ${val[@]}
do
    for i in {1..5}
    do 
        slim -d Rep=${i} -d Mu=${Mu} -d AnalysisName="\"$AnalysisName\""  Inversion_model.slim &
    done
done