#!/bin/bash

val=( '0.0' '0.1' '0.25' '0.5' )
AnalysisName = 'MutationDominanceSensitivityAnalysis'

for H in ${val[@]}
do
    for i in {1..5}
    do 
        slim -d Rep=${i} -d H=${H} -d AnalysisName = "\"$AnalysisName\""  Inversion_model.slim &
    done
done