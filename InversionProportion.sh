#!/bin/bash

val=( '0.02' '0.05' '0.1' '0.25' '0.5' '0.75' )
AnalysisName="InversionProportionSensitivityAnalysis"

for I in ${val[@]}
do
    for i in {1..5}
    do 
        slim -d Rep=${i} -d I=${I} -d AnalysisName="\"$AnalysisName\""  Inversion_model.slim &
    done
done