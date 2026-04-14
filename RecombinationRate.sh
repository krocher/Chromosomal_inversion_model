#!/bin/bash

val=('1e-8' '1e-7' '1e-6' '1e-5')
AnalysisName="RecombinationRateSensitivityAnalysis"

for R in ${val[@]}
do
    for i in {1..2}
    do 
        slim -d Rep=${i} -d R=${R} -d AnalysisName="\"$AnalysisName\""  Inversion_model.slim &
    done
done