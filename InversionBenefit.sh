#!/bin/bash

val=( '0.001' '0.005' '0.01' '0.05' '0.1')
AnalysisName="InversionBenefitSensitivityAnalysis"

for SINV in ${val[@]}
do
    for i in {1..5}
    do 
        slim -d Rep=${i} -d S_INV=${SINV} -d AnalysisName="\"$AnalysisName\""  Inversion_model.slim &
    done
done