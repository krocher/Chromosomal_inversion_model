#!/bin/bash

AnalysisName = 'NoAnalysis'

for i in {1..5}
do 
    slim -d Rep=${i} -d AnalysisName = AnalysisName  Inversion_model.slim
done
