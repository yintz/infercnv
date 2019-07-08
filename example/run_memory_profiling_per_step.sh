#!/bin/bash

for i in `seq 1 21`; do
    gtime -v Rscript run_test.R $i > profiling/up_to_step_${i}_1.log 2> profiling/up_to_step_${i}_1.times
    gtime -v Rscript run_test.R $i > profiling/up_to_step_${i}_2.log 2> profiling/up_to_step_${i}_2.times
done
