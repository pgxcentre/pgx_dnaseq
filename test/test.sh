#!/usr/bin/env bash

execute_pipeline.py \
    --input input_files.txt \
    --tool-config tool.conf \
    -n 2 \
    --pipeline-config pipeline.conf

rc=$?
if [[ $rc == 0 ]]
then
    automatic_report.py \
        --input input_files.txt \
        --pipeline-config pipeline.conf \
        --run-name "Dummy_Run"
fi

