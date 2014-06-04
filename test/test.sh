#!/usr/bin/env bash

PYTHONPATH=..:$PYTHONPATH \
../scripts/execute_pipeline.py \
    --input input_files.txt \
    --tool-config tool.conf \
    --pipeline-config pipeline.conf
