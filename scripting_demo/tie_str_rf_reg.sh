#!/bin/bash

# AutoML random forest runs for tie strength score prediction
python run_automl.py \
       ../data/final_features/all_tie_str_baseline \   # input features
       final_results/tie_str/tie_str_baseline_rf_reg \ # output name
       tie_str_score \                                 # outcome variable (regression)
       --run_time 1440 --task_time 21600 \             # training time
       --rand_forest \                                 # only train random forest estimators
       > tie_str_baseline_rf_reg.out;


