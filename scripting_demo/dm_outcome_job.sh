#!/bin/bash
# Executes BSUB job for pulling diabetes outcomes

# schedule processes for data in each quarter of every year
for yr in `seq 2001 2016`;
do
    for q in `seq 1 4`;
    do
        { time python data_extract/med_extract.py /origdata/Optum/2017update/SES/ $yr $q m_data/outcome_data/ m 100000 --m_dm_outcome ; } &
    done
done
