#!/bin/bash
EPSILON_LIST=$(awk 'BEGIN{for(i=0.0;i<=2.0;i+=0.2)print i}')
for EPSILON in $EPSILON_LIST
do
    echo "Running: epsilon=$EPSILON"
    ./build/main --batch-mode false --epsilon $EPSILON --fit-subrange false --energy-cut false --output-file "of_data_testing.txt" > "log_testing_EPS_"$EPSILON".txt"
done
