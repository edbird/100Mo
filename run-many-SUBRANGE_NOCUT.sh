#!/bin/bash
EPSILON_LIST=$(awk 'BEGIN{for(i=0.0;i<=2.0;i+=0.01)print i}')
for EPSILON in $EPSILON_LIST
do
    echo "Running: epsilon=$EPSILON"
    ./build/main --batch-mode true --epsilon $EPSILON --fit-subrange true --energy-cut false --output-file "of_data_SUBRANGE_NOCUT.txt" > "log_SUBRANGE_NOCUT_EPS_"$EPSILON".txt"
done
