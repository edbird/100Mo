#!/bin/bash
EPSILON_LIST=$(awk 'BEGIN{for(i=0.0;i<=0.08;i+=0.001)print i}')
for EPSILON in $EPSILON_LIST
do
    ./build/main --batch-mode true --epsilon $EPSILON
done
