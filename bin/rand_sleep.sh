#!/bin/bash
r=$(( $RANDOM % 1000 ))
ms=$(echo $r / 1000 | bc -l)
ms_rounded=$(printf "%.2f" $ms)
echo "Sleeping for $ms_rounded seconds!"
sleep $ms_rounded
