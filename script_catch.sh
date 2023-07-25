#!/bin/bash

for i in $1
do
    grep -f $i -A 3 $2
    
done