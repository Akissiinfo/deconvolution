#!/bin/bash

for i in $1 
do
    awk -F',' '{print $1}' $i;
done