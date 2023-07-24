#!/bin/bash
$x_paramètre = $1
for i in $paramètre
do
    grep -f $i -A 3 $2
    
done