#!/bin/bash

for file in $*
do
	echo "altering file $file"
	sed -i 's/\/\/C/\/\//g' $file
	for word in float0 int0 char0 short0 double0
	do
		sed -i "s/fem::$word/0/g" $file
	done
	for key in cosh sinh asin acos asinh acosh log10 pow cos sin exp log sqrt
	do
		sed -i "s/fem::$key/$key/g" $file
	done
done
