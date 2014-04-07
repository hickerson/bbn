#!/bin/sh

for file in bbn.cpp cmn.hpp
do
	sed -i 's/\/\/C/\/\//g' $file
	sed -i 's/fem::float0/0/g' $file
	sed -i 's/fem::int0/0/g' $file
done
