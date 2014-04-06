#!/bin/sh

sed -i 's/\/\/C/\/\//g' bbn.cpp
sed -i 's/fem::float0/0/g' bbn.cpp
sed -i 's/fem::int0/0/g' bbn.cpp
