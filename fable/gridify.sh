#!/bin/sh

COL=1
ROW=1

printf "{\n";

for WORD in `cat $1`
do
	if [ $COL -eq 1 ]
	then
		printf "\t{";
	fi
	printf "%s" "$WORD";
	if [ $COL -eq 8 ]
	then
		printf "},\n";
		COL=1;
		(( ROW += 1 ))
	else
		printf ",\t";
		(( COL += 1 ))
	fi
done

printf "}\n";
