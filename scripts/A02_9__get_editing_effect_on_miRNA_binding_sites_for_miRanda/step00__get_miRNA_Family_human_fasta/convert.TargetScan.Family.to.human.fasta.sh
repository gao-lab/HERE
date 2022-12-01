#!/bin/bash


cat $1 | 
    awk '
OFS="\t"{
	if($6==2 && $3==9606) {
		 print ">"$4"__"$1
		 print $5
	}
}'
