#!/bin/bash


awk '{
    if ($2 == 9606) {
       print ">"$1"__3UTR"
       print $3
    }
}' $1
