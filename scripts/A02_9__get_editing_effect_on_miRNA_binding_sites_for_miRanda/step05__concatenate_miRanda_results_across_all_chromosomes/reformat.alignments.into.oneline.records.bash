#!/bin/bash

## pick informative blocks
grep -B 7 -P "^>[^>]+" |
    ## pick informative lines from each block
    grep -P "(Query:|^[ :|]+$|Ref:|^>)" |
    ## remove decorative precedings
    sed -E -e "s@(^[ ]{3}Query:[ ]{4}3'[ ]{1}|[ ]{16}|[ ]{3}Ref:[ ]{6}5'[ ]{1}|>)@@" |
    ## remove decorative trailings
    sed -E -e "s@( [35]+'$)@@" |
    ## print each four line into one line, separated by comma
    paste -d, - - - - |
    ## replace tabs with commas
    tr '\t' ','
