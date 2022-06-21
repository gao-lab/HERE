#!/bin/bash

grep -P "\tNH:i:1\t" $1 | sed -E -e "s/.*\tMD:Z:([^\t]+)\t.*/\1/"
