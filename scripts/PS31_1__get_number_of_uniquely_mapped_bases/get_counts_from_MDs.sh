#!/bin/bash

sed -E -e 's/[A-Z^]+/\n/g' $1 | awk '{sum += $0} END {print sum}'
