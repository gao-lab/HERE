#!/bin/sh

awk '{if ($5 !~ /^[.,*]+$/) print}'
