#!/bin/bash

echo This is a list of people who have contributed to the development of BAMM
echo and/or part of BAMMtools. It is sorted by the number of commits made to
echo to the master branch of BAMM\'s Git repository. This file was generated
echo using the script tools/generate-authors.sh.
echo

git log --format='%an' |
    python tools/fix-authors.py | sort | uniq -c | sort -nr |
    awk '{ for (i = 2; i <= NF; i++) printf("%s ", $i); printf("\n") }'
