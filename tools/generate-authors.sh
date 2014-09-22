#!/bin/bash

echo This is a list of people who have contributed to the development of BAMM
echo and/or part of BAMMtools. It is sorted by the number of commits made to
echo to the master branch of BAMM\'s Git repository. This file was generated
echo using the script tools/generate-authors.sh.
echo

git log --format='%an' |
    awk '{ if ($0 == "blueraleigh")
               print "Mike Grundler"
           else if ($0 == "Dan Rabosky")
               print "Daniel Rabosky"
           else if ($0 == "pascaltitle")
               print "Pascal Title"
           else if ($0 == "josephwb")
               print "Joseph W. Brown"
           else if ($0 == "SimonGreenhill")
               print "Simon Greenhill"
           else
               print
         }' | sort | uniq -c | sort -nr |
    awk '{ $1 = ""; print }' | cut -d' ' -f2-
