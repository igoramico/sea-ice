#!/bin/bash

echo -ne "MYFLAGS =" >> myflags.tmp
for var in "$@"
 do
  echo -ne -D$var >> myflags.tmp
  echo -ne " " >> myflags.tmp
 done
 echo >> myflags.tmp

cat myflags.tmp Makefile-template > Makefile

rm -f myflags.tmp 



