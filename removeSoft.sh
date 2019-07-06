#!/bin/bash

if [ -e $1 ] 
 then 
  awk '{if($19!=3) print $0}' $1 > $1.h
 fi
