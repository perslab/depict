#! /bin/bash 
 
source /broad/software/scripts/useuse 

reuse -q use Python-2.7
reuse -q use Java-1.8

"$@"
