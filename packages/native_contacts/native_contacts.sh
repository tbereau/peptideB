#!/bin/bash
#
# Author: David Stone
# Date  : Feb, 2009

# Read a PDB file and compare the structure to a number of reference
# contacts. These contacts need to be explicitly written in a 'contact file'. The file will be of the form:
#
# 'numbonds part_ID_pair1_1 part_ID_pair1_2 part_ID_pair2_1 part_ID_pair2_2 ... dist1 dist2 ... '
# where numbonds is the number of contacts, part_ID_pairX_k is the ID of particle k of the pair X.
# dist1 is the distance between part_ID_pair1_1 and part_ID_pair1_2.
#
#
#

i=1
k=7
m=0
q=0
x0=0
x1=0
x2=0
y0=0
y1=0
y2=0
q0=0 #for incrementing q by exponential
dist=0
scale=0.2

# Command line arguments
if [ -z  "$2" ]; then
	echo "Usage: $0 pdb_file contacts_file"
	echo ""
	echo "contacts_file should be of the form:"
	echo "'numbonds part_ID_pair1_1 part_ID_pair1_2 part_ID_pair2_1 part_ID_pair2_2 ... dist1 dist2 ... '"
	echo "where numbonds is the number of contacts, part_ID_pairX_k is the ID of particle k of the pair X."
	echo "dist1 is the distance between part_ID_pair1_1 and part_ID_pair1_2."
	exit 1
fi


#assume incoming second command line file of form:
#numbonds ID1 ID2 ID3 ID4...dist1 dist2 dist3... 
#this assumes ID bonds are unique (one ID per bond)
bonds=`cat ${2} | sed 's/ * / /g' | cut -d ' ' -f1`;
h=$(echo "2*${bonds}"|bc)

while [ $m -lt $h ] 
do
	n=$(echo "${m}+2"|bc);
	id[$m]=`cat ${2}| cut -d ' ' -f${n}`
	m=$(($m+1));
done

#distance threshold loop
m=0
while [ $m -lt $bonds ]
do
	n=$(echo "2+${m}+${h}"|bc)
	thresh[$m]=`cat ${2}|cut -d ' ' -f${n}`
	m=$(($m+1))
done

date=`date +"%m-%d"`;
pdb=$1;
#output=`echo "${1}" | cut -d '.' -f1`;

h=$(echo "${h}-1"|bc)
m=0
#had to run the numbers backwards so the thresh vars correctly lined up with m
for m in $(seq 0 2 $h)
do
	k=7; #for selecting distances
	dist=0;
	n=$(echo "${m}+1"|bc)
	p=$(echo "${m}/2"|bc) #for thresh value
	
	while [ $k -lt 10 ]
	do
		#need EACH coordinate for distance
		x=`cat ${1} | sed 's/ * / /g' | grep "ATOM ${id[$m]} " | cut -d ' ' -f${k}`
		y=`cat ${1} | sed 's/ * / /g' | grep "ATOM ${id[$n]} " | cut -d ' ' -f${k}`
		change=$(echo "(${x}-(${y}))^2"|bc);
		dist=$(echo "${dist}+${change}"|bc);
		k=$(($k+1));
	done
	#calculating distance for THIS contact
	dist=$(echo "sqrt(${dist})"|bc);
	thresh=${thresh[$p]}
	bool=$(echo "${dist} <= ${thresh}"|bc);
	if [ $bool -eq 1 ]
	then
		q=$(echo "${q}+1"|bc);
	else 
		q0=$(echo "e (-${scale}*(${dist}-${thresh}))"|bc -l); #INSERT THRESHOLD VALUES HERE
		# Make sure we never go above 1.
		if [[ $q0 > 1. ]]; then
			q0=1.
		fi
		q=$(echo "${q}+${q0}"|bc);
		#for incrementing when the distance is greater than the threshold bond length
	fi
	i=$(($i+1));
done

#now compute the entire Q parameter
q=$(echo "${q} / ${bonds}" |bc -l);

# print result
printf "%1.3f\n" $q;


exit
