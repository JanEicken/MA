#!/bin/bash

cp template/headertrajcomp.png html

for f in $(find -name *.md |grep -E "^./html"); do 
	echo Processing $f
	HTMLNAME=$(echo $f | rev |cut -d"." -f2- |rev).html
#	echo $HTMLNAME
	# Find out the relative prefix of asset links
	NUMSLASH=$(echo $HTMLNAME | tr -dc "/" |wc -c)
	DEPTH=$(echo $NUMSLASH-2 |bc)
#	echo Depth: $DEPTH
	PREFIX=$(eval "printf '../'%.0s {0..$DEPTH}" |cut -b2-)
#	echo $PREFIX
	HEADERURI="$PREFIX"headertrajcomp.png
#	echo $HEADERURI
	cat template/template-head-beforeheader.html > $HTMLNAME
	echo -n $HEADERURI >> $HTMLNAME
	cat template/template-head-afterheader.html >> $HTMLNAME
	
	markdown $f >> $(echo $f | rev |cut -d"." -f2- |rev).html
	cat template/template-foot.html >> $HTMLNAME

done
