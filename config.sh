#!/bin/bash

VERSION="0.1"
OK="[OK]"
ERR="[ERR]"
NOTE="      Note>"
#
# This script tries to build / install trajcomp only for
# those modules, that have the needed dependencies.
#

cat <<EOF 
trajcomp $VERSION
-----------------

This is the trajcomp build system. Unlike many other libraries, you
don't have to compile trajcomp in order to use it: It is an header-only
library. By using 

	sudo make install 

the needed header is installed in /usr/local/include and can be used
in your C++ projects on this machine.

However, we provide bindings for the following languages, which need
compilation and installation:
EOF
echo
echo "Do you want to check dependencies for those language bindings?"
echo "Press CTRL-C to stop here." 
read YESNO


# check and use python
echo ""
echo "P Y T H O N     M O D U L E"
echo "---------------------------"
echo 
echo "Looking for Python.h somewhere below /usr/include"
echo 

PYTHON_FOUND=$(find /usr/include -name Python.h |wc -l)
if [ "$PYTHON_FOUND" -eq 0 ]; then 
  PYTHON_ENABLED=0;
  echo "$ERR Did not find python development, please install it."
  echo "$NOTE Ubuntu: sudo apt-get install python-dev"
  echo "$NOTE Debian: sudo apt-get install python-all-dev"
  echo 
else
  PYTHON_ENABLED=1;
  echo "$OK Found $PYTHON_FOUND different versions of python development."
fi
# check octave
echo ""
echo "O C T A V E    M O D U L E"
echo "---------------------------"
echo 
echo "Looking for mkoctfile in the path, calling mkoctfile --version"
echo 

MKOCTFILE_FOUND=$(mkoctfile --version 2>&1)
if [ $? -ne 0 ]; then
  OCTAVE_ENABLED=0;
  echo "$ERR Did not find octave development files."
  echo "$NOTE Debian: sudo apt-get install octave-dev"
  echo 
else
  OCTAVE_ENABLED=1;
  echo "$OK Found $MKOCTFILE_FOUND."
fi

# Results block
echo "R E S U L T S"
echo "Python :" $([ $PYTHON_ENABLED -eq 1 ] && echo "enabled" || echo "disabled" )
echo "Octave :" $([ $OCTAVE_ENABLED -eq 1 ] && echo "enabled" || echo "disabled" )

# Write to 
cat << EOF > config-data.mk
PYTHON_ENABLED=$OCTAVE_ENABLED
OCTAVE_ENABLED=$OCTAVE_ENABLED
EOF
