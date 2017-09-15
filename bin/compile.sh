#!/bin/bash

# target name
target=int2lm
host=`hostname | sed 's/[0-9]//g'`

# get number of threads for compilation
if [ $# -eq 1 ] ; then
  nthread=$1
else
  nthread=1
fi

# generation function
genexe() {
  # check arguments
  if [ -z "$1" -o -z "$2" ] ; then
    echo "ERROR: genexe calle with no parameters"
    exit 1
  fi
  exename=${target}_${1}
  exeopts=$2
  # generate executable
  echo "> Generating ${exename} (make -j${nthread} $exeopts)"
  \rm -f ${target} ${exename} ${exename}_*.log 2>/dev/null 1>&2
  log=${exename}_`date "+%y%m%d"`.log
  touch ${log}
  echo "*** date: `date`" >> ${log}
  echo "*** command: make clean" >> ${log}
  make clean 1>>${log} 2>&1
  echo "*** command: make -j${nthread} $exeopts" >> ${log}
  make -j${nthread} VERBOSE=1 $exeopts 1>${log} 2>&1
  if [ ! -f ${target} ] ; then
    echo "ERROR: Could not generate executable ${exename}"
    echo "       See ${log} for details"
    exit 1
  else
    mv -f ${target} ${exename}
  fi
}

# clean directory
\rm ${target}_* *.log 2>/dev/null 1>&2

# no libraries, debugging options, no parallelization, for interactive nodes
#genexe "seqdebug_nolib"  "DEBUG=debug ${target}"

# with libraries, debugging options, no parallelization, for interactive nodes
#genexe "seqdebug"        "seqdebug"

# with libraries, optimization options, no parallelization, for interactive nodes
#genexe "seqopt"          "seqopt"

# no libraries, debugging options, parallelization, for compute nodes
#genexe "pardebug_nolib"  "DEBUG=debug MPI=mpi ${target}"

# with libraries, debugging options, parallelization, for compute nodes
genexe "pardebug"        "pardebug"

# with libraries, optimization options, parallelization, for compute nodes
genexe "paropt"          "paropt"

# goodbye
