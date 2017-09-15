#!/bin/bash

#  Script module updating info_int2lm.f90 (which prints version information
#  on binary call to the stdout)
#
# Author  : Oliver Fuhrer (oliver.fuhrer@meteoswiss.ch)
#
# History :
#   25.10.2010  fuo     First release
#   30.11.2013  fuo     Modification to work with long lines
#
# Notes   :
# 1. This script must receive a fresh .config file as argument 1
# 2. It will update the info_int2lm.f90 file in the directory received as argument 2

# set some globals (need to be manually updated)
code_name="int2lm"
code_version="2.0.4"

# parse command line arguments
fconfig="$1"
srcdir="$2"

# check if file is available
if [ ! -f "${srcdir}/info_int2lm.f90" ] ; then
  echo "ERROR: could not find ${srcdir}/info_int2lm.f90 file" 1>&2
  exit 1
fi

# get host and date
creation_host=`hostname | sed 's/[0-9]*$//g'`
creation_date=`date`
creation_user=`whoami`

# get source code checksum (excluding info_int2lm.f90 file)
code_files=`/bin/ls -1d ${srcdir}/*.f90 | grep -v 'info_int2lm_f90\.f90'`
code_checksum=`cksum ${code_files} | cksum | awk '{print $1}'`

# get compiler information
function get_config {
  cat $fconfig | \
  grep "$1" | \
  cut -d ":" -f2 | \
  sed 's/^ *//g' | \
  sed 's/  */ /g' | \
  sed 's/ *$//g' | \
  tr -d "\n" | \
  tr "\." "%" | \
  sed 's/^/___________________________________________________/g' | \
  fold -w 100 -s | \
  sed 's/^__*//g' | \
  tr "%" "." | \
  tr "\n" "@" | \
  sed 's/@ */@@@/g' | \
  sed "s|@@@|' // \\\\\&\\\n    '|g"
}
compiler=`get_config 'Compiler command'`
cversion=`get_config 'Compiler version'`
cincludes=`get_config 'Compiler includes'`
cflags=`get_config 'Compiler flags'`
linker=`get_config 'Linker command'`
lversion=`get_config 'Linker version'`
lflags=`get_config 'Linker flags'`
llibraries=`get_config 'Linker libraries'`
string=`get_config 'Compiler *:'`

# information related to version control system
svn info . 1>/dev/null 2>/dev/null
if [ "$?" -ne 0 ] ; then
  
  # this directory is not under version control
  INFO_LibraryName="${code_name}${code_version}"
  INFO_RevisionTag="(missing)"
  INFO_CheckinDate="(missing)"
  INFO_RevisionNumber="(missing)"
  INFO_CheckoutDate="(missing)"
  INFO_ProductionDate="(missing)"
  INFO_CodeIsModified="unknown"

else
  
  # this is a subversion working copy
  INFO_LibraryName="${code_name}${code_version}"
  INFO_RevisionTag=`cd ..; svn info | grep '^URL: ' | sed 's/^URL: //g'`
  INFO_CheckinDate=`cd ..; svn info | grep '^Last Changed Date: ' | sed 's/^Last Changed Date: //g'`
  INFO_RevisionNumber=`cd ..; svnversion`
  INFO_CheckoutDate="(missing)"
  INFO_ProductionDate="(missing)"
  echo "$INFO_RevisionNumber" | grep 'M' 2>/dev/null 1>/dev/null
  if [ "$?" -eq 0 ] ; then
    INFO_CodeIsModified="true"
  else
    INFO_CodeIsModified="false"
  fi
  
fi

# information related to compilation
INFO_CodeChecksum="${code_checksum}"
INFO_CompilerCall="${compiler} ${cflags}"
INFO_CompilerVersion="${cversion}"
INFO_DefinedMacros="${cincludes}"
INFO_UndefinedMacros="(missing)"
INFO_DebugOptions="(missing)"
INFO_LinkOptions="${linker} ${llibraries}"
INFO_CompiledBy="${creation_user}"
INFO_CompileTime="${creation_date}"
INFO_CompileMachine="${creation_host}"

# information which will have to be defined runtime
INFO_StartTime=""
INFO_BinaryName="${code_name}"
INFO_RunMachine=""
INFO_Nodes=""
INFO_Domain=""
INFO_Options=""

# replace the placeholders in info_int2lm.f90 with the correct strings
ftmp="/tmp/info_int2lm_f90.`whoami`.$$"
ierr=0
rm -f ${ftmp} 2>/dev/null
cp ${srcdir}/info_int2lm.f90 ${ftmp}
function replace_placeholder {
  name=$1
  eval value=\$$name
  cat ${ftmp} | perl -0777 -pe 's|'"${name} *= *'"'.*?([^\&] *\n)|'"${name} = '${value}'"'\n|igs' > ${ftmp}.1
  ((ierr=ierr+$?))
  mv -f ${ftmp}.1 ${ftmp}
  ((ierr=ierr+$?))
}
replace_placeholder INFO_LibraryName
replace_placeholder INFO_LibraryName
replace_placeholder INFO_RevisionTag
replace_placeholder INFO_CheckinDate
replace_placeholder INFO_RevisionNumber
replace_placeholder INFO_CheckoutDate
replace_placeholder INFO_ProductionDate
replace_placeholder INFO_CodeIsModified
replace_placeholder INFO_CodeChecksum
replace_placeholder INFO_CompilerCall
replace_placeholder INFO_CompilerVersion
replace_placeholder INFO_DefinedMacros
replace_placeholder INFO_UndefinedMacros
replace_placeholder INFO_DebugOptions
replace_placeholder INFO_LinkOptions
replace_placeholder INFO_StartTime
replace_placeholder INFO_CompiledBy
replace_placeholder INFO_CompileTime
replace_placeholder INFO_CompileMachine
replace_placeholder INFO_BinaryName
replace_placeholder INFO_RunMachine
replace_placeholder INFO_Nodes
replace_placeholder INFO_Domain
replace_placeholder INFO_Options
  
# check if everything ok
if [ "$ierr" -eq 0 ] ; then
  mv ${ftmp} ${srcdir}/info_int2lm.f90
fi
/bin/rm -f ${ftmp}*

exit 0



