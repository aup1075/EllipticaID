##! /bin/bash
## TODO? this file should be checked.
#
#################################################################################
## Prepare
#################################################################################
#
## Set up shell
#if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
#    set -x                      # Output commands
#fi
#set -e                          # Abort on errors
#
#
#
#################################################################################
## Search
#################################################################################
#
#if [ -z "${EIDR_DIR}" ]; then
#    echo "BEGIN MESSAGE"
#    echo "EIDR selected, but EIDR_DIR not set."
#    echo "END MESSAGE"
#else
#    echo "BEGIN MESSAGE"
#    echo "Using EIDR in ${EIDR_DIR}"
#    echo "END MESSAGE"
#fi
#
#THORN=EIDR
#
#
#
#################################################################################
## Build
#################################################################################
#
##if [ -z "${EIDR_DIR}"                                                 \
##     -o "$(echo "${EIDR_DIR}" | tr '[a-z]' '[A-Z]')" = 'BUILD' ]
##then
##    echo "BEGIN MESSAGE"
##    echo "Using bundled EIDR..."
##    echo "END MESSAGE"
##    
##    # Check for required tools. Do this here so that we don't require
##    # them when using the system library.
##    if [ "x$TAR" = x ] ; then
##        echo 'BEGIN ERROR'
##        echo 'Could not find tar command.'
##        echo 'Please make sure that the (GNU) tar command is present,'
##        echo 'and that the TAR variable is set to its location.'
##        echo 'END ERROR'
##        exit 1
##    fi
##    if [ "x$PATCH" = x ] ; then
##        echo 'BEGIN ERROR'
##        echo 'Could not find patch command.'
##        echo 'Please make sure that the patch command is present,'
##        echo 'and that the PATCH variable is set to its location.'
##        echo 'END ERROR'
##        exit 1
##    fi
##
##    # Set locations
##    NAME=EIDR
##    SRCDIR="$(dirname $0)"
##    BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
##    if [ -z "${EIDR_INSTALL_DIR}" ]; then
##        INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
##    else
##        echo "BEGIN MESSAGE"
##        echo "Installing EIDR into ${EIDR_INSTALL_DIR}"
##        echo "END MESSAGE"
##        INSTALL_DIR=${EIDR_INSTALL_DIR}
##    fi
##    EIDR_BUILD=1
##    EIDR_DIR=${INSTALL_DIR}
##else
##    EIDR_BUILD=
##    DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
##    if [ ! -e ${DONE_FILE} ]; then
##        mkdir ${SCRATCH_BUILD}/done 2> /dev/null || true
##        date > ${DONE_FILE}
##    fi
##fi
#
#
#
################################################################################
# Configure Cactus
################################################################################

# Pass configuration options to build script
echo "BEGIN MAKE_DEFINITION"
echo "EIDR_BUILD          = ${EIDR_BUILD}"
echo "EIDRBHNS_EXTRA_LIB_DIRS = ${EIDR_EXTRA_LIB_DIRS}"
echo "EIDRBHNS_EXTRA_LIBS     = ${EIDR_EXTRA_LIBS}"
echo "EIDRBHNS_INSTALL_DIR    = ${EIDR_INSTALL_DIR}"
echo "END MAKE_DEFINITION"

# Set options
EIDRBHNS_INC_DIRS="${EIDR_DIR}/include"
EIDRBHNS_LIB_DIRS="${EIDR_DIR}/lib"
EIDRBHNS_LIBS="${EIDR_EXTRA_LIBS}"

EIDRBHNS_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${EIDR_INC_DIRS})"
EIDRBHNS_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${EIDR_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "EIDR_DIR      = ${EIDR_DIR}"
echo "EIDRBHNS_INC_DIRS = ${EIDR_INC_DIRS}"
echo "EIDRBHNS_LIB_DIRS = ${EIDR_LIB_DIRS}"
echo "EIDRBHNS_LIBS     = ${EIDR_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(EIDRBHNS_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(EIDRBHNS_LIB_DIRS)'
echo 'LIBRARY           $(EIDRBHNS_LIBS)'
