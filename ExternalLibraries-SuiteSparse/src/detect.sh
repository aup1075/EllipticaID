#! /bin/bash

################################################################################
# Prepare
################################################################################

# Set up shell
if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
    set -x                      # Output commands
fi
set -e                          # Abort on errors



################################################################################
# Search
################################################################################

if [ -z "${SUITESPARSE_DIR}" ]; then
    echo "BEGIN MESSAGE"
    echo "SuiteSparse selected, but SuiteSparse_DIR not set."
    echo "END MESSAGE"
else
    echo "BEGIN MESSAGE"
    echo "Using SuiteSparse in ${SUITESPARSE_DIR}"
    echo "END MESSAGE"
fi

THORN=SuiteSparse



################################################################################
# Build
################################################################################

if [ -z "${SUITESPARSE_DIR}"                                                 \
     -o "$(echo "${SUITESPARSE_DIR}" | tr '[a-z]' '[A-Z]')" = 'BUILD' ]
then
    echo "BEGIN MESSAGE"
    echo "Using bundled SuiteSparse..."
    echo "END MESSAGE"
    
    # Check for required tools. Do this here so that we don't require
    # them when using the system library.
    if [ "x$TAR" = x ] ; then
        echo 'BEGIN ERROR'
        echo 'Could not find tar command.'
        echo 'Please make sure that the (GNU) tar command is present,'
        echo 'and that the TAR variable is set to its location.'
        echo 'END ERROR'
        exit 1
    fi
    if [ "x$PATCH" = x ] ; then
        echo 'BEGIN ERROR'
        echo 'Could not find patch command.'
        echo 'Please make sure that the patch command is present,'
        echo 'and that the PATCH variable is set to its location.'
        echo 'END ERROR'
        exit 1
    fi

    # Set locations
    NAME=SuiteSparse
    SRCDIR="$(dirname $0)"
    BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
    if [ -z "${SUITESPARSE_INSTALL_DIR}" ]; then
        INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
    else
        echo "BEGIN MESSAGE"
        echo "Installing SuiteSparse into ${SUITESPARSE_INSTALL_DIR}"
        echo "END MESSAGE"
        INSTALL_DIR=${SUITESPARSE_INSTALL_DIR}
    fi
    SUITESPARSE_BUILD=1
    SUITESPARSE_DIR=${INSTALL_DIR}
else
    SUITESPARSE_BUILD=
    DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
    if [ ! -e ${DONE_FILE} ]; then
        mkdir ${SCRATCH_BUILD}/done 2> /dev/null || true
        date > ${DONE_FILE}
    fi
fi



################################################################################
# Configure Cactus
################################################################################

# Pass configuration options to build script
echo "BEGIN MAKE_DEFINITION"
echo "SUITESPARSE_BUILD          = ${SUITESPARSE_BUILD}"
echo "SUITESPARSE_EXTRA_LIB_DIRS = ${SUITESPARSE_EXTRA_LIB_DIRS}"
echo "SUITESPARSE_EXTRA_LIBS     = ${SUITESPARSE_EXTRA_LIBS}"
echo "SUITESPARSE_INSTALL_DIR    = ${SUITESPARSE_INSTALL_DIR}"
echo "END MAKE_DEFINITION"

# Set options -- TODO: CAN WE PUT UMFPACK HERE?
SUITESPARSE_INC_DIRS="${SUITESPARSE_DIR}/Export/Include ${SUITESPARSE_DIR}/Include"
SUITESPARSE_LIB_DIRS="${SUITESPARSE_DIR}/Lib ${SUITESPARSE_EXTRA_LIB_DIRS}"
SUITESPARSE_LIBS="${SUITESPARSE_EXTRA_LIBS}"

SUITESPARSE_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${SUITESPARSE_INC_DIRS})"
SUITESPARSE_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${SUITESPARSE_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "SUITESPARSE_DIR      = ${SUITESPARSE_DIR}"
echo "SUITESPARSE_INC_DIRS = ${SUITESPARSE_INC_DIRS}"
echo "SUITESPARSE_LIB_DIRS = ${SUITESPARSE_LIB_DIRS}"
echo "SUITESPARSE_LIBS     = ${SUITESPARSE_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(SUITESPARSE_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(SUITESPARSE_LIB_DIRS)'
echo 'LIBRARY           $(SUITESPARSE_LIBS)'
