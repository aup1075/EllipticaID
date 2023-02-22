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

if [ -z "${ELLIPTICA_DIR}" ]; then
    echo "BEGIN MESSAGE"
    echo "Elliptica selected, but Elliptica_DIR not set."
    echo "END MESSAGE"
else
    echo "BEGIN MESSAGE"
    echo "Using Elliptica in ${ELLIPTICA_DIR}"
    echo "END MESSAGE"
fi

THORN=Elliptica



################################################################################
# Build
################################################################################

if [ -z "${ELLIPTICA_DIR}"                                                 \
     -o "$(echo "${ELLIPTICA_DIR}" | tr '[a-z]' '[A-Z]')" = 'BUILD' ]
then
    echo "BEGIN MESSAGE"
    echo "Using bundled Elliptica..."
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
    NAME=Lorene
    SRCDIR="$(dirname $0)"
    BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
    if [ -z "${ELLIPTICA_INSTALL_DIR}" ]; then
        INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
    else
        echo "BEGIN MESSAGE"
        echo "Installing Elliptica into ${ELLIPTICA_INSTALL_DIR}"
        echo "END MESSAGE"
        INSTALL_DIR=${ELLIPTICA_INSTALL_DIR}
    fi
    ELLIPTICA_BUILD=1
    ELLIPTICA_DIR=${INSTALL_DIR}
else
    ELLIPTICA_BUILD=
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
echo "ELLIPTICA_BUILD          = ${ELLIPTICA_BUILD}"
echo "ELLIPTICA_EXTRA_LIB_DIRS = ${ELLIPTICA_EXTRA_LIB_DIRS}"
echo "ELLIPTICA_EXTRA_LIBS     = ${ELLIPTICA_EXTRA_LIBS}"
echo "ELLIPTICA_INSTALL_DIR    = ${ELLIPTICA_INSTALL_DIR}"
echo "END MAKE_DEFINITION"

# Set options -- TODO: CAN WE PUT UMFPACK HERE?
ELLIPTICA_INC_DIRS="${ELLIPTICA_DIR}/Export/Include ${ELLIPTICA_DIR}/Include"
ELLIPTICA_LIB_DIRS="${ELLIPTICA_DIR}/Lib ${ELLIPTICA_EXTRA_LIB_DIRS}"
ELLIPTICA_LIBS="${ELLIPTICA_EXTRA_LIBS}"

ELLIPTICA_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${ELLIPTICA_INC_DIRS})"
ELLIPTICA_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${ELLIPTICA_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "ELLIPTICA_DIR      = ${ELLIPTICA_DIR}"
echo "ELLIPTICA_INC_DIRS = ${ELLIPTICA_INC_DIRS}"
echo "ELLIPTICA_LIB_DIRS = ${ELLIPTICA_LIB_DIRS}"
echo "ELLIPTICA_LIBS     = ${ELLIPTICA_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(ELLIPTICA_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(ELLIPTICA_LIB_DIRS)'
echo 'LIBRARY           $(ELLIPTICA_LIBS)'
