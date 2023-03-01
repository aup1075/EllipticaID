#! /bin/bash

################################################################################
# Prepare 
################################################################################

# Set up shell
if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
    set -x                      # Output commands
fi
set -e                          # Abort on errors

function printtime() {
  if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
    perl -e 'eval { require("Time/HiRes.pm"); print(&Time::HiRes::gettimeofday()+0," "); }; if($@) { print time(); }'
  fi
}

# Set locations
THORN=SUITESPARSE
NAME=SuiteSparse
SRCDIR="$(dirname $0)"

#==========================================================
#===================== HERE ===============================
#==========================================================

BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
if [ -z "${LORENE_INSTALL_DIR}" ]; then
    INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
else
    echo "BEGIN MESSAGE"
    echo "Installing LORENE into ${LORENE_INSTALL_DIR}"
    echo "END MESSAGE"
    INSTALL_DIR=${LORENE_INSTALL_DIR}
fi
DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
LORENE_DIR=${INSTALL_DIR}
    
# Set up environment
unset LIBS
if echo '' ${ARFLAGS} | grep 64 > /dev/null 2>&1; then
    export OBJECT_MODE=64
fi

echo "LORENE: Preparing directory structure..."
cd ${SCRATCH_BUILD}
mkdir build external done 2> /dev/null || true
rm -rf ${BUILD_DIR} ${INSTALL_DIR}
mkdir ${BUILD_DIR} ${INSTALL_DIR}

echo "LORENE: Unpacking archive..."
pushd ${BUILD_DIR}
${TAR?} xzf ${SRCDIR}/../dist/${NAME}.tar.gz
# Some (ancient but still used) versions of patch don't support the
# patch format used here but also don't report an error using the exit
# code. So we use this patch to test for this
${PATCH?} -p0 < ${SRCDIR}/../dist/patchtest.patch
if [ ! -e Lorene/.patch_tmp ]; then
    echo 'BEGIN ERROR'
    echo 'The version of patch is too old to understand this patch format.'
    echo 'Please set the PATCH environment variable to a more recent '
    echo 'version of the patch command.'
    echo 'END ERROR'
    exit 1
fi

echo "LORENE: Configuring..."
cd ${NAME}

export HOME_LORENE=${BUILD_DIR}/${NAME}
cat > local_settings <<EOF
C = ${CXX}
CFLAGS = ${CFLAGS} $(addprefix -I,${SYS_INC_DIRS}) ${LDFLAGS} -fPIC
INC = -I\$(HOME_LORENE)/src \$(addprefix -I,${GSL_INC_DIRS})
RANLIB = ${RANLIB}
# We don't need dependencies since we always build from scratch
#MAKEDEPEND = ${CXX_DEPEND} \$(INC) \$< ${CXX_DEPEND_OUT} && mv \$@ \$(df).d
MAKEDEPEND = : > \$(df).d
DEPDIR = .deps
FFT_DIR = FFT991
LIB_CXX = ${LIBS}
LIB_LAPACK = ${LAPACK_LIBS} ${BLAS_LIBS}
LIB_PGPLOT =
LIB_GSL = ${GSL_LIBS}
DONTBUILDDEBUGLIB = yes
EOF
if [ -n "$XARGS" ]; then echo "XARGS = $XARGS" >> local_settings; fi
if [ -n "$FIND"  ]; then echo "FIND = $FIND"   >> local_settings; fi

printtime
echo "LORENE: Building..."
# Note that this builds two versions of the library, a "regular"
# version and a "debug" version. Both are identical (since we
# specified identical build options above), and we ignore the "debug"
# version.
${MAKE} cpp fortran export

printtime
echo "LORENE: Installing..."
mv ${BUILD_DIR}/${NAME}/Lib                ${INSTALL_DIR}
mkdir ${INSTALL_DIR}/C++
mv ${BUILD_DIR}/${NAME}/C++/Include        ${INSTALL_DIR}/C++
mkdir ${INSTALL_DIR}/Export
mkdir ${INSTALL_DIR}/Export/C++
mv ${BUILD_DIR}/${NAME}/Export/C++/Include ${INSTALL_DIR}/Export/C++
popd

echo "LORENE: Cleaning up..."
rm -rf ${BUILD_DIR}

date > ${DONE_FILE}
echo "LORENE: Done."
