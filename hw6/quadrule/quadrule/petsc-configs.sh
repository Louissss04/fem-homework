#!/bin/sh
# Shell script to find out PETSc configurations
# $Id: petsc-configs.sh,v 1.5 2020/10/20 12:36:55 zlb Exp $

if test -z "$PETSC_DIR"; then
    if test -r ${MPI_LIB}/petsc/conf/variables; then
	petsc_conf_base=${MPI_LIB}/petsc/conf/variables
    else
	echo 1>&2 "$0: will not use PETSc since \$PETSC_DIR is not set!"
	exit 1
    fi
elif test -r "${PETSC_DIR}/conf/base"; then
    petsc_conf_base="${PETSC_DIR}/conf/base"
elif test -r "${PETSC_DIR}/bmake/common/base"; then
    petsc_conf_base="${PETSC_DIR}/bmake/common/base"
elif test -r "${PETSC_DIR}/conf/variables"; then
    petsc_conf_base="${PETSC_DIR}/conf/variables"
elif test -r "${PETSC_DIR}/lib/petsc/conf/variables"; then
    petsc_conf_base="${PETSC_DIR}/lib/petsc/conf/variables"
else
    echo 1>&2 "$0: error: can't find PETSc configuration files!"
    exit 1
fi

# figure out PETSC_ARCH
if test -z "$PETSC_ARCH" -a -n "$PETSC_dir" -a -x ${PETSC_DIR}/bin/configarch
then
    PETSC_ARCH=`${PETSC_DIR}/bin/configarch 2>/dev/null`
    export PETSC_ARCH
fi

# Build a makefile to retrieve PETSc compiler/linker opts from $petsc_conf_base
trap "/bin/rm -f Makefile.tmp.$$" 0 1 2 3 15
cat <<END >Makefile.tmp.$$
BOPT=O
# Note: PETSc >= 3.2 no longer provides the macro ${PETSC_INCLUDE}
default:
	true
CPPFLAGS:
	@if test -n "\${PETSC_CC_INCLUDES}"; then \\
	    echo \${PETSC_CC_INCLUDES}; \\
	else \\
	    echo \${PETSC_INCLUDE}; \\
	fi
LDFLAGS:
	@echo "-L\${PETSC_LIB_DIR}"
RPATHS:
	@echo \${PETSC_C_SH_LIB_PATH}
LIBS:
	@echo \${PETSC_KSP_LIB}
###	@echo \${PETSC_KSP_LIB_BASIC} \${BLASLAPACK_LIB}
.PHONY: default CPPFLAGS LDFLAGS LIBS RPATHS
include $petsc_conf_base
END

make -s -f Makefile.tmp.$$ "$@"
exit 0
