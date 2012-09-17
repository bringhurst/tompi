AC_INIT([libtompi],[0.1],[github-hpc@lanl.gov],[libtompi],[http://github.com/hpc/libtompi])

# Keep artifacts in directories as much as possible.
AC_CONFIG_AUX_DIR(auxdir)
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall no-define])
AC_CONFIG_HEADERS([config.h])

# Check for compilers.
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_LEX
AC_PROG_MAKE_SET

# Pull in libtool. We need to support ancient autotools, so don't use LT_INIT.
#LT_INIT
AC_PROG_LIBTOOL

# Checks for typedefs, structs, and compiler options.
AC_TYPE_SIZE_T
AC_C_INLINE

# Checks for library functs.
AC_FUNC_MALLOC

# Define the library version.
AC_SUBST([LIBTOMPI_SO_VERSION], [0:0:0])
AC_SUBST([LIBTOMPI_API_VERSION], [0.1])

AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 examples/mpich/Makefile
                 include/Makefile
                 mpicc/Makefile
                 src/Makefile
                 src/coll/Makefile
                 src/comm/Makefile
                 src/context/Makefile
                 src/error/Makefile
                 src/group/Makefile
                 src/misc/Makefile
                 src/pt2pt/Makefile
                 src/topo/Makefile
                 src/types/Makefile])
                 libtompi.pc:libtompi.pc.in])

AC_OUTPUT
