AC_INIT(CompactHash, v1.1)

AC_CONFIG_AUX_DIR(./config)

AM_INIT_AUTOMAKE([foreign])

dnl AC_CONFIG_HEADERS([config.h])
dnl AC_CONFIG_HEADERS([hash/config.h])

dnl
dnl Caching is usually WRONG for systems with cross-mounted file systems
dnl  the cache file may correspond to a different system.
dnl
define([AC_CACHE_LOAD], )dnl
define([AC_CACHE_SAVE], )dnl
dnl

#
echo "Configuring hash with $ac_configure_args"
#

AC_ARG_ENABLE( [debug],
  [ AS_HELP_STRING([--enable-debug],
    [compiles code with debug flags @<:@default=no@:>@])],
  [],
  [enable_debug=no]
)

AC_ARG_ENABLE( [bounds-check],
  [ AS_HELP_STRING([--enable-bounds-check],
    [compiles code with bounds check flags @<:@default=no@:>@])],
  [],
  [enable_bounds_check=no]
)

AC_ARG_ENABLE( [prof],
  [ AS_HELP_STRING([--enable-prof],
    [compiles code with profiling flags @<:@default=no@:>@])],
  [],
  [enable_prof=no]
)

AC_ARG_ENABLE( [opencl], 
  [ AS_HELP_STRING([--enable-opencl],
    [compiles code with opencl library @<:@default=yes@:>@])],
  [],
  [enable_opencl=yes]
)

AC_ARG_ENABLE( [detailed-timing], 
  [ AS_HELP_STRING([--enable-detailed-timing],
    [compiles code with detailed timing on @<:@default=no@:>@])],
  [],
  [enable_detailed_timing=no]
)

AC_ARG_ENABLE( [write-mem-usage], 
  [ AS_HELP_STRING([--enable-write-mem-usage],
    [compiles code with write memory usage on @<:@default=no@:>@])],
  [],
  [enable_write_mem_usage=no]
)

#AM_WITH_DMALLOC

dnl so we don't have shared libraries
dnl AC_ENABLE_STATIC(yes)

# Checks for programs
AC_PROG_CC
AC_PROG_CXX
AM_PROG_CC_C_O
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL

AC_CONFIG_MACRO_DIR([m4])

INDENT="indent"
CPREPROCESSOR="cpp -P"

AC_CHECK_PROG([HAVE_ADDR2LINE], [addr2line], [yes])
if test "x$HAVE_ADDR2LINE" == "xyes"; then
   AC_DEFINE([HAVE_ADDR2LINE],[1],[Has addr2line program])
fi

# Checks for header files
AC_HEADER_STDC

# Checks for typedefs, structures, and compiler characteristics
AC_SYS_LARGEFILE

# Checks for library functions

AC_CANONICAL_HOST
case "$host" in
    *linux*)
        CFLAGS="-std=gnu99 -static"
        if test "${enable_debug}" = "yes"; then
          CFLAGS="-g -DDEBUG=1 ${CFLAGS}"
        else
          CFLAGS="-g -O3 ${CFLAGS}"
        fi
        if test "${enable_prof}" = "yes"; then
          CFLAGS="-pg ${CFLAGS}"
          LFLAGS="-pg ${LFLAGS}"
        fi

        if test "${enable_bounds_check}" = "yes"; then
          CFLAGS="-fbounds-check ${CFLAGS}"
        fi

        RPATH_PREFIX="-Wl,"

        ;;
    *darwin*)
        CFLAGS="-std=gnu99 -pedantic -Wno-overlength-strings"
        INDENT="gindent"
        CPREPROCESSOR="gcc -E -P -U__APPLE_CC__"
        if test "${enable_debug}" = "yes"; then
          CFLAGS="-g -DDEBUG=1 ${CFLAGS}"
        else
          CFLAGS="-g -O3 ${CFLAGS}"
        fi

        if test "${enable_bounds_check}" = "yes"; then
          CFLAGS="-fbounds-check ${CFLAGS}"
        fi

        RPATH_PREFIX="-Wl,"

        ;;
    *)
        echo "Unknown machine"
        exit 1
        ;;
esac

if test "${enable_detailed_timing}" = "yes"; then
   CFLAGS="-DDETAILED_TIMING=1 ${CFLAGS}"
fi

if test "${enable_write_mem_usage}" = "no"; then
   CFLAGS="-DWRITE_MEM_USAGE=0 ${CFLAGS}"
fi

if test "${enable_opencl}" != "no"; then
  case "$host" in
    *linux*)
      AC_MSG_CHECKING([OPENCL_INCLUDE environment variable to set OpenCL includes])
      opencl_inc_found="no"
      if test ! -z "${OPENCL_INCLUDE}" ; then 
        AC_MSG_RESULT(-I${OPENCL_INCLUDE})
        opencl_inc_found="yes"
        AM_CPPFLAGS="${AM_CPPFLAGS} -I${OPENCL_INCLUDE}"
      else
        AC_MSG_RESULT(no)
      fi

      opencl_lib_found="no"
      AC_MSG_CHECKING([OPENCL_LIBS environment variable to set OpenCL libs])
      if test ! -z "${OPENCL_LIBS}" ; then
        AC_MSG_RESULT(-L${OPENCL_LIBS} -lOpenCL)
        opencl_lib_found="yes"
        LIBS="${LIBS} -L${OPENCL_LIBS} -lOpenCL"
        if test "${RPATH_LIBS}x" = "x"; then
          RPATH_LIBS="${OPENCL_LIBS}"
        else
          RPATH_LIBS="${RPATH_LIBS}:${OPENCL_LIBS}"
        fi
      else
        AC_MSG_RESULT(no)
      fi

      if test "$opencl_inc_found" = "no"; then
         AC_CHECK_PROG(lspci_found, lspci, yes, no)
         if test "$lspci_found" = "no"; then
            AC_CHECK_FILE(/projects/opt/cuda/cuda-5.5/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/projects/opt/cuda/cuda-5.5/include" && opencl_inc_found="yes")
            if test "$opencl_inc_found" = "no"; then
               AC_CHECK_FILE(/opt/cudatoolkit-6.0/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/opt/cudatoolkit-6.0/include" && opencl_inc_found="yes")
            fi
            if test "$opencl_inc_found" = "no"; then
               AC_CHECK_FILE( /projects/opt/opencl1d2_additional_files/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/projects/opt/opencl1d2_additional_files/include" && opencl_inc_found="yes")
            fi
            if test "$opencl_inc_found" = "no"; then
               AC_CHECK_FILE( /projects/opt/AMD_APP_SDK/AMDAPP-2.8/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/projects/opt/AMD_APP_SDK/AMDAPP-2.8/include" && opencl_inc_found="yes")
            fi
            if test "$opencl_inc_found" = "no"; then
               AC_CHECK_FILE( /opt/intel/opencl-1.2-3.0.67279/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/opt/intel/opencl-1.2-3.0.67279/include" && opencl_inc_found="yes")
            fi

            if test "$opencl_lib_found" = "no"; then
              AC_CHECK_FILE(/opt/cudatoolkit-6.0/lib64/libOpenCL.so, LIBS="${LIBS} -L/opt/cudatoolkit-6.0/lib64 -lOpenCL" && opencl_lib_found="yes")
              if test "${opencl_lib_found}" = "yes"; then
                if test "${RPATH_LIBS}x" = "x"; then
                  RPATH_LIBS="/opt/cudatoolkit-6.0/lib64"
                else
                  RPATH_LIBS="${RPATH_LIBS}:/opt/cudatoolkit-6.0/lib64"
                fi
              fi
            fi
            if test "$opencl_lib_found" = "no"; then
              AC_CHECK_FILE(/usr/lib/libOpenCL.so, LIBS="${LIBS} -lOpenCL" && opencl_lib_found="yes")
              if test "${opencl_lib_found}" = "yes"; then
                if test "${RPATH_LIBS}x" = "x"; then
                  RPATH_LIBS="/usr/lib"
                else
                  RPATH_LIBS="${RPATH_LIBS}:/usr/lib"
                fi
              fi
            fi
            if test "$opencl_lib_found" = "no"; then
              AC_CHECK_FILE(/usr/lib64/libOpenCL.so, LIBS="${LIBS} -lOpenCL" && opencl_lib_found="yes")
              if test "${opencl_lib_found}" = "yes"; then
                if test "${RPATH_LIBS}x" = "x"; then
                  RPATH_LIBS="/usr/lib64"
                else
                  RPATH_LIBS="${RPATH_LIBS}:/usr/lib64"
                fi
              fi
            fi
            if test "$opencl_lib_found" = "no"; then
              AC_CHECK_FILE(/projects/opt/AMD_APP_SDK/AMDAPP-2.8/lib/x86_64/libOpenCL.so, LIBS="${LIBS} -L/projects/opt/AMD_APP_SDK/AMDAPP-2.8/lib/x86_64 -lOpenCL" && opencl_lib_found="yes")
              if test "${opencl_lib_found}" = "yes"; then
                if test "${RPATH_LIBS}x" = "x"; then
                  RPATH_LIBS="/projects/opt/AMD_APP_SDK/AMDAPP-2.8/lib/x86_64"
                else
                  RPATH_LIBS="${RPATH_LIBS}:/projects/opt/AMD_APP_SDK/AMDAPP-2.8/lib/x86_64"
                fi
              fi
            fi
            if test "$opencl_lib_found" = "no"; then
              AC_CHECK_FILE(/opt/intel/opencl-1.2-3.0.67279/lib64/libOpenCL.so, LIBS="${LIBS} -L/opt/intel/opencl-1.2-3.0.67279/lib64 -L/opt/intel/opencl-1.2-3.0.67279/libmic" && opencl_lib_found="yes")
              if test "${opencl_lib_found}" = "yes"; then
                if test "${RPATH_LIBS}x" = "x"; then
                  RPATH_LIBS="/opt/intel/opencl-1.2-3.0.67279/lib64"
                else
                  RPATH_LIBS="${RPATH_LIBS}:/opt/intel/opencl-1.2-3.0.67279/lib64"
                fi
              fi
            fi
         elif test "`lspci | grep "VGA compatible controller" | grep -i nVidia | wc -l`" != "0"; then
            AC_CHECK_FILE(/usr/local/NVIDIA_GPU_Computing_SDK/OpenCL/common/inc/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/usr/local/NVIDIA_GPU_Computing_SDK/OpenCL/common/inc" && opencl_inc_found="yes")
            if test "$opencl_inc_found" = "no"; then
              AC_CHECK_FILE(/usr/local/cuda/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/usr/local/cuda/include" && opencl_inc_found="yes")
            fi
            if test "$opencl_inc_found" = "no"; then
              AC_CHECK_FILE(/projects/opt/cuda/cuda-5.5/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/projects/opt/cuda/cuda-5.5/include" && opencl_inc_found="yes")
	   fi
	   if test "$opencl_inc_found" = "no"; then
              AC_CHECK_FILE( /projects/opt/opencl1d2_additional_files/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/projects/opt/opencl1d2_additional_files/include" && opencl_inc_found="yes")
            fi
            if test "$opencl_inc_found" = "no"; then
              AC_CHECK_FILE(/home/cuda/sdk4.0.13/OpenCL/common/inc/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/home/cuda/sdk4.0.13/OpenCL/common/inc" && opencl_inc_found="yes")
            fi

            if test "$opencl_lib_found" = "no"; then
              AC_CHECK_FILE(/usr/lib/libOpenCL.so, LIBS="${LIBS} -lOpenCL" && opencl_lib_found="yes")
              if test "${opencl_lib_found}" = "yes"; then
                if test "${RPATH_LIBS}x" = "x"; then
                  RPATH_LIBS="/usr/lib"
                else
                  RPATH_LIBS="${RPATH_LIBS}:/usr/lib"
                fi
              fi
            fi
            if test "$opencl_lib_found" = "no"; then
              AC_CHECK_FILE(/usr/lib64/libOpenCL.so, LIBS="${LIBS} -lOpenCL" && opencl_lib_found="yes")
              if test "${opencl_lib_found}" = "yes"; then
                if test "${RPATH_LIBS}x" = "x"; then
                  RPATH_LIBS="/usr/lib64"
                else
                  RPATH_LIBS="${RPATH_LIBS}:/usr/lib64"
                fi
              fi
            fi
         elif test "`lspci | grep "VGA compatible controller" | grep ATI | wc -l`" != "0"; then 
            if test "$opencl_inc_found" = "no"; then
              AC_CHECK_FILE(/home/opt/AMDAPP_v2.6/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/home/opt/AMDAPP_v2.6/include" && opencl_inc_found="yes")
            fi
            if test "$opencl_inc_found" = "no"; then
              AC_CHECK_FILE(/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/include/CL/cl.h, AM_CPPFLAGS="${AM_CPPFLAGS} -I/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/include" && opencl_inc_found="yes")
            fi

            if test "$opencl_lib_found" = "no"; then
              AC_CHECK_FILE(/home/opt/AMDAPP_v2.6/lib/x86_64/libOpenCL.so, LIBS="${LIBS} -L/home/opt/AMDAPP_v2.6/lib/x86_64 -lOpenCL" && opencl_lib_found="yes")
              if test "${opencl_lib_found}" = "yes"; then
                if test "${RPATH_LIBS}x" = "x"; then
                  RPATH_LIBS="/home/opt/AMDAPP_v2.6/lib/x86_64"
                else
                  RPATH_LIBS="${RPATH_LIBS}:/home/opt/AMDAPP_v2.6/lib/x86_64"
                fi
              fi
            fi
            if test "$opencl_lib_found" = "no"; then
              AC_CHECK_FILE(/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/lib/x86_64/libOpenCL.so, LIBS="${LIBS} -L/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/lib/x86_64 -lOpenCL" && opencl_lib_found="yes")
              if test "${opencl_lib_found}" = "yes"; then
                if test "${RPATH_LIBS}x" = "x"; then
                  RPATH_LIBS="/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/lib/x86_64"
                else
                  RPATH_LIBS="${RPATH_LIBS}:/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/lib/x86_64"
                fi
              fi
            fi
         fi
      fi
      dnl For debugging fail
      opencl_inc_found=no
      if test "$opencl_inc_found" = "yes"; then
         AC_DEFINE(HAVE_OPENCL,1,[Has OpenCL])
      else
         AC_MSG_WARN(=================================================================)
         AC_MSG_WARN(OpenCL was not found, but will try to build without it)
         AC_MSG_WARN(If it fails, try setting OPENCL_INCLUDE and OPENCL_LIBS)
         AC_MSG_WARN(to location of include directory and the OpenCL library)
         AC_MSG_WARN(or build without OpenCL with "./configure --disable-opencl".)
         AC_MSG_WARN(If you disable opencl, you may need to go into the hash directory)
         AC_MSG_WARN(and force the rebuild of hash.c and hash.h "make hash.c hash.h".)
         AC_MSG_WARN(You may need to remove the hash.c and hash.h files if make gets)
         AC_MSG_WARN(confused by timestamps)
         AC_MSG_WARN(=================================================================)
      fi
      ;;
    *darwin*)
      LIBS="${LIBS} -framework OpenCL"
      ;;
  esac


dnl ---------------------------------------------------------------
   AC_MSG_CHECKING(for OpenCL double type)

   SAVE_C_FLAGS=$CFLAGS
   CFLAGS="$CFLAGS $AM_CPPFLAGS"

   AC_LANG([C])

   AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
      [
       [#ifdef __APPLE_CC__]
       [#include <OpenCL/OpenCL.h>]
       [#else]
       [#include <CL/cl.h>]
       [#endif]
      ],
      [[cl_double x;]])],
      ac_cv_opencl_double=yes, ac_cv_opencl_double=no
   )
   if test "x$ac_cv_opencl_double" = "xyes"; then
      AC_DEFINE(HAVE_CL_DOUBLE,1,[Has OpenCL double type])
   else
      AC_DEFINE(NO_CL_DOUBLE,1,[Does not have OpenCL double type])
   fi

   CFLAGS=$SAVE_C_FLAGS

   AC_MSG_RESULT($ac_cv_opencl_double)
dnl ---------------------------------------------------------------
fi

if test ! -z "${RPATH_LIBS}"; then
   LIBS="${LIBS} ${RPATH_PREFIX}-rpath ${RPATH_PREFIX}${RPATH_LIBS}"
fi

LIBS="${LIBS} -lm"

AC_SUBST(CC)
AC_SUBST(CFLAGS)
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(LIBS)
AC_SUBST(INDENT)
AC_SUBST(CPREPROCESSOR)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([kdtree/Makefile])
AC_CONFIG_FILES([hash/Makefile])
AC_CONFIG_FILES([oldHash/Makefile])
AC_CONFIG_FILES([timer/Makefile])

AC_OUTPUT

