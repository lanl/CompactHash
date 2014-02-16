AC_INIT(hash, v1.0)

AC_CONFIG_AUX_DIR(./config)

AM_INIT_AUTOMAKE([foreign])

AC_CONFIG_HEADERS([config.h hash/config.h])

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

AC_ARG_ENABLE( debug,
 [--enable-debug                     - Compiles code with debug flags ],
 )

AC_ARG_ENABLE( bounds-check,
 [--enable-bounds-check              - Compiles code with bounds check flags ],
 )

AC_ARG_ENABLE( prof,
 [--enable-prof                      - Compiles code with profiling flags ],
 )

AC_ARG_ENABLE( mpi, 
 [--enable-mpi                       - Compiles code with mpi library ],
 enable_mpi=no)

AC_ARG_ENABLE( mpe, 
 [--enable-mpe                       - Compiles code with mpe library ],
 enable_mpe=yes)

AC_ARG_ENABLE( opencl, 
 [--enable-opencl                    - Compiles code with opencl library ],
 enable_opencl=yes)

#AM_WITH_DMALLOC

dnl so we don't have shared libraries
dnl AC_ENABLE_STATIC(yes)

# Checks for programs
AC_PROG_CC
AC_PROG_CXX
AM_PROG_CC_C_O
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_LIBTOOL

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

if test "${enable_mpe}" != "no"; then
  AC_PATH_PROG(MPE_BIN, clog2alog)

  if test ! -z "${MPE_BIN}"; then
    MPE_BIN=`dirname ${MPE_BIN}`
    MPEHOME=`dirname ${MPE_BIN}`
  fi

  if test -z "${MPE_BIN}"; then
    if test -z "${MPEHOME}"; then
      FOUND="no"
      AC_CHECK_FILE(${HOME}/mpe-1.9.1/include/mpe.h, MPEHOME="${HOME}/mpe-1.9.1" && FOUND="yes")
      if test "${FOUND}" = "no"; then
        AC_CHECK_FILE(${HOME}/mpe/include/mpe.h, MPEHOME="${HOME}/mpe" && FOUND="yes")
      fi
    else
      AC_CHECK_FILE(${MPEHOME}/include/mpe.h, INCLUDES="${INCLUDES} -I${MPEHOME}/include")
      LIBS="${LIBS} -L${MPEHOME}/lib -lmpe_nompi -L/usr/X11/lib -lX11 -lm"
      AC_DEFINE([USE_MPE],"yes","Use MPE for graphics")
      if test "${RPATH_LIBS}x" = "x"; then
        RPATH_LIBS="${MPEHOME}/lib"
      else
        RPATH_LIBS="${RPATH_LIBS}:${MPEHOME}/lib"
      fi
    fi
  fi
fi

if test "${enable_opencl}" != "no"; then
  case "$host" in
    *linux*)
      found="no"
      AC_CHECK_PROG(found, lspci, yes)
      if test ! -z "${OPENCL_INCLUDE}" ; then 
        INCLUDES="${INCLUDES} -I/${OPENCL_INCLUDE}"
        if test ! -z "${OPENCL_LIBS}" ; then
          LIBS="${LIBS} -L${OPENCL_LIBS} -lOpenCL"
          if test "${RPATH_LIBS}x" = "x"; then
            RPATH_LIBS="${OPENCL_LIBS}"
          else
            RPATH_LIBS="${RPATH_LIBS}:${OPENCL_LIBS}"
          fi
        fi
      elif test "$found" = "no"; then
         found="no";
         AC_CHECK_FILE(/projects/opt/cuda/cuda-5.5/include/CL/cl.h, INCLUDES="${INCLUDES} -I/projects/opt/cuda/cuda-5.5/include" && found="yes")
         if test "$found" = "no"; then
            AC_CHECK_FILE( /projects/opt/opencl1d2_additional_files/include/CL/cl.h, INCLUDES="${INCLUDES} -I/projects/opt/opencl1d2_additional_files/include" && found="yes")
         fi
         if test "$found" = "no"; then
            AC_CHECK_FILE( /projects/opt/AMD_APP_SDK/AMDAPP-2.8/include/CL/cl.h, INCLUDES="${INCLUDES} -I/projects/opt/AMD_APP_SDK/AMDAPP-2.8/include" && found="yes")
         fi
         if test "$found" = "no"; then
            AC_CHECK_FILE( /opt/intel/opencl-1.2-3.0.67279/include/CL/cl.h, INCLUDES="${INCLUDES} -I/opt/intel/opencl-1.2-3.0.67279/include" && found="yes")
         fi
         found="no"
         if test "$found" = "no"; then
           AC_CHECK_FILE(/usr/lib/libOpenCL.so, LIBS="${LIBS} -lOpenCL" && found="yes")
           if test "${found}" = "yes"; then
             if test "${RPATH_LIBS}x" = "x"; then
               RPATH_LIBS="/usr/lib"
             else
               RPATH_LIBS="${RPATH_LIBS}:/usr/lib"
             fi
           fi
         fi
         if test "$found" = "no"; then
           AC_CHECK_FILE(/usr/lib64/libOpenCL.so, LIBS="${LIBS} -lOpenCL" && found="yes")
           if test "${found}" = "yes"; then
             if test "${RPATH_LIBS}x" = "x"; then
               RPATH_LIBS="/usr/lib64"
             else
               RPATH_LIBS="${RPATH_LIBS}:/usr/lib64"
             fi
           fi
         fi
         if test "$found" = "no"; then
           AC_CHECK_FILE(/projects/opt/AMD_APP_SDK/AMDAPP-2.8/lib/x86_64/libOpenCL.so, LIBS="${LIBS} -L/projects/opt/AMD_APP_SDK/AMDAPP-2.8/lib/x86_64 -lOpenCL" && found="yes")
           if test "${found}" = "yes"; then
             if test "${RPATH_LIBS}x" = "x"; then
               RPATH_LIBS="/projects/opt/AMD_APP_SDK/AMDAPP-2.8/lib/x86_64"
             else
               RPATH_LIBS="${RPATH_LIBS}:/projects/opt/AMD_APP_SDK/AMDAPP-2.8/lib/x86_64"
             fi
           fi
         fi
         if test "$found" = "no"; then
           AC_CHECK_FILE(/opt/intel/opencl-1.2-3.0.67279/lib64/libOpenCL.so, LIBS="${LIBS} -L/opt/intel/opencl-1.2-3.0.67279/lib64 -L/opt/intel/opencl-1.2-3.0.67279/libmic" && found="yes")
           if test "${found}" = "yes"; then
             if test "${RPATH_LIBS}x" = "x"; then
               RPATH_LIBS="/opt/intel/opencl-1.2-3.0.67279/lib64"
             else
               RPATH_LIBS="${RPATH_LIBS}:/opt/intel/opencl-1.2-3.0.67279/lib64"
             fi
           fi
         fi
      elif test "`lspci | grep "VGA compatible controller" | grep -i nVidia | wc -l`" != "0"; then
         found="no"
         AC_CHECK_FILE(/usr/local/NVIDIA_GPU_Computing_SDK/OpenCL/common/inc/CL/cl.h, INCLUDES="${INCLUDES} -I/usr/local/NVIDIA_GPU_Computing_SDK/OpenCL/common/inc" && found="yes")
         if test "$found" = "no"; then
           AC_CHECK_FILE(/usr/local/cuda/include/CL/cl.h, INCLUDES="${INCLUDES} -I/usr/local/cuda/include" && found="yes")
         fi
         if test "$found" = "no"; then
           AC_CHECK_FILE(/projects/opt/cuda/cuda-5.5/include/CL/cl.h, INCLUDES="${INCLUDES} -I/projects/opt/cuda/cuda-5.5/include
" && found="yes")
	fi
	if test "$found" = "no"; then
           AC_CHECK_FILE( /projects/opt/opencl1d2_additional_files/include/CL/cl.h, INCLUDES="${INCLUDES} -I/projects/opt/opencl1d2_additional_files/include" && found="yes")
         fi
         if test "$found" = "no"; then
           AC_CHECK_FILE(/home/cuda/sdk4.0.13/OpenCL/common/inc/CL/cl.h, INCLUDES="${INCLUDES} -I/home/cuda/sdk4.0.13/OpenCL/common/inc" && found="yes")
         fi
         found="no"
         if test "$found" = "no"; then
           AC_CHECK_FILE(/usr/lib/libOpenCL.so, LIBS="${LIBS} -lOpenCL" && found="yes")
           if test "${found}" = "yes"; then
             if test "${RPATH_LIBS}x" = "x"; then
               RPATH_LIBS="/usr/lib"
             else
               RPATH_LIBS="${RPATH_LIBS}:/usr/lib"
             fi
           fi
         fi
         if test "$found" = "no"; then
           AC_CHECK_FILE(/usr/lib64/libOpenCL.so, LIBS="${LIBS} -lOpenCL" && found="yes")
           if test "${found}" = "yes"; then
             if test "${RPATH_LIBS}x" = "x"; then
               RPATH_LIBS="/usr/lib64"
             else
               RPATH_LIBS="${RPATH_LIBS}:/usr/lib64"
             fi
           fi
         fi
      elif test "`lspci | grep "VGA compatible controller" | grep ATI | wc -l`" != "0"; then 
         found="no"
         if test "$found" = "no"; then
           AC_CHECK_FILE(/home/opt/AMDAPP_v2.6/include/CL/cl.h, INCLUDES="${INCLUDES} -I/home/opt/AMDAPP_v2.6/include" && found="yes")
         fi
         if test "$found" = "no"; then
           AC_CHECK_FILE(/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/include/CL/cl.h, INCLUDES="${INCLUDES} -I/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/include" && found="yes")
         fi
         found="no"
         if test "$found" = "no"; then
           AC_CHECK_FILE(/home/opt/AMDAPP_v2.6/lib/x86_64/libOpenCL.so, LIBS="${LIBS} -L/home/opt/AMDAPP_v2.6/lib/x86_64 -lOpenCL" && found="yes")
           if test "${found}" = "yes"; then
             if test "${RPATH_LIBS}x" = "x"; then
               RPATH_LIBS="/home/opt/AMDAPP_v2.6/lib/x86_64"
             else
               RPATH_LIBS="${RPATH_LIBS}:/home/opt/AMDAPP_v2.6/lib/x86_64"
             fi
           fi
         fi
         if test "$found" = "no"; then
           AC_CHECK_FILE(/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/lib/x86_64/libOpenCL.so, LIBS="${LIBS} -L/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/lib/x86_64 -lOpenCL" && found="yes")
           if test "${found}" = "yes"; then
             if test "${RPATH_LIBS}x" = "x"; then
               RPATH_LIBS="/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/lib/x86_64"
             else
               RPATH_LIBS="${RPATH_LIBS}:/projects/opt/AMD_APP_SDK/AMDAPP_v2.6/lib/x86_64"
             fi
           fi
         fi
      fi
      ;;
    *darwin*)
      LIBS="${LIBS} -framework OpenCL"
      ;;
  esac


dnl ---------------------------------------------------------------
   AC_MSG_CHECKING(for OpenCL double type)

   SAVE_C_FLAGS=$CFLAGS
   CFLAGS="$CFLAGS $INCLUDES"

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

