dnl aux.m4 -- Auxiliary macros for WHIZARD's configure.ac
dnl

dnl Horizontal line for readability:
AC_DEFUN([WO_HLINE],
[AC_MSG_NOTICE([--------------------------------------------------------------])])

dnl Message at the beginning of a configure section
AC_DEFUN([WO_CONFIGURE_SECTION],
[WO_HLINE()
AC_MSG_NOTICE([--- ]$1[ ---])
AC_MSG_NOTICE([])
])


dnl Define a variable and export it
dnl WO_SET(variable, value)
AC_DEFUN([WO_SET],
[$1=$2
AC_SUBST($1)])

dnl Add a list of names to the list of output or executable files
dnl WO_OUTFILES(subdir, files)
AC_DEFUN([WO_OUTFILES],
[for file in $2; do
  OUTFILES="$OUTFILES $1/$file";
done])

dnl WO_EXECUTABLES(subdir, files)
AC_DEFUN([WO_EXECUTABLES],
[for file in $2; do
  OUTFILES="$OUTFILES $1/$file";
  EXECUTABLES="$EXECUTABLES $1/$file"
done])

dnl Add a list of names to the list of automake-controlled files
AC_DEFUN([WO_AUTOMAKE],
[for file in $2; do
  AUTOMAKEFILES="$AUTOMAKEFILES $1/$file";
done])

dnl This adds a package which resides in a subdirectory of 'src'
dnl The `variable' is inserted into AC_SUBST; it will refer to 
enl the package subdirectory in Makefiles.  The package
dnl resides in src/XXX and  and optionally has a Makefile (or similar).
dnl WO_PACKAGE(variable, identifier [,Makefiles [,Executables]])
AC_DEFUN([WO_PACKAGE],
[WO_SET($1,src/$2)
ifelse($#, 3, 
[WO_OUTFILES($$1, $3)
WO_AUTOMAKE($$1, $3)])
ifelse($#, 4, 
[WO_OUTFILES($$1, $3)
WO_EXECUTABLES($$1, $4)])
])

dnl This is like WO_PACKAGE, but it calls AC_ARG_ENABLE in addition.
dnl If the package `id' is disabled or the 'file' is not found,
dnl the variable `var' is set to "no".
dnl WO_PACKAGE_ENABLE(var, id, help, file [,Makefiles [,Executables]])
AC_DEFUN([WO_PACKAGE_ENABLE],
[AC_ARG_ENABLE($2,[$3])
if test "$enable_$2" = "no"; then
AC_MSG_CHECKING([for src/$2/$4])
WO_SET($1, no)
AC_MSG_RESULT([(disabled)])
else
AC_CHECK_FILE(src/$2/$4, enable_$2=yes, enable_$2=no)
if test "$enable_$2" = "no"; then
WO_SET($1, no)
else
ifelse($#, 4, [WO_PACKAGE($1, $2)])
ifelse($#, 5, [WO_PACKAGE($1, $2, $5)])
ifelse($#, 6, [WO_PACKAGE($1, $2, $5, $6)])
fi
fi
])

dnl The same, but disabled by default.
dnl If the package `id' is disabled or the 'file' is not found,
dnl the variable `var' is set to "no".
dnl WO_PACKAGE_DISABLE(var, id, help, file [,Makefiles [,Executables]])
AC_DEFUN([WO_PACKAGE_DISABLE],
[AC_ARG_ENABLE($2,[$3])
if test "$enable_$2" = "yes"; then
AC_CHECK_FILE(src/$2/$4, enable_$2=yes, enable_$2=no)
if test "$enable_$2" = "no"; then
WO_SET($1, no)
else
ifelse($#, 4, [WO_PACKAGE($1, $2)])
ifelse($#, 5, [WO_PACKAGE($1, $2, $5)])
ifelse($#, 6, [WO_PACKAGE($1, $2, $5, $6)])
fi
else
enable_$2="no"
AC_MSG_CHECKING([for src/$2/$4])
WO_SET($1, no)
AC_MSG_RESULT([(disabled)])
fi
])


dnl Extension of AC_PATH_PROG: Search for libraries for which the name
dnl is not exactly known (because it may have the version number in it)
dnl Set output variables $var, $var_DIR, $var_LIB accordingly
dnl WO_PATH_LIB(var, id, name-list, path)
AC_DEFUN([WO_PATH_LIB],
[AC_CACHE_CHECK(for $3, wo_cv_path_$1,
[case "$$1" in
  /*)
  wo_cv_path_$1="$$1" # User-supplied path
  ;;
  *)
  IFS="${IFS=   }"; ac_save_ifs="$IFS"; IFS=":"
  ac_dummy=$4
  for ac_dir in $ac_dummy; do
    test -z "$ac_dir" && ac_dir=.
    unset wo_cv_path_$1
    ac_pwd=`pwd`
    if test -d "$ac_dir"; then
      cd $ac_dir
      for ac_word in $3; do
        test -f "$ac_word" && wo_cv_path_$1="$ac_dir/$ac_word"
      done
      cd $ac_pwd
    fi
    if test -n "$wo_cv_path_$1"; then
      break
    fi
  done
  IFS="$ac_save_ifs"
  if test -z "$wo_cv_path_$1"; then
    wo_cv_path_$1="no"
  fi
  ;;
esac
])
$1=$wo_cv_path_$1
if test "$$1" != "no"; then
$1_DIR=`echo $$1 | sed -e 's|/lib$2.*\.a$||'`
$1_LIB=`echo $$1 | sed -e 's|^.*/lib\($2.*\)\.a$|\1|'`
fi
AC_SUBST($1)
AC_SUBST($1_DIR)
AC_SUBST($1_LIB)
])



