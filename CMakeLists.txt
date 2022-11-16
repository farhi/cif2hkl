cmake_minimum_required(VERSION 3.12...3.24)

#First a bit of magic to extract the version number from cif2hkl.F90, where we
#expect there will be a single line like >> VERSION="1.2.1"<< somewhere. This is
#a bit backwards from how one might do it in a purely cmake-based project, but
#hopefully it will be good enough for our purposes. One could also have used a
#command like >>git describe --tags --dirty --match "v[0-9]*"<< to extract it
#from the last git tag of the form "vx.y.z", however that would not work when
#the code was downloaded via an archive rather than git cloned. The present
#scheme should hopefully always make it compatible with what the installed
#"cif2hkl --version" will report.

set( cif2hkl_src_name "cif2hkl.F90" )
file( READ "${CMAKE_CURRENT_SOURCE_DIR}/${cif2hkl_src_name}" tmp )
string( REGEX REPLACE "\r?\n" ";" tmp "${tmp}" )
set( cif2hkl_extracted_version "" )
foreach( line ${tmp} )
  if ( line MATCHES "VERSION" )
    string(REGEX MATCH "^[ \\t]*VERSION[ \\t]*=[ \\t]*\"([0-9]*)\\.([0-9]*)\\.([0-9]*)\"" line "${line}")
    if ( line )
      set ( line "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}.${CMAKE_MATCH_3}" )
      if ( line MATCHES "^[1-9][0-9]*\\.[0-9]+\\.[0-9]+$" )
        set( cif2hkl_extracted_version "${line}" )
        break()
      endif()
    endif()
  endif()
endforeach()
if ( NOT cif2hkl_extracted_version )# MATCHES "^[1-9][0-9]*\\.[0-9]+\\.[0-9]+$" )
  message(FATAL_ERROR "Could not extract version number from ${cif2hkl_src_name}")
endif()
message( STATUS "Extracted cif2hkl VERSION from ${cif2hkl_src_name}: ${cif2hkl_extracted_version}")

#Now declare the project:
project( cif2hkl VERSION ${cif2hkl_extracted_version} LANGUAGES Fortran )

#Pick installation directory below CMAKE_INSTALL_PREFIX in standard manner,
#allowing both the standard CMAKE_INSTALL_BINDIR and the project specific
#cif2hkl_BINDIR as means to override:
include( GNUInstallDirs )
if ( NOT cif2hkl_BINDIR )
  set( cif2hkl_BINDIR "${CMAKE_INSTALL_BINDIR}" )
endif()

#Now configure the build and installation of the cif2hkl executable:

file( GLOB CFMLSRC CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/CFML/*.f90" )
add_executable( cif2hkl "${PROJECT_SOURCE_DIR}/${cif2hkl_src_name}" ${CFMLSRC} )

#For each compiler flag, only add it if it is supported by the compiler, in the
#hope that this will make us slightly more portable:
include(CheckFortranCompilerFlag)
foreach( compilerflag -ffree-line-length-512 -lm )
  check_fortran_compiler_flag( "${compilerflag}" is_supported )
  if ( is_supported )
    target_compile_options( cif2hkl PRIVATE "${compilerflag}" )
  endif()
endforeach()

install( TARGETS cif2hkl DESTINATION "${cif2hkl_BINDIR}" )

#Tests (this makes a "test" target):

include(CTest)
add_test( NAME cif2hkl-get-version COMMAND cif2hkl --version )
file( GLOB EXAMPLE_TEST_FILES CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/examples/*.cfl"  "${PROJECT_SOURCE_DIR}/examples/*.cif" )
foreach( fn ${EXAMPLE_TEST_FILES} )
  get_filename_component( bn ${fn} NAME )
  add_test( NAME "cif2hkl-convert-${bn}" COMMAND cif2hkl --verbose --xtal --no-outout-files "${fn}" )
endforeach()
