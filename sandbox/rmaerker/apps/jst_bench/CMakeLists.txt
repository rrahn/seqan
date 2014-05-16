# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/rmaerker/apps/jst_bench/CMakeLists.txt
#
# CMakeLists.txt file for jst_bench.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_rmaerker_apps_jst_bench)
message (STATUS "Configuring sandbox/rmaerker/apps/jst_bench")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
#set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (jst_bench jst_bench.cpp
						  jst_bench_base.h
						  jst_bench_io.h
						  jst_bench_find.h)


# Add dependencies found by find_package (SeqAn).
target_link_libraries (jst_bench ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (jst_bench)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install jst_bench in ${PREFIX}/bin directory
install (TARGETS jst_bench
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/jst_bench for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (jst_bench)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:jst_bench")
  set (CPACK_PACKAGE_NAME "jst_bench")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "jst_bench")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Rene Rahn <rene.rahn@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Rene Rahn <rene.rahn@fu-berlin.de>")

  seqan_configure_cpack_app (jst_bench "jst_bench")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:jst_bench")

