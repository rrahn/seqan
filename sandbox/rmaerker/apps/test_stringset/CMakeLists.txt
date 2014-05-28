# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/rmaerker/apps/test_stringset/CMakeLists.txt
#
# CMakeLists.txt file for test_stringset.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_rmaerker_apps_test_stringset)
message (STATUS "Configuring sandbox/rmaerker/apps/test_stringset")

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
add_executable (test_stringset test_stringset.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_stringset ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (test_stringset)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install test_stringset in ${PREFIX}/bin directory
install (TARGETS test_stringset
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/test_stringset for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (test_stringset)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:test_stringset")
  set (CPACK_PACKAGE_NAME "test_stringset")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "test_stringset")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Rene Rahn <rene.rahn@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Rene Rahn <rene.rahn@fu-berlin.de>")

  seqan_configure_cpack_app (test_stringset "test_stringset")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:test_stringset")

