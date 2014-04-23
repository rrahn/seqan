# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/apps/jseq_tools/CMakeLists.txt
#
# CMakeLists.txt file for gdf_tools.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_extras_apps_gdf_tools)
message (STATUS "Configuring extras/apps/gdf_tools")

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
								  
add_executable (gdf_converter gdf_converter.cpp
						      gdf_converter_vcf.h
							  gdf_tools_io.h)
									  
# Add dependencies found by find_package (SeqAn).
target_link_libraries (gdf_converter ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (gdf_tools)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install jseq_tools in ${PREFIX}/bin directory
install (TARGETS gdf_converter
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/gdf_tools for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (gdf_tools)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:gdf_tools")
  set (CPACK_PACKAGE_NAME "jseq_tools")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "gdf_tools")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Rene Rahn <rene.rahn@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Rene Rahn <rene.rahn@fu-berlin.de>")

  seqan_configure_cpack_app (jseq_tools "gdf_tools")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:gdf_tools")

