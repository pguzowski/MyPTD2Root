# - Basic CMake setup
# Check version meets ou requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.3)
project(MyPTD2Root VERSION 0.1.0)

# Will install module and support program
include(GNUInstallDirs)

# Module uses Falaise, so we need to locate this or fail
find_package(Falaise REQUIRED)

# Ensure our code can see the Falaise etc headers
include_directories(${Falaise_INCLUDE_DIRS})

# Build a dynamic library from our sources
add_library(MyPTD2Root SHARED myptd2root.h myptd2root.cpp)

# - Don't link directly as Bayeux/Falaise will export the
#   needed symbols.
# Ensure dynamic lookup on OS X
if(APPLE)
  set_target_properties(MyPTD2Root
    PROPERTIES LINK_FLAGS "-undefined dynamic_lookup"
    )
endif()

# - configure the example pipeline file
configure_file(myp2r.conf myp2r.conf @ONLY)

install(TARGETS MyPTD2Root DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(PROGRAMS ${PROJECT_SOURCE_DIR}/flmyptd2root.py DESTINATION ${CMAKE_INSTALL_BINDIR})
install(PROGRAMS ${PROJECT_SOURCE_DIR}/myana.py DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES ${PROJECT_BINARY_DIR}/myp2r.conf DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES readme.pdf DESTINATION ${CMAKE_INSTALL_DOCDIR})
