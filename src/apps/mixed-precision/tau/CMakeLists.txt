cmake_minimum_required(VERSION 3.22)
set(app_name tau)
project(${app_name} CXX)

# Universal is a C++ header-only library, so we do not need to build anything.
# The balanced tau number system ships in this repository's include tree, so we
# add both include roots: the project-local one and the Universal submodule.
include_directories(${MPADAO_NUMBER_SYSTEMS_INCLUDE_DIR})
include_directories(${MPADAO_ROOT_DIR}/ext/stillwater-sc/universal/include/sw)

# source files that make up the command
set(SOURCE_FILES
	main.cpp
)

add_executable(${app_name} ${SOURCE_FILES})
set(folder "Applications/Mixed-Precision")
set_target_properties(${app_name} PROPERTIES FOLDER ${folder})

install(TARGETS ${app_name} DESTINATION ${MPADAO_INSTALL_BIN_DIR})
