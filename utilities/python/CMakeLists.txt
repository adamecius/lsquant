cmake_minimum_required(VERSION 2.4.0)

# CMake useful variables
project("lsquant_utilities")

set (CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR} )
set( ENV_PYTHON_ROOT $ENV{CONDA_PREFIX})


install(DIRECTORY lib/ DESTINATION ${ENV_PYTHON_ROOT}/lib/python${PYTHON_VERSION}
        PATTERN "lib/*.py"
        PERMISSIONS OWNER_EXECUTE OWNER_READ
        GROUP_EXECUTE GROUP_READ)

#install(DIRECTORY bin/ DESTINATION test/bin/
#        PATTERN "bin/*.py"
#        PERMISSIONS OWNER_EXECUTE OWNER_READ
#        GROUP_EXECUTE GROUP_READ)
