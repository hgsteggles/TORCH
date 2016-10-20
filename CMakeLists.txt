cmake_minimum_required(VERSION 2.6)
project(TORCH)

if( NOT CMAKE_BUILD_TYPE )
    set( CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: None Debug Release RelWithDebInfo
     MinSizeRel."
            FORCE )
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
   set(ENABLE_CXX11 "-std=c++11")

   execute_process(COMMAND "${CMAKE_CXX_COMPILER} -dumpversion" OUTPUT_VARIABLE GCC_VERSION)
   if (GCC_VERSION VERSION_LESS 4.7)
      set(ENABLE_CXX11 "-std=c++0x")
   endif()
endif()

if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g ${ENABLE_CXX11} -stdlib=libc++")
elseif(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g ${ENABLE_CXX11}")
elseif(WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g ${ENABLE_CXX11}")
endif(APPLE)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

find_package(MPI REQUIRED)
find_package(ZLIB REQUIRED)

#message(WARNING "MPI_CXX_INCLUDE_PATH = " ${MPI_CXX_INCLUDE_PATH})
#message(WARNING "MPI_CXX_LIBRARIES = " ${MPI_CXX_LIBRARIES})
#message(WARNING "MPI_CXX_COMPILER = " ${MPI_CXX_COMPILER})

#message(WARNING "ZLIB_INCLUDE_DIRS = " ${ZLIB_INCLUDE_DIRS})
#message(WARNING "ZLIB_LIBRARIES = " ${ZLIB_LIBRARIES})

include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories(${ZLIB_INCLUDE_DIRS})

function(copy_torch_config filename)
	configure_file(${filename} 
					${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${filename} 
					COPYONLY)
endfunction()

copy_torch_config(config/torch-config.lua)
copy_torch_config(config/torch-setup.lua)
copy_torch_config(config/data2D_025.txt)

add_subdirectory(src)
