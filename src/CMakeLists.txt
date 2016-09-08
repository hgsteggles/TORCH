cmake_minimum_required(VERSION 2.6)

set(TORCH_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Torch/Torch.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/MPI/MPI_Wrapper.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/IO/DataPrinter.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/IO/DataReader.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/IO/Logger.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/IO/ProgressBar.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/IO/StreamGZ.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/IO/FileManagement.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Torch/Constants.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Torch/Converter.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Torch/Parameters.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Fluid/Fluid.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Fluid/GridCellCollection.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Fluid/Grid.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Fluid/GridCell.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Fluid/Star.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Fluid/PartitionManager.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Integrators/Hydro.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Integrators/Riemann.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Integrators/SlopeLimiter.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Integrators/Radiation.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Integrators/Thermodynamics.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Integrators/SplineData.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Misc/Timer.cpp)

include_directories("${TORCH_SOURCE_DIR}/lib")
include_directories("${TORCH_SOURCE_DIR}/lib/lua-5.2.3")
include_directories("${TORCH_SOURCE_DIR}/src")

add_executable(torch ${TORCH_SRCS})
#target_link_libraries(radio ${LUA_LIBRARIES} cfitsio)
target_link_libraries(torch ${TORCH_SOURCE_DIR}/lib/liblua.a dl 
						${MPI_CXX_LIBRARIES} ${ZLIB_LIBRARIES})
