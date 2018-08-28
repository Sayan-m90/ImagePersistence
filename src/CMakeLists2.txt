

#cmake_minimum_required( VERSION 2.8 FATAL_ERROR )
#project( Persistence of Images )
#cmake_policy( SET CMP0017 NEW )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#file(GLOB ANN RELATIVE "ANN" "*.h")
#file(GLOB GIComplex RELATIVE "GIComplex" "*.cpp" "*.h")
#file(GLOB Graphs RELATIVE "Graphs" "*.cpp" "*.hpp")
#file(GLOB SimPers RELATIVE "SimPers" "*.cpp" "*.h")
#file(GLOB Wrappers RELATIVE "Wrappers" "*.cpp" "*.hpp")

#file(GLOB Wrappersh RELATIVE "Wrappers" "*.hpp")
#file(GLOB GIComplexh RELATIVE "GIComplex" "*.h")
#file(GLOB SimPersh RELATIVE "SimPers" "*.h")
#file(GLOB Graphsh RELATIVE "Graphs" "*.hpp")

#include_directories( ${CMAKE_SOURCE_DIR}/src "${CMAKE_SOURCE_DIR}/GIComplex" "${CMAKE_SOURCE_DIR}/ANN" "${CMAKE_SOURCE_DIR}/Graphs" "${CMAKE_SOURCE_DIR}/SimPers" "${CMAKE_SOURCE_DIR}/Wrappes")
link_directories(${CMAKE_SOURCE_DIR}/lib)

include_directories( ${CMAKE_SOURCE_DIR}/src "${CMAKE_SOURCE_DIR}/src/include" )
link_directories(${CMAKE_SOURCE_DIR}/src/lib)
# Boost
#set(Boost_INCLUDE_DIR /usr/include/boost)
find_package( Boost 1.5 COMPONENTS filesystem program_options system REQUIRED )
include_directories( ${Boost_INCLUDE_DIRS} )
link_directories(${Boost_LIBRARY_DIRS})

#OpenCV
find_package(OpenCV REQUIRED COMPONENTS core imgproc highgui imgcodecs REQUIRED)
#set(OpenCV_DIR "/usr/local/lib")
#set(OpenCV_LIBS "/usr/local/lib")
#highgui imgcodecs imgproc core)

#COMPONENTS core imgproc highgui imgcodecs)
include_directories(${OpenCV_INCLUDE_DIRS})

include_directories( "." ".." ) 

add_executable( pers
	Wrappers/GICWrapper.cpp
	Wrappers/GICWrapper.hpp
	Wrappers/SimPersWrapper.cpp
	Wrappers/SimPersWrapper.hpp
	SimPers/AnnotationMatrix.cpp
	SimPers/SimplexNode.h
	SimPers/SimplicialComplex.cpp
	SimPers/SimplicialComplex.h
	SimPers/UnionFindDeletion.cpp
	Graphs/Graph.hpp
	Graphs/KNN_Graph.cpp
	Graphs/KNN_Graph.hpp
	Graphs/NormalGraph.cpp
	Graphs/NormalGraph.hpp
	GIComplex/ANNSearchSampling.cpp
	GIComplex/ANNSearchSampling.h
	GIComplex/DijkstraAlgorithm.cpp
	GIComplex/DijkstraAlgorithm.h
	GIComplex/FibonacciHeap.cpp
	GIComplex/FibonacciHeap.h
	GIComplex/FlagComplex.h
	GIComplex/FlagComplex.cpp
	GIComplex/gic.cpp
	GIComplex/GIComplex.cpp
	GIComplex/GIComplex.h
	GIComplex/LongVector.h
	GIComplex/PointSet.cpp
	GIComplex/PointSet.h
	GIComplex/SimpleGraph.cpp
	GIComplex/SimpleGraph.h
	GIComplex/SimplexNodeGIC.h
	GIComplex/SimplicialComplexGIC.h
	#ANN
	#GIComplex
	#Graphs
	#SimPers
	#Wrappers
	Collapse.hpp
	Collapse.cpp
	Constants.cpp
	Constants.hpp
	GIC.hpp
	GIC.cpp
	main.cpp
	MortonCode.cpp
	MortonCode.hpp
	RandomHomology.cpp
	RandomHomology.hpp
	Utilities.cpp
	Utilities.hpp
	) 
#target_link_libraries( pers ann.lib )
target_link_libraries( pers ${Boost_LIBRARIES} ${OpenCV_LIBS} ann.lib)
