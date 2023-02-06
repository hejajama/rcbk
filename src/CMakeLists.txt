configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gitsha1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp" @ONLY)


add_executable(
	rcbk 	
	main.cpp
	amplitude.cpp
	gbw.cpp
	ic.cpp
	ic_datafile.cpp
	ic_special.cpp
	mv.cpp
	solver.cpp	
	interpolation.cpp
	config.cpp
	${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp
	)

include_directories(.)

target_link_libraries(
	rcbk
	PRIVATE
	GSL::gsl
	GSL::gslcblas
)


