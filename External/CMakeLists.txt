if (${CMAKE_VERSION} VERSION_LESS "3.14") 
    message ("Unable to clone Eigen and STB repos")
else()
	include(FetchContent)

	FetchContent_Declare(
	  Eigen
	  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
	  GIT_TAG master
	  GIT_SHALLOW TRUE
	  GIT_PROGRESS TRUE)

	set(EIGEN_BUILD_DOC OFF)
	set(BUILD_TESTING OFF)
	set(EIGEN_BUILD_PKGCONFIG OFF)
	set( OFF)
	FetchContent_MakeAvailable(Eigen)

	FetchContent_Declare(
	stb
	GIT_REPOSITORY https://github.com/nothings/stb.git
	GIT_TAG master
	)

	FetchContent_MakeAvailable(eigen-master stb)
endif()