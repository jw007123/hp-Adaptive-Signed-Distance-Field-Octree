if (${CMAKE_VERSION} VERSION_LESS "3.14") 
    message ("Unable to clone Eigen and STB repos")
else()
	include (FetchContent)

	FetchContent_Declare(
	eigen
	GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
	GIT_TAG f913a4067806312f007cc86d65f74ada26460b20
	)

	FetchContent_Declare(
	stb
	GIT_REPOSITORY https://github.com/nothings/stb.git
	GIT_TAG 8b5f1f37b5b75829fc72d38e7b5d4bcbf8a26d55
	)

	FetchContent_MakeAvailable(eigen-master stb)
endif()