# Switch between Intel and GCC:
set(USEINTEL FALSE)

# Set compiler flags / options:
if(USECUDA)
    set(USER_CXX_FLAGS "-std=c++17 -fopenmp")
    set(USER_CXX_FLAGS_RELEASE "-O3")
    add_definitions(-DRESTRICTKEYWORD=__restrict__)
else()
    set(USER_CXX_FLAGS "-std=c++17")
    set(USER_CXX_FLAGS_RELEASE "-O3 -march=native")
    set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

    set(USER_FC_FLAGS "-fdefault-real-8 -fdefault-double-8 -fPIC -ffixed-line-length-none -fno-range-check")
    set(USER_FC_FLAGS_RELEASE "-DNDEBUG -O3 -march=native")

    add_definitions(-DRESTRICTKEYWORD=__restrict__)
endif()

set(LIBS -L/user-environment/env/default/lib -L/user-environment/env/default/lib64 netcdf fftw3 fftw3f hdf5)
set(INCLUDE_DIRS /user-environment/env/default/include)

if(USECUDA)
    set(CMAKE_CUDA_ARCHITECTURES "90a")
    set(USER_CUDA_NVCC_FLAGS "--expt-relaxed-constexpr -lineinfo")
    set(USER_CUDA_NVCC_FLAGS_RELEASE "-DNDEBUG")
    set(USER_CUDA_NVCC_FLAGS_DEBUG "-O0 -g -DCUDACHECKS")
    add_definitions(-DRTE_RRTMGP_GPU_MEMPOOL_CUDA)
endif()
