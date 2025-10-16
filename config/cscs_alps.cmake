# Generic CMake file that should work on many systems, provided that the environmental variables are set properly
# Provides pathways for Intel and GCC compilers
# Make sure to set the following variables: CPLUS_INCLUDE_PATH and LIBRARY_PATH (depending on the system configuriation, they may be called differently)
# The following libraries should be included: netcdf fftw3 fftw3f hdf5
#
# Tested on the following systems:
# ARM Cumulus with GCC CUDA/MPICH (RedHat)
# Cleveland State with GCC CUDA/MPICH (Ubuntu 22.04)

# Switch between Intel and GCC:
set(USEINTEL FALSE)

# GPU builds are always with GCC:
if(USECUDA)
    set(USEINTEL FALSE)
endif()
# Select correct compilers for Intel/GCC + parallel/serial:
if(USEMPI)
    if(USEINTEL)
        set(ENV{CC} mpiicc )
        set(ENV{CXX} mpiicpc)
        set(ENV{FC} mpiifort)
    else()
        set(ENV{CC} mpicc )
        set(ENV{CXX} mpicxx)
        set(ENV{FC} mpif90)
    endif()
else()
    if(USEINTEL)
        set(ENV{CC} icc )
        set(ENV{CXX} icpc)
        set(ENV{FC} ifort)
    else()
        set(ENV{CC} mpicc)
        set(ENV{CXX} mpic++)
        set(ENV{FC} mpif90)
    endif()
endif()

# Set library paths using environment variable USER
if(USECUDA)
    set(BOOST_INCLUDE_DIR  "/user-environment/env/default/include")
    set(NETCDF_INCLUDE_DIR "/capstor/scratch/cscs/$ENV{USER}/microHH-env/view/include")
    set(INCLUDE_DIRS ${BOOST_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})
    set(FFTW_LIB "/user-environment/env/default/lib/libfftw3.so")
    set(FFTWF_LIB "/user-environment/env/default/lib/libfftw3f.so")
    set(HDF5_LIB "/user-environment/env/default/lib/libhdf5.so")
    set(NETCDF_LIB_C "/capstor/scratch/cscs/$ENV{USER}/microHH-env/view/lib64/libnetcdf.so")
    set(LIBS ${FFTW_LIB} ${FFTWF_LIB} ${NETCDF_LIB_C} ${HDF5_LIB}) #It may be necessary to add m z curl sz if necessary
else()
    set(NETCDF_LIB_C "netcdf")
    set(FFTW_LIB "fftw3")
    set(FFTWF_LIB "fftw3f")
    set(HDF5_LIB "hdf5")
    set(LIBS -L/user-environment/env/default/lib -L/user-environment/env/default/lib64 ${FFTW_LIB} ${FFTWF_LIB} ${NETCDF_LIB_C} ${HDF5_LIB}) #It may be necessary to add m z curl sz if necessary
endif()
    
# Set compiler flags / options:
if(USECUDA)
    set(USER_CXX_FLAGS "-std=c++17 -fopenmp")
    set(USER_CXX_FLAGS_RELEASE "-O3")
    add_definitions(-DRESTRICTKEYWORD=__restrict__)
else()
    if(USEINTEL)
        set(USER_CXX_FLAGS "-std=c++17 -restrict")
        set(USER_CXX_FLAGS_RELEASE "-O3 -march=native")
        set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")
        add_definitions(-DRESTRICTKEYWORD=restrict)
    else()
        set(USER_CXX_FLAGS "-std=c++17 -I/user-environment/env/default/include")
        set(USER_CXX_FLAGS_RELEASE "-O3 -march=native")
        set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

        set(USER_FC_FLAGS "-fdefault-real-8 -fdefault-double-8 -fPIC -ffixed-line-length-none -fno-range-check")
        set(USER_FC_FLAGS_RELEASE "-DNDEBUG -O3 -march=native")

        add_definitions(-DRESTRICTKEYWORD=__restrict__)
    endif()
endif()


if(USECUDA)
    set(CMAKE_CUDA_ARCHITECTURES "70;80;90")
    set(USER_CUDA_NVCC_FLAGS "--expt-relaxed-constexpr -lineinfo")
    set(USER_CUDA_NVCC_FLAGS_RELEASE "-DNDEBUG")
    set(USER_CUDA_NVCC_FLAGS_DEBUG "-O0 -g -DCUDACHECKS")
    add_definitions(-DRTE_RRTMGP_GPU_MEMPOOL_CUDA)
endif()

# Disable MPI-IO for cross-sections on GPFS file systems. This may or may not be necessary, depending on the system
add_definitions(-DDISABLE_2D_MPIIO=1)
add_definitions(-DRTE_USE_CBOOL)
