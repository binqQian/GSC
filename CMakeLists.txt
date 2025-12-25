cmake_minimum_required(VERSION 3.16)
project(GSC LANGUAGES CXX)

# C++ standard and common flags
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CheckCXXCompilerFlag)
set(COMMON_COMPILE_FLAGS -Wall -O3 -g -m64 -pthread -mpopcnt)
check_cxx_compiler_flag(-msse4.2 HAS_SSE42)
check_cxx_compiler_flag(-msse2 HAS_SSE2)
if(HAS_SSE42)
    list(APPEND COMMON_COMPILE_FLAGS -msse4.2)
elseif(HAS_SSE2)
    list(APPEND COMMON_COMPILE_FLAGS -msse2)
endif()
check_cxx_compiler_flag(-mavx HAS_AVX)
if(HAS_AVX)
    list(APPEND COMMON_COMPILE_FLAGS -mavx)
endif()

# Allow using an existing vcpkg install without a toolchain.
if(NOT DEFINED VCPKG_ROOT AND DEFINED ENV{VCPKG_ROOT})
    set(VCPKG_ROOT "$ENV{VCPKG_ROOT}")
endif()
if(VCPKG_ROOT AND EXISTS "${VCPKG_ROOT}/installed/x64-linux")
    list(APPEND CMAKE_PREFIX_PATH "${VCPKG_ROOT}/installed/x64-linux")
elseif(EXISTS "/home/binq/vcpkg/installed/x64-linux")
    list(APPEND CMAKE_PREFIX_PATH "/home/binq/vcpkg/installed/x64-linux")
endif()

find_package(Threads REQUIRED)
find_package(spdlog CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(Brotli CONFIG QUIET)
if(NOT Brotli_FOUND)
    find_package(unofficial-brotli CONFIG REQUIRED)
    set(BROTLi_LIBS unofficial::brotli::brotlienc unofficial::brotli::brotlidec unofficial::brotli::brotlicommon)
else()
    set(BROTLi_LIBS Brotli::brotlienc Brotli::brotlidec Brotli::brotlicommon)
endif()
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)

# Local static deps shipped with the repo.
set(LIB_DIR "${CMAKE_SOURCE_DIR}/lib")
foreach(lib_name hts sdsl bsc zstd)
    find_library(LIB_${lib_name} NAMES ${lib_name} lib${lib_name}
        PATHS ${LIB_DIR}
        NO_DEFAULT_PATH)
    if(NOT LIB_${lib_name})
        message(FATAL_ERROR "Missing ${lib_name} in ${LIB_DIR}")
    endif()
endforeach()

# Sources
set(SOURCES
    src/bit_memory.cpp
    src/block_processing.cpp
    src/vint_code.cpp
    src/decompression_reader.cpp
    src/decompressor.cpp
    src/compressor.cpp
    src/main.cpp
    src/samples.cpp
    src/compression_reader.cpp
    src/parallel_vcf_reader.cpp
    src/parallel_vcf_writer.cpp
    src/file_handle.cpp
    src/utils.cpp
    src/bsc.cpp
    src/zstd_compress.cpp
    src/logger.cpp
    src/compression_strategy.cpp
    include/cpp-mmf/memory_mapped_file.cpp
)

add_executable(gsc ${SOURCES})

target_include_directories(gsc
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/cpp-mmf
)

target_compile_options(gsc PRIVATE ${COMMON_COMPILE_FLAGS})

target_link_libraries(gsc
    PRIVATE
        Threads::Threads
        spdlog::spdlog
        fmt::fmt
        ${BROTLi_LIBS}
        ZLIB::ZLIB
        BZip2::BZip2
        LibLZMA::LibLZMA
        ${LIB_hts}
        ${LIB_sdsl}
        ${LIB_bsc}
        ${LIB_zstd}
)
