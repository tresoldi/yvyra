cmake_minimum_required(VERSION 3.16)
project(yvyra VERSION 0.2.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Compiler flags
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-g)
else()
    add_compile_options(-O3 -DNDEBUG)
endif()
add_compile_options(-Wall -pedantic)

# --- Optional: SSE ---
option(ENABLE_SSE "Use SSE SIMD extensions if available" ON)
if(ENABLE_SSE)
    include(CheckCCompilerFlag)
    check_c_compiler_flag(-msse HAS_SSE_FLAG)
    include(CheckCSourceCompiles)
    set(CMAKE_REQUIRED_FLAGS "-msse")
    check_c_source_compiles("
        #include <xmmintrin.h>
        int main() { __m128 v = _mm_setzero_ps(); (void)v; return 0; }
    " HAVE_SSE)
    unset(CMAKE_REQUIRED_FLAGS)
    if(HAVE_SSE)
        add_compile_options(-msse)
    else()
        set(DISABLE_SSE 1)
    endif()
else()
    set(DISABLE_SSE 1)
endif()

# --- Optional: AVX ---
option(ENABLE_AVX "Use AVX SIMD extensions if available" ON)
if(ENABLE_AVX)
    include(CheckCSourceCompiles)
    set(CMAKE_REQUIRED_FLAGS "-mavx")
    check_c_source_compiles("
        #include <immintrin.h>
        int main() { __m256 v = _mm256_setzero_ps(); (void)v; return 0; }
    " HAVE_AVX)
    unset(CMAKE_REQUIRED_FLAGS)
    if(HAVE_AVX)
        add_compile_options(-mavx)
    else()
        set(DISABLE_AVX 1)
    endif()
else()
    set(DISABLE_AVX 1)
endif()

# --- Optional: FMA ---
option(ENABLE_FMA "Use FMA SIMD extensions if available" ON)
if(ENABLE_FMA)
    include(CheckCSourceCompiles)
    set(CMAKE_REQUIRED_FLAGS "-mfma")
    check_c_source_compiles("
        #include <immintrin.h>
        int main() { __m256 a = _mm256_setzero_ps(); __m256 r = _mm256_fmadd_ps(a,a,a); (void)r; return 0; }
    " HAVE_FMA3)
    unset(CMAKE_REQUIRED_FLAGS)
    if(HAVE_FMA3)
        add_compile_options(-mfma)
    else()
        set(DISABLE_FMA 1)
    endif()
else()
    set(DISABLE_FMA 1)
endif()

# --- Optional: MPI ---
option(ENABLE_MPI "Build with MPI support" OFF)
if(ENABLE_MPI)
    find_package(MPI REQUIRED)
    set(MPI_ENABLED 1)
endif()

# --- Optional: readline/libedit ---
if(NOT ENABLE_MPI)
    find_library(READLINE_LIB NAMES edit readline)
    if(READLINE_LIB)
        set(HAVE_LIBREADLINE 1)
        # Check for readline/readline.h vs readline.h
        include(CheckIncludeFile)
        check_include_file("readline/readline.h" HAVE_READLINE_READLINE_H)
        check_include_file("readline.h" HAVE_READLINE_H)
        check_include_file("readline/history.h" HAVE_READLINE_HISTORY_H)
        check_include_file("history.h" HAVE_HISTORY_H)
        # Check for add_history
        include(CheckSymbolExists)
        set(CMAKE_REQUIRED_LIBRARIES ${READLINE_LIB})
        check_symbol_exists(add_history "readline/history.h" HAVE_READLINE_HISTORY)
        unset(CMAKE_REQUIRED_LIBRARIES)
    endif()
endif()

# --- Math library ---
find_library(MATH_LIB m)

# --- Host info ---
cmake_host_system_information(RESULT HOST_CPU QUERY OS_PLATFORM)
set(HOST_TYPE "${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_NAME}")

# --- Compiler info ---
set(COMPILER_VENDOR "${CMAKE_C_COMPILER_ID}")
set(COMPILER_VERSION "${CMAKE_C_COMPILER_VERSION}")

# --- Standard headers (always present on modern systems, but keep for compat) ---
include(CheckIncludeFile)
check_include_file("float.h" HAVE_FLOAT_H)
check_include_file("limits.h" HAVE_LIMITS_H)
check_include_file("unistd.h" HAVE_UNISTD_H)

# --- Generate config.h ---
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/config.h.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/src/config.h
    @ONLY
)

add_subdirectory(src)

# --- Testing ---
enable_testing()

find_package(Python3 COMPONENTS Interpreter)

add_test(
    NAME smoke_test
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/testing/smoke_test.sh $<TARGET_FILE:yvyra>
)

if(Python3_FOUND)
    add_test(
        NAME mk_tests
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/testing/run_tests.py
            --binary $<TARGET_FILE:yvyra>
    )
    set_tests_properties(mk_tests PROPERTIES TIMEOUT 300)
endif()
