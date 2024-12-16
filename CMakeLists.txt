cmake_minimum_required(VERSION 3.22)
project(NART VERSION 0.1.0 LANGUAGES CXX)

include(FetchContent)

# Fetch Zlib
FetchContent_Declare(
  zlib
  GIT_REPOSITORY https://github.com/madler/zlib.git
  GIT_TAG v1.2.13
)
FetchContent_MakeAvailable(zlib)

# We know zlib and its targets exist now.
if(TARGET example)
  target_include_directories(example PRIVATE
    ${zlib_SOURCE_DIR}
    ${zlib_BINARY_DIR}
  )
endif()

if(TARGET minigzip)
  target_include_directories(minigzip PRIVATE
    ${zlib_SOURCE_DIR}
    ${zlib_BINARY_DIR}
  )
endif()

# If 64-bit variants are built:
if(TARGET example64)
  target_include_directories(example64 PRIVATE
    ${zlib_SOURCE_DIR}
    ${zlib_BINARY_DIR}
  )
endif()

if(TARGET minigzip64)
  target_include_directories(minigzip64 PRIVATE
    ${zlib_SOURCE_DIR}
    ${zlib_BINARY_DIR}
  )
endif()

# Create an alias for consistency
add_library(ZLIB::ZLIB ALIAS zlib)

# Mark Zlib as found and point to correct include dir
set(ZLIB_FOUND TRUE)
set(ZLIB_INCLUDE_DIR "${zlib_BINARY_DIR}")
# You can try leaving ZLIB_LIBRARIES unset since we have a target. If needed:
set(ZLIB_LIBRARIES ZLIB::ZLIB)

include_directories(${zlib_SOURCE_DIR} ${zlib_BINARY_DIR})

# Fetch Imath
FetchContent_Declare(
    imath
    GIT_REPOSITORY https://github.com/AcademySoftwareFoundation/Imath.git
    GIT_TAG v3.1.12
)
FetchContent_MakeAvailable(imath)

# Configure OpenEXR to use system zlib
set(OPENEXR_USE_SYSTEM_ZLIB ON CACHE BOOL "Use system zlib" FORCE)
set(OPENEXR_INSTALL OFF CACHE BOOL "Disable OpenEXR installation" FORCE)
set(OPENEXR_BUILD_UTILS ON CACHE BOOL "Enable OpenEXR utilities" FORCE)
set(OPENEXR_ENABLE_TESTS OFF CACHE BOOL "No tests" FORCE)
set(OPENEXR_INSTALL_TOOLS OFF CACHE BOOL "No tools" FORCE)
set(OPENEXR_BUILD_EXAMPLES OFF CACHE BOOL "No examples" FORCE)
set(OPENEXR_BUILD_STATIC OFF CACHE BOOL "Build shared libraries" FORCE)
set(OPENEXR_BUILD_RGBA ON CACHE BOOL "Enable legacy Rgba API" FORCE)
set(OPENEXR_BUILD_LEGACY ON CACHE BOOL "Enable legacy Rgba API" FORCE)

# Fetch OpenEXR
message(STATUS "Fetching OpenEXR via FetchContent")
FetchContent_Declare(
    openexr
    GIT_REPOSITORY https://github.com/AcademySoftwareFoundation/openexr.git
    GIT_TAG v3.2.3
)
FetchContent_MakeAvailable(openexr)

include_directories(
	${zlib_SOURCE_DIR}
    ${openexr_SOURCE_DIR}/include
    ${openexr_SOURCE_DIR}/src/lib
    ${openexr_BINARY_DIR}/src/lib/OpenEXR
    ${openexr_BINARY_DIR}/src/lib/Half
    ${imath_SOURCE_DIR}/include
)

# Fetch GLM
FetchContent_Declare(
    glm
    GIT_REPOSITORY https://github.com/g-truc/glm.git
    GIT_TAG 0.9.9.8
)
FetchContent_MakeAvailable(glm)

# Fetch nlohmann_json
FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.3
)
FetchContent_MakeAvailable(nlohmann_json)

# Fetch TBB
FetchContent_Declare(
    tbb
    GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
    GIT_TAG v2021.9.0
)
set(TBB_BUILD_SHARED ON CACHE BOOL "Build TBB as shared library" FORCE)
set(TBB_TEST OFF CACHE BOOL "Disable TBB tests" FORCE)
set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL "Do not build TBB allocator" FORCE)
FetchContent_MakeAvailable(tbb)

set(SOURCES
    NART/src/core/main.cpp
    NART/src/core/bvh.cpp
    NART/src/core/geometry.cpp
    NART/src/core/reflection.cpp
    NART/src/core/render.cpp
    NART/src/core/sampling.cpp
    NART/src/core/scene.cpp
    NART/src/cameras/pinholecamera.cpp
    NART/src/lights/disklight.cpp
    NART/src/lights/distantlight.cpp
    NART/src/lights/ringlight.cpp
    NART/src/materials/diffusematerial.cpp
    NART/src/materials/glassmaterial.cpp
    NART/src/materials/glossydielectricmaterial.cpp
    NART/src/materials/plasticmaterial.cpp
    NART/src/materials/specularmaterial.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
    ${glm_SOURCE_DIR}
    ${openexr_SOURCE_DIR}/include
)

target_link_libraries(${PROJECT_NAME} PRIVATE
    Imath::Imath
    OpenEXR::OpenEXRUtil
    OpenEXR::OpenEXR
    nlohmann_json::nlohmann_json
    TBB::tbb
    ZLIB::ZLIB
)


