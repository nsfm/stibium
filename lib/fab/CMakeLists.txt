find_package(PNG REQUIRED)
find_package(ZLIB REQUIRED)
find_package(FLEX REQUIRED)
find_package(Threads REQUIRED)

################################################################################

# Custom command to generate the parser sources with lemon
if(WIN32)
    find_program(LEMON_EXECUTABLE lemon.exe DOC "path to the lemon executable")
else()
    find_program(LEMON_EXECUTABLE lemon DOC "path to the lemon executable")
endif()

mark_as_advanced(LEMON_EXECUTABLE)
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/v2syntax.lemon.hpp
           ${CMAKE_CURRENT_BINARY_DIR}/v2syntax.lemon.cpp
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/tree/v2syntax.y
    COMMAND ${LEMON_EXECUTABLE} -q -c -s ${CMAKE_CURRENT_SOURCE_DIR}/src/tree/v2syntax.y
    COMMAND mv ${CMAKE_CURRENT_SOURCE_DIR}/src/tree/v2syntax.c ${CMAKE_CURRENT_BINARY_DIR}/v2syntax.lemon.cpp
    COMMAND mv ${CMAKE_CURRENT_SOURCE_DIR}/src/tree/v2syntax.h ${CMAKE_CURRENT_BINARY_DIR}/v2syntax.lemon.hpp)

# Custom command to generate the lexer sources with flex
if(WIN32)
    find_program(FLEX_EXECUTABLE flex.exe DOC "path to the flex executable")
else()
    find_program(FLEX_EXECUTABLE flex DOC "path to the flex executable")
endif()
mark_as_advanced(FLEX_EXECUTABLE)
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/v2syntax.yy.hpp
           ${CMAKE_CURRENT_BINARY_DIR}/v2syntax.yy.cpp
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/tree/v2syntax.l
    COMMAND ${FLEX_EXECUTABLE} --outfile=${CMAKE_CURRENT_BINARY_DIR}/v2syntax.yy.cpp
                               --header-file=${CMAKE_CURRENT_BINARY_DIR}/v2syntax.yy.hpp
                               ${CMAKE_CURRENT_SOURCE_DIR}/src/tree/v2syntax.l)

add_library(SbFab STATIC
    src/fab.cpp
    src/formats/dxf.cpp
    src/formats/png.c
    src/mesh/mesh_import.cpp
    src/mesh/mesh_query.cpp
    src/mesh/trimesh.cpp
    src/tree/grid.cpp
    src/formats/stl.c
    src/formats/svg.cpp
    src/formats/threemf.cpp
    src/tree/contour.cpp
    src/tree/analytics.cpp
    src/tree/parallel_eval.cpp
    src/tree/render_mt.cpp
    src/tree/eval.c
    src/tree/math/math_f.c
    src/tree/math/math_g.c
    src/tree/math/math_i.c
    src/tree/math/math_r.c
    src/tree/node/node_c.c
    src/tree/node/opcodes.c
    src/tree/node/printers.c
    src/tree/node/printers_ss.cpp
    src/tree/node/results.c
    src/tree/parser.c
    src/tree/render.c
    src/tree/gpu.cpp
    src/tree/tape.cpp
    src/tree/tree.c
    src/tree/v2parser.cpp
    src/tree/triangulate/delaunay.cpp
    src/tree/triangulate/mesher.cpp
    src/tree/triangulate.cpp
    src/types/bounds.cpp
    src/types/shape.cpp
    src/types/transform.cpp
    src/util/region.c

    # Generated files
    ${CMAKE_CURRENT_BINARY_DIR}/v2syntax.lemon.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/v2syntax.yy.cpp
)

target_include_directories(SbFab SYSTEM PRIVATE
    ${Python_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${PNG_INCLUDE_DIR}
)

# The kernel never reads FP exception flags, and min_f/max_f's
# NaN-aware selects only if-convert (and hence vectorize) when the
# compiler may ignore them.  Computed values are unchanged — renders
# and golden meshes stay bit-identical (verified 2026-07-13).
target_compile_options(SbFab PRIVATE -fno-trapping-math)

# GPU depth renderer (src/tree/gpu.cpp): needs only EGL/GL headers at
# compile time - libEGL is dlopen'd at runtime.  Without the headers
# the stubs compile and gpu_render16 reports unavailable.
find_package(OpenGL COMPONENTS EGL)
if (OpenGL_EGL_FOUND)
    target_compile_definitions(SbFab PRIVATE STIBIUM_HAS_EGL)
    target_link_libraries(SbFab ${CMAKE_DL_LIBS})
endif()

# Adaptive-Delaunay mesher (src/tree/triangulate/delaunay.cpp):
# stages B+C need CGAL's robust incremental Delaunay (header-only,
# exact predicates).  Without it, stage A still builds and
# delaunay_mesh reports unavailable.
# (CGAL's config runs its own find_package(Boost), which clobbers
# the parent's component results - save and restore.)
set(_stibium_boost_libs ${Boost_LIBRARIES})
find_package(CGAL QUIET)
set(Boost_LIBRARIES ${_stibium_boost_libs})
if (CGAL_FOUND)
    target_compile_definitions(SbFab PRIVATE STIBIUM_HAS_CGAL)
    target_link_libraries(SbFab CGAL::CGAL)
endif()

# Add the binary dir to include generated files
target_include_directories(SbFab PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(SbFab
    ${Boost_LIBRARIES}
    ${Python_LIBRARY_RELEASE}
    ${PNG_LIBRARIES}
    ZLIB::ZLIB
    Threads::Threads
)
target_include_directories(SbFab PUBLIC inc)
target_include_directories(SbFab SYSTEM PRIVATE
    vendor)

################################################################################

add_executable(SbFabTest
    tests/analytics.cpp
    tests/contour.cpp
    tests/grid.cpp
    tests/main.cpp
    tests/mesh.cpp
    tests/mesher.cpp
    tests/parser.cpp
    tests/aa.cpp
    tests/dmesh.cpp
    tests/shape.cpp
    tests/tape.cpp
)
target_link_libraries(SbFabTest SbFab)

# GPU tape-interpreter prototype (hidden "[.gpu]" tag): headless EGL
# compute.  Only needs libEGL - GL entry points load at runtime.
find_package(OpenGL COMPONENTS EGL)
if (OpenGL_EGL_FOUND)
    target_sources(SbFabTest PRIVATE tests/gpu.cpp)
    target_link_libraries(SbFabTest OpenGL::EGL)
endif()

target_include_directories(SbFabTest SYSTEM PRIVATE
    ../../vendor
    ${Boost_INCLUDE_DIRS}
    ${Python_INCLUDE_DIRS}
)

################################################################################

set_property(TARGET SbFab PROPERTY CXX_STANDARD 17)
set_property(TARGET SbFab PROPERTY C_STANDARD 99)
set_property(TARGET SbFabTest PROPERTY CXX_STANDARD 17)
set_property(TARGET SbFabTest PROPERTY C_STANDARD 99)

add_test(NAME SbFabTest COMMAND SbFabTest)
