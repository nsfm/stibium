# Instruct CMake to run moc, uic, and rrc automatically when needed.
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

# .ui headers are wrapped explicitly (qt6_wrap_ui) rather than via
# AUTOUIC: with Ninja, AUTOUIC's timestamp handling ran one build
# behind on .ui edits, leaving dependents compiled against a stale
# Ui:: layout (shifted member offsets = startup segfaults). Explicit
# wrapping makes the generated headers first-class build nodes.

# Make the build directory an include target (for generated moc_ files)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Find the Qt libraries
find_package(Qt6 REQUIRED COMPONENTS
    Core Gui Widgets OpenGL OpenGLWidgets Concurrent)

set(SRC_FILES
    app/main.cpp
    app/app.cpp
    app/colors.cpp
    app/theme.cpp
    app/settings.cpp
    canvas/scene.cpp
    canvas/canvas_view.cpp
    canvas/add_node_popup.cpp
    canvas/info.cpp
    canvas/inspector/frame.cpp
    canvas/inspector/export.cpp
    canvas/inspector/title.cpp
    canvas/inspector/buttons.cpp
    canvas/inspector/util.cpp
    canvas/datum_row.cpp
    canvas/datum_editor.cpp
    canvas/datum_port.cpp
    canvas/subdatum/subdatum_frame.cpp
    canvas/subdatum/subdatum_editor.cpp
    canvas/subdatum/subdatum_row.cpp
    canvas/connection/base.cpp
    canvas/connection/connection.cpp
    canvas/connection/dummy.cpp
    window/base.cpp
    window/canvas.cpp
    window/viewport.cpp
    window/quad.cpp
    window/base_viewport_window.cpp
    window/script_window.cpp
    graph/constructor/populate.cpp
    graph/proxy/graph.cpp
    graph/proxy/node.cpp
    graph/proxy/script.cpp
    graph/proxy/datum.cpp
    graph/proxy/subdatum.cpp
    graph/proxy/base_datum.cpp
    graph/hooks/hooks.cpp
    graph/hooks/export.cpp
    graph/hooks/ui.cpp
    graph/hooks/title.cpp
    graph/serialize/serializer.cpp
    graph/serialize/deserializer.cpp
    script/syntax.cpp
    script/editor.cpp
    script/frame.cpp
    undo/undo_command.cpp
    undo/undo_add_node.cpp
    undo/undo_add_datum.cpp
    undo/undo_delete_datum.cpp
    undo/undo_delete_multi.cpp
    undo/undo_delete_node.cpp
    undo/undo_delete_link.cpp
    undo/undo_move_node.cpp
    undo/undo_move_datum.cpp
    undo/undo_change_script.cpp
    undo/undo_change_expr.cpp
    undo/undo_stack.cpp
    dialog/exporting.cpp
    dialog/resolution.cpp
    dialog/render_image_dialog.cpp
    export/export_mesh.cpp
    export/export_heightmap.cpp
    export/export_svg.cpp
    export/export_image.cpp
    export/export_worker.cpp
    viewport/scene.cpp
    viewport/gl.cpp
    viewport/view.cpp
    viewport/image.cpp
    viewport/control/control.cpp
    viewport/control/control_instance.cpp
    viewport/control/point.cpp
    viewport/control/wireframe.cpp
    viewport/render/instance.cpp
    viewport/render/task.cpp

    vendor/meshoptimizer/allocator.cpp
    vendor/meshoptimizer/indexgenerator.cpp
    vendor/meshoptimizer/simplifier.cpp

    gl/gl.qrc
)

qt6_wrap_ui(UI_HEADERS
    window/base_window.ui
    dialog/resolution_dialog.ui
    dialog/exporting_dialog.ui
)
list(APPEND SRC_FILES ${UI_HEADERS})

set(PY_FAB_FILES
    ../py/fab/__init__.py
    ../py/fab/shapes.py
    ../py/fab/types.py
)

if (APPLE)
    set(STIBIUM_APP Stibium)
    set(CMAKE_MACOSX_RPATH ON)

    # Add the icon file and deploy it to Resources
    set(SB_ICON ../deploy/mac/sb.icns)
    set_source_files_properties(${SB_ICON}
        PROPERTIES MACOSX_PACKAGE_LOCATION "Resources"
    )

    # Deploy fab module's files into Resources/fab
    set_source_files_properties(${PY_FAB_FILES}
        PROPERTIES MACOSX_PACKAGE_LOCATION "Resources/fab"
    )

    add_executable(${STIBIUM_APP} MACOSX_BUNDLE
        ${SRC_FILES}
        ${SB_ICON}
        ${PY_FAB_FILES}
    )

    # Set application Info.plist file
    set_target_properties(${STIBIUM_APP} PROPERTIES
        MACOSX_BUNDLE_INFO_PLIST ../deploy/mac/Info.plist
    )

    # Copy node directory into application bundle
    add_custom_command(TARGET ${STIBIUM_APP} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_SOURCE_DIR}/py/nodes
            ${CMAKE_CURRENT_BINARY_DIR}/${STIBIUM_APP}.app/Contents/Resources/nodes)
else()
    set(STIBIUM_APP stibium)
    add_executable(${STIBIUM_APP} ${SRC_FILES})

    # Copy node and fab directories into the build directory whenever
    # any py source changes (a POST_BUILD copy only fires on relink,
    # leaving stale python after py-only edits). CONFIGURE_DEPENDS
    # re-globs when files are added or removed.
    file(GLOB_RECURSE SB_PY_SOURCES CONFIGURE_DEPENDS
         ${CMAKE_SOURCE_DIR}/py/*)
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/sb_py.stamp
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_SOURCE_DIR}/py/nodes
            ${CMAKE_CURRENT_BINARY_DIR}/sb/nodes
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_SOURCE_DIR}/py/fab
            ${CMAKE_CURRENT_BINARY_DIR}/sb/fab
        COMMAND ${CMAKE_COMMAND} -E touch
            ${CMAKE_CURRENT_BINARY_DIR}/sb_py.stamp
        DEPENDS ${SB_PY_SOURCES})
    add_custom_target(sb_py_assets ALL
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/sb_py.stamp)
    add_dependencies(${STIBIUM_APP} sb_py_assets)

    install(TARGETS stibium
            DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
    install(DIRECTORY ../py/nodes
            DESTINATION ${CMAKE_INSTALL_PREFIX}/share/stibium)
    install(DIRECTORY ../py/fab
            DESTINATION ${CMAKE_INSTALL_PREFIX}/share/stibium)
endif()

target_link_libraries(${STIBIUM_APP}
    Qt6::Widgets
    Qt6::Gui
    Qt6::OpenGL
    Qt6::OpenGLWidgets
    Qt6::Concurrent
    SbGraph
    SbFab)

################################################################################

execute_process(COMMAND git log --pretty=format:'%h' -n 1
                OUTPUT_VARIABLE GITREV)
execute_process(COMMAND bash -c "git diff --quiet --exit-code || echo +"
                OUTPUT_VARIABLE GITDIFF)
execute_process(COMMAND git describe --exact-match --tags
                OUTPUT_VARIABLE GITTAG
                ERROR_QUIET)
execute_process(COMMAND git rev-parse --abbrev-ref HEAD
                OUTPUT_VARIABLE GITBRANCH)

add_definitions(-D'GITREV="${GITREV}${GITDIFF}"'
                -D'GITTAG="${GITTAG}"'
                -D'GITBRANCH="${GITBRANCH}"')

################################################################################

target_include_directories(${STIBIUM_APP} SYSTEM PRIVATE
    ${Python_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${AUTOGEN_TARGETS_FOLDER}
)

################################################################################

################################################################################
# App-level integration tests: every example must load, evaluate,
# render, and round-trip through the serializer deterministically.
if (NOT APPLE)
    foreach(example csg gear ring hi label showcase_gear)
        add_test(NAME app_validate_${example}
                 COMMAND ${STIBIUM_APP} --validate
                         ${CMAKE_SOURCE_DIR}/examples/${example}.sb)
        add_test(NAME app_render_${example}
                 COMMAND ${STIBIUM_APP}
                         --render ${CMAKE_CURRENT_BINARY_DIR}/ctest_${example}.png
                         --size 128
                         ${CMAKE_SOURCE_DIR}/examples/${example}.sb)
        add_test(NAME app_resave_roundtrip_${example}
                 COMMAND bash -c "\
$<TARGET_FILE:${STIBIUM_APP}> --resave ${CMAKE_CURRENT_BINARY_DIR}/ctest_rt1_${example}.sb ${CMAKE_SOURCE_DIR}/examples/${example}.sb && \
$<TARGET_FILE:${STIBIUM_APP}> --resave ${CMAKE_CURRENT_BINARY_DIR}/ctest_rt2_${example}.sb ${CMAKE_CURRENT_BINARY_DIR}/ctest_rt1_${example}.sb && \
$<TARGET_FILE:${STIBIUM_APP}> --validate ${CMAKE_CURRENT_BINARY_DIR}/ctest_rt1_${example}.sb && \
cmp ${CMAKE_CURRENT_BINARY_DIR}/ctest_rt1_${example}.sb ${CMAKE_CURRENT_BINARY_DIR}/ctest_rt2_${example}.sb")
        set_tests_properties(
                app_validate_${example}
                app_render_${example}
                app_resave_roundtrip_${example}
                PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen")
    endforeach()

    # Mesh import example: its relative asset path resolves against
    # the project file, so the roundtrip copy needs meshes/ beside it
    add_test(NAME app_validate_import_bead
             COMMAND ${STIBIUM_APP} --validate
                     ${CMAKE_SOURCE_DIR}/examples/import_bead.sb)
    add_test(NAME app_render_import_bead
             COMMAND ${STIBIUM_APP}
                     --render ${CMAKE_CURRENT_BINARY_DIR}/ctest_import_bead.png
                     --size 128
                     ${CMAKE_SOURCE_DIR}/examples/import_bead.sb)
    add_test(NAME app_resave_roundtrip_import_bead
             COMMAND bash -c "\
mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/meshes && \
cp ${CMAKE_SOURCE_DIR}/examples/meshes/pebble.stl ${CMAKE_CURRENT_BINARY_DIR}/meshes/ && \
$<TARGET_FILE:${STIBIUM_APP}> --resave ${CMAKE_CURRENT_BINARY_DIR}/ctest_rt1_import_bead.sb ${CMAKE_SOURCE_DIR}/examples/import_bead.sb && \
$<TARGET_FILE:${STIBIUM_APP}> --resave ${CMAKE_CURRENT_BINARY_DIR}/ctest_rt2_import_bead.sb ${CMAKE_CURRENT_BINARY_DIR}/ctest_rt1_import_bead.sb && \
cmp ${CMAKE_CURRENT_BINARY_DIR}/ctest_rt1_import_bead.sb ${CMAKE_CURRENT_BINARY_DIR}/ctest_rt2_import_bead.sb")
    set_tests_properties(
            app_validate_import_bead
            app_render_import_bead
            app_resave_roundtrip_import_bead
            PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen")
endif()

set_property(TARGET ${STIBIUM_APP} PROPERTY CXX_STANDARD 17)
set_property(TARGET ${STIBIUM_APP} PROPERTY C_STANDARD 99)

