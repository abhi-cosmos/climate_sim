cmake_minimum_required(VERSION 3.20)
project(climate_lab LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 1) Find system-installed libraries
find_package(glfw3 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(NETCDF REQUIRED netcdf)


# 2) Point to the ImGui + ImPlot source folders inside this repo
set(IMGUI_DIR ${CMAKE_SOURCE_DIR}/external/imgui)
set(IMPLOT_DIR ${CMAKE_SOURCE_DIR}/external/implot)

# 3) Build an executable and include all needed .cpp files
add_executable(climate_lab
    src/main.cpp

    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_demo.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl2.cpp

    ${IMPLOT_DIR}/implot.cpp
    ${IMPLOT_DIR}/implot_demo.cpp
    ${IMPLOT_DIR}/implot_items.cpp
)

# 4) Tell the compiler where to find header files (.h)
target_include_directories(climate_lab PRIVATE
    ${IMGUI_DIR}
    ${IMGUI_DIR}/backends
    ${IMPLOT_DIR}
    ${NETCDF_INCLUDE_DIRS}
)

target_link_directories(climate_lab PRIVATE ${NETCDF_LIBRARY_DIRS})


# 5) Link the final executable with GLFW + OpenGL
target_link_libraries(climate_lab PRIVATE
    glfw
    OpenGL::GL
    ${NETCDF_LIBRARIES}
)
