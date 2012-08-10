#Define the project we are working on
project(proteinvis)

# Look for various software we need
cmake_minimum_required(VERSION 2.6)

set(QT_MIN_VERSION "4.5.0")
find_package(Qt4 REQUIRED)


set(QT_USE_QTGUI TRUE)
set(QT_USE_QTOPENGL TRUE)
set(QT_USE_QTXML TRUE)

include(${QT_USE_FILE})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(QGLViewer REQUIRED)

find_package(GLEW REQUIRED)

find_package(Boost COMPONENTS regex REQUIRED)

set(SECONDARY_NUM_SPLINESEGS 20 CACHE STRING
 "number of spline segments for each residue in secondary structure renderings")

subdirs(utls)

file(GLOB proteinvis_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/src/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/include/*.h)

file(GLOB qtcolorpicker_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/qtcolorpicker/src/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/qtcolorpicker/src/*.h
  ${CMAKE_CURRENT_SOURCE_DIR}/qtcolorpicker/include/*.h)

QT4_WRAP_CPP(  proteinvis_MOC_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/include/proteinModel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/include/qtglmaterialeditor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/include/qtgllighteditor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/include/mainwindow.h
  )

QT4_WRAP_CPP(qtcolorpicker_MOC_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/qtcolorpicker/include/qtcolorpicker.h
  ${CMAKE_CURRENT_SOURCE_DIR}/qtcolorpicker/src/qtcolorpicker_widgets.h
  )

file(GLOB proteinvis_UIS ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/forms/*.ui)

QT4_WRAP_UI(proteinvis_UI_SRCS ${proteinvis_UIS})

file(GLOB proteinvis_RC_SRCS 
  ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/resources/*.png
  ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/resources/*.glsl
  ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/resources/*.vert
  ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/resources/*.geom
  ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/resources/*.frag)


QT4_ADD_RESOURCES(proteinvis_RC_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/resources/proteinvis.qrc)

include_directories(
  ${QT_ADDITIONAL_INCLUDE_PATHS}
  ${GLEW_INCLUDE_PATH}
  ${Boost_INCLUDE_DIRS}
  ${utls_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/qtcolorpicker/include
  ${CMAKE_CURRENT_SOURCE_DIR}/qtcolorpicker/src
  ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/include
  ${CMAKE_CURRENT_BINARY_DIR}
  ${QGLVIEWER_INCLUDE_DIR}
  )

configure_file(${PROJECT_SOURCE_DIR}/pv_config.h.in ${PROJECT_BINARY_DIR}/pv_config.h)

add_executable(proteinvis
  ${proteinvis_SRCS}
  ${proteinvis_MOC_SRCS}
  ${proteinvis_UI_SRCS}
  ${proteinvis_RC_SRCS}
  ${qtcolorpicker_SRCS}
  ${qtcolorpicker_MOC_SRCS}

)

target_link_libraries(proteinvis utls ${Boost_LIBRARIES} ${QGLVIEWER_LIBRARY})

