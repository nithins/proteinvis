#Define the project we are working on
project(proteinvis)

# Look for various software we need 
cmake_minimum_required(VERSION 2.6)

set(QT_MIN_VERSION "4.5.0")
find_package(Qt4 REQUIRED)

if(NOT QT4_FOUND)
  message(FATAL_ERROR "Failed to find Qt 4.5 or greater.")
endif(NOT QT4_FOUND)

set(QT_USE_QTGUI TRUE)
set(QT_USE_QTOPENGL TRUE)

include(${QT_USE_FILE})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(GLEW)

if(NOT GLEW_FOUND)
  message(FATAL_ERROR "Failed to find GLEW .. quitting")
endif(NOT GLEW_FOUND)

find_package(Boost COMPONENTS regex)

if(NOT Boost_FOUND)
  message(FATAL_ERROR "Cannot find boost")
endif(NOT Boost_FOUND)

subdirs(glapp_framework)

file(GLOB proteinvis_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/include/*.h)

QT4_WRAP_CPP(  proteinvis_MOC_SRCS   ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/include/proteinModel.h )

file(GLOB proteinvis_UIS ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/forms/*.ui)

QT4_WRAP_UI(proteinvis_UI_SRCS ${proteinvis_UIS})

QT4_ADD_RESOURCES(proteinvis_RC_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/resources/proteinvis.qrc)

include_directories(
  ${QT_ADDITIONAL_INCLUDE_PATHS}
  ${GLEW_INCLUDE_PATH}
  ${Boost_INCLUDE_DIRS}
  ${glapp_framework_SOURCE_DIR}/utls/include 
  ${glapp_framework_SOURCE_DIR}/qtcolorpicker/include
  ${glapp_framework_SOURCE_DIR}/framework/include
  ${CMAKE_CURRENT_SOURCE_DIR}/pv_app/include
  ${CMAKE_CURRENT_BINARY_DIR}
  )

add_executable(proteinvis ${proteinvis_SRCS} ${proteinvis_MOC_SRCS} ${proteinvis_UI_SRCS} ${proteinvis_RC_SRCS})

target_link_libraries(proteinvis glapp_framework ${Boost_LIBRARIES})

