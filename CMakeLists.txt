CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(CUDA-PF)

FIND_PACKAGE(CUDA REQUIRED)

INCLUDE(FindCUDA)

INCLUDE_DIRECTORIES(/usr/local/cuda/include)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
FILE(GLOB SOURCES "*.cu" "*.cpp" "*.c" "*.h")
CUDA_ADD_EXECUTABLE(CUDA-PF ${SOURCES})
 
if(VTK_LIBRARIES)
  target_link_libraries(CUDA-PF /usr/local/cuda/lib64/libcudart.so ${VTK_LIBRARIES})
else()
  target_link_libraries(CUDA-PF /usr/local/cuda/lib64/libcudart.so vtkHybrid vtkWidgets)
endif()
