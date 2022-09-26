cmake_minimum_required(VERSION 3.3)

PROJECT(test)

find_package(VTKm REQUIRED)

add_executable(test test.cpp Recenter.cpp Filter.cpp DataSet.cpp vtkm_dataset_info.cpp)
target_link_libraries(test vtkm_io vtkm_filter)
