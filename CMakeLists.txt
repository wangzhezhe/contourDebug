cmake_minimum_required(VERSION 3.3)

PROJECT(test)

find_package(VTKm REQUIRED)

add_executable(test1 test1.cpp Recenter.cpp Filter.cpp DataSet.cpp vtkm_dataset_info.cpp)
target_link_libraries(test1 vtkm_io vtkm_filter)


add_executable(test2 test2.cpp)
target_link_libraries(test2 vtkm_io vtkm_filter)