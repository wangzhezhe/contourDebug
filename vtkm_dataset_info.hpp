//-----------------------------------------------------------------------------
///
/// file: vtkm_dataset_info.hpp
///
//-----------------------------------------------------------------------------

#ifndef VTKM_DATASET_INFO_HPP
#define VTKM_DATASET_INFO_HPP

#include <vtkm/cont/DataSet.h>

namespace vtkh {

class VTKMDataSetInfo
{
public:
 typedef typename vtkm::cont::ArrayHandleUniformPointCoordinates UniformArrayHandle;
 typedef vtkm::cont::ArrayHandle<vtkm::FloatDefault> DefaultHandle;
 typedef typename vtkm::cont::ArrayHandleCartesianProduct<DefaultHandle,
                                                          DefaultHandle,
                                                          DefaultHandle> CartesianArrayHandle;
//-----------------------------------------------------------------------------
  static bool IsStructured(const vtkm::cont::DataSet &data_set,
                           int &topo_dims);


  static bool IsStructured(const vtkm::cont::DynamicCellSet &cell_set, int &topo_dims);

  static bool IsRectilinear(const vtkm::cont::DataSet &data_set);


  static bool IsRectilinear(const vtkm::cont::CoordinateSystem &coords);

  static bool IsUniform(const vtkm::cont::DataSet &data_set);


  static bool IsUniform(const vtkm::cont::CoordinateSystem &coords);

  static bool GetPointDims(const vtkm::cont::DataSet &data_set, int *dims);


  static bool GetPointDims(const vtkm::cont::DynamicCellSet &cell_set, int *dims);

  static bool GetCellDims(const vtkm::cont::DataSet &data_set, int *dims);


  static bool GetCellDims(const vtkm::cont::DynamicCellSet &cell_set, int *dims);

  static bool IsSingleCellShape(const vtkm::cont::DynamicCellSet &cell_set, vtkm::UInt8 &shape_id);

};

} // namespace vtkh

#endif
//-----------------------------------------------------------------------------
// -- end header ifdef guard
//-----------------------------------------------------------------------------
