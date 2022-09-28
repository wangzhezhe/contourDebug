#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <vtkm/filter/Contour.h>
#include <vtkm/filter/FieldSelection.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include "Recenter.hpp"

int main(int argc, char **argv) {

  // load the vtk data
  // uniform rectilinear grid
  std::string fileName = "../dataset/test.vtk";

  std::vector<double> iso_values = {300};
  //std::vector<double> iso_values = {167.5};
  std::string fieldToOperateOn = "density";

  vtkm::io::VTKDataSetReader reader(fileName);
  vtkm::cont::DataSet ds = reader.ReadDataSet();

  vtkm::filter::FieldSelection sel;
  sel.AddField(fieldToOperateOn);

  vtkm::filter::contour::Contour marcher;

  marcher.SetFieldsToPass(sel);
  marcher.SetIsoValues(iso_values);
  marcher.SetMergeDuplicatePoints(false);
  marcher.SetActiveField(fieldToOperateOn);

  auto dataset = marcher.Execute(ds);


   vtkm::io::VTKDataSetWriter writer("../dataset/contour_output_test.vtk");
    writer.WriteDataSet(dataset);

  return 0;
}
