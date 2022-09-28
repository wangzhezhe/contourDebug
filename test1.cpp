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
  std::string filePrefix = "../dataset/fileReader_output_step000_partition00";

  vtkh::DataSet data_set;
  data_set.SetCycle(6);

  for (int i = 0; i < 4; i++) {
    std::string fileName = filePrefix + std::to_string(i) + ".vtk";
    std::cout << "process " << fileName << std::endl;
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet ds = reader.ReadDataSet();
    data_set.AddDomain(ds, i);
    // ds.PrintSummary(std::cout);
  }

  std::vector<double> iso_values = {167.5, 334, 500.5, 667, 833.5};
  //std::vector<double> iso_values = {167.5};
  std::string fieldToOperateOn = "density";

  vtkm::filter::FieldSelection sel;
  sel.AddField(fieldToOperateOn);

  // recenter
  vtkh::Recenter recenter;

  recenter.SetInput(&data_set);
  recenter.SetField(fieldToOperateOn);
  recenter.SetResultAssoc(vtkm::cont::Field::Association::Points);
  recenter.Update();

  vtkh::DataSet *m_input = recenter.GetOutput();

  vtkh::DataSet temp_data;
  const int num_domains = m_input->GetNumberOfDomains();
  for (int i = 0; i < num_domains; ++i) {

    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    m_input->GetDomain(i, dom, domain_id);

    //vtkm::io::VTKDataSetWriter writerA(filePrefix + "recenter_output_" +
    //                                  std::to_string(i) + ".vtk");
    //writerA.WriteDataSet(dom);

    //if (!dom.HasField(fieldToOperateOn)) {
    //  std::cout << "error, no specified field" << std::endl;
    //}
    
    vtkm::filter::contour::Contour marcher;

    marcher.SetFieldsToPass(sel);
    marcher.SetIsoValues(iso_values);
    marcher.SetMergeDuplicatePoints(false);
    marcher.SetActiveField(fieldToOperateOn);

    auto dataset = marcher.Execute(dom);

    // output dataset for checking
    //vtkm::io::VTKDataSetWriter writer(filePrefix + "contour_output_" +
    //                                  std::to_string(i) + ".vtk");
    //writer.WriteDataSet(dataset);

    temp_data.AddDomain(dataset, domain_id);
    
  }

  return 0;
}
