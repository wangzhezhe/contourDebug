#ifndef VTK_H_RECENTER_HPP
#define VTK_H_RECENTER_HPP

#include "Filter.hpp"
#include "DataSet.hpp"

namespace vtkh
{

class Recenter : public Filter
{
public:
  Recenter();
  virtual ~Recenter();
  std::string GetName() const override;
  void SetField(const std::string &field_name);
  void SetResultAssoc(vtkm::cont::Field::Association assoc);

protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;

  std::string m_field_name;
  vtkm::cont::Field::Association m_assoc;
};

} //namespace vtkh
#endif
