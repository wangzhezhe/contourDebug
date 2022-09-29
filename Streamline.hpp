#ifndef VTK_H_STREAMLINE_HPP
#define VTK_H_STREAMLINE_HPP

#include "Filter.hpp"
#include "DataSet.hpp"

#include <vtkm/Particle.h>

namespace vtkh
{

class Streamline : public Filter
{
public:
  Streamline();
  virtual ~Streamline();
  std::string GetName() const override { return "vtkh::Streamline";}
  void SetField(const std::string &field_name) {  m_field_name = field_name; }
  void SetStepSize(const double &step_size) {   m_step_size = step_size; }
  void SetSeeds(const std::vector<vtkm::Particle>& seeds) { m_seeds = seeds; }
  void SetNumberOfSteps(int numSteps) { m_num_steps = numSteps; }

protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;

  std::string m_field_name;
  double m_step_size;
  int m_num_steps;
  std::vector<vtkm::Particle> m_seeds;
};

} //namespace vtkh
#endif
