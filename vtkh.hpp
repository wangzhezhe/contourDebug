#ifndef VTK_H_HPP
#define VTK_H_HPP


#include <string>

namespace vtkh
{

   std::string AboutVTKH();
   void        Initialize();

  // is backend support compiled in
   bool        IsSerialAvailable();
   bool        IsOpenMPAvailable();
   bool        IsCUDAAvailable();
   bool        IsKokkosAvailable();

  // is backend enabled (e.g., ForceX)
   bool        IsSerialEnabled();
   bool        IsOpenMPEnabled();
   bool        IsCUDAEnabled();
   bool        IsKokkosEnabled();

   bool        IsMPIEnabled();

   int         CUDADeviceCount();
   void        SelectCUDADevice(int device_index);
   void        SelectKokkosDevice(int device_index);

   void        ForceSerial();
   void        ForceOpenMP();
   void        ForceCUDA();
   void        ForceKokkos();
   void        ResetDevices();
   std::string GetCurrentDevice();

   int         GetMPIRank();
   int         GetMPISize();

   void        SetMPICommHandle(int mpi_comm_id);
   int         GetMPICommHandle();
}
#endif
