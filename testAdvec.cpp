#include <vector>
#include "Streamline.hpp"
#include "Recenter.hpp"
#include "mpi.h"
#include <random>
#include <vtkm/io/VTKDataSetReader.h>
#include "vtkh.hpp"
static vtkm::FloatDefault random01()
{
    return (vtkm::FloatDefault)rand()/(vtkm::FloatDefault)RAND_MAX;
}

void createBoxOfSeeds(std::vector<vtkm::Particle> &seeds,
                      vtkm::FloatDefault xMin, vtkm::FloatDefault xMax,
                      vtkm::FloatDefault yMin, vtkm::FloatDefault yMax,
                      vtkm::FloatDefault zMin, vtkm::FloatDefault zMax,
                      int numSeeds, int rank, int numRanks) {
  vtkm::Particle diff, startPoint, endPoint;

  startPoint.Pos = {xMin, yMin, zMin};
  endPoint.Pos = {xMax, yMax, zMax};
  diff.Pos = endPoint.Pos - startPoint.Pos;

  for (int i = 0; i < numSeeds; i++) {
    vtkm::Particle p;
    p.Pos = {startPoint.Pos[0] + (diff.Pos[0] * random01()),
             startPoint.Pos[1] + (diff.Pos[1] * random01()),
             startPoint.Pos[2] + (diff.Pos[2] * random01())};
    p.ID = static_cast<vtkm::Id>(i);
    seeds.push_back(p);
  }

  // Make the paricle ID's unique
  std::vector<int> particlesPerRank(numRanks, 0);
  particlesPerRank[rank] = numSeeds;
  MPI_Allreduce(MPI_IN_PLACE, particlesPerRank.data(), numRanks, MPI_INT,
                MPI_SUM, MPI_COMM_WORLD);

  int offset = 0;
  for (int i = 0; i < rank; i++)
    offset += particlesPerRank[i];

  if (offset > 0) {
    for (auto &p : seeds)
      p.ID += offset;
  }
}

void runAdvection(vtkh::DataSet *data_set, int rank, int numRanks) {

  vtkh::Streamline streamline;

  std::vector<vtkm::Particle> seeds;
  // createLineOfSeeds(seeds, startPoint, endPoint, GLOBAL_ADVECT_NUM_SEEDS,
  // rank);
  int GLOBAL_ADVECT_NUM_SEEDS = 2048;
 
  createBoxOfSeeds(seeds, 0, 0.28, 0, 0.055, 0, 0.085,
                     GLOBAL_ADVECT_NUM_SEEDS, rank, numRanks);

  std::string fieldToOperateOn = "density";
  int GLOBAL_ADVECT_STEP_SIZE = 0.01;
  int GLOBAL_ADVECT_NUM_STEPS = 512;
  streamline.SetSeeds(seeds);
  streamline.SetInput(data_set);
  streamline.SetField(fieldToOperateOn);
  streamline.SetStepSize(GLOBAL_ADVECT_STEP_SIZE);
  streamline.SetNumberOfSteps(GLOBAL_ADVECT_NUM_STEPS);

  if (rank == 0) {
    std::cerr << "\nstep size " << GLOBAL_ADVECT_STEP_SIZE << std::endl;
    std::cerr << "maxSteps " << GLOBAL_ADVECT_NUM_STEPS << std::endl;
  }

  vtkh::DataSet *streamline_output = NULL;

  streamline.Update();

  streamline_output = streamline.GetOutput();
}

int main(int argc, char **argv) {

  int rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS) {
    printf("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  int numTasks,mpiRank;
  MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
 
  if (numTasks != 4) {
    std::cout << "there are only 4 blocks, use mpi 4" << std::endl;
    exit(0);
  }

    vtkh::SetMPICommHandle(MPI_Comm_c2f(MPI_COMM_WORLD));

  // load data
  vtkh::DataSet data_set;
  data_set.SetCycle(6);

  std::string filePrefix = "../dataset/fileReader_output_step000_partition00";

  std::string fileName = filePrefix + std::to_string(mpiRank) + ".vtk";
  std::cout << "process " << fileName << std::endl;
  vtkm::io::VTKDataSetReader vtkreader(fileName);
  vtkm::cont::DataSet ds = vtkreader.ReadDataSet();
  data_set.AddDomain(ds, 0);
  // ds.PrintSummary(std::cout);

  // call streamline
  for (vtkm::Id i = 0; i < 1; i++) {
    if (data_set.GetDomain(i).HasField("topo_ghosts")) {
      auto temp = data_set.GetDomain(i).GetField("topo_ghosts");

      if (temp.GetNumberOfValues() >= 1) {
        auto ghostArr =
            temp.GetData()
                .AsArrayHandle<
                    vtkm::cont::ArrayHandleBasic<vtkm::FloatDefault>>();
        const vtkm::FloatDefault *buff = ghostArr.GetReadPointer();
        vtkm::cont::ArrayHandle<vtkm::UInt8> ghosts;
        ghosts.Allocate(temp.GetNumberOfValues());
        for (vtkm::Id z = 0; z < temp.GetNumberOfValues(); z++) {
          ghosts.WritePortal().Set(z, static_cast<vtkm::UInt8>(buff[z]));
        }
        data_set.GetDomain(i).AddCellField("vtkmGhostCells", ghosts);
      }
    } else {
      std::cout << "topo_ghosts does not exist in data_set" << std::endl;
      exit(0);
    }
  }

  // create velocity field from x y z velocity
  for (int currentPartition = 0;
       currentPartition < data_set.GetNumberOfDomains(); currentPartition++) {
    auto xarrayV =
        data_set.GetDomain(currentPartition)
            .GetField("velocityx")
            .GetData()
            .AsArrayHandle<vtkm::cont::ArrayHandle<vtkm::FloatDefault>>();
    auto yarrayV =
        data_set.GetDomain(currentPartition)
            .GetField("velocityy")
            .GetData()
            .AsArrayHandle<vtkm::cont::ArrayHandle<vtkm::FloatDefault>>();
    auto zarrayV =
        data_set.GetDomain(currentPartition)
            .GetField("velocityz")
            .GetData()
            .AsArrayHandle<vtkm::cont::ArrayHandle<vtkm::FloatDefault>>();
    data_set.GetDomain(currentPartition)
        .AddCellField("velocityVector",
                      vtkm::cont::make_ArrayHandleSOA<vtkm::Vec3f_64>(
                          {xarrayV, yarrayV, zarrayV}));
  }

  vtkh::Recenter recenter;
  recenter.SetInput(&data_set);
  recenter.SetField("velocityVector");
  recenter.SetResultAssoc(vtkm::cont::Field::Association::POINTS);
  recenter.Update();
  vtkh::DataSet *recenteredData = recenter.GetOutput();

  runAdvection(recenteredData, mpiRank, numTasks);

  return 0;
}