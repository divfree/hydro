/*******************************************************************
******************       CFD SOLVER 2D        **********************
********************************************************************/
#include "test_module.hpp"


namespace test_module
{

TModule* CreateTestDouble(TExperiment* ex) {
  return new hydro<geom::geom2d::MeshStructured<double>>(ex); 
}

TModule* CreateTestSingle(TExperiment* ex) {
  return new hydro<geom::geom2d::MeshStructured<float>>(ex);
}

} // namespace test_module

