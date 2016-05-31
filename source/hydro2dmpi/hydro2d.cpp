/*******************************************************************
******************       CFD SOLVER 2D        **********************
********************************************************************/
#include "hydro2d.hpp"


namespace hydro2D_uniform_MPI
{

TModule* CreateHydro2dDouble(TExperiment* ex) {
  return new hydro<geom::geom2d::MeshStructured<double>>(ex);
}

//TModule* CreateHydro2dSingle(TExperiment* ex) {
//  return new hydro<geom::geom2d::MeshStructured<float>>(ex);
//}

TModule* CreateHydro3dDouble(TExperiment* ex) {
  return new hydro<geom::geom3d::MeshStructured<double>>(ex);
}

//TModule* CreateHydro3dSingle(TExperiment* ex) {
//  return new hydro<geom::geom3d::MeshStructured<float>>(ex);
//}

// namespace end
}

