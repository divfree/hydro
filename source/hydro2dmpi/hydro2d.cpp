/*******************************************************************
******************       CFD SOLVER 2D        **********************
********************************************************************/
#include "hydro2d.hpp"


namespace hydro2D_uniform_MPI
{

namespace registrators {

ModuleRegistrator<hydro<geom::geom2d::MeshStructured<double>>>
reg_2d_double({"hydro2D_uniform_MPI", "hydro2d"});

ModuleRegistrator<hydro<geom::geom3d::MeshStructured<double>>>
reg_3d_double({"hydro3D_uniform_MPI", "hydro3d"});

} // namespace registrators

// namespace end
}

