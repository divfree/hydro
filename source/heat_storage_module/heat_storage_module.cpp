/*******************************************************************
******************       CFD SOLVER 2D        **********************
********************************************************************/
#include "heat_storage_module.hpp"

namespace heat_storage_module
{

namespace registrators {

#ifdef MODULE_HEAT_STORAGE
ModuleRegistrator<hydro<geom::geom2d::MeshStructured<double>>>
reg_double({"test"});
#endif

} // namespace registrators

} // namespace heat_storage_module

