/*******************************************************************
******************       CFD SOLVER 2D        **********************
********************************************************************/
#include "heat_storage_module.hpp"

namespace heat_storage_module
{

namespace registrators {

#ifdef MODULE_HEAT_STORAGE
ModuleRegistrator<hydro<geom::geom1d::MeshStructured<double>>>
reg_double({"heat_storage"});
#endif

} // namespace registrators

} // namespace heat_storage_module

