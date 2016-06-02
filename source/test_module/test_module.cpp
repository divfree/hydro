/*******************************************************************
******************       CFD SOLVER 2D        **********************
********************************************************************/
#include "test_module.hpp"

namespace test_module
{

namespace registrators {

#ifdef MODULE_TEST
ModuleRegistrator<hydro<geom::geom2d::MeshStructured<double>>>
reg_double({"test"});
#endif

} // namespace registrators

} // namespace test_module

