#include "hydro2d.hpp"

namespace hydro2D_uniform_MPI {

template <>
geom::Vect<double, 3> GetVect<geom::Vect<double, 3>>(const column<double>& v) {
  return geom::Vect<double, 3>(v[0], v[1], v[2]);
}

template <>
geom::Vect<float, 3> GetVect<geom::Vect<float, 3>>(const column<double>& v) {
  return geom::Vect<float, 3>(float(v[0]), float(v[1]), float(v[2]));
}

namespace registrators {

#ifdef MODULE_HYDRO_3D
ModuleRegistrator<hydro<geom::geom3d::MeshStructured<double>>>
reg_3d_double({"hydro3D_uniform_MPI", "hydro3d"});
#endif

} // namespace registrators

} // namespace hydro2D_uniform_MPI
