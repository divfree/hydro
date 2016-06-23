#include "hydro2d.hpp"

namespace hydro2D_uniform_MPI {

template <>
geom::Vect<double, 2> GetVect<geom::Vect<double, 2>>(const column<double>& v) {
  return geom::Vect<double, 2>(v[0], v[1]);
}

template <>
geom::Vect<float, 2> GetVect<geom::Vect<float, 2>>(const column<double>& v) {
  return geom::Vect<float, 2>(float(v[0]), float(v[1]));
}

namespace registrators {

#ifdef MODULE_HYDRO_2D
ModuleRegistrator<hydro<geom::geom2d::MeshStructured<double>>>
reg_2d_double({"hydro2D_uniform_MPI", "hydro2d"});
#endif

} // namespace registrators

} // namespace hydro2D_uniform_MPI
