/*******************************************************************
******************       CFD SOLVER 2D        **********************
********************************************************************/

#pragma once

#include "../control/experiment.hpp"
#include "../common/vect.hpp"
#include "../hydro2dmpi/mesh.hpp"
#include "../hydro2dmpi/mesh1d.hpp"
#include "../hydro2dmpi/output.hpp"
#include "../hydro2dmpi/output_tecplot.hpp"
#include "../hydro2dmpi/linear.hpp"
#include "../hydro2dmpi/solver.hpp"
#include <memory>

// TODO: Change parameters on events (e.g. certain time moments)

namespace heat_storage_module
{

template <class Vect>
Vect GetVect(const column<double>& v);

template <>
geom::Vect<double, 1> GetVect<geom::Vect<double, 1>>(const column<double>& v) {
  return geom::Vect<double, 1>(v[0]);
}


template <>
geom::Vect<float, 1> GetVect<geom::Vect<float, 1>>(const column<double>& v) {
  return geom::Vect<float, 1>(float(v[0]));
}

template <class Mesh>
class hydro : public TModule
{
  static constexpr size_t dim = Mesh::dim;
  using Scal = typename Mesh::Scal;
  using MIdx = typename Mesh::MIdx;
  using Direction = typename Mesh::Direction;
  using Vect = typename Mesh::Vect;

  using IntIdx = geom::IntIdx;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using IdxNode = geom::IdxNode;

  template <class T>
  using FieldCell = geom::FieldCell<T>;
  template <class T>
  using FieldFace = geom::FieldFace<T>;
  template <class T>
  using FieldNode = geom::FieldNode<T>;
 public:
  hydro(TExperiment* _ex);
  ~hydro() {}
  void step();
  void write_results(bool force=false);
  Mesh mesh;
  output::Content content, content_scalar;
  std::shared_ptr<output::Session> session, session_scalar;

  double last_frame_time_;
  double last_frame_scalar_time_;
  void StepMultigrid();
};

template <class Mesh>
hydro<Mesh>::hydro(TExperiment* _ex)
    : TExperiment_ref(_ex), TModule(_ex)
{
  P_int.set("last_s", 0);
  P_double.set("last_R", 0);
  P_double.set("last_Rn", 0);

  P_int.set("s_sum", 0);
  P_int.set("s_max", 0);
  P_int.set("s", 0);
  P_int.set("current_frame", 0);
  P_int.set("current_frame_scalar", 0);

  // Prepare mesh nodes
  MIdx mesh_size;
  for (size_t i = 0; i < dim; ++i) {
    mesh_size[i] = P_int[std::string("N") + Direction(i).GetLetter()];
  }
  geom::Rect<Vect> domain(GetVect<Vect>(P_vect["A"]),
                          GetVect<Vect>(P_vect["B"]));
  auto domain_size = domain.GetDimensions();
  // Create mesh
  geom::InitUniformMesh(mesh, domain, mesh_size);

  P_int.set("cells_number", static_cast<int>(mesh.GetNumCells()));

  content = {
      std::make_shared<output::EntryFunction<Scal, IdxNode, Mesh>>(
          "x", mesh,
          [this](IdxNode idx) { return mesh.GetNode(idx)[0]; })
  };

  auto P = [this](std::string entry, std::string parameter) {
    return std::make_shared<output::EntryScalarFunction<Scal>>(
        entry, [this, parameter](){ return static_cast<Scal>(P_double[parameter]); });
  };
  content_scalar = { P("time", "t") };

  if (!P_string.exist(_plt_title)) {
    P_string.set(_plt_title, P_string[_exp_name]);
  }
  if (!P_string.exist("filename_field")) {
    P_string.set("filename_field", P_string[_exp_name] + ".field.plt");
  }
  if (!P_string.exist("filename_scalar")) {
    P_string.set("filename_scalar", P_string[_exp_name] + ".scalar.plt");
  }

//  session = std::make_shared<output::SessionTecplotBinaryStructured<Mesh>>(
//      content, P_string[_plt_title], mesh, P_string["filename_field"]);

  session_scalar = std::make_shared<output::SessionTecplotAsciiScalar<Scal>>(
      content_scalar, P_string[_plt_title],
      "zone_title", P_string["filename_scalar"]);

  last_frame_time_ = 0;
  last_frame_scalar_time_ = 0;
//  session->Write(0., "initial");
  session_scalar->Write();


  // Adhoc: Write mesh nodes
  std::ofstream f("mesh.dat");
  for (auto idxcell : mesh.Cells()) {
    f << idxcell.GetRaw() << " " << mesh.GetCenter(idxcell) << std::endl;
  }
}

template <class Mesh>
void hydro<Mesh>::step() {
  ex->timer_.Push("step");

  ex->timer_.Pop();
}

template <class Mesh>
void hydro<Mesh>::write_results(bool force) {
  if (ecast(P_bool("no_output"))) {
    return;
  }

  const double time = P_int["n"];
  const double total_time = P_double["T"];
  const size_t max_frame_index = P_int["max_frame_index"];
  const double frame_duration = total_time / max_frame_index;

  if (force || (!ecast(P_bool("no_mesh_output")) &&
      time >= last_frame_time_ + frame_duration)) {
    last_frame_time_ = time;
//    session->Write(time, "step");
    logger() << "Frame " << (P_int["current_frame"]++) << ": t=" << time;
  }

  const size_t max_frame_scalar_index = P_int["max_frame_scalar_index"];
  const double frame_scalar_duration = total_time / max_frame_scalar_index;

  if (force || time >= last_frame_scalar_time_ + frame_scalar_duration) {
    last_frame_scalar_time_ = time;
    session_scalar->Write();
    logger() << "Frame_scalar " << (P_int["current_frame_scalar"]++)
        << ": t=" << time;
  }
}

} // namespace heat_storage_module
