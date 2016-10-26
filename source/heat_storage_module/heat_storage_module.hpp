/*******************************************************************
******************       CFD SOLVER 2D        **********************
********************************************************************/

#pragma once

#include "../control/experiment.hpp"
#include "../common/vect.hpp"
#include "../hydro2dmpi/mesh.hpp"
#include "../hydro2dmpi/mesh1d.hpp"
#include "../hydro2dmpi/output.hpp"
#include "../hydro2dmpi/linear.hpp"
#include "../hydro2dmpi/solver.hpp"
#include <memory>
#include <utility>

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
class HeatStorage : public solver::UnsteadySolver {
 public:
  static constexpr size_t dim = Mesh::dim;
  using Scal = typename Mesh::Scal;
  using MIdx = typename Mesh::MIdx;
  using Direction = typename Mesh::Direction;
  using Vect = typename Mesh::Vect;

  using IntIdx = geom::IntIdx;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using IdxNode = geom::IdxNode;

  using Layers = solver::Layers;

  template <class T>
  using FieldCell = geom::FieldCell<T>;
  template <class T>
  using FieldFace = geom::FieldFace<T>;
  template <class T>
  using FieldNode = geom::FieldNode<T>;

  using FuncRhs = std::function<Scal(Scal, Scal)>;

  static Scal RhsCos(Scal t, Scal x) {
    return std::cos(x);
  }

  class TesterMms {
   public:
    struct Diff {
      size_t num_cells;
      Scal diff;
      Scal h;
    };
    TesterMms(size_t num_cells_initial, size_t num_stages, size_t factor,
              Scal domain_length, size_t num_steps, double time_step,
              Scal fluid_velocity, Scal conductivity, Scal Tleft) {
      for (size_t num_cells = num_cells_initial, i = 0;
          i < num_stages; ++i, num_cells *= factor) {
        solvers_.emplace_back(num_cells, domain_length, time_step, fluid_velocity, conductivity, Tleft);
        auto& solver = solvers_.back();
        for (size_t n = 0; n < num_steps; ++n) {
          solver.StartStep();
          solver.CalcStep();
          solver.FinishStep();
          std::cout << num_cells << std::endl;
        }
      }
    }
    const std::vector<Diff>& GetDiffSeries() const {
      return diff_series_;
    }

   private:
    std::vector<HeatStorage> solvers_;
    std::vector<Diff> diff_series_;
  };

  class Scheduler {
   public:
    enum class State { Charging, Idle, Discharging };

    Scheduler(double d1, double d2, double d3, double d4)
        : d1_(d1), d2_(d2), d3_(d3), d4_(d4) {}

    State GetState(double t) {
      double cycle_duration = d1_ + d2_ + d3_ + d4_;
      size_t cycle = static_cast<size_t>(t / cycle_duration);
      double offset = t - cycle * cycle_duration;
      return offset < d1_ ? State::Charging :
          offset < d1_ + d2_ ? State::Idle :
          offset < d1_ + d2_ + d3_ ? State::Discharging : State::Idle;
    }
    size_t GetStateIdx(double t) {
      switch (GetState(t)) {
        case State::Charging:
          return 1;
        case State::Discharging:
          return 2;
        case State::Idle:
          return 3;
        default:
          assert(false);
      }
    }
   private:
    double d1_, d2_, d3_, d4_;
  };

  HeatStorage(size_t num_cells, Scal domain_length, double time_step,
              double fluid_velocity, double conductivity, double Tleft)
      : UnsteadySolver(0., time_step),
        fluid_velocity_(fluid_velocity), conductivity_(conductivity), Tleft_(Tleft) {

    //,
    //      scheduler_(P_double["duration_1"], P_double["duration_2"],
    //                 P_double["duration_3"], P_double["duration_4"])

    // Prepare mesh nodes
    MIdx mesh_size(num_cells);

    geom::Rect<Vect> domain(Vect(0.), Vect(domain_length));

    //auto domain_size = domain.GetDimensions();
    // Create mesh
    geom::InitUniformMesh(mesh, domain, mesh_size);

    // Init fields
    fc_temperature_fluid_.time_curr.Reinit(mesh, 0.);
    fc_temperature_fluid_.time_prev.Reinit(mesh, 0.);
    fc_temperature_solid_.time_curr.Reinit(mesh, 0.);
    fc_temperature_solid_.time_prev.Reinit(mesh, 0.);


    std::ofstream f("mesh" + IntToStr(num_cells));
    for (auto idxcell : mesh.Cells()) {
      f << mesh.GetVolume(idxcell) << " ";
    }
  }
  const FieldCell<Scal>& GetFluidTemperature(Layers layer) {
    return fc_temperature_fluid_.Get(layer);
  }
  const FieldCell<Scal>& GetFluidTemperature() {
    return fc_temperature_fluid_.time_curr;
  }
  void CalcStep() {
    auto& tf = fc_temperature_fluid_.time_prev;
    auto& tf_new = fc_temperature_fluid_.time_curr;
    std::swap(tf, tf_new);
    auto& ts = fc_temperature_solid_.time_prev;
    auto& ts_new = fc_temperature_solid_.time_curr;
    std::swap(ts, ts_new);

    const Scal h = mesh.GetVolume(IdxCell(0));
    const Scal dt = this->GetTimeStep();
    const Scal uf = fluid_velocity_;
    const Scal alpha = conductivity_;
    const Scal Tleft = Tleft_;

    for (IdxCell c : mesh.Cells()) {
      IdxCell cm = mesh.GetNeighbourCell(c, 0);
      IdxCell cp = mesh.GetNeighbourCell(c, 1);
      if (cm.IsNone()) { // left boundary
        tf_new[c] = tf[c] - dt * uf * (tf[c] - Tleft) / h
            + dt * alpha * (-tf[c] + tf[cp]) / sqr(h);
      } else if (cp.IsNone()) { // right boundary
        tf_new[c] = tf[c] - dt * uf * (tf[c] - tf[cm]) / h
            + dt * alpha * (tf[cm] - tf[c]) / sqr(h);
      } else {
        tf_new[c] = tf[c] - dt * uf * (tf[c] - tf[cm]) / h
            + dt * alpha * (tf[cm] - 2. * tf[c] + tf[cp]) / sqr(h);
      }
    }
  }
  const Mesh& GetMesh() const { return mesh; }

 private:
  Mesh mesh;
  solver::LayersData<FieldCell<Scal>>
  fc_temperature_fluid_, fc_temperature_solid_;
  Scal fluid_velocity_, conductivity_, Tleft_;
  //Scheduler scheduler_;
};

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

  using Layers = solver::Layers;

  template <class T>
  using FieldCell = geom::FieldCell<T>;
  template <class T>
  using FieldFace = geom::FieldFace<T>;
  template <class T>
  using FieldNode = geom::FieldNode<T>;

  // TODO: pack common declarations somewhere (Scal, Vect, FieldCell, etc)
 public:
  hydro(TExperiment* _ex);
  ~hydro() {}
  void step();
  void write_results(bool force=false);

  output::Content content, content_scalar;
  std::shared_ptr<output::Session> session, session_scalar;

  double last_frame_time_;
  double last_frame_scalar_time_;

  std::shared_ptr<HeatStorage<Mesh>> solver_;
  const Mesh * pmesh;
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

  geom::Rect<Vect> domain(GetVect<Vect>(P_vect["A"]),
                          GetVect<Vect>(P_vect["B"]));

  solver_ = std::make_shared<HeatStorage<Mesh>>(
      P_int["Nx"], domain.GetDimensions()[0], dt,
      P_double["uf"], P_double["alpha"], P_double["T_left"]);

  pmesh = &solver_->GetMesh();

  P_int.set("cells_number", static_cast<int>(pmesh->GetNumCells()));

  content = {
      std::make_shared<output::EntryFunction<Scal, IdxNode, Mesh>>(
          "x", *pmesh,
          [this](IdxNode idx) { return pmesh->GetNode(idx)[0]; })
//      , std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
//          "tf", *pmesh,
//          [this](IdxCell idx) { return fc_temperature_fluid_.time_curr[idx]; })
  };

  auto P = [this](std::string entry, std::string parameter) {
    return std::make_shared<output::EntryScalarFunction<Scal>>(
        entry, [this, parameter](){ return static_cast<Scal>(P_double[parameter]); });
  };
  auto Pint = [this](std::string entry, std::string parameter) {
    return std::make_shared<output::EntryScalarFunction<Scal>>(
        entry, [this, parameter](){ return static_cast<Scal>(P_int[parameter]); });
  };
  content_scalar = { P("time", "t"), Pint("n", "n")
      //, std::make_shared<output::EntryScalarFunction<Scal>>(
      //    "status", [this](){ return static_cast<Scal>(scheduler_.GetStateIdx(P_double["t"])); })
  };

  if (!P_string.exist(_plt_title)) {
    P_string.set(_plt_title, P_string[_exp_name]);
  }
  if (!P_string.exist("filename_field")) {
    P_string.set("filename_field", P_string[_exp_name] + ".field.dat");
  }
  if (!P_string.exist("filename_scalar")) {
    P_string.set("filename_scalar", P_string[_exp_name] + ".scalar.dat");
  }

  session = std::make_shared<output::SessionPlain<Mesh>>(
      content, P_string["filename_field"], *pmesh);

  session_scalar = std::make_shared<output::SessionPlainScalar<Scal>>(
      content_scalar, P_string["filename_scalar"]);

  last_frame_time_ = 0;
  last_frame_scalar_time_ = 0;
  session->Write(0., "initial");
  session_scalar->Write();

  typename HeatStorage<Mesh>::TesterMms tester(
      4, 4, 2, 1., 1000, dt,
      P_double["uf"], P_double["alpha"], P_double["T_left"]);
}

template <class Mesh>
void hydro<Mesh>::step() {
  ex->timer_.Push("step");


  ++P_int["n"];
  ex->timer_.Pop();
}

template <class Mesh>
void hydro<Mesh>::write_results(bool force) {
  if (ecast(P_bool("no_output"))) {
    return;
  }

  const double time = P_double["t"];
  const double total_time = P_double["T"];
  const size_t max_frame_index = P_int["max_frame_index"];
  const double frame_duration = total_time / max_frame_index;

  if (force || (!ecast(P_bool("no_mesh_output")) &&
      time >= last_frame_time_ + frame_duration)) {
    last_frame_time_ = time;
    session->Write(time, "step");
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
