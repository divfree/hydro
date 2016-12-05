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
#include <map>

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

template <class T>
std::map<std::string, T> ConvertToMap(const BinSearchSet2<std::string, T>& p) {
  std::map<std::string, T> res;
  for (size_t i = 0; i < p.get_length(); ++i) {
    res[p.get_key(i)] = p.get_data(i);
  }
  return res;
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

  using FuncTX = std::function<Scal(Scal, Scal)>;

  static FieldCell<Scal> Evaluate(FuncTX func, double t, const Mesh& mesh) {
    FieldCell<Scal> res(mesh);
    for (auto idxcell : mesh.Cells()) {
      Scal x = mesh.GetCenter(idxcell)[0];
      res[idxcell] = func(t, x);
    }
    return res;
  }

  class TesterMms {
   public:
    struct Entry {
      Scal diff_prev_fluid;
      Scal diff_prev_solid;
      Scal error_fluid;
      Scal error_solid;
      Scal h;
      Mesh mesh;
      FieldCell<Scal> fc_fluid_temperature;
      FieldCell<Scal> fc_solid_temperature;
      FieldCell<Scal> fc_exact_fluid_temperature;
      FieldCell<Scal> fc_exact_solid_temperature;
      std::shared_ptr<HeatStorage> solver;
    };
    TesterMms(size_t num_cells_initial, size_t num_stages, size_t factor,
              Scal domain_length, size_t num_steps, double time_step,
              Scal step_threshold,
              Scal fluid_velocity, Scal conductivity,
              Scal T_left, Scal T_right,
              FuncTX func_rhs_fluid, FuncTX func_exact_fluid_temperature,
              FuncTX func_rhs_solid, FuncTX func_exact_solid_temperature) {
      std::ofstream stat("mms_statistics.dat");
      stat.precision(20);
      stat << "num_cells error_fluid error_solid diff_fluid diff_solid dt num_steps step_diff_fluid step_diff_solid" << std::endl;
      
      std::string field_name_prefix = "field_T_";

      for (size_t num_cells = num_cells_initial, i = 0;
          i < num_stages; ++i, num_cells *= factor) {
        series_.emplace_back();
        Entry& entry = series_.back();

        // Create mesh
        Mesh& mesh = entry.mesh;
        MIdx mesh_size(num_cells);
        geom::Rect<Vect> domain(Vect(0.), Vect(domain_length));
        geom::InitUniformMesh(mesh, domain, mesh_size);

        FieldCell<Scal> fc_rhs_fluid(mesh, 0.), fc_rhs_solid(mesh, 0.);
        fc_rhs_fluid = Evaluate(func_rhs_fluid, 0., mesh);
        fc_rhs_solid = Evaluate(func_rhs_solid, 0., mesh);

        std::map<std::string, double> p_double = {
            {"fluid_velocity", fluid_velocity},
            {"density_fluid", 2.},
            {"density_solid", 2.},
            {"specific_heat_fluid", 1.},
            {"specific_heat_solid", 1.},
            {"conductivity_fluid", conductivity},
            {"conductivity_solid", conductivity},
            {"temperature_hot", T_left},
            {"temperature_cold", T_right},
            {"temperature_reference", 1.},
            {"heat_exchange", 0.},
            {"porosity", 0.5},
        };

        entry.solver = std::make_shared<HeatStorage>(
            mesh, time_step, p_double,
            &fc_rhs_fluid, &fc_rhs_solid,
            Scheduler(1e9, 0., 0., 0.),
            true, func_exact_fluid_temperature(0., 0.),
            func_exact_fluid_temperature(0., domain_length) // solid assumed equal
        );
        auto& solver = *entry.solver;
        size_t actual_num_steps = num_steps;
        Scal step_diff_fluid = 0., step_diff_solid = 0.;
        for (size_t n = 0; n < num_steps; ++n) {
          solver.StartStep();
          solver.CalcStep();
          solver.FinishStep();
          step_diff_fluid = solver::CalcDiff(solver.GetFluidTemperature(Layers::time_curr),
                                             solver.GetFluidTemperature(Layers::time_prev),
                                             mesh);
          step_diff_solid = solver::CalcDiff(solver.GetSolidTemperature(Layers::time_curr),
                                             solver.GetSolidTemperature(Layers::time_prev),
                                             mesh);
          if (std::max(step_diff_fluid, step_diff_solid) < step_threshold) {
            actual_num_steps = n;
            break;
          }
        }

        entry.fc_fluid_temperature = solver.GetFluidTemperature();
        entry.diff_prev_fluid = series_.empty() ? 0. :
            solver::CalcDiff(entry.fc_fluid_temperature,
                             GetInterpolated(series_.back().fc_fluid_temperature, series_.back().mesh, mesh),
                             mesh);

        entry.fc_solid_temperature = solver.GetSolidTemperature();
        entry.diff_prev_solid = series_.empty() ? 0. :
            solver::CalcDiff(entry.fc_solid_temperature,
                             GetInterpolated(series_.back().fc_solid_temperature, series_.back().mesh, mesh),
                             mesh);

        entry.fc_exact_fluid_temperature = Evaluate(func_exact_fluid_temperature, 0., mesh);
        entry.error_fluid = solver::CalcDiff(entry.fc_exact_fluid_temperature,
                                       entry.fc_fluid_temperature,
                                       mesh);
        
        entry.fc_exact_solid_temperature = Evaluate(func_exact_solid_temperature, 0., mesh);
        entry.error_solid = solver::CalcDiff(entry.fc_exact_solid_temperature,
                                       entry.fc_solid_temperature,
                                       mesh);

        stat << num_cells << ' '
            << entry.error_fluid << ' ' << entry.error_solid << ' '
            << entry.diff_prev_fluid << ' ' << entry.diff_prev_solid << ' '
            << time_step << ' ' << actual_num_steps << ' '
            << step_diff_fluid << ' ' << step_diff_solid << std::endl;

        solver.WriteField(solver.GetFluidTemperature(),
                          solver.GetSolidTemperature(),
                          field_name_prefix + IntToStr(num_cells) + ".dat");
      }

      // Write the exact solution on the finest mesh
      if (!series_.empty()) {
        auto& entry = series_.back();
        entry.solver->WriteField(entry.fc_exact_fluid_temperature,
                                 entry.fc_exact_solid_temperature,
                                 field_name_prefix + "exact" + ".dat");
      }
    }
    const std::vector<Entry>& GetSeries() const {
      return series_;
    }

   private:
    std::vector<Entry> series_;
  };

  class Scheduler {
   public:
    enum class State { Charging, Idle, Discharging };

    Scheduler(double d1, double d2, double d3, double d4)
        : d1_(d1), d2_(d2), d3_(d3), d4_(d4) {}

    State GetState(double t) const {
      double cycle_duration = d1_ + d2_ + d3_ + d4_;
      size_t cycle = static_cast<size_t>(t / cycle_duration);
      double offset = t - cycle * cycle_duration;
      return offset < d1_ ? State::Charging :
          offset < d1_ + d2_ ? State::Idle :
          offset < d1_ + d2_ + d3_ ? State::Discharging : State::Idle;
    }
    size_t GetStateIdx(double t) const {
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
    Scal GetStateFactor(double t) const {
      switch (GetState(t)) {
        case State::Charging:
          return 1.;
        case State::Discharging:
          return -1.;
        case State::Idle:
          return 0.;
        default:
          assert(false);
      }
    }
   private:
    double d1_, d2_, d3_, d4_;
  };

  class ExergyCalculator {
   public:
    // requires: constant mass flux, constant fluid specific heat
    // ensures: exergy is calculated up to a factor and an additive constant
    ExergyCalculator(const HeatStorage* parent) : parent(parent) {}
    void Update() {
      Scal curr_time = parent->GetTime();
      State curr_state = parent->scheduler_.GetState(curr_time);
      if (fresh_) {
        fresh_ = false;
        prev_state_ = curr_state;
        prev_time_ = curr_time;
        return;
      }
      Scal T0 = parent->p_double.at("temperature_reference");
      IdxCell cleft(0);
      IdxCell cright(parent->mesh.GetNumCells() - 1);
      if (prev_state_ == curr_state) {
        const Scal dt = curr_time - prev_time_;
        auto f = [=](Scal T) { return (T - T0 * std::log(T)) * dt; };
        integral_left_ += f(parent->GetFluidTemperature()[cleft]);
        integral_right_ += f(parent->GetFluidTemperature()[cright]);
      } else {
        if (prev_state_ == State::Charging) {
          exergy_c_in = integral_left_;
          exergy_c_out = integral_right_;
          temperature_increase =
              parent->GetFluidTemperature()[cright] -
              parent->temperature_cold_;
        } else if (prev_state_ == State::Discharging) {
          exergy_d_in = integral_right_;
          exergy_d_out = integral_left_;
        }
        // update efficieny value if the denominator is non-zero
        if (exergy_c_in != exergy_c_out) {
          efficiency = (exergy_d_out - exergy_d_in) /
              (exergy_c_in - exergy_c_out);
        }
        // update thermal energy values
        CalcThermalEnergy();
        if (curr_state == State::Charging) {
          thermal_energy_before_charging = thermal_energy;
        } else if (curr_state == State::Discharging) {
          thermal_energy_before_discharging = thermal_energy;
          capacity_factor =
              (thermal_energy_before_discharging - thermal_energy_before_charging) /
              thermal_energy_limit;
        }

        integral_left_ = 0.;
        integral_right_ = 0.;
      }
      prev_state_ = curr_state;
      prev_time_ = curr_time;
    }
    Scal exergy_c_in{0.}, exergy_c_out{0.}, exergy_d_in{0.}, exergy_d_out{0.};
    Scal efficiency{0.};
    Scal thermal_energy{0.}; // without D^2 * pi / 4 factor
    Scal thermal_energy_before_charging{0.};
    Scal thermal_energy_before_discharging{0.};
    Scal thermal_energy_limit{0.};
    Scal capacity_factor{0.};
    Scal temperature_increase{0.}; // temp. incr. at end of charging

   private:
    using State = typename Scheduler::State;
    const HeatStorage* parent;
    bool fresh_{0.};
    Scal integral_left_{0.}, integral_right_{0.};
    State prev_state_;
    Scal prev_time_;
    void CalcThermalEnergy() {
      auto& fc_fluid = parent->GetFluidTemperature();
      auto& fc_solid = parent->GetSolidTemperature();
      const Mesh& mesh = parent->GetMesh();
      Scal int_fluid = 0., int_solid = 0.;
      for (auto idxcell : mesh.Cells()) {
        int_fluid += mesh.GetVolume(idxcell) *
            (fc_fluid[idxcell] - parent->temperature_cold_);
        int_solid += mesh.GetVolume(idxcell) *
            (fc_solid[idxcell] - parent->temperature_cold_);
      }

      thermal_energy =
          parent->porosity_ * parent->density_fluid_ *
          parent->specific_heat_fluid_ * int_fluid +
          (1. - parent->porosity_) * parent->density_solid_ *
          parent->specific_heat_solid_ * int_solid;

      Scal height = 0.;
      for (auto idxcell : mesh.Cells()) {
        height += mesh.GetVolume(idxcell);
      }

      thermal_energy_limit =
          (parent->porosity_ * parent->density_fluid_ *
          parent->specific_heat_fluid_ +
          (1. - parent->porosity_) * parent->density_solid_ *
          parent->specific_heat_solid_) *
          (parent->temperature_hot_ - parent->temperature_cold_) *
          height;
    }
  };

  HeatStorage(const Mesh& mesh,
              double time_step,
              const std::map<std::string, double>& p_double,
              const FieldCell<Scal>* p_fc_rhs_fluid,
              const FieldCell<Scal>* p_fc_rhs_solid,
              const Scheduler& scheduler,
              bool given_temperature,
              Scal given_temperature_left, Scal given_temperature_right)
      : UnsteadySolver(0., time_step),
        mesh(mesh), p_double(p_double),
        fluid_velocity_(p_double.at("fluid_velocity")),
        density_fluid_(p_double.at("density_fluid")),
        density_solid_(p_double.at("density_solid")),
        conductivity_fluid_(p_double.at("conductivity_fluid")),
        conductivity_solid_(p_double.at("conductivity_solid")),
        specific_heat_fluid_(p_double.at("specific_heat_fluid")),
        specific_heat_solid_(p_double.at("specific_heat_solid")),
        temperature_hot_(p_double.at("temperature_hot")),
        temperature_cold_(p_double.at("temperature_cold")),
        heat_exchange_(p_double.at("heat_exchange")),
        porosity_(p_double.at("porosity")),
        p_fc_rhs_fluid_(p_fc_rhs_fluid),
        p_fc_rhs_solid_(p_fc_rhs_solid),
        scheduler_(scheduler),
        given_temperature_(given_temperature),
        given_temperature_left_(given_temperature_left),
        given_temperature_right_(given_temperature_right)
  {

    // Init fields
    fc_temperature_fluid_.time_curr.Reinit(mesh, temperature_cold_);
    fc_temperature_fluid_.time_prev.Reinit(mesh, temperature_cold_);
    fc_temperature_solid_.time_curr.Reinit(mesh, temperature_cold_);
    fc_temperature_solid_.time_prev.Reinit(mesh, temperature_cold_);
  }
  const FieldCell<Scal>& GetFluidTemperature(Layers layer) const {
    return fc_temperature_fluid_.Get(layer);
  }
  const FieldCell<Scal>& GetFluidTemperature() const {
    return fc_temperature_fluid_.time_curr;
  }
  const FieldCell<Scal>& GetSolidTemperature(Layers layer) const {
    return fc_temperature_solid_.Get(layer);
  }
  const FieldCell<Scal>& GetSolidTemperature() const {
    return fc_temperature_solid_.time_curr;
  }
  static Scal GetInterpolated(Scal x, Scal x_left, Scal x_right, Scal u_left, Scal u_right) {
    return ((x - x_left) * u_right + (x_right - x) * u_left) / (x_right - x_left);
  }
  // Applies linear interpolation on 1D mesh
  static FieldCell<Scal> GetInterpolated(const FieldCell<Scal>& fc_src, const Mesh& mesh_src, const Mesh& mesh_dest) {
	  static_assert(Mesh::dim == 1, "HeatStorage::GetInterpolated() requires a 1D mesh");
	  FieldCell<Scal> res(mesh_dest);
	  IdxCell idx_left_src = IdxCell(0);
	  IdxCell idx_right_src = idx_left_src;
	  for (IdxCell idx_dest : mesh_dest.Cells()) {
	    Scal x_dest = mesh_dest.GetCenter(idx_dest)[0];
	    while (mesh_src.GetCenter(idx_right_src)[0] < x_dest) {
	      IdxCell idx_new_src = mesh_src.GetNeighbourCell(idx_right_src, 1);
	      if (idx_new_src.IsNone()) {
	        break;
	      }
	      idx_left_src = idx_right_src;
	      idx_right_src = idx_new_src;
	    }
	    res[idx_dest] = GetInterpolated(x_dest, mesh_src.GetCenter(idx_left_src)[0], mesh_src.GetCenter(idx_right_src)[0],
	                                    fc_src[idx_left_src], fc_src[idx_right_src]);
	  }
	  return res;
  }
  void WriteField(const FieldCell<Scal>& fc_u, std::string filename) const {
    output::Content content = {
        std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
            "x", mesh,
            [this](IdxCell idx) { return mesh.GetCenter(idx)[0]; })
        , std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
            "u", mesh,
            [&](IdxCell idx) { return fc_u[idx]; })
    };
    output::SessionPlain<Mesh> session(content, filename, mesh);
    session.Write(0., "field");
  }
  void WriteField(const FieldCell<Scal>& fc_t_fluid,
                  const FieldCell<Scal>& fc_t_solid,
                  std::string filename) const {
    output::Content content = {
        std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
            "x", mesh,
            [this](IdxCell idx) { return mesh.GetCenter(idx)[0]; })
        , std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
            "T_f", mesh,
            [&](IdxCell idx) { return fc_t_fluid[idx]; })
        , std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
            "T_s", mesh,
            [&](IdxCell idx) { return fc_t_solid[idx]; })
    };
    output::SessionPlain<Mesh> session(content, filename, mesh);
    session.Write(0., "field");
  }
  void CalcStep() {
    auto& Tf = fc_temperature_fluid_.time_prev;
    auto& Tf_new = fc_temperature_fluid_.time_curr;
    std::swap(Tf, Tf_new);
    auto& Ts = fc_temperature_solid_.time_prev;
    auto& Ts_new = fc_temperature_solid_.time_curr;
    std::swap(Ts, Ts_new);

    const Scal h = mesh.GetVolume(IdxCell(0)); // uniform mesh assumed
    const Scal dt = this->GetTimeStep();
    const Scal uf = fluid_velocity_ * scheduler_.GetStateFactor(GetTime());
    const Scal rho_f = density_fluid_;
    const Scal rho_s = density_solid_;
    const Scal C_f = specific_heat_fluid_;
    const Scal C_s = specific_heat_solid_;
    const Scal k_f = conductivity_fluid_;
    const Scal k_s = conductivity_solid_;
    const Scal eps = porosity_;
    const Scal alpha_f = k_f / (eps * rho_f * C_f);
    const Scal alpha_s = k_s / ((1. - eps) * rho_s * C_s);
    const Scal T_in = uf > 0. ? temperature_hot_ : temperature_cold_;

    // Equation: dT/dt + div(fluxes) = 0
    FieldFace<Scal> ff_flux_fluid(mesh, 0.);
    FieldFace<Scal> ff_flux_solid(mesh, 0.);
    for (IdxFace idxface : mesh.Faces()) {
      IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
      IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
      auto& flux_fluid = ff_flux_fluid[idxface];
      auto& flux_solid = ff_flux_solid[idxface];
      if (cm.IsNone()) { // left boundary
        flux_fluid += uf * (uf > 0. ? T_in : Tf[cp]);
        if (given_temperature_) {
          auto val = given_temperature_left_;
          flux_fluid += -alpha_f * (Tf[cp] - val) / (0.5 * h);
          flux_solid += -alpha_s * (Ts[cp] - val) / (0.5 * h);
        }
      } else if (cp.IsNone()) { // right boundary
        flux_fluid += uf * (uf < 0. ? T_in : Tf[cm]);
        if (given_temperature_) {
          auto val = given_temperature_right_;
          flux_fluid += -alpha_f * (val - Tf[cm]) / (0.5 * h);
          flux_solid += -alpha_s * (val - Ts[cm]) / (0.5 * h);
        }
      } else {
        // convection: first order upwind
        flux_fluid += uf * (uf > 0. ? Tf[cm] : Tf[cp]);
        // diffusion: central second order
        flux_fluid += -alpha_f * (Tf[cp] - Tf[cm]) / h;
        flux_solid += -alpha_s * (Ts[cp] - Ts[cm]) / h;
      }
    }

    // Time integration of flux terms
    for (IdxCell idxcell : mesh.Cells()) {
      IdxFace fm = mesh.GetNeighbourFace(idxcell, 0);
      IdxFace fp = mesh.GetNeighbourFace(idxcell, 1);
      const Scal fluxsum_fluid = ff_flux_fluid[fp] - ff_flux_fluid[fm];
      Tf_new[idxcell] = Tf[idxcell] - dt / h * fluxsum_fluid;
      const Scal fluxsum_solid = ff_flux_solid[fp] - ff_flux_solid[fm];
      Ts_new[idxcell] = Ts[idxcell] - dt / h * fluxsum_solid;
    }

    // Time integration of source terms
    if (p_fc_rhs_fluid_) {
      for (IdxCell idxcell : mesh.Cells()) {
        Tf_new[idxcell] += dt * (*p_fc_rhs_fluid_)[idxcell];
      }
    }
    if (p_fc_rhs_solid_) {
      for (IdxCell idxcell : mesh.Cells()) {
        Ts_new[idxcell] += dt * (*p_fc_rhs_solid_)[idxcell];
      }
    }

    // Implicit heat exchange
    const Scal h_v = heat_exchange_;
    const Scal hf = h_v / (eps * rho_f * C_f);
    const Scal hs = h_v / ((1. - eps) * rho_s * C_s);
    const Scal a = 1. + hf * dt;
    const Scal b = -hf * dt;
    const Scal c = -hs * dt;
    const Scal d = 1. + hs * dt;
    const Scal det = a * d - b * c;
    const Scal ai = d / det;
    const Scal bi = -b / det;
    const Scal ci = -c / det;
    const Scal di = a / det;
    for (IdxCell idxcell : mesh.Cells()) {
      const Scal tf = Tf_new[idxcell];
      const Scal ts = Ts_new[idxcell];
      Tf_new[idxcell] = ai * tf + bi * ts;
      Ts_new[idxcell] = ci * tf + di * ts;
    }

    // Update exergy
    exergy_.Update();
  }
  const Mesh& GetMesh() const { return mesh; }
  Scal GetState() const { 
    return scheduler_.GetStateFactor(GetTime());
  }
  const ExergyCalculator& GetExergy() const {
    return exergy_;
  }

 private:
  const Mesh& mesh;
  std::map<std::string, double> p_double;
  solver::LayersData<FieldCell<Scal>>
  fc_temperature_fluid_, fc_temperature_solid_;
  Scal fluid_velocity_;
  Scal density_fluid_, density_solid_;
  Scal conductivity_fluid_, conductivity_solid_;
  Scal specific_heat_fluid_, specific_heat_solid_;
  Scal temperature_hot_, temperature_cold_;
  Scal heat_exchange_;
  Scal porosity_;
  const FieldCell<Scal>* p_fc_rhs_fluid_;
  const FieldCell<Scal>* p_fc_rhs_solid_;
  Scheduler scheduler_;
  ExergyCalculator exergy_{this};
  bool given_temperature_;
  Scal given_temperature_left_, given_temperature_right_;
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
  Mesh mesh;
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

  // Create mesh
  MIdx mesh_size(P_int["Nx"]);
  geom::InitUniformMesh(mesh, domain, mesh_size);

  auto scheduler = typename HeatStorage<Mesh>::Scheduler(
      P_double["duration_charging"], P_double["duration_idle_charging"],
      P_double["duration_discharging"], P_double["duration_idle_discharging"]);

  std::map<std::string, double> p_double = {
      {"fluid_velocity", P_double["uf"]},
      {"density_fluid", P_double["density_fluid"]},
      {"density_solid", P_double["density_solid"]},
      {"specific_heat_fluid", P_double["specific_heat_fluid"]},
      {"specific_heat_solid", P_double["specific_heat_solid"]},
      {"conductivity_fluid", P_double["conductivity_fluid"]},
      {"conductivity_solid", P_double["conductivity_solid"]},
      {"temperature_hot", P_double["T_hot"]},
      {"temperature_cold", P_double["T_cold"]},
      {"temperature_reference", P_double["T0"]},
      {"heat_exchange", P_double["heat_exchange"]},
      {"porosity", P_double["porosity"]},
  };

  solver_ = std::make_shared<HeatStorage<Mesh>>(
      mesh, dt,
      p_double,
      nullptr, nullptr, scheduler,
      false, 0., 0. // given temperature
  );

  P_int.set("cells_number", static_cast<int>(mesh.GetNumCells()));

  content = {
      std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
          "x", mesh,
          [this](IdxCell idx) { return mesh.GetCenter(idx)[0]; })
      , std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
          "Tf", mesh,
          [this](IdxCell idx) { return solver_->GetFluidTemperature()[idx]; })
      , std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
          "Ts", mesh,
          [this](IdxCell idx) { return solver_->GetSolidTemperature()[idx]; })

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
      , std::make_shared<output::EntryScalarFunction<Scal>>(
          "state", [this](){ return static_cast<Scal>(solver_->GetState()); })
      , std::make_shared<output::EntryScalarFunction<Scal>>(
          "exergy_c_in", [this](){ return static_cast<Scal>(solver_->GetExergy().exergy_c_in); })
      , std::make_shared<output::EntryScalarFunction<Scal>>(
          "exergy_c_out", [this](){ return static_cast<Scal>(solver_->GetExergy().exergy_c_out); })
      , std::make_shared<output::EntryScalarFunction<Scal>>(
          "exergy_d_in", [this](){ return static_cast<Scal>(solver_->GetExergy().exergy_d_in); })
      , std::make_shared<output::EntryScalarFunction<Scal>>(
          "exergy_d_out", [this](){ return static_cast<Scal>(solver_->GetExergy().exergy_d_out); })
      , std::make_shared<output::EntryScalarFunction<Scal>>(
          "exergy_efficiency", [this](){ return static_cast<Scal>(solver_->GetExergy().efficiency); })
      , std::make_shared<output::EntryScalarFunction<Scal>>(
          "capacity_factor", [this](){ return static_cast<Scal>(solver_->GetExergy().capacity_factor); })
      , std::make_shared<output::EntryScalarFunction<Scal>>(
          "temperature_increase", [this](){ return static_cast<Scal>(solver_->GetExergy().temperature_increase); })
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
      content, P_string["filename_field"], mesh);

  session_scalar = std::make_shared<output::SessionPlainScalar<Scal>>(
      content_scalar, P_string["filename_scalar"]);

  last_frame_time_ = 0;
  last_frame_scalar_time_ = 0;
  session->Write(0., "initial");
  session_scalar->Write();


  if (flag("MMS")) {
    const Scal uf = P_double["MMS_fluid_velocity"];
    const Scal alpha = P_double["MMS_alpha"];
    const Scal wavenumber = P_double["MMS_wavenumber"];
    const Scal domain_length = P_double["MMS_domain_length"];

    using FuncTX = typename HeatStorage<Mesh>::FuncTX;
    FuncTX func_exact, func_rhs_fluid, func_rhs_solid;

    const std::string mmsfunc = P_string["MMS_exact_solution"];
    if (mmsfunc == "cos(kx)") {
      func_exact = [=](Scal, Scal x) { return std::cos(x * wavenumber); };
      func_rhs_fluid = [=](Scal, Scal x) {
        return -uf * wavenumber * std::sin(x * wavenumber) +
            alpha * sqr(wavenumber) * std::cos(x * wavenumber);
      };
      func_rhs_solid = [=](Scal, Scal x) {
        return alpha * sqr(wavenumber) * std::cos(x * wavenumber);
      };
    } else if (mmsfunc == "cos(kx^2)") {
      Scal modwavenumber = wavenumber / domain_length;
      func_exact = [=](Scal, Scal x) {
        return std::cos(sqr(x) * modwavenumber);
      };
      func_rhs_fluid = [=](Scal, Scal x) {
        return -uf * modwavenumber * 2. * x * std::sin(sqr(x) * modwavenumber) +
            alpha * (sqr(modwavenumber) * 4. * sqr(x) * std::cos(sqr(x) * modwavenumber) +
                    modwavenumber * 2. * std::sin(sqr(x) * modwavenumber));
      };
      func_rhs_solid = [=](Scal, Scal x) {
        return alpha * (sqr(modwavenumber) * 4. * sqr(x) * std::cos(sqr(x) * modwavenumber) +
                    modwavenumber * 2. * std::sin(sqr(x) * modwavenumber));
      };
    } else if (mmsfunc == "test") {
      func_exact = [=](Scal, Scal x) { return x; };
      func_rhs_fluid = [=](Scal, Scal) {
        return uf;
      };
      func_rhs_solid= [=](Scal, Scal) {
        return 0.;
      };
    } else {
      throw std::string("Unknown MMS_exact_solution: " + mmsfunc);
    }

    typename HeatStorage<Mesh>::TesterMms tester(
        P_int["MMS_mesh_initial"],
        P_int["MMS_num_stages"],
        P_int["MMS_factor"],
        domain_length,
        P_int["MMS_num_steps"],
        P_double["MMS_time_step"],
        P_double["MMS_step_threshold"],
        P_double["MMS_fluid_velocity"],
        P_double["MMS_alpha"],
        func_exact(0., 0.), // T_left is left boundary value from the exact solution
        0., // T_right acts as initial value (zero)
        func_rhs_fluid, func_exact,
        func_rhs_solid, func_exact // func_exact for solid expected to have integral T_right * length
    );
  }
}

template <class Mesh>
void hydro<Mesh>::step() {
  ex->timer_.Push("step");

  solver_->StartStep();
  solver_->CalcStep();
  solver_->FinishStep();

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
  const double first_frame_time = ecast(P_double("T_first_frame"));
  const double frame_duration = (total_time - first_frame_time) / max_frame_index;

  if (force || (!ecast(P_bool("no_mesh_output")) &&
      time >= first_frame_time &&
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
