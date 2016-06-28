/*******************************************************************
******************       CFD SOLVER 2D        **********************
********************************************************************/

#pragma once

#include "../control/experiment.hpp"
#include "../common/vect.hpp"
#include "mesh.hpp"
#include "mesh2d.hpp"
#include "mesh3d.hpp"
#include "output.hpp"
#include "output_tecplot.hpp"
#include "output_paraview.hpp"
#include "heat.hpp"
#include "advection.hpp"
#include <memory>
#include "fluid.hpp"
#include "chemistry.hpp"
#include "parallel.hpp"

#ifdef MPI_ENABLE
#include <mpi.h>
#endif


// TODO: Change parameters on events (e.g. certain time moments)

namespace hydro2D_uniform_MPI
{

template <class Vect>
Vect GetVect(const column<double>& v);

template <>
geom::Vect<double, 2> GetVect<geom::Vect<double, 2>>(const column<double>& v);

template <>
geom::Vect<double, 3> GetVect<geom::Vect<double, 3>>(const column<double>& v);

template <>
geom::Vect<float, 2> GetVect<geom::Vect<float, 2>>(const column<double>& v);

template <>
geom::Vect<float, 3> GetVect<geom::Vect<float, 3>>(const column<double>& v);


template <class Mesh>
class hydro_core : public TModule
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

  struct Flags {
    bool no_output = false;
  };
  Flags flags;

  Mesh mesh;
  Mesh outmesh;
  FieldCell<IdxCell> out_to_mesh_;
  std::shared_ptr<solver::
  FluidSolver<Mesh>> fluid_solver;
  std::shared_ptr<solver::
  AdvectionSolverMulti<Mesh, FieldFace<Scal>>> advection_solver;
  std::shared_ptr<solver::
  HeatSolver<Mesh>> heat_solver;
  output::Content content, content_scalar;
  std::shared_ptr<output::Session> session, session_scalar;
  FieldCell<Scal> fc_density, fc_density_smooth, fc_viscosity_smooth;
  FieldCell<Scal> fc_volume_source, fc_mass_source;
  std::vector<Scal> v_true_density, v_viscosity, v_conductivity;
  std::vector<FieldCell<Scal>> v_fc_mass_source;
  // Velocity slip is phase velocity relative to mixture velocity
  std::vector<FieldCell<Vect>> v_fc_velocity_slip;
  std::vector<FieldFace<Scal>> v_ff_volume_flux_slip;
  std::vector<FieldCell<Scal>> v_fc_volume_fraction;
  std::vector<FieldCell<Scal>> v_fc_true_density_target, v_fc_true_density;

  std::vector<const FieldFace<Scal>*> v_p_ff_volume_flux_slip;
  FieldCell<Vect> fc_force;
  FieldCell<Scal> fc_c1_smooth, fc_c1_laplacian;
  FieldCell<Scal> fc_radiation;
  geom::MapFace<std::shared_ptr<solver::ConditionFace>> mf_cond_radiation_shared;
  std::vector<geom::MapFace<std::shared_ptr<solver::ConditionFace>>>
  v_mf_partial_density_cond_;
  // TODO: GetBoundaryConditions() in solvers
  geom::MapFace<solver::ConditionFace*> mf_cond_radiation;
  geom::Rect<Vect> domain;
  geom::MapFace<std::shared_ptr<solver::ConditionFaceFluid>> fluid_cond;
  FieldCell<Scal> fc_conductivity, fc_temperature_source;

  double last_frame_time_;
  double last_frame_scalar_time_;

  // Parallel (MPI)
  std::shared_ptr<solver::ParallelTools<Mesh>> parallel;
  bool is_master, is_worker;

  const size_t num_phases;
  const geom::Range<size_t> phases;
  void UpdateFluidProperties();
  void CalcStat();
  std::shared_ptr<const solver::LinearSolverFactory>
  GetLinearSolverFactory(std::string linear_name,
                         std::string first_prefix = "");
  void InitMesh();
  void InitParallel();
  void InitFluidSolver();
  void InitAdvectionSolver();
  void InitRadiation();
  void InitHeatSolver();
  void InitOutput();
  void CalcPhasesVolumeFraction();
  std::vector<Scal> GetPhaseProperty(std::string name);
  std::vector<Scal> GetPhaseProperty(std::string name,
                                     const std::vector<bool>& v_enabled);
  std::vector<bool> GetPhasePropertyFlag(std::string name);
  void CalcPhaseVelocitySlip();
  void CalcPhasesTrueDensity();
  void CalcPhasesMassSource();

  FieldCell<Scal> GetVolumeAveraged(const std::vector<Scal>& v_value);
  void AppendAntidiffusionVolumeFlux();
  void CalcRadiation();
  void CalcForce();
  void CalcMixtureVolumeSource();
  template <class T>
  FieldCell<T> GetVolumeAveraged(const std::vector<FieldCell<T>>& v_fc_field);
  template <class T>
  FieldCell<T> GetSum(const std::vector<FieldCell<T>>& v_fc_field);

 public:
  hydro_core(TExperiment* _ex);
  ~hydro_core() {}
  void step();
  void write_results(bool force=false);
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

  template <class T>
  using FieldCell = geom::FieldCell<T>;
  template <class T>
  using FieldFace = geom::FieldFace<T>;
  template <class T>
  using FieldNode = geom::FieldNode<T>;

  std::shared_ptr<hydro_core<Mesh>> hydro_core_;

 public:
  hydro(TExperiment* _ex)
      : TExperiment_ref(_ex), TModule(_ex) {
    hydro_core_ = std::make_shared<hydro_core<Mesh>>(_ex);
  }
  ~hydro() {}
  void step() {
    hydro_core_->step();
  }
  void write_results(bool force=false) {
    hydro_core_->write_results(force);
  }
};

template <class Mesh>
std::shared_ptr<const solver::LinearSolverFactory>
hydro_core<Mesh>::GetLinearSolverFactory(std::string linear_name,
                                    std::string /*first_prefix*/) {
  if (linear_name == "lu") {
    return std::make_shared<const solver::LinearSolverFactory>(
        std::make_shared<const solver::LuDecompositionFactory>());
  } else if (linear_name == "lu_relaxed") {
    return std::make_shared<const solver::LinearSolverFactory>(
        std::make_shared<const solver::LuDecompositionRelaxedFactory>(
            P_double["lu_relaxed_tolerance"],
            P_int["lu_relaxed_num_iters_limit"],
            P_double["lu_relaxed_relaxation_factor"]));
  } else if (linear_name == "gauss_seidel") {
    return std::make_shared<const solver::LinearSolverFactory>(
        std::make_shared<const solver::GaussSeidelFactory>(
            P_double["lu_relaxed_tolerance"],
            P_int["lu_relaxed_num_iters_limit"],
            P_double["lu_relaxed_relaxation_factor"]));
  } /*else if (linear_name == "pardiso") {
    std::string second_prefix = "pardiso_";

    auto try_parameter = [this, &first_prefix, &second_prefix](
        std::string parameter, MKL_INT& value) {
      if (auto ptr = P_int(first_prefix + second_prefix + parameter)) {
        value = *ptr;
        return true;
      } else if (auto ptr = P_int(second_prefix + parameter)) {
        value = *ptr;
        return true;
      }
      return false;
    };

    MKL_INT mtype = 0;
    try_parameter("mtype", mtype);

    std::map<size_t, MKL_INT> iparm_map;
    for (size_t i = 0; i < 64; ++i) {
      MKL_INT value = 0;
      if (try_parameter(std::string("iparm_") + IntToStr(i), value)) {
        iparm_map[i] = value;
      }
    }
    return std::make_shared<const solver::LinearSolverFactory>(
        std::make_shared<const solver::PardisoFactory>(iparm_map, mtype));
  }*/
  throw std::runtime_error("Unknown linear solver '" + linear_name + "'");
}

template <class Mesh>
void hydro_core<Mesh>::InitMesh() {
  MIdx mesh_size;
  for (size_t i = 0; i < dim; ++i) {
    mesh_size[i] = P_int[std::string("N") + Direction(i).GetLetter()];
  }

  domain = geom::Rect<Vect>(GetVect<Vect>(P_vect["A"]),
                            GetVect<Vect>(P_vect["B"]));

  geom::InitUniformMesh(mesh, domain, mesh_size);

  P_int.set("cells_number", static_cast<int>(mesh.GetNumCells()));
}

template <class Mesh>
void hydro_core<Mesh>::InitParallel() {
  parallel = std::make_shared<solver::ParallelTools<Mesh>>(mesh);

  is_master = (parallel->GetRank() == 0);
  is_worker = !is_master;
  if (is_worker) {
    flags.no_output = true;
  }
}

template <class Mesh>
void hydro_core<Mesh>::InitFluidSolver() {
  // Parse linear solver parameters
  std::shared_ptr<const solver::LinearSolverFactory>
  p_linear_factory_velocity =
      GetLinearSolverFactory(P_string["linear_solver_velocity"], "velocity_");

  std::shared_ptr<const solver::LinearSolverFactory>
  p_linear_factory_pressure =
      GetLinearSolverFactory(P_string["linear_solver_pressure"], "pressure_");

  auto is_left_boundary = [this](IdxFace idxface) {
    return mesh.GetDirection(idxface) == Direction::i &&
        mesh.GetBlockFaces().GetMIdx(idxface)[0] == 0;
  };
  auto is_right_boundary = [this](IdxFace idxface) {
    return mesh.GetDirection(idxface) == Direction::i &&
        mesh.GetBlockFaces().GetMIdx(idxface)[0] ==
            mesh.GetBlockCells().GetDimensions()[0];
  };
  auto is_bottom_boundary = [this](IdxFace idxface) {
    return mesh.GetDirection(idxface) == Direction::j &&
        mesh.GetBlockFaces().GetMIdx(idxface)[1] == 0;
  };
  auto is_top_boundary = [this](IdxFace idxface) {
    return mesh.GetDirection(idxface) == Direction::j &&
        mesh.GetBlockFaces().GetMIdx(idxface)[1] ==
            mesh.GetBlockCells().GetDimensions()[1];
  };
  auto is_close_boundary = [this](IdxFace idxface) {
    return dim == 3 && mesh.GetDirection(idxface) == Direction::k &&
        mesh.GetBlockFaces().GetMIdx(idxface)[2] == 0;
  };
  auto is_far_boundary = [this](IdxFace idxface) {
    return dim == 3 && mesh.GetDirection(idxface) == Direction::k &&
        mesh.GetBlockFaces().GetMIdx(idxface)[2] ==
            mesh.GetBlockCells().GetDimensions()[2];
  };

  // Rigid box
  geom::Rect<Vect> rigid_box = geom::Rect<Vect>(
      GetVect<Vect>(P_vect["box_A"]), GetVect<Vect>(P_vect["box_B"]));

  // Initial conditions for velocity
  Vect initial_velocity(Vect::kZero);
  if (P_vect.exist("initial_velocity")) {
    initial_velocity = GetVect<Vect>(P_vect["initial_velocity"]);
  }
  FieldCell<Vect> fc_velocity_initial(mesh, initial_velocity);
  if (P_bool["deforming_velocity"]) {
    fc_velocity_initial = solver::GetDeformingVelocity(mesh);
  }

  // Boundary conditions for fluid
  for (auto idxface : mesh.Faces()) {
    if (!mesh.IsInner(idxface)) {
      if (is_top_boundary(idxface)) {
        fluid_cond[idxface] =
            solver::Parse(P_string["condition_top"], idxface, mesh);
      } else if (is_bottom_boundary(idxface)) {
        fluid_cond[idxface] =
            solver::Parse(P_string["condition_bottom"], idxface, mesh);
      } else if (is_left_boundary(idxface)) {
        fluid_cond[idxface] =
            solver::Parse(P_string["condition_left"], idxface, mesh);
      } else if (is_right_boundary(idxface)) {
        fluid_cond[idxface] =
            solver::Parse(P_string["condition_right"], idxface, mesh);
      } else if (is_close_boundary(idxface)) {
        fluid_cond[idxface] =
            solver::Parse(P_string["condition_close"], idxface, mesh);
      } else if (is_far_boundary(idxface)) {
        fluid_cond[idxface] =
            solver::Parse(P_string["condition_far"], idxface, mesh);
      } else {
        throw std::runtime_error("Unhandled boundary");
      }
    } else {
      // Rigid box boundary conditions
      IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
      IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
      auto is_inside = [this, &rigid_box](IdxCell idxcell) {
        return  !idxcell.IsNone() &&
            rigid_box.IsInside(mesh.GetCenter(idxcell));
      };
      if (is_inside(cm) != is_inside(cp)) {
        fluid_cond[idxface] =
            solver::Parse("wall 0 0 0", idxface, mesh);
      }
    }
  }

  // Exclude rigid box cells
  std::vector<IdxCell> exclusion_list;
  for (auto idxcell : mesh.Cells()) {
    if (rigid_box.IsInside(mesh.GetCenter(idxcell))) {
      exclusion_list.push_back(idxcell);
    }
  }
  mesh.ExcludeCells(exclusion_list);

  // Remove conditions on excluded boundaries
  for (auto idxface : mesh.Faces()) {
    if (mesh.IsExcluded(idxface)) {
      fluid_cond.erase(idxface);
    }
  }

  // Checkout for unhandled boundaries
  for (auto idxface : mesh.Faces()) {
    if (!mesh.IsInner(idxface) && !mesh.IsExcluded(idxface)) {
      if (!fluid_cond.find(idxface)) {
        throw std::runtime_error("Unhandled boundary");
      }
    }
  }

  // Cell conditions for fluid
  geom::MapCell<std::shared_ptr<solver::ConditionCellFluid>> mc_cond_fluid;

  // Fixed pressure point
  if (auto* ptr_point = P_vect("pressure_fixed_point")) {
     Vect point = GetVect<Vect>(*ptr_point);
     Scal value = 0.;
     if (double* ptr_value = P_double("pressure_fixed_value")) {
       value = *ptr_value;
     }
     mc_cond_fluid[mesh.FindNearestCell(point)] =
         std::make_shared<solver::fluid_condition::
         GivenPressureFixed<Mesh>>(value);
  }

  fluid_solver = std::make_shared<solver::FluidSimple<Mesh>>(
    mesh, fc_velocity_initial, fluid_cond, mc_cond_fluid,
    P_double["velocity_relaxation_factor"],
    P_double["pressure_relaxation_factor"],
    P_double["rhie_chow_factor"],
    &fc_density_smooth, &fc_viscosity_smooth, &fc_force,
    &fc_volume_source, &fc_mass_source,
    0., dt, *p_linear_factory_velocity, *p_linear_factory_pressure,
    P_double["convergence_tolerance"], P_int["num_iterations_limit"],
    &ex->timer_, P_bool["time_second_order"], P_bool["simpler"],
    P_bool["force_geometric_average"], P_double["guess_extrapolation"]);

  /*fluid_solver = std::make_shared<solver::FluidSimpleParallel<Mesh>>(
      mesh, fc_velocity_initial, fluid_cond, mc_cond_fluid,
      P_double["velocity_relaxation_factor"],
      P_double["pressure_relaxation_factor"],
      &fc_density, &fc_viscosity, &fc_force,
      &fc_volume_source, &fc_mass_source,
      0., dt, linear_id_velocity, linear_id_pressure,
      P_double["convergence_tolerance"], P_int["num_iterations_limit"],
      &ex->timer_, P_int["num_local_mesh"], P_int["overlap_width"]);*/
}

template <class Mesh>
void hydro_core<Mesh>::InitAdvectionSolver() {
  // Properties of phases
  v_true_density = GetPhaseProperty("density_");
  v_viscosity = GetPhaseProperty("viscosity_");
  v_conductivity = GetPhaseProperty("conductivity_");

  v_fc_mass_source.resize(num_phases);
  v_ff_volume_flux_slip.resize(num_phases);
  v_fc_volume_fraction.resize(num_phases);
  v_fc_velocity_slip.resize(num_phases);
  v_p_ff_volume_flux_slip.resize(num_phases);
  v_fc_true_density_target.resize(num_phases);
  v_fc_true_density.resize(num_phases);

  std::vector<FieldCell<Scal>> v_fc_partial_density_initial(num_phases);

  for (auto i : phases) {
    v_fc_mass_source[i].Reinit(mesh, 0.);
    v_fc_velocity_slip[i].Reinit(mesh);
    v_ff_volume_flux_slip[i].Reinit(mesh);
    v_fc_volume_fraction[i].Reinit(mesh);
    v_fc_true_density_target[i].Reinit(mesh, v_true_density[i]);
    v_fc_true_density[i].Reinit(mesh, v_true_density[i]);
    v_fc_partial_density_initial[i].Reinit(
        mesh, v_true_density[i] *
        ecast(P_double("initial_volume_fraction_" + IntToStr(i))));
    v_p_ff_volume_flux_slip[i] = &v_ff_volume_flux_slip[i];
  }

  // Initial partial density for phases [1, num_phases)
  geom::Rect<Vect> block1(GetVect<Vect>(P_vect["A1"]),
                    GetVect<Vect>(P_vect["B1"]));
  geom::Rect<Vect> block2(GetVect<Vect>(P_vect["A2"]),
                    GetVect<Vect>(P_vect["B2"]));
  for (auto idx : mesh.Cells()) {
    Vect x = mesh.GetCenter(idx);
    if (block2.IsInside(x)) {
      v_fc_partial_density_initial[2][idx] = v_fc_true_density[2][idx];
    } else if (block1.IsInside(x)) {
      v_fc_partial_density_initial[1][idx] = v_fc_true_density[1][idx];
    }
  }

  for (size_t i = 1; i < num_phases; ++i) {
    v_fc_partial_density_initial[i] =
        solver::GetSmoothField(v_fc_partial_density_initial[i], mesh,
                               P_int["initial_volume_fraction_smooth_times"]);
  }

  // Initial partial density of 0-phase compensates the others
  for (auto idx : mesh.Cells()) {
    Scal volume_sum = 0.;
    for (size_t i = 1; i < num_phases; ++i) {
      volume_sum +=
          v_fc_partial_density_initial[i][idx] / v_fc_true_density[i][idx];
    }
    v_fc_partial_density_initial[0][idx] =
        (1. - volume_sum) * v_fc_true_density[0][idx];
  }

  // Boundary conditions for concentration
  v_mf_partial_density_cond_.resize(num_phases);
  for (auto i : phases) {
    auto& mf_cond = v_mf_partial_density_cond_[i];
    for (auto it = fluid_cond.cbegin(); it != fluid_cond.cend(); ++it) {
      IdxFace idxface = it->GetIdx();
      if (dynamic_cast<solver::fluid_condition::
          Inlet<Mesh>*>(it->GetValue().get())) {
        IdxCell idxcell = mesh.GetNeighbourCell(idxface,
            mesh.GetValidNeighbourCellId(idxface));
        mf_cond[idxface] =
            std::make_shared<
            solver::ConditionFaceValueFixed<Scal>>(
                v_fc_partial_density_initial[i][idxcell]);
      } else {
        mf_cond[idxface] =
            std::make_shared<
            solver::ConditionFaceDerivativeFixed<Scal>>(0.);
      }
    }
  }

  // Make vector of pointers to mass sources
  std::vector<const geom::FieldCell<Scal>*> v_p_fc_mass_source(num_phases);
  for (auto i : phases) {
    v_p_fc_mass_source[i] = &v_fc_mass_source[i];
  }

  std::string advection_solver_name = P_string["advection_solver"];
  if (advection_solver_name == "tvd") {
    advection_solver = std::make_shared<solver::
        AdvectionSolverMultiExplicit<Mesh, FieldFace<Scal>>>(
        mesh, v_fc_partial_density_initial, v_mf_partial_density_cond_,
        &fluid_solver->GetVolumeFlux(solver::Layers::iter_curr),
        v_p_ff_volume_flux_slip,
        v_p_fc_mass_source,
        0., dt * P_double["advection_dt_factor"], P_bool["tvd_split"]);
  } else if (advection_solver_name == "pic") {
    advection_solver = std::make_shared<solver::
        AdvectionSolverMultiParticles<Mesh, FieldFace<Scal>>>(
        mesh, v_fc_partial_density_initial, v_mf_partial_density_cond_,
        &fluid_solver->GetVolumeFlux(solver::Layers::iter_curr),
        v_p_fc_mass_source,
        0., dt * P_double["advection_dt_factor"],
        P_double["spawning_gap"], P_double["particle_radius"],
        P_int["min_num_particles"], P_int["max_num_particles"],
        P_double["back_relaxation_factor"]);
  } else {
    std::runtime_error(
        "Unknown advection solver = '" + advection_solver_name + "'");
  }
}

template <class Mesh>
void hydro_core<Mesh>::InitRadiation() {
  // Boundary conditions for radiation
  Vect radiation_box_lb = GetVect<Vect>(P_vect["radiation_box_lb"]);
  Vect radiation_box_rt = GetVect<Vect>(P_vect["radiation_box_rt"]);
  Scal radiation_intensity = P_double["radiation_intensity"];
  Vect radiation_direction = GetVect<Vect>(P_vect["radiation_direction"]);
  for (auto idxface : mesh.Faces()) {
    Vect xface = mesh.GetCenter(idxface);
    if (!mesh.IsInner(idxface)) {
      if (radiation_box_lb < xface && xface < radiation_box_rt) {
        mf_cond_radiation_shared[idxface] =
            std::make_shared<
            solver::ConditionFaceValueFixed<Scal>>(radiation_intensity);
      } else {
        if (radiation_direction.dot(
                mesh.GetNormal(idxface) *
                (mesh.GetValidNeighbourCellId(idxface) == 0 ? 1. : -1.)) > 0.) {

          mf_cond_radiation_shared[idxface] =
              std::make_shared<
              solver::ConditionFaceDerivativeFixed<Scal>>(0.);
        } else {
          mf_cond_radiation_shared[idxface] =
              std::make_shared<
              solver::ConditionFaceValueFixed<Scal>>(0.);

        }
      }
    }
  }
  fc_radiation.Reinit(mesh, 0.);
}

template <class Mesh>
void hydro_core<Mesh>::InitHeatSolver() {
  std::shared_ptr<const solver::LinearSolverFactory>
  p_linear_factory_heat =
      GetLinearSolverFactory(P_string["linear_solver_heat"], "heat_");

  geom::FieldCell<Scal> fc_temperature_initial(
      mesh, P_double["temperature_initial"]);

  // Boundary conditions for temperature
  geom::MapFace<std::shared_ptr<solver::ConditionFace>>
  mf_temperature_cond;
  geom::Rect<Vect> heat_box(GetVect<Vect>(P_vect["heat_box_lb"]),
                            GetVect<Vect>(P_vect["heat_box_rt"]));
  for (auto idxface : mesh.Faces()) {
    if (!mesh.IsInner(idxface)) {
      if (heat_box.IsInside(mesh.GetCenter(idxface))) {
        mf_temperature_cond[idxface] =
            std::make_shared<solver::
            ConditionFaceValueFixed<Scal>>(P_double["heat_box_temperature"]);
      } else {
        mf_temperature_cond[idxface] =
            std::make_shared<solver::
            ConditionFaceDerivativeFixed<Scal>>(0.);
      }
    }
  }

  heat_solver = std::make_shared<solver::HeatSolver<Mesh>>(
      mesh, fc_temperature_initial,
      mf_temperature_cond,
      P_double["heat_relaxation_factor"],
      &fc_conductivity, &fc_temperature_source,
      &fluid_solver->GetVolumeFlux(solver::Layers::iter_curr),
      0., dt, *p_linear_factory_heat,
      0.5, 1, P_bool["time_second_order_heat"]);
}

template <class Mesh>
void hydro_core<Mesh>::InitOutput() {
  if (flags.no_output) {
    return;
  }

  // Create output mesh
  MIdx output_factor;
  for (size_t i = 0; i < dim; ++i) {
    output_factor[i] = P_int[std::string("output_factor_") +
                         Direction(i).GetLetter()];
  }
  MIdx mesh_size = mesh.GetBlockCells().GetDimensions();
  geom::InitUniformMesh(outmesh, domain, mesh_size * output_factor);
  out_to_mesh_.Reinit(outmesh);
  for (auto outidx : outmesh.Cells()) {
    MIdx outmidx = outmesh.GetBlockCells().GetMIdx(outidx);
    out_to_mesh_[outidx] =
        mesh.GetBlockCells().GetIdx(outmidx / output_factor);
  }

  // TODO: EntryFunctionField
  output::Content content_pool = {
      std::make_shared<output::EntryFunction<Scal, IdxNode, Mesh>>(
          "x", outmesh, [this](IdxNode idx) { return outmesh.GetNode(idx)[0]; })
      , std::make_shared<output::EntryFunction<Scal, IdxNode, Mesh>>(
          "y", outmesh, [this](IdxNode idx) { return outmesh.GetNode(idx)[1]; })
      , std::make_shared<output::EntryFunction<Scal, IdxNode, Mesh>>(
          "z", outmesh, [this](IdxNode idx) {
              return dim == 2 ? 0. : outmesh.GetNode(idx)[2]; })
      , std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
          "velocity_x", outmesh, [this](IdxCell idx) {
              return fluid_solver->GetVelocity()[out_to_mesh_[idx]][0]; })
      , std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
          "velocity_y", outmesh, [this](IdxCell idx) {
              return fluid_solver->GetVelocity()[out_to_mesh_[idx]][1]; })
      , std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
          "velocity_z", outmesh, [this](IdxCell idx) {
              return dim == 2 ? 0. :
                  fluid_solver->GetVelocity()[out_to_mesh_[idx]][2]; })
      , std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
          "pressure", outmesh, [this](IdxCell idx) {
              return fluid_solver->GetPressure()[out_to_mesh_[idx]]; })
      , std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
          "density", outmesh, [this](IdxCell idx) {
              return fc_density[out_to_mesh_[idx]]; })
      , std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
          "viscosity", outmesh, [this](IdxCell idx) {
              return fc_viscosity_smooth[out_to_mesh_[idx]]; })
      , std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
          "radiation", outmesh, [this](IdxCell idx) {
              return fc_radiation[out_to_mesh_[idx]]; })
      , std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
          "volume_source", outmesh, [this](IdxCell idx) {
              return fc_volume_source[out_to_mesh_[idx]]; })
      , std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
          "temperature", outmesh, [this](IdxCell idx) {
              return heat_solver->GetTemperature()[out_to_mesh_[idx]]; })
      , std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
          "excluded", outmesh, [this](IdxCell idx) {
              return mesh.IsExcluded(out_to_mesh_[idx]) ? 1. : 0.; })
      , std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
          "divergence", outmesh, [this](IdxCell idx) {
              auto idxcell = out_to_mesh_[idx];
              Scal sum = 0.;
              for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
                sum += fluid_solver->
                    GetVolumeFlux()[mesh.GetNeighbourFace(idxcell, i)] *
                    mesh.GetOutwardFactor(idxcell, i);
              }
              return sum / mesh.GetVolume(idxcell);})
  };

  // Mass source output field
  for (auto i : phases) {
    std::string phase = IntToStr(i);
    content_pool.push_back(
        std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
                "mass_source_" + phase, outmesh, [this, i](IdxCell idx) {
                    return v_fc_mass_source[i][out_to_mesh_[idx]]; }));
  }

  // Mass fraction output field
  for (auto i : phases) {
    std::string phase = IntToStr(i);
    content_pool.push_back(
        std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
            "mass_fraction_" + phase, outmesh, [this, i](IdxCell idx) {
                return advection_solver->GetField(i)[out_to_mesh_[idx]] /
                    fc_density[out_to_mesh_[idx]]; }));
  }

  // Density output field
  for (auto i : phases) {
    std::string phase = IntToStr(i);
    content_pool.push_back(
        std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
            "density_" + phase, outmesh, [this, i](IdxCell idx) {
                return v_fc_true_density[i][out_to_mesh_[idx]]; }));
  }

  // Target density output field
  for (auto i : phases) {
    std::string phase = IntToStr(i);
    content_pool.push_back(
        std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
            "target_density_" + phase, outmesh, [this, i](IdxCell idx) {
                return v_fc_true_density_target[i][out_to_mesh_[idx]]; }));
  }

  // Partial density output field
  for (auto i : phases) {
    std::string phase = IntToStr(i);
    content_pool.push_back(
        std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
            "partial_density_" + phase, outmesh, [this, i](IdxCell idx) {
                return advection_solver->GetField(i)[out_to_mesh_[idx]]; }));
  }

  // Volume fraction output field
  for (auto i : phases) {
    std::string phase = IntToStr(i);
    content_pool.push_back(
        std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
            "volume_fraction_" + phase, outmesh, [this, i](IdxCell idx) {
                return v_fc_volume_fraction[i][out_to_mesh_[idx]]; }));
  }

  for (auto& entry : content_pool) {
    if (P_bool[std::string("output_") + entry->GetName()]) {
      content.push_back(entry);
    }
  }

  // Scalar output
  auto P = [this](std::string entry, std::string parameter) {
    return std::make_shared<output::EntryScalarFunction<Scal>>(
        entry, [this, parameter](){ return P_double[parameter]; });
  };
  content_scalar = {
      P("time", "t"),
      std::make_shared<output::EntryScalarFunction<Scal>>(
          "num_iters", [this](){
              return P_int["s"];
          }),
      std::make_shared<output::EntryScalarFunction<Scal>>(
          "iter_diff_velocity", [this](){
              return fluid_solver->GetConvergenceIndicator();
          })
  };


  // stat_mass scalar output
  for (auto i : phases) {
    std::string phase = IntToStr(i);
    content_scalar.push_back(P("mass_" + phase,
                               "stat_mass_" + phase));
    content_scalar.push_back(P("mass_in_" + phase,
                               "stat_mass_in_" + phase));
    content_scalar.push_back(P("mass_out_" + phase,
                               "stat_mass_out_" + phase));
  }

  // stat_volume scalar output
  for (auto i : phases) {
    std::string phase = IntToStr(i);
    content_scalar.push_back(P("volume_" + phase,
                               "stat_volume_" + phase));
    content_scalar.push_back(P("volume_in_" + phase,
                               "stat_volume_in_" + phase));
    content_scalar.push_back(P("volume_out_" + phase,
                               "stat_volume_out_" + phase));
  }

  if (!P_string.exist(_plt_title)) {
    P_string.set(_plt_title, P_string[_exp_name]);
  }
  if (!P_string.exist("filename_field")) {
    P_string.set("filename_field", P_string[_exp_name] + ".field");
  }
  if (!P_string.exist("filename_scalar")) {
    P_string.set("filename_scalar", P_string[_exp_name] + ".scalar");
  }

  std::string field_output_format = P_string["field_output_format"];
  if (field_output_format == "tecplot_binary") {
    session = std::make_shared<output::SessionTecplotBinaryStructured<Mesh>>(
        content, P_string[_plt_title], outmesh,
        P_string["filename_field"] + ".plt");
  } else if (field_output_format == "paraview") {
    session = std::make_shared<output::SessionParaviewStructured<Mesh>>(
        content, P_string[_plt_title],
        P_string["filename_field"], outmesh);
  }

  session_scalar = std::make_shared<output::SessionTecplotAsciiScalar<Scal>>(
      content_scalar, P_string[_plt_title],
      P_string[_plt_title], P_string["filename_scalar"] + ".dat");

  last_frame_time_ = 0;
  last_frame_scalar_time_ = 0;
}

template <class Mesh>
hydro_core<Mesh>::hydro_core(TExperiment* _ex)
    : TExperiment_ref(_ex), TModule(_ex)
    , num_phases(P_int["num_phases"])
    , phases(0, num_phases)
{
  flags.no_output = P_bool["no_output"];

  P_int.set("last_s", 0);
  P_double.set("last_R", 0);
  P_double.set("last_Rn", 0);

  P_int.set("s_sum", 0);
  P_int.set("s_max", 0);
  P_int.set("s", 0);
  P_int.set("current_frame", 0);
  P_int.set("current_frame_scalar", 0);

  InitMesh();
  InitParallel();
  InitFluidSolver();
  InitAdvectionSolver();
  InitRadiation();
  InitHeatSolver();

  UpdateFluidProperties();
  CalcStat();

  InitOutput();

  if (!flags.no_output) {
    if (P_int["max_frame_index"] > 0) {
      session->Write(0., P_string[_plt_title] + ":0");
    }
    session_scalar->Write();
  }
}

template <class Scal>
void Limit(Scal& source, Scal concentration, Scal density, Scal dt) {
  source = std::max(source, -concentration * density / dt);
  source = std::min(source, (1 - concentration) * density / dt);
}

template <class Mesh>
void hydro_core<Mesh>::CalcPhasesVolumeFraction() {
  // Calc and normalize volume fraction for phases
  for (auto idxcell : mesh.Cells()) {
    Scal sum = 0.;
    for (auto i : phases) {
      v_fc_volume_fraction[i][idxcell] =
          advection_solver->GetField(i)[idxcell] /
          v_fc_true_density[i][idxcell];
      sum += v_fc_volume_fraction[i][idxcell];
    }
  }
}

template <class Mesh>
auto hydro_core<Mesh>::GetPhaseProperty(std::string name) -> std::vector<Scal> {
  std::vector<Scal> v_property(num_phases);
  for (auto i : phases) {
    v_property[i] = P_double[name + IntToStr(i)];
  }
  return v_property;
}

template <class Mesh>
auto hydro_core<Mesh>::GetPhaseProperty(std::string name,
                                   const std::vector<bool>& v_enabled)
    -> std::vector<Scal> {
  std::vector<Scal> v_property(num_phases, 0);
  for (auto i : phases) {
    if (v_enabled[i]) {
      v_property[i] = P_double[name + IntToStr(i)];
    }
  }
  return v_property;
}

template <class Mesh>
auto hydro_core<Mesh>::GetPhasePropertyFlag(std::string name) -> std::vector<bool> {
  std::vector<bool> v_property(num_phases);
  for (auto i : phases) {
    v_property[i] = flag(name + IntToStr(i));
  }
  return v_property;
}

template <class Mesh>
void hydro_core<Mesh>::CalcPhaseVelocitySlip() {
  Vect gravity = GetVect<Vect>(P_vect["gravity"]);
  auto v_enable_settling = GetPhasePropertyFlag("enable_settling_");
  auto v_bubble_radius = GetPhaseProperty("bubble_radius_", v_enable_settling);

  bool velocity_is_carrier = flag("velocity_is_carrier");

  for (auto idxcell : mesh.Cells()) {
    // Calc phase velocities relative to carrier velocity
    // using the Stokes law
    std::vector<Vect> v_velocity_relative_to_carrier(num_phases, Vect::kZero);
    for (auto i : phases) {
      if (v_enable_settling[i]) {
        Scal phase_c = v_fc_volume_fraction[i][idxcell];
        Scal mixture_density = fc_density[idxcell];
        Scal mixture_viscosity = fc_viscosity_smooth[idxcell];
        Scal phase_density = v_true_density[i];
        v_velocity_relative_to_carrier[i] =
            phase_c < 0.01 || phase_c > 0.99 ?
                Vect::kZero :
                gravity * (phase_density - mixture_density) *
                sqr(v_bubble_radius[i]) / (18. * mixture_viscosity);
      }
    }

    if (!velocity_is_carrier) {
      // Calc carrier velocity relative to mixture velocity
      Vect carrier_velocity_relative_to_mixture = Vect::kZero;
      for (auto i : phases) {
        carrier_velocity_relative_to_mixture +=
            v_velocity_relative_to_carrier[i] *
            v_fc_volume_fraction[i][idxcell];
      }

      // Calc phase velocities relative to mixture velocity
      for (auto i : phases) {
        v_fc_velocity_slip[i][idxcell] =
            v_velocity_relative_to_carrier[i] -
            carrier_velocity_relative_to_mixture;
      }
    } else {
      for (auto i : phases) {
        v_fc_velocity_slip[i][idxcell] = v_velocity_relative_to_carrier[i];
      }
    }
  }

  // Calc volume flux slip on faces
  for (auto idxface : mesh.Faces()) {
    if (mesh.IsInner(idxface)) {
      for (auto i : phases) {
        v_ff_volume_flux_slip[i][idxface] =
            solver::GetInterpolatedInner(
                v_fc_velocity_slip[i], idxface, mesh).dot(
                    mesh.GetSurface(idxface));
      }
    } else {
      for (auto i : phases) {
        v_ff_volume_flux_slip[i][idxface] = 0.;
      }
    }
  }

  if (!velocity_is_carrier) {
    // Correct volume flux slip to make average flux zero
    for (auto idxface : mesh.Faces()) {
      if (mesh.IsInner(idxface)) {
        Scal aver = 0.;
        for (auto i : phases) {
          Scal c = solver::GetInterpolatedInner(
              v_fc_volume_fraction[i], idxface, mesh);
          aver += c * v_ff_volume_flux_slip[i][idxface];
        }
        for (auto i : phases) {
          v_ff_volume_flux_slip[i][idxface] -= aver;
        }
      }
    }
  } else {
    // Add carrier volume flux
    for (auto idxcell : mesh.Cells()) {
      Scal volume_fraction_defect = 1.;
      for (auto i : phases) {
        volume_fraction_defect -=
            advection_solver->GetField(i)[idxcell] /
            v_fc_true_density[i][idxcell];
      }
      fc_volume_source[idxcell] = -volume_fraction_defect /
          fluid_solver->GetTimeStep();
    }
  }
}

// TODO: Resource manager
// Allow one to specify attributes like REQUIRE(pressure)

template <class Mesh>
void hydro_core<Mesh>::CalcPhasesTrueDensity() {
  if (P_bool["compressible_enable"]) {
    // Calc target density
    for (auto i : phases) {
      std::string phase = IntToStr(i);
      Scal rate = P_double["temperature_expansion_rate_" + phase];
      Scal base = P_double["temperature_expansion_base_" + phase];
      for (auto idxcell : mesh.Cells()) {
        v_fc_true_density_target[i][idxcell] = v_true_density[i]
            * (1. - rate * (heat_solver->GetTemperature()[idxcell] - base));
      }
    }
    // Calc true density from target density
    // and normalize the value to preserve correct volume fraction
    for (auto idxcell : mesh.Cells()) {
      Scal alpha = 0.;
      for (auto i : phases) {
        alpha += advection_solver->GetField(i)[idxcell]
            / v_fc_true_density_target[i][idxcell];
      }
      for (auto i : phases) {
        v_fc_true_density[i][idxcell] = v_fc_true_density_target[i][idxcell]
            * alpha;
      }
    }
  } else {
    // Just keep initial true density
  }
}

template <class Mesh>
void hydro_core<Mesh>::CalcPhasesMassSource() {
  std::vector<Scal> v_molar_mass = GetPhaseProperty("molar_");
  std::shared_ptr<solver::Kinetics<Mesh>> chem_solver;

  std::string chemistry = P_string["chemistry"];
  auto P = [this](std::string name) {
    return P_double["chem_" + name];
  };
  if (chemistry == "steady") {
    chem_solver = std::make_shared<solver::
        KineticsSteady<Mesh>>(
            P("intensity"), v_molar_mass, mesh);
  } else if (chemistry == "Nadirov") {
    chem_solver = std::make_shared<solver::
        KineticsNadirov<Mesh>>(
            v_molar_mass,
            P("k_P"), P("k"), P("I_0"), P("E_T"),
            P("k_B"), P("T"), P("G"), P("D"), mesh);
  } else if (chemistry == "steady_radiation") {
    chem_solver = std::make_shared<solver::
        KineticsSteadyRadiation<Mesh>>(
            P("intensity"), v_molar_mass, fc_radiation, mesh);
  } else {
    // TODO: Quotemarks string wrapper function
    std::runtime_error(
        "UpdateFluidProperties: Unknown chemistry = '" + chemistry + "'");
  }
  Vect block2_lb = GetVect<Vect>(P_vect["reaction_zone_lb"]);
  Vect block2_rt = GetVect<Vect>(P_vect["reaction_zone_rt"]);
  for (auto idxcell : mesh.Cells()) {
    Vect x = mesh.GetCenter(idxcell);
    bool in_reaction_zone = (x <= block2_rt && block2_lb <= x);
    Scal factor = in_reaction_zone ? 1. : 0.;

    std::vector<Scal> molar_concentration(num_phases);
    for (auto i : phases) {
      molar_concentration[i] =
          v_fc_volume_fraction[i][idxcell] *
          v_true_density[i] / v_molar_mass[i];
    }

    auto reaction_rate =
        chem_solver->GetReactionRate(idxcell, molar_concentration);

    for (auto i : phases) {
      v_fc_mass_source[i][idxcell] =
          factor * v_molar_mass[i] * reaction_rate[i];
    }
  }

  // Limit mass sources to prevent negative concentration values
//  for (auto idxcell : mesh.Cells()) {
//    for (auto i : phases) {
//      Limit(v_fc_mass_source[i][idxcell],
//            v_fc_volume_fraction[i][idxcell],
//            v_true_density[i], Scal(fluid_solver->GetTimeStep()));
//    }
//  }
}

// TODO: Add feature "requre smth" to indicate which fields should be known

template <class Mesh>
template <class T>
auto hydro_core<Mesh>::GetVolumeAveraged(
    const std::vector<FieldCell<T>>& v_fc_field) -> FieldCell<T> {
  FieldCell<T> res(mesh, T(0));
  for (auto i : phases) {
    for (auto idxcell : mesh.Cells()) {
      res[idxcell] +=
          v_fc_field[i][idxcell] * v_fc_volume_fraction[i][idxcell];
    }
  }
  return res;
}

template <class Mesh>
auto hydro_core<Mesh>::GetVolumeAveraged(
    const std::vector<Scal>& v_value) -> FieldCell<Scal> {
  FieldCell<Scal> res(mesh, 0);
  for (auto i : phases) {
    for (auto idxcell : mesh.Cells()) {
      res[idxcell] +=
          v_value[i] * v_fc_volume_fraction[i][idxcell];
    }
  }
  return res;
}

template <class Mesh>
template <class T>
auto hydro_core<Mesh>::GetSum(
    const std::vector<FieldCell<T>>& v_fc_field) -> FieldCell<T> {
  FieldCell<T> res(mesh, T(0));
  for (auto i : phases) {
    for (auto idxcell : mesh.Cells()) {
      res[idxcell] += v_fc_field[i][idxcell] ;
    }
  }
  return res;
}

// TODO: Split templates into files

template<class Mesh>
void hydro_core<Mesh>::AppendAntidiffusionVolumeFlux() {
  if (double* factor = P_double("antidiffusion_factor")) {
    std::vector<geom::FieldCell<Vect> > v_fc_antidiffusion(num_phases);
    Scal antidiffusion_factor = *factor;
    for (auto i : phases) {
      v_fc_antidiffusion[i] = solver::Gradient(
          solver::Interpolate(advection_solver->GetField(i),
                              v_mf_partial_density_cond_[i], mesh),
          mesh);
      for (auto idxcell : mesh.Cells()) {
        Scal c = v_fc_volume_fraction[i][idxcell];
        v_fc_antidiffusion[i][idxcell] *= -antidiffusion_factor * c * (1. - c);
      }
    }
    // Add anti-diffusion terms to velocity split fields
    for (auto idxface : mesh.Faces()) {
      if (mesh.IsInner(idxface)) {
        for (auto i : phases) {
          v_ff_volume_flux_slip[i][idxface] += solver::GetInterpolatedInner(
              v_fc_antidiffusion[i], idxface, mesh).dot(
              mesh.GetSurface(idxface));
        }
      }
    }
  }
}

template<class Mesh>
void hydro_core<Mesh>::CalcRadiation() {
  // Calc radiation field
  if (flag("radiation_enable")) {
    Vect radiation_direction = GetVect<Vect>(P_vect["radiation_direction"]);
    radiation_direction /= radiation_direction.norm();
    FieldCell<Scal> fc_absorption_rate = GetVolumeAveraged(
        GetPhaseProperty("absorption_rate_"));
    ex->timer_.Push("fluid_properties.radiation");
    fc_radiation = solver::CalcRadiationField(mesh, mf_cond_radiation_shared,
                                              fc_absorption_rate,
                                              radiation_direction,
                                              fc_radiation);
    ex->timer_.Pop();
  }
}

template<class Mesh>
void hydro_core<Mesh>::CalcForce() {
  Vect gravity = GetVect<Vect>(P_vect["gravity"]);
  fc_force.Reinit(mesh);
  for (auto idxcell : mesh.Cells()) {
    fc_force[idxcell] = gravity * fc_density[idxcell];
  }
  fc_force = solver::GetSmoothField(fc_force, mesh,
                                    P_int["force_smooth_times"]);
}

template<class Mesh>
void hydro_core<Mesh>::CalcMixtureVolumeSource() {
  fc_volume_source.Reinit(mesh, 0.);
  for (auto idxcell : mesh.Cells()) {
    for (auto i : phases) {
      fc_volume_source[idxcell] += v_fc_mass_source[i][idxcell]
          / v_fc_true_density[i][idxcell];
    }
  }
  if (P_bool["compressible_enable"]) {
    // Correct the volume source to account for target density
    Scal compressible_relaxation = P_double["compressible_relaxation"];
    for (auto idxcell : mesh.Cells()) {
      Scal alpha = 0.;
      for (auto i : phases) {
        alpha += advection_solver->GetField(i)[idxcell]
            / v_fc_true_density_target[i][idxcell];
      }
      fc_volume_source[idxcell] += (alpha - 1.) / fluid_solver->GetTimeStep()
          * compressible_relaxation;
    }
  }
}

// TODO: Make alias v_fc_density for advection_solver->GetField(i)

template <class Mesh>
void hydro_core<Mesh>::UpdateFluidProperties() {
  CalcPhasesTrueDensity();
  CalcPhasesVolumeFraction();
  CalcPhasesMassSource();

  fc_density = GetVolumeAveraged(v_fc_true_density);
  fc_density_smooth = solver::GetSmoothField(
      fc_density, mesh, P_int["density_smooth_times"]);

  fc_viscosity_smooth = solver::GetSmoothField(
      GetVolumeAveraged(v_viscosity), mesh, P_int["viscosity_smooth_times"]);

  fc_mass_source = GetSum(v_fc_mass_source);
  CalcMixtureVolumeSource();

  CalcForce();

  CalcPhaseVelocitySlip();

  AppendAntidiffusionVolumeFlux();

  CalcRadiation();

  fc_conductivity = GetVolumeAveraged(v_conductivity);
  fc_temperature_source.Reinit(mesh, 0.);
}

template <class Mesh>
void hydro_core<Mesh>::CalcStat() {
  for (auto i : phases) {
    Scal density = v_true_density[i];
    Scal volume = 0.;
    for (auto idxcell : mesh.Cells()) {
      Scal c = v_fc_volume_fraction[i][idxcell];
      volume += c * mesh.GetVolume(idxcell);
    }
    Scal mass = volume * density;
    P_double.set("stat_volume_" + IntToStr(i), volume);
    P_double.set("stat_mass_" + IntToStr(i), mass);

    Scal volume_flux_in = 0.;
    Scal volume_flux_out = 0.;
    for (auto idxface : mesh.Faces()) {
      if (!mesh.IsInner(idxface)) {
        auto cellid = mesh.GetValidNeighbourCellId(idxface);
        auto idxcell = mesh.GetNeighbourCell(idxface, cellid);
        Scal factor = (cellid == 1 ? 1. : -1.);
        Scal flux = fluid_solver->GetVolumeFlux()[idxface];
        Scal flux_to_cell = flux * factor;
        Scal volume_fraction = v_fc_volume_fraction[i][idxcell];
        if (flux_to_cell >=0 ) {
          volume_flux_in += flux_to_cell * volume_fraction;
        } else {
          volume_flux_out -= flux_to_cell * volume_fraction;
        }
      }
    }

    Scal mass_flux_in = volume_flux_in * density;
    Scal mass_flux_out = volume_flux_out * density;

    auto increase = [this](std::string name, Scal flux) {
      if (!P_double.exist(name)) {
        P_double.set(name, 0.);
      }
      P_double[name] += flux * fluid_solver->GetTimeStep();
      return;
    };

    std::string phase = IntToStr(i);

    increase("stat_volume_in_" + phase, volume_flux_in);
    increase("stat_volume_out_" + phase, volume_flux_out);
    increase("stat_mass_in_" + phase, mass_flux_in);
    increase("stat_mass_out_" + phase, mass_flux_out);
  }
}

template <class Mesh>
void hydro_core<Mesh>::step() {
  ex->timer_.Push("step");

  ex->timer_.Push("step.fluid");
  fluid_solver->StartStep();

  bool iter_history =
      P_bool["iter_history_enable"] && P_int["iter_history_n"] == P_int["n"]
      && !flags.no_output;
  output::Content content_iter;
  std::shared_ptr<output::Session> session_iter;
  if (!P_string.exist("filename_iter_history")) {
    P_string.set("filename_iter_history",
                 P_string[_exp_name] + ".iter_history.plt");
  }
  if (iter_history) {
    content_iter = {
        std::make_shared<output::EntryScalarFunction<Scal>>(
            "iter", [this](){
                return static_cast<Scal>(fluid_solver->GetIterationCount());
            }),
        std::make_shared<output::EntryScalarFunction<Scal>>(
            "iter_diff_velocity", [this](){
                return fluid_solver->GetConvergenceIndicator();
            })
    };
    session_iter = std::make_shared<output::SessionTecplotAsciiScalar<Scal>>(
        content_iter, P_string[_plt_title],
        P_string[_plt_title], P_string["filename_iter_history"]);
  }


  if (P_bool["fluid_enable"]) {
    size_t iter_history_sfixed = P_int["iter_history_sfixed"];
    while (!fluid_solver->IsConverged() ||
        (iter_history && P_int["iter_history_sfixed"] > 0 &&
        fluid_solver->GetIterationCount() < iter_history_sfixed)) {
      fluid_solver->MakeIteration();

      ++P_int["s_sum"];
      P_int["s"] = static_cast<int>(fluid_solver->GetIterationCount());

      logger() << ".....s=" << fluid_solver->GetIterationCount()
          << ", Rs=" << fluid_solver->GetConvergenceIndicator();

      if (iter_history) {
        session_iter->Write();
      }
    }
  }
  fluid_solver->FinishStep();
  ex->timer_.Pop();

  if (P_bool["advection_enable"]) {
    ex->timer_.Push("step.advection");
    while (advection_solver->GetTime() <
        fluid_solver->GetTime() - 0.5 * advection_solver->GetTimeStep())  {
      advection_solver->StartStep();
      advection_solver->CalcStep();
      advection_solver->FinishStep();
    };
    ex->timer_.Pop();
  }

  if (P_bool["heat_enable"]) {
    ex->timer_.Push("step.heat");
    heat_solver->StartStep();
    heat_solver->CalcStep();
    heat_solver->FinishStep();
    ex->timer_.Pop();
  }

  ex->timer_.Push("step.fluid_properties");
  UpdateFluidProperties();
  CalcStat();
  ex->timer_.Pop();

  ex->timer_.Pop();

  int rank = parallel->GetRank();

  // Parallel test:
  FieldCell<Scal> l_fc_radiation(parallel->GetLocalMesh(), rank);

  if (rank == 1) {
    parallel->SendOverlap(l_fc_radiation, 2);
  } else if (rank == 2) {
    parallel->RecvOverlap(l_fc_radiation, 1);
  }

  if (rank == 0) {
    parallel->RecvExtLocal(fc_radiation, 2);
  } else if (rank == 2) {
    parallel->SendLocal(l_fc_radiation, 0);
  }
}

template <class Mesh>
void hydro_core<Mesh>::write_results(bool force)
{
  if (flags.no_output) {
    return;
  }

  const double time = fluid_solver->GetTime();
  const double total_time = P_double["T"];
  const size_t max_frame_index = P_int["max_frame_index"];
  const double frame_duration = total_time / max_frame_index;

  if (force || (!ecast(P_bool("no_mesh_output")) &&
      time >= last_frame_time_ + frame_duration)) {
    last_frame_time_ = time;
    session->Write(time, P_string[_plt_title] + ":" + IntToStr(P_int["n"]));
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

}
