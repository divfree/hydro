/*
 * simple.hpp
 *
 *  Created on: Feb 14, 2016
 *      Author: Petr Karnakov (petr.karnakov@gmail.com)
 */

#pragma once

#include "mesh.hpp"
#include "linear.hpp"
#include <exception>
#include "solver.hpp"
#include "../control/metrics.hpp"
#include <cmath>
#include <sstream>
#include "conv_diff.hpp"
#include <memory>

namespace solver {

// TODO: Pass parameters as constructor arguments

template <class Mesh>
class ConvectionDiffusion : public UnsteadyIterativeSolver {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  static constexpr size_t dim = Mesh::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

 protected:
  geom::FieldCell<Scal>* p_fc_density_;
  geom::FieldFace<Scal>* p_ff_kinematic_viscosity_; // adhoc: dynamic viscosity
  geom::FieldCell<Vect>* p_fc_force_;
  geom::FieldFace<Scal>* p_ff_vol_flux_;

 public:
  ConvectionDiffusion(double time, double time_step,
                      geom::FieldCell<Scal>* p_fc_density,
                      geom::FieldFace<Scal>* p_ff_kinematic_viscosity,
                      geom::FieldCell<Vect>* p_fc_force,
                      geom::FieldFace<Scal>* p_ff_vol_flux,
                      double convergence_tolerance,
                      size_t num_iterations_limit)
      : UnsteadyIterativeSolver(time, time_step, convergence_tolerance,
                                num_iterations_limit)
      , p_fc_density_(p_fc_density)
      , p_ff_kinematic_viscosity_(p_ff_kinematic_viscosity)
      , p_fc_force_(p_fc_force)
      , p_ff_vol_flux_(p_ff_vol_flux)
  {

  }
  virtual const geom::FieldCell<Vect>& GetVelocity() = 0;
  virtual const geom::FieldCell<Vect>& GetVelocity(Layers layer) = 0;
  virtual void CorrectVelocity(Layers layer,
                               const geom::FieldCell<Vect>& fc_corr) = 0;
  virtual const geom::FieldCell<Expr>& GetVelocityEquations(size_t comp) = 0;
};

template <class Mesh>
class ConvectionDiffusionImplicit : public ConvectionDiffusion<Mesh> {
  const Mesh& mesh;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using Solver = ConvectionDiffusionScalarImplicit<Mesh>;

  static constexpr size_t dim = Mesh::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;
  template <class T>
  using VectGeneric = std::array<T, dim>;
  LayersData<geom::FieldCell<Vect>> fc_velocity_;

  geom::MapFace<std::shared_ptr<ConditionFace>> mf_velocity_cond_;
  geom::MapCell<std::shared_ptr<ConditionCell>> mc_velocity_cond_;
  geom::FieldCell<Vect>* p_fc_force_;

  VectGeneric<geom::MapFace<std::shared_ptr<ConditionFace>>>
  v_mf_velocity_cond_;
  // TODO: Extract scalar CellCondition
  VectGeneric<std::shared_ptr<Solver>>
  v_solver_;
  VectGeneric<geom::FieldCell<Scal>>
  v_fc_force_;

 public:
  void CopyToVector(Layers layer) {
    fc_velocity_.Get(layer).Reinit(mesh);
    for (size_t n = 0; n < dim; ++n) {
      SetComponent(fc_velocity_.Get(layer), n, v_solver_[n]->GetField(layer));
    }
  }
  ConvectionDiffusionImplicit(
      const Mesh& mesh,
      const geom::FieldCell<Vect>& fc_velocity_initial,
      const geom::MapFace<std::shared_ptr<ConditionFace>>&
      mf_velocity_cond,
      const geom::MapCell<std::shared_ptr<ConditionCell>>&
      mc_velocity_cond,
      Scal relaxation_factor,
      geom::FieldCell<Scal>* p_fc_density,
      geom::FieldFace<Scal>* p_ff_kinematic_viscosity,
      geom::FieldCell<Vect>* p_fc_force,
    geom::FieldFace<Scal>* p_ff_vol_flux,
    double time, double time_step,
    const LinearSolverFactory& linear_factory,
    double convergence_tolerance,
    size_t num_iterations_limit,
    bool time_second_order = true,
    Scal guess_extrapolation = 0.)
    : ConvectionDiffusion<Mesh>(
        time, time_step, p_fc_density, p_ff_kinematic_viscosity, p_fc_force,
        p_ff_vol_flux,
          convergence_tolerance, num_iterations_limit)
      , mesh(mesh)
      , mf_velocity_cond_(mf_velocity_cond)
      , mc_velocity_cond_(mc_velocity_cond)
      , p_fc_force_(p_fc_force)
  {
    for (size_t n = 0; n < dim; ++n) {
      // Boundary conditions for each velocity component
      // (copied from given vector conditions)
      for (auto it = mf_velocity_cond_.cbegin();
          it != mf_velocity_cond_.cend(); ++it) {
        IdxFace idxface = it->GetIdx();
        if (auto cond = dynamic_cast<ConditionFaceValue<Vect>*>(
            mf_velocity_cond_[idxface].get())) {
          v_mf_velocity_cond_[n][idxface] =
              std::make_shared<ConditionFaceValueExtractComponent<Vect>>(
                  cond, n);
        } else {
          throw std::runtime_error("Unknown boudnary condition type");
        }
      }

      // Initialize solver
      v_solver_[n] = std::make_shared<Solver>(
          mesh, GetComponent(fc_velocity_initial, n),
          v_mf_velocity_cond_[n],
          geom::MapCell<std::shared_ptr<ConditionCell>>() /*empty*/,
          relaxation_factor, p_fc_density, p_ff_kinematic_viscosity,
          &(v_fc_force_[n]), p_ff_vol_flux, time, time_step,
          linear_factory, convergence_tolerance,
          num_iterations_limit, time_second_order, guess_extrapolation);
    }
    CopyToVector(Layers::time_curr);
    CopyToVector(Layers::time_prev);
  }
  void StartStep() override {
    for (size_t n = 0; n < dim; ++n) {
      v_solver_[n]->SetTimeStep(this->GetTimeStep());
      v_solver_[n]->StartStep();
    }
    CopyToVector(Layers::iter_curr);
    this->ClearIterationCount();
  }
  void MakeIteration() override {
    for (size_t n = 0; n < dim; ++n) {
      v_fc_force_[n] = GetComponent(*p_fc_force_, n);
      v_solver_[n]->MakeIteration();
    }
    CopyToVector(Layers::iter_prev);
    CopyToVector(Layers::iter_curr);
    this->IncIterationCount();
  }
  void FinishStep() override {
    for (size_t n = 0; n < dim; ++n) {
      v_solver_[n]->FinishStep();
    }
    CopyToVector(Layers::time_prev);
    CopyToVector(Layers::time_curr);
    this->IncTime();
  }
  double GetConvergenceIndicator() const override {
    if (this->GetIterationCount() == 0) {
      return 1.;
    }
    return CalcDiff(fc_velocity_.iter_curr, fc_velocity_.iter_prev, mesh);
  }
  const geom::FieldCell<Vect>& GetVelocity() override {
    return fc_velocity_.time_curr;
  }
  const geom::FieldCell<Vect>& GetVelocity(Layers layer) override {
    return fc_velocity_.Get(layer);
  }
  void CorrectVelocity(Layers layer,
                       const geom::FieldCell<Vect>& fc_corr) override {
    for (size_t n = 0; n < dim; ++n) {
      v_solver_[n]->CorrectField(layer, GetComponent(fc_corr, n));
    }
    CopyToVector(layer);
  }
  const geom::FieldCell<Expr>& GetVelocityEquations(size_t comp) override {
    return v_solver_[comp]->GetEquations();
  }
  geom::MapFace<std::shared_ptr<ConditionFace>>&
  GetVelocityCond(size_t comp) {
    return v_mf_velocity_cond_[comp];
  }
};

template <class Mesh>
class FluidSolver : public UnsteadyIterativeSolver {
  static constexpr size_t dim = Mesh::dim;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

 protected:
  geom::FieldCell<Scal>* p_fc_density_;
  geom::FieldCell<Scal>* p_fc_viscosity_;
  geom::FieldCell<Vect>* p_fc_force_;
  geom::FieldCell<Vect>* p_fc_stforce_;
  geom::FieldFace<Vect>* p_ff_stforce_;
  geom::FieldCell<Scal>* p_fc_volume_source_;
  geom::FieldCell<Scal>* p_fc_mass_source_;
  Vect meshvel_;

 public:
  FluidSolver(double time, double time_step,
                      geom::FieldCell<Scal>* p_fc_density,
                      geom::FieldCell<Scal>* p_fc_viscosity,
                      geom::FieldCell<Vect>* p_fc_force,
                      geom::FieldCell<Vect>* p_fc_stforce,
                      geom::FieldFace<Vect>* p_ff_stforce,
                      geom::FieldCell<Scal>* p_fc_volume_source,
                      geom::FieldCell<Scal>* p_fc_mass_source,
                      double convergence_tolerance,
                      size_t num_iterations_limit,
                      Vect meshvel)
      : UnsteadyIterativeSolver(time, time_step,
                                convergence_tolerance, num_iterations_limit)
      , p_fc_density_(p_fc_density)
      , p_fc_viscosity_(p_fc_viscosity)
      , p_fc_force_(p_fc_force)
      , p_fc_stforce_(p_fc_stforce)
      , p_ff_stforce_(p_ff_stforce)
      , p_fc_volume_source_(p_fc_volume_source)
      , p_fc_mass_source_(p_fc_mass_source)
      , meshvel_(meshvel)
  {

  }
  virtual const geom::FieldCell<Vect>& GetVelocity() = 0;
  virtual const geom::FieldCell<Vect>& GetVelocity(Layers layer) = 0;
  virtual const geom::FieldCell<Scal>& GetPressure() = 0;
  virtual const geom::FieldCell<Scal>& GetPressure(Layers layer) = 0;
  virtual const geom::FieldFace<Scal>& GetVolumeFlux() = 0;
  virtual const geom::FieldFace<Scal>& GetVolumeFlux(Layers layer) = 0;
  virtual Vect GetMeshVel() { return meshvel_; }
  virtual void SetMeshVel(Vect meshvel) { meshvel_ = meshvel; }
  virtual double GetAutoTimeStep() { return GetTimeStep(); }
};

class ConditionFaceFluid : public ConditionFace {};

class ConditionCellFluid : public ConditionCell {};

namespace fluid_condition {

template <class Mesh>
class NoSlipWall : public ConditionFaceFluid {
  using Vect = typename Mesh::Vect;
 public:
  virtual Vect GetVelocity() const = 0;
};

template <class Mesh>
class NoSlipWallFixed : public NoSlipWall<Mesh> {
  using Vect = typename Mesh::Vect;
  Vect velocity_;
 public:
  NoSlipWallFixed(Vect velocity)
      : velocity_(velocity)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
};

template <class Mesh>
class Inlet : public ConditionFaceFluid {
  using Vect = typename Mesh::Vect;
 public:
  virtual Vect GetVelocity() const = 0;
  virtual size_t GetNeighbourCellId() const = 0;
};

template <class Mesh>
class InletFixed : public Inlet<Mesh> {
  using Vect = typename Mesh::Vect;
  Vect velocity_;
  size_t neighbour_cell_id_;

 public:
  InletFixed(Vect velocity, size_t neighbour_cell_id)
      : velocity_(velocity)
      , neighbour_cell_id_(neighbour_cell_id)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
  size_t GetNeighbourCellId() const override {
    return neighbour_cell_id_;
  }
};

template <class Mesh>
class Outlet : public ConditionFaceFluid {
  using Vect = typename Mesh::Vect;
 public:
  virtual Vect GetVelocity() const = 0;
  virtual void SetVelocity(Vect velocity) = 0;
  virtual size_t GetNeighbourCellId() const = 0;
};

template <class Mesh>
class OutletAuto : public Outlet<Mesh> {
  using Vect = typename Mesh::Vect;
  Vect velocity_;
  size_t neighbour_cell_id_;

 public:
  OutletAuto(size_t neighbour_cell_id)
      : velocity_(Vect::kZero)
      , neighbour_cell_id_(neighbour_cell_id)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
  void SetVelocity(Vect velocity) override {
    velocity_ = velocity;
  }
  size_t GetNeighbourCellId() const override {
    return neighbour_cell_id_;
  }
};

template <class Mesh>
class GivenVelocityAndPressure : public ConditionCellFluid {
  using Vect = typename Mesh::Vect;
  using Scal = typename Mesh::Scal;
 public:
  virtual Vect GetVelocity() const = 0;
  virtual Scal GetPressure() const = 0;
};

template <class Mesh>
class GivenVelocityAndPressureFixed : public GivenVelocityAndPressure<Mesh> {
  using Vect = typename Mesh::Vect;
  using Scal = typename Mesh::Scal;
  Vect velocity_;
  Scal pressure_;

 public:
  GivenVelocityAndPressureFixed(Vect velocity, Scal pressure)
      : velocity_(velocity)
      , pressure_(pressure)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
  Scal GetPressure() const override {
    return pressure_;
  }
};

template <class Mesh>
class GivenPressure : public ConditionCellFluid {
  using Scal = typename Mesh::Scal;
 public:
  virtual Scal GetPressure() const = 0;
};

template <class Mesh>
class GivenPressureFixed : public GivenPressure<Mesh> {
  using Scal = typename Mesh::Scal;
  Scal pressure_;

 public:
  GivenPressureFixed(Scal pressure)
      : pressure_(pressure)
  {}
  Scal GetPressure() const override {
    return pressure_;
  }
};

} // namespace fluid_condition


template <class Mesh>
std::shared_ptr<ConditionFaceFluid> Parse(std::string argstr,
                                          geom::IdxFace idxface,
                                          const Mesh& mesh) {
  using namespace fluid_condition;
  std::stringstream arg(argstr);

  std::string name;
  arg >> name;

  if (name == "wall") {
    typename Mesh::Vect vel;
    arg >> vel;
    return std::make_shared<NoSlipWallFixed<Mesh>>(vel);
  } else if (name == "inlet") {
    typename Mesh::Vect vel;
    arg >> vel;
    return std::make_shared<InletFixed<Mesh>>(
        vel, mesh.GetValidNeighbourCellId(idxface));
  } else if (name == "outlet") {
    return std::make_shared<OutletAuto<Mesh>>(
        mesh.GetValidNeighbourCellId(idxface));
  } else {
    throw std::runtime_error("Parse: Unknown boundary condition type");
  }
}

// TODO: Second order time step
// TODO: Extrapolation for first iteration
template <class Mesh>
class FluidSimple : public FluidSolver<Mesh> {
  const Mesh& mesh;
  static constexpr size_t dim = Mesh::dim;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using IdxNode = geom::IdxNode;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

  geom::FieldCell<Vect> fc_force_;
  Scal velocity_relaxation_factor_;
  Scal pressure_relaxation_factor_;
  Scal rhie_chow_factor_;
  LayersData<geom::FieldFace<Scal>> ff_vol_flux_;
  std::shared_ptr<ConvectionDiffusionImplicit<Mesh>> conv_diff_solver_;

  LayersData<geom::FieldCell<Scal>> fc_pressure_;
  geom::FieldCell<Scal> fc_kinematic_viscosity_;
  geom::FieldFace<Scal> ff_kinematic_viscosity_;

  geom::MapFace<std::shared_ptr<ConditionFaceFluid>> mf_cond_;

  geom::MapFace<std::shared_ptr<ConditionFace>> mf_velocity_cond_;
  // TODO: Const specifier for ConditionFace*

  geom::MapFace<std::shared_ptr<ConditionFace>> mf_pressure_cond_;

  geom::MapFace<std::shared_ptr<ConditionFace>> mf_pressure_grad_cond_;

  geom::MapFace<std::shared_ptr<ConditionFace>> mf_force_cond_;

  geom::MapFace<std::shared_ptr<ConditionFace>> mf_pressure_corr_cond_;

  geom::MapFace<std::shared_ptr<ConditionFace>> mf_viscosity_cond_;

  geom::MapCell<std::shared_ptr<ConditionCellFluid>> mc_cond_;
  geom::MapCell<std::shared_ptr<ConditionCell>> mc_pressure_cond_;
  geom::MapCell<std::shared_ptr<ConditionCell>> mc_velocity_cond_;

  std::shared_ptr<LinearSolver<Scal, IdxCell, Expr>> linear_;

  geom::FieldFace<bool> is_boundary_;

  // common buffers
  geom::FieldFace<Scal> ff_pressure_;
  geom::FieldCell<Vect> fc_pressure_grad_;
  geom::FieldFace<Vect> ff_pressure_grad_;
  geom::FieldCell<Vect> fc_velocity_asterisk_;
  geom::FieldFace<Vect> ff_velocity_asterisk_;
  geom::FieldFace<Scal> ff_volume_flux_asterisk_;
  geom::FieldFace<Scal> ff_volume_flux_interpolated_;
  geom::FieldCell<Scal> fc_diag_coeff_;
  geom::FieldFace<Scal> ff_diag_coeff_;
  geom::FieldFace<Expr> ff_volume_flux_corr_;
  geom::FieldCell<Expr> fc_pressure_corr_system_;
  geom::FieldCell<Scal> fc_pressure_corr_;
  geom::FieldCell<Vect> fc_pressure_corr_grad_;
  geom::FieldFace<Vect> ff_ext_force_;
  geom::FieldCell<Vect> fc_ext_force_restored_;
  geom::FieldFace<Vect> ff_ext_force_restored_;
  geom::FieldFace<Vect> ff_stforce_restored_;
  // / needed for MMIM, now disabled
  //geom::FieldFace<Vect> ff_velocity_iter_prev_;

  MultiTimer<std::string>* timer_;
  bool time_second_order_;
  bool simpler_;
  bool force_geometric_average_;
  Scal guess_extrapolation_;

  void UpdateDerivedConditions() {
    using namespace fluid_condition;

    for (auto it = mf_cond_.cbegin();
        it != mf_cond_.cend(); ++it) {
      IdxFace idxface = it->GetIdx();
      ConditionFaceFluid* cond_generic = it->GetValue().get();

      if (auto cond = dynamic_cast<NoSlipWall<Mesh>*>(cond_generic)) {
        *dynamic_cast<ConditionFaceValueFixed<Vect>*>(
            mf_velocity_cond_[idxface].get()) =
            ConditionFaceValueFixed<Vect>(cond->GetVelocity());
      } else if (auto cond = dynamic_cast<Inlet<Mesh>*>(cond_generic)) {
        *dynamic_cast<ConditionFaceValueFixed<Vect>*>(
            mf_velocity_cond_[idxface].get()) =
            ConditionFaceValueFixed<Vect>(cond->GetVelocity());
      } else if (auto cond = dynamic_cast<Outlet<Mesh>*>(cond_generic)) {
        *dynamic_cast<ConditionFaceValueFixed<Vect>*>(
            mf_velocity_cond_[idxface].get()) =
            ConditionFaceValueFixed<Vect>(cond->GetVelocity());
      } else {
        throw std::runtime_error("Unknown fluid condition");
      }
    }


    for (auto it = mc_cond_.cbegin();
        it != mc_cond_.cend(); ++it) {
      IdxCell idxcell = it->GetIdx();
      ConditionCellFluid* cond_generic = it->GetValue().get();

      if (auto cond = dynamic_cast<GivenPressure<Mesh>*>(cond_generic)) {
        *dynamic_cast<ConditionCellValueFixed<Scal>*>(
            mc_pressure_cond_[idxcell].get()) =
                ConditionCellValueFixed<Scal>(cond->GetPressure());
      } else if (auto cond =
          dynamic_cast<GivenVelocityAndPressure<Mesh>*>(cond_generic)) {
        *dynamic_cast<ConditionCellValueFixed<Vect>*>(
            mc_velocity_cond_[idxcell].get()) =
                ConditionCellValueFixed<Vect>(cond->GetVelocity());
        *dynamic_cast<ConditionCellValueFixed<Scal>*>(
            mc_pressure_cond_[idxcell].get()) =
                ConditionCellValueFixed<Scal>(cond->GetPressure());
      } else {
        throw std::runtime_error("Unknown fluid cell condition");
      }
    }
  }
  // TODO: Think about seperate canals in one domain
  void UpdateOutletBaseConditions() {
    using namespace fluid_condition;
    // Extrapolate velocity on outlet faces from cell centers
    // and calculate total inlet and outlet volumetric fluxes
    Scal inlet_volume_flux = 0.;   // Both should be positive
    Scal outlet_volume_flux = 0.;
    Scal outlet_area = 0.;
    for (auto it = mf_cond_.cbegin();
        it != mf_cond_.cend(); ++it) {
      IdxFace idxface = it->GetIdx();
      ConditionFaceFluid* cond_generic = it->GetValue().get();

      if (auto cond = dynamic_cast<Outlet<Mesh>*>(cond_generic)) {
        size_t id = cond->GetNeighbourCellId();
        IdxCell idxcell = mesh.GetNeighbourCell(idxface, id);
        Scal factor = (id == 0 ? 1. : -1.);

        cond->SetVelocity(
            this->GetVelocity(Layers::iter_curr)[idxcell]);

        outlet_volume_flux +=
            cond->GetVelocity().dot(mesh.GetSurface(idxface)) * factor;

        outlet_area += mesh.GetArea(idxface);
      } else if (auto cond = dynamic_cast<Inlet<Mesh>*>(cond_generic)) {
        size_t id = cond->GetNeighbourCellId();
        Scal factor = (id == 0 ? -1. : 1.);

        inlet_volume_flux +=
            cond->GetVelocity().dot(mesh.GetSurface(idxface)) * factor;
      }
    }

    for (auto idxcell : mesh.Cells()) {
      inlet_volume_flux +=
          (*this->p_fc_volume_source_)[idxcell] * mesh.GetVolume(idxcell);
    }

    Scal average_outlet_velocity_correction = // Additive correction
        (inlet_volume_flux - outlet_volume_flux) / outlet_area;

    // Apply correction on outlet faces
    for (auto it = mf_cond_.cbegin();
        it != mf_cond_.cend(); ++it) {
      IdxFace idxface = it->GetIdx();
      ConditionFaceFluid* cond_generic = it->GetValue().get();

      if (auto cond = dynamic_cast<Outlet<Mesh>*>(cond_generic)) {
        size_t id = cond->GetNeighbourCellId();
        Scal factor = (id == 0 ? 1. : -1.);

        Vect normal = mesh.GetSurface(idxface) / mesh.GetArea(idxface);

        cond->SetVelocity(
            cond->GetVelocity() +
            normal * average_outlet_velocity_correction * factor);
      }
    }
  }

  void CalcExtForce() {
    timer_->Push("fluid.1.force-correction");
    // Interpolate force to faces (considered as given force)
    ff_ext_force_ = Interpolate(
        *this->p_fc_force_, mf_force_cond_, mesh,
        force_geometric_average_);
    ff_stforce_restored_ = Interpolate(
        *this->p_fc_stforce_, mf_force_cond_, mesh,
        force_geometric_average_);
    fc_ext_force_restored_.Reinit(mesh);

#pragma omp parallel for
    for (IntIdx i = 0; i < static_cast<IntIdx>(mesh.Cells().size()); ++i) {
      IdxCell idxcell(i);
    //for (auto idxcell : mesh.Cells()) {
      Vect sum = Vect::kZero;
      for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
        IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
        sum += mesh.GetSurface(idxface) *
            ff_ext_force_[idxface].dot(mesh.GetNormal(idxface)) *
            mesh.GetCenter(idxcell).dist(mesh.GetCenter(idxface));
      }
      fc_ext_force_restored_[idxcell] = sum / mesh.GetVolume(idxcell);
    }
    // Interpolated restored force to faces (needed later)
    ff_ext_force_restored_ = Interpolate(
        fc_ext_force_restored_, mf_force_cond_, mesh,
        force_geometric_average_);
    timer_->Pop();
  }
  void CalcKinematicViscosity() {
    fc_kinematic_viscosity_.Reinit(mesh);
    for (auto idxcell : mesh.Cells()) {
      fc_kinematic_viscosity_[idxcell] =
          (*this->p_fc_viscosity_)[idxcell];
    }
    ff_kinematic_viscosity_ = Interpolate(
        fc_kinematic_viscosity_, mf_viscosity_cond_, mesh,
        force_geometric_average_);
  }

 public:
  FluidSimple(const Mesh& mesh,
              const geom::FieldCell<Vect>& fc_velocity_initial,
              const geom::MapFace<std::shared_ptr<ConditionFaceFluid>>&
              mf_cond,
              const geom::MapCell<std::shared_ptr<ConditionCellFluid>>&
              mc_cond,
              Scal velocity_relaxation_factor,
              Scal pressure_relaxation_factor,
              Scal rhie_chow_factor,
              geom::FieldCell<Scal>* p_fc_density,
              geom::FieldCell<Scal>* p_fc_viscosity,
              geom::FieldCell<Vect>* p_fc_force,
              geom::FieldCell<Vect>* p_fc_stforce,
              geom::FieldFace<Vect>* p_ff_stforce,
              geom::FieldCell<Scal>* p_fc_volume_source,
              geom::FieldCell<Scal>* p_fc_mass_source,
              double time, double time_step,
              const LinearSolverFactory& linear_factory_velocity,
              const LinearSolverFactory& linear_factory_pressure,
              double convergence_tolerance,
              size_t num_iterations_limit,
              MultiTimer<std::string>* timer,
              bool time_second_order,
              bool simpler,
              bool force_geometric_average,
              Scal guess_extrapolation = 0.,
              Vect meshvel=0)
      : FluidSolver<Mesh>(time, time_step, p_fc_density, p_fc_viscosity,
                    p_fc_force, p_fc_stforce, p_ff_stforce,
                    p_fc_volume_source, p_fc_mass_source,
                    convergence_tolerance, num_iterations_limit, meshvel)
      , mesh(mesh)
      , fc_force_(mesh)
      , velocity_relaxation_factor_(velocity_relaxation_factor)
      , pressure_relaxation_factor_(pressure_relaxation_factor)
      , rhie_chow_factor_(rhie_chow_factor)
      , mf_cond_(mf_cond)
      , mc_cond_(mc_cond)
      , ff_volume_flux_corr_(mesh)
      , fc_pressure_corr_system_(mesh)
      , timer_(timer)
      , time_second_order_(time_second_order)
      , simpler_(simpler)
      , force_geometric_average_(force_geometric_average)
      , guess_extrapolation_(guess_extrapolation)
  {
    linear_ = linear_factory_pressure.Create<Scal, IdxCell, Expr>();

    using namespace fluid_condition;

    is_boundary_.Reinit(mesh, false);
    for (auto it = mf_cond_.cbegin();
        it != mf_cond_.cend(); ++it) {
      IdxFace idxface = it->GetIdx();
      is_boundary_[idxface] = true;
      ConditionFaceFluid* cond_generic = it->GetValue().get();

      if (auto cond = dynamic_cast<NoSlipWall<Mesh>*>(cond_generic)) {
        mf_velocity_cond_[idxface] =
            std::make_shared<
            ConditionFaceValueFixed<Vect>>(cond->GetVelocity());
        mf_pressure_cond_[idxface] =
            std::make_shared<ConditionFaceExtrapolation>();
      } else if (auto cond = dynamic_cast<Inlet<Mesh>*>(cond_generic)) {
        mf_velocity_cond_[idxface] =
            std::make_shared<
            ConditionFaceValueFixed<Vect>>(cond->GetVelocity());
        mf_pressure_cond_[idxface] =
            std::make_shared<ConditionFaceExtrapolation>();
      } else if (auto cond = dynamic_cast<Outlet<Mesh>*>(cond_generic)) {
        mf_velocity_cond_[idxface] =
            std::make_shared<
            ConditionFaceValueFixed<Vect>>(cond->GetVelocity());
        mf_pressure_cond_[idxface] =
            std::make_shared<ConditionFaceExtrapolation>();
      } else {
        throw std::runtime_error("Unknown fluid condition");
      }

      mf_pressure_grad_cond_[idxface] =
          std::make_shared<ConditionFaceDerivativeFixed<Vect>>(Vect::kZero);
      mf_force_cond_[idxface] =
          std::make_shared<ConditionFaceDerivativeFixed<Vect>>(Vect::kZero);
      mf_pressure_corr_cond_[idxface] =
          std::make_shared<ConditionFaceExtrapolation>();
      mf_viscosity_cond_[idxface] =
          std::make_shared<ConditionFaceDerivativeFixed<Scal>>(0.);
    }

    for (auto it = mc_cond_.cbegin();
        it != mc_cond_.cend(); ++it) {
      IdxCell idxcell = it->GetIdx();
      ConditionCellFluid* cond_generic = it->GetValue().get();

      if (auto cond = dynamic_cast<GivenPressure<Mesh>*>(cond_generic)) {
        mc_pressure_cond_[idxcell] =
            std::make_shared<
            ConditionCellValueFixed<Scal>>(cond->GetPressure());
      } else if (auto cond =
          dynamic_cast<GivenVelocityAndPressure<Mesh>*>(cond_generic)) {
        mc_pressure_cond_[idxcell] =
            std::make_shared<
            ConditionCellValueFixed<Scal>>(cond->GetPressure());
        mc_velocity_cond_[idxcell] =
            std::make_shared<
            ConditionCellValueFixed<Vect>>(cond->GetVelocity());
      } else {
        throw std::runtime_error("Unknown fluid cell condition");
      }
    }

    conv_diff_solver_ = std::make_shared<
        ConvectionDiffusionImplicit<Mesh>>(
            mesh, fc_velocity_initial,
            mf_velocity_cond_, mc_velocity_cond_,
            velocity_relaxation_factor_,
            p_fc_density, &ff_kinematic_viscosity_, &fc_force_,
            &ff_vol_flux_.iter_prev,
            time, time_step,
            linear_factory_velocity,
            convergence_tolerance, num_iterations_limit,
            time_second_order_, guess_extrapolation_);

    fc_pressure_.time_curr.Reinit(mesh, 0.);
    fc_pressure_.time_prev = fc_pressure_.time_curr;

    // Calc initial volume fluxes
    fc_velocity_asterisk_ = conv_diff_solver_->GetVelocity();
    ff_velocity_asterisk_ = Interpolate(
        fc_velocity_asterisk_, mf_velocity_cond_, mesh);
    ff_vol_flux_.time_curr.Reinit(mesh, 0.);
    for (auto idxface : mesh.Faces()) {
      ff_vol_flux_.time_curr[idxface] =
          ff_velocity_asterisk_[idxface].dot(mesh.GetSurface(idxface));
    }
    // Apply meshvel
    for (auto idxface : mesh.Faces()) {
      ff_vol_flux_.time_curr[idxface] -= 
          this->meshvel_.dot(mesh.GetSurface(idxface));
    }

    ff_vol_flux_.time_prev = ff_vol_flux_.time_curr;
    ff_volume_flux_interpolated_.Reinit(mesh, 0.);

    ff_pressure_ = Interpolate(fc_pressure_.time_curr, mf_pressure_cond_, mesh);
    fc_pressure_grad_ = Gradient(ff_pressure_, mesh);
    ff_pressure_grad_ = Interpolate(
        fc_pressure_grad_, mf_pressure_grad_cond_, mesh);
  }
  void StartStep() override {
    this->ClearIterationCount();
    if (IsNan(fc_pressure_.time_curr)) {
      throw std::string("NaN initial pressure");
    }
    conv_diff_solver_->SetTimeStep(this->GetTimeStep());
    conv_diff_solver_->StartStep();
    fc_pressure_.iter_curr = fc_pressure_.time_curr;
    ff_vol_flux_.iter_curr = ff_vol_flux_.time_curr;
    for (auto idxcell : mesh.Cells()) {
      fc_pressure_.iter_curr[idxcell] +=
          (fc_pressure_.time_curr[idxcell] - fc_pressure_.time_prev[idxcell]) *
          guess_extrapolation_;
    }
    for (auto idxface : mesh.Faces()) {
      ff_vol_flux_.iter_curr[idxface] +=
          (ff_vol_flux_.time_curr[idxface] - ff_vol_flux_.time_prev[idxface]) *
          guess_extrapolation_;
    }
  }
  // TODO: rewrite norm() using dist() where needed
  void MakeIteration() override {
    auto& fc_pressure_prev = fc_pressure_.iter_prev;
    auto& fc_pressure_curr = fc_pressure_.iter_curr;
    fc_pressure_prev = fc_pressure_curr;
    ff_vol_flux_.iter_prev = ff_vol_flux_.iter_curr;

    UpdateOutletBaseConditions();
    UpdateDerivedConditions();

    CalcExtForce();

    CalcKinematicViscosity();

    timer_->Push("fluid.0.pressure-gradient");
    ff_pressure_ = Interpolate(fc_pressure_prev, mf_pressure_cond_, mesh);
    fc_pressure_grad_ = Gradient(ff_pressure_, mesh);
    ff_pressure_grad_ = Interpolate(
        fc_pressure_grad_, mf_pressure_grad_cond_, mesh);
    timer_->Pop();

    // initialize force with zero
    fc_force_.Reinit(mesh, Vect(0));
    // append viscous term
    timer_->Push("fluid.1a.explicit-viscosity");
    for (size_t n = 0; n < dim; ++n) {
      geom::FieldCell<Scal> fc = GetComponent(
          conv_diff_solver_->GetVelocity(Layers::iter_curr), n);
      auto ff = Interpolate(fc, conv_diff_solver_->GetVelocityCond(n), mesh);
      auto gc = Gradient(ff, mesh);
      auto gf = Interpolate(gc, mf_force_cond_, mesh); // adhoc: zero-der cond
      for (auto idxcell : mesh.Cells()) {
        Vect sum = Vect::kZero;
        for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
          IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
          sum += gf[idxface] * (ff_kinematic_viscosity_[idxface] * 
              mesh.GetOutwardSurface(idxcell, i)[n]);
        }
        fc_force_[idxcell] += sum / mesh.GetVolume(idxcell);
      }
    }
    timer_->Pop();

    // append to force
#pragma omp parallel for
    for (IntIdx i = 0; i < static_cast<IntIdx>(mesh.Cells().size()); ++i) {
      IdxCell idxcell(i);
    //for (auto idxcell : mesh.Cells()) {
      fc_force_[idxcell] +=
          fc_pressure_grad_[idxcell] * (-1.) +
          fc_ext_force_restored_[idxcell] +
          (*this->p_fc_stforce_)[idxcell] +
          // Volume source momentum compensation:
          conv_diff_solver_->GetVelocity(Layers::iter_curr)[idxcell] *
          ((*this->p_fc_density_)[idxcell] *
          (*this->p_fc_volume_source_)[idxcell] -
          (*this->p_fc_mass_source_)[idxcell]);
    }

    timer_->Push("fluid.2.convection-diffusion");
    conv_diff_solver_->MakeIteration();
    timer_->Pop();

    fc_diag_coeff_.Reinit(mesh);
    for (auto idxcell : mesh.Cells()) {
      Scal sum = 0.;
      for (size_t n = 0; n < dim; ++n) {
        sum += conv_diff_solver_->GetVelocityEquations(n)[idxcell].CoeffSum();
      }
      fc_diag_coeff_[idxcell] = sum / dim;
    }

    // Define ff_diag_coeff_ on inner faces only
    ff_diag_coeff_ = Interpolate(
        fc_diag_coeff_, geom::MapFace<std::shared_ptr<ConditionFace>>(), mesh);

    fc_velocity_asterisk_ = conv_diff_solver_->GetVelocity(Layers::iter_curr);

    ff_velocity_asterisk_ = Interpolate(
        fc_velocity_asterisk_, mf_velocity_cond_, mesh);

    // // needed for MMIM, now disabled
    //ff_velocity_iter_prev_ = Interpolate(
    //    conv_diff_solver_->GetVelocity(Layers::iter_prev),
    //    mf_velocity_cond_, mesh);

    // Calc volumetric flux (asterisk)
    // using momentum interpolation (Rhie-Chow)
    // including hydrostatics-correction
    // TODO: Extend hydrostatics-correction on a non-uniform mesh
    timer_->Push("fluid.3.momentum-interpolation");
    ff_volume_flux_asterisk_.Reinit(mesh);
    for (auto idxface : mesh.Faces()) {
      const auto volume_flux_interpolated =
          ff_velocity_asterisk_[idxface].dot(mesh.GetSurface(idxface));
      if (!is_boundary_[idxface] && !mesh.IsExcluded(idxface)) {
        IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
        IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
        Vect dm = mesh.GetVectToCell(idxface, 0);
        Vect dp = mesh.GetVectToCell(idxface, 1);
        const auto pressure_surface_derivative_wide =
            (ff_pressure_grad_[idxface] -
            //ff_stforce_restored_[idxface] -
            ff_ext_force_restored_[idxface]).dot(mesh.GetSurface(idxface));
        const auto pressure_surface_derivative_compact =
            (fc_pressure_prev[cp] - fc_pressure_prev[cm]) /
            (dp - dm).norm() * mesh.GetArea(idxface) -
            //(*this->p_ff_stforce_)[idxface].dot(mesh.GetNormal(idxface)) - 
            ff_ext_force_[idxface].dot(mesh.GetSurface(idxface));
        //const auto mmim = (1. - velocity_relaxation_factor_) *
        //    (ff_vol_flux_.iter_prev[idxface] -
        //     ff_velocity_iter_prev_[idxface].dot(mesh.GetSurface(idxface)));
        ff_volume_flux_asterisk_[idxface] =
            volume_flux_interpolated +
            rhie_chow_factor_ * (pressure_surface_derivative_wide -
            pressure_surface_derivative_compact) / ff_diag_coeff_[idxface] +
            0; //mmim; // TODO: Test MMIM
      } else {
        ff_volume_flux_asterisk_[idxface] =
            ff_velocity_asterisk_[idxface].dot(mesh.GetSurface(idxface));
      }
    }

    // Apply meshvel
    for (auto idxface : mesh.Faces()) {
      ff_volume_flux_asterisk_[idxface] -= 
          this->meshvel_.dot(mesh.GetSurface(idxface));
    }

    timer_->Pop();

    // TODO: Rename SurfaceVelocity to MassFlux or VolumeFlux

    timer_->Push("fluid.4.volume-flux");
    // TODO: Rename velocity_corr to smth
    // (it's actually full velocity not just correction)
    // (same for ff_volume_flux_corr_)
#pragma omp parallel for
    for (IntIdx i = 0; i < static_cast<IntIdx>(mesh.Faces().size()); ++i) {
      IdxFace idxface(i);
    //for (auto idxface : mesh.Faces()) {
      auto& expr = ff_volume_flux_corr_[idxface];
      expr.Clear();
      if (!is_boundary_[idxface] && !mesh.IsExcluded(idxface)) {
        IdxCell cm = mesh.GetNeighbourCell(idxface, 0);
        IdxCell cp = mesh.GetNeighbourCell(idxface, 1);
        Vect dm = mesh.GetVectToCell(idxface, 0);
        Vect dp = mesh.GetVectToCell(idxface, 1);
        auto coeff = - mesh.GetArea(idxface) /
            ((dp - dm).norm() * ff_diag_coeff_[idxface]);
        expr.InsertTerm(-coeff, cm);
        expr.InsertTerm(coeff, cp);
        // adhoc for periodic
        expr.SortTerms(true);
      }
      expr.SetConstant(ff_volume_flux_asterisk_[idxface]);
    }
    timer_->Pop();

    timer_->Push("fluid.5.pressure-system");
#pragma omp parallel for
    for (IntIdx i = 0; i < static_cast<IntIdx>(mesh.Cells().size()); ++i) {
      IdxCell idxcell(i);
      auto& eqn = fc_pressure_corr_system_[idxcell];
      if (!mesh.IsExcluded(idxcell)) {
        Expr flux_sum;
        for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
          IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
          flux_sum +=
              ff_volume_flux_corr_[idxface] *
              mesh.GetOutwardFactor(idxcell, i);
        }
        eqn =
            flux_sum -
            Expr((*this->p_fc_volume_source_)[idxcell] *
                 mesh.GetVolume(idxcell));
      } else {
        eqn.Clear();
        eqn.InsertTerm(1., idxcell);
        eqn.SetConstant(0.);
      }
    }

    // Account for cell conditions for pressure
    for (auto it = mc_pressure_cond_.cbegin();
        it != mc_pressure_cond_.cend(); ++it) {
      IdxCell idxcell(it->GetIdx());
      ConditionCell* cond = it->GetValue().get();
      if (auto cond_value = dynamic_cast<ConditionCellValue<Scal>*>(cond)) {
        for (auto idxlocal : mesh.Cells()) {
          auto& eqn = fc_pressure_corr_system_[idxlocal];
          if (idxlocal == idxcell) {
            eqn.Clear();
            eqn.InsertTerm(1., idxcell);
            eqn.SetConstant(-cond_value->GetValue());
          } else {
            // Substitute value to obtain symmetrix matrix
            eqn.SetKnownValue(idxcell, cond_value->GetValue());
          }
        }
      }
    }
/*
    if (dynamic_cast<Pardiso<Scal, IdxCell, Expr>*>(linear_.get())) {
      static bool first = true;
      if (first) {
        std::cout << "pressure_system restricted" << std::endl;
        first = false;
      }
      for (auto idxcell : mesh.Cells()) {
        fc_pressure_corr_system_[idxcell].RestrictTerms(
            idxcell, IdxCell(mesh.Cells().size()));
      }
    }*/
    timer_->Pop();

    timer_->Push("fluid.6.pressure-solve");
    fc_pressure_corr_ = linear_->Solve(fc_pressure_corr_system_);
    timer_->Pop();

    timer_->Push("fluid.7.correction");
    // Correct pressure
    for (auto idxcell : mesh.Cells()) {
      fc_pressure_curr[idxcell] = fc_pressure_prev[idxcell] +
          pressure_relaxation_factor_ * fc_pressure_corr_[idxcell];
    }

    fc_pressure_corr_grad_ = Gradient(
        Interpolate(fc_pressure_corr_, mf_pressure_corr_cond_, mesh),
        mesh);

    // Correct the velocity
    geom::FieldCell<Vect> fc_velocity_corr(mesh);
    for (auto idxcell : mesh.Cells()) {
      fc_velocity_corr[idxcell] =
          fc_pressure_corr_grad_[idxcell] / (-fc_diag_coeff_[idxcell]);
    }
    conv_diff_solver_->CorrectVelocity(Layers::iter_curr, fc_velocity_corr);

    // Calc divergence-free volume fluxes
    for (auto idxface : mesh.Faces()) {
      ff_vol_flux_.iter_curr[idxface] =
          ff_volume_flux_corr_[idxface].Evaluate(fc_pressure_corr_);
    }
    timer_->Pop();

    // Apply SIMPLER pressure correction
    if (simpler_) {
      timer_->Push("fluid.8.simpler");
      const auto& fc_velocity =
          conv_diff_solver_->GetVelocity(Layers::iter_curr);
          //fc_velocity_asterisk_;
      const auto& ff_volume_flux =
          ff_vol_flux_.iter_curr;
          //ff_volume_flux_asterisk_;
      // Evaluate momentum equations using new velocity
      geom::FieldFace<Vect> ff_velocity(mesh);
      ff_velocity = Interpolate(
          fc_velocity, mf_velocity_cond_, mesh);

      auto fc_velocity_delta = fc_velocity;
#pragma omp parallel for
      for (IntIdx rawcell = 0; rawcell < static_cast<IntIdx>(mesh.Cells().size()); ++rawcell) {
        IdxCell idxcell(rawcell);
        fc_velocity_delta[idxcell] -=
            conv_diff_solver_->GetVelocity(Layers::iter_prev)[idxcell];
      }

      geom::FieldCell<Vect> fc_evaluated(mesh);
#pragma omp parallel for
      for (IntIdx rawcell = 0; rawcell < static_cast<IntIdx>(mesh.Cells().size()); ++rawcell) {
        IdxCell idxcell(rawcell);
        for (size_t n = 0; n < dim; ++n) {
          fc_evaluated[idxcell][n] =
              conv_diff_solver_->GetVelocityEquations(n)[idxcell].Evaluate(
                  fc_velocity_delta)[n];
        }

        fc_evaluated[idxcell] +=
            fc_ext_force_restored_[idxcell] - fc_pressure_grad_[idxcell];
      }

      auto ff_evaluated = Interpolate(
          fc_evaluated, geom::MapFace<std::shared_ptr<ConditionFace>>(), mesh);

      geom::FieldFace<Scal> ff_rhs(mesh, 0);
#pragma omp parallel for
      for (IntIdx rawface = 0; rawface < static_cast<IntIdx>(mesh.Faces().size()); ++rawface) {
        IdxFace idxface(rawface);
        if (!is_boundary_[idxface] && !mesh.IsExcluded(idxface)) {
          Vect eval = ff_evaluated[idxface] - ff_ext_force_[idxface];
          ff_rhs[idxface] = eval.dot(mesh.GetSurface(idxface)) /
              ff_diag_coeff_[idxface] +
              (ff_volume_flux[idxface] -
                  ff_velocity[idxface].dot(mesh.GetSurface(idxface))) /
                  rhie_chow_factor_;
        }
      }

#pragma omp parallel for
      for (IntIdx rawcell = 0; rawcell < static_cast<IntIdx>(mesh.Cells().size()); ++rawcell) {
        IdxCell idxcell(rawcell);
        Scal sum = 0.;
        for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
          IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
          sum += ff_rhs[idxface] * mesh.GetOutwardFactor(idxcell, i);
        }
        fc_pressure_corr_system_[idxcell].SetConstant(-sum);
      }

      // Account for cell conditions for pressure
      for (auto it = mc_pressure_cond_.cbegin();
          it != mc_pressure_cond_.cend(); ++it) {
        IdxCell idxcell(it->GetIdx());
        ConditionCell* cond = it->GetValue().get();
        if (auto cond_value = dynamic_cast<ConditionCellValue<Scal>*>(cond)) {
          for (auto idxlocal : mesh.Cells()) {
            auto& eqn = fc_pressure_corr_system_[idxlocal];
            if (idxlocal == idxcell) {
              eqn.Clear();
              eqn.InsertTerm(1., idxcell);
              eqn.SetConstant(-cond_value->GetValue());
            } else {
              // Substitute value to obtain symmetrix matrix
              eqn.SetKnownValue(idxcell, cond_value->GetValue());
            }
          }
        }
      }

      for (auto idxcell : mesh.Cells()) {
        auto& eqn = fc_pressure_corr_system_[idxcell];
        eqn.SetConstant(eqn.Evaluate(fc_pressure_curr));
      }

      fc_pressure_corr_ = linear_->Solve(fc_pressure_corr_system_);

      for (auto idxcell : mesh.Cells()) {
        fc_pressure_curr[idxcell] += fc_pressure_corr_[idxcell];
      }

      timer_->Pop();
    }

    this->IncIterationCount();
  }
  void FinishStep() override {
    fc_pressure_.time_prev = fc_pressure_.time_curr;
    ff_vol_flux_.time_prev = ff_vol_flux_.time_curr;
    fc_pressure_.time_curr = fc_pressure_.iter_curr;
    ff_vol_flux_.time_curr = ff_vol_flux_.iter_curr;
    if (IsNan(fc_pressure_.time_curr)) {
      throw std::string("NaN pressure");
    }
    conv_diff_solver_->FinishStep();
    this->IncTime();
  }
  double GetConvergenceIndicator() const override {
    return conv_diff_solver_->GetConvergenceIndicator();
  }
  const geom::FieldCell<Vect>& GetVelocity() override {
    return conv_diff_solver_->GetVelocity();
  }
  const geom::FieldCell<Vect>& GetVelocity(Layers layer) override {
    return conv_diff_solver_->GetVelocity(layer);
  }
  const geom::FieldCell<Scal>& GetPressure() override {
    return fc_pressure_.time_curr;
  }
  const geom::FieldCell<Scal>& GetPressure(Layers layer) override {
    return fc_pressure_.Get(layer);
  }
  const geom::FieldFace<Scal>& GetVolumeFlux() override {
    return ff_vol_flux_.time_curr;
  }
  const geom::FieldFace<Scal>& GetVolumeFlux(Layers layer) override {
    return ff_vol_flux_.Get(layer);
  }
  double GetAutoTimeStep() override { 
    double dt = 1e10;
    auto& flux = ff_vol_flux_.time_curr;
    for (auto idxcell : mesh.Cells()) {
      if (!mesh.IsExcluded(idxcell)) {
        for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
          IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
          if (flux[idxface] != 0.) {
            dt = std::min(dt, 
                std::abs(mesh.GetVolume(idxcell) / flux[idxface]));
          }
        }
      }
    }
    return dt; 
  }
};



template <class Mesh>
class FluidSimpleParallel : public FluidSolver<Mesh> {
  static constexpr size_t dim = Mesh::dim;
  const Mesh& global_mesh;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using MIdx = typename Mesh::MIdx;
  using Direction = typename Mesh::Direction;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using IdxNode = geom::IdxNode;
  using BlockNodes = geom::BlockNodes<dim>;
  using BlockCells = geom::BlockCells<dim>;
  using BlockFaces = geom::BlockFaces<dim>;
  template <class T>
  using FieldCell = geom::FieldCell<T>;
  template <class T>
  using FieldFace = geom::FieldFace<T>;
  template <class T>
  using MapFace = geom::MapFace<T>;
  template <class T>
  using MapCell = geom::MapCell<T>;

  MapFace<std::shared_ptr<ConditionFaceFluid>> mf_cond_;
  MapCell<std::shared_ptr<ConditionCellFluid>> mc_cond_;
  Scal velocity_relaxation_factor_;
  Scal pressure_relaxation_factor_;
  LayersData<FieldFace<Scal>> ff_vol_flux_;
  LayersData<FieldCell<Scal>> fc_pressure_;
  LayersData<FieldCell<Vect>> fc_velocity_;
  MultiTimer<std::string>* timer_;

  struct LocalData {
    Mesh mesh;
    FieldCell<IdxCell> fc_global;
    FieldCell<bool> fc_active;
    FieldFace<IdxFace> ff_global;
    FieldFace<bool> ff_active;
    std::shared_ptr<FluidSolver<Mesh>> fluid_solver;
    MapFace<std::shared_ptr<ConditionFaceFluid>> mf_cond;
    MapCell<std::shared_ptr<ConditionCellFluid>> mc_cond;
    FieldCell<Scal> fc_density;
    FieldCell<Scal> fc_viscosity;
    FieldCell<Vect> fc_force;
    FieldCell<Scal> fc_volume_source;
    FieldCell<Scal> fc_mass_source;
  };

  class ConditionCellReferenceToGlobal :
      public fluid_condition::GivenVelocityAndPressureFixed<Mesh> {
   public:
    ConditionCellReferenceToGlobal(Vect velocity, Scal pressure)
        : fluid_condition::GivenVelocityAndPressureFixed<Mesh>(
            velocity, pressure)
    {}
  };

  class ConditionFaceReferenceToGlobal :
      public fluid_condition::NoSlipWallFixed<Mesh> {
   public:
    ConditionFaceReferenceToGlobal(Vect velocity)
        : fluid_condition::NoSlipWallFixed<Mesh>(velocity)
    {}
  };

  std::vector<LocalData> local_;

 public:
  template <class T>
  void CopyToLocal(const FieldCell<T>& global_field,
                   FieldCell<T>& local_field,
                   const LocalData& local_data) {
    local_field.Reinit(local_data.mesh);
    for (auto lidx : local_data.mesh.Cells()) {
      local_field[lidx] = global_field[local_data.fc_global[lidx]];
    }
  }
  template <class T>
  void CopyToGlobal(const FieldCell<T>& local_field,
                    const LocalData& local_data,
                    FieldCell<T>& global_field) {
    for (auto idxcell : local_data.mesh.Cells()) {
      if (local_data.fc_active[idxcell]) {
        global_field[local_data.fc_global[idxcell]] = local_field[idxcell];
      }
    }
  }
  template <class T>
  void CopyToGlobal(const FieldFace<T>& local_field,
                    const LocalData& local_data,
                    FieldFace<T>& global_field) {
    for (auto idxface : local_data.mesh.Faces()) {
      if (local_data.ff_active[idxface]) {
        global_field[local_data.ff_global[idxface]] = local_field[idxface];
      }
    }
  }
  void CopyFields(Layers layer_src, Layers layer_dest) {
    fc_velocity_.Get(layer_dest) = fc_velocity_.Get(layer_src);
    fc_pressure_.Get(layer_dest) = fc_pressure_.Get(layer_src);
    ff_vol_flux_.Get(layer_dest) = ff_vol_flux_.Get(layer_src);
  }
  void CopyFieldsToGlobal(Layers layer) {
    for (auto& local : local_) {
      CopyToGlobal(local.fluid_solver->GetVelocity(layer), local,
                   fc_velocity_.Get(layer));
      CopyToGlobal(local.fluid_solver->GetPressure(layer), local,
                   fc_pressure_.Get(layer));
      CopyToGlobal(local.fluid_solver->GetVolumeFlux(layer), local,
                   ff_vol_flux_.Get(layer));
    }
  }
  void CopyFluidPropertiesToLocal() {
    for (auto& local : local_) {
      CopyToLocal(*this->p_fc_density_, local.fc_density, local);
      CopyToLocal(*this->p_fc_viscosity_, local.fc_viscosity, local);
      CopyToLocal(*this->p_fc_force_, local.fc_force, local);
      CopyToLocal(*this->p_fc_volume_source_, local.fc_volume_source, local);
      CopyToLocal(*this->p_fc_mass_source_, local.fc_mass_source, local);
    }
  }
  void InitLocalMeshes(size_t num_local_mesh, size_t overlap_width) {
    using namespace geom::geom2d;
    using namespace geom;

    local_.resize(num_local_mesh);

    MIdx global_bcells_size = global_mesh.GetBlockCells().GetDimensions();
    for (IntIdx part = 0; part < IntIdx(num_local_mesh); ++part) {
      MIdx local_active_cells_begin(
          part * global_bcells_size[0] / IntIdx(num_local_mesh), 0);
      MIdx local_active_cells_end(
          (part + 1) * global_bcells_size[0] / IntIdx(num_local_mesh),
          global_bcells_size[1]);

      std::cout << "\npart = " << part << std::endl;
      std::cout << "begin = " << local_active_cells_begin
          << " end = " << local_active_cells_end << std::endl;

      MIdx local_cells_begin(
          std::max<IntIdx>(
              0, local_active_cells_begin[0] - overlap_width),
          local_active_cells_begin[1]);
      MIdx local_cells_end(
          std::min<IntIdx>(
              global_bcells_size[0], local_active_cells_end[0] + overlap_width),
          local_active_cells_end[1]);

      // Initialize local nodes and mesh
      BlockNodes local_bnodes(local_cells_end - local_cells_begin + MIdx(1, 1));
      BlockNodes global_bnodes = global_mesh.GetBlockNodes();
      FieldNode<Vect> local_fn_node(local_bnodes);
      for (auto local_midx : local_bnodes) {
        IdxNode local_idxnode = local_bnodes.GetIdx(local_midx);
        MIdx global_midx = local_cells_begin + local_midx;
        IdxNode global_idxnode = global_bnodes.GetIdx(global_midx);
        local_fn_node[local_idxnode] = global_mesh.GetNode(global_idxnode);
      }
      local_[part].mesh = Mesh(local_bnodes, local_fn_node);
      auto& local = local_[part];

      // Fill references to global cells
      BlockCells local_bcells = local.mesh.GetBlockCells();
      BlockCells global_bcells = global_mesh.GetBlockCells();
      local.fc_global.Reinit(local.mesh);
      local.fc_active.Reinit(local.mesh);
      for (IdxCell local_idxcell : local.mesh.Cells()) {
        MIdx local_midx = local_bcells.GetMIdx(local_idxcell);
        MIdx global_midx = local_cells_begin + local_midx;
        IdxCell global_idxcell = global_bcells.GetIdx(global_midx);
        local.fc_global[local_idxcell] = global_idxcell;
        local.fc_active[local_idxcell] =
            (local_active_cells_begin <= global_midx &&
                global_midx < local_active_cells_end);
//        std::cout << "cell:"
//            << "  local: " << local_midx << " " << local_idxcell.GetRaw()
//            << "; global: " << global_midx << " " << global_idxcell.GetRaw()
//            << "; " << (local.fc_active[local_idxcell] ? "active" : "overlap")
//            << std::endl;
      }

      // Fill references to global faces
      BlockFaces local_bfaces = local.mesh.GetBlockFaces();
      BlockFaces global_bfaces = global_mesh.GetBlockFaces();
      local.ff_global.Reinit(local.mesh);
      for (IdxFace local_idxface : local.mesh.Faces()) {
        MIdx local_midx = local_bfaces.GetMIdx(local_idxface);
        Direction dir = local_bfaces.GetDirection(local_idxface);
        MIdx global_midx = local_cells_begin + local_midx;
        IdxFace global_idxface = global_bfaces.GetIdx(global_midx, dir);
        local.ff_global[local_idxface] = global_idxface;
//        std::cout << "face: " << GetMIdx(dir)
//            << "  local: " << local_midx << " " << local_idxface.GetRaw()
//            << "; global: " << global_midx << " " << global_idxface.GetRaw()
//            << std::endl;
      }

      // Traverse all active cells and mark their neighbour faces active
      local.ff_active.Reinit(local.mesh, false);
      for (IdxCell local_idxcell : local.mesh.Cells()) {
        if (local.fc_active[local_idxcell]) {
          for (size_t k = 0;
              k < local.mesh.GetNumNeighbourFaces(local_idxcell); ++k) {
            local.ff_active[local.mesh.GetNeighbourFace(local_idxcell, k)] =
                true;
          }
        }
      }
    }
  }
  void UpdateReferencesToGlobal() {
    IdxCell idxcell_pressure_fixed(0);

    for (auto& local : local_) {
      for (auto it = local.mc_cond.begin();
          it != local.mc_cond.end(); ++it) {
        auto idxcell = it->GetIdx();
        auto cond = it->GetValue().get();
        if (auto cond_ref =
            dynamic_cast<ConditionCellReferenceToGlobal*>(cond)) {
          (*cond_ref) = ConditionCellReferenceToGlobal(
              fc_velocity_.iter_curr[local.fc_global[idxcell]],
              fc_pressure_.iter_curr[local.fc_global[idxcell]] -
              fc_pressure_.iter_curr[idxcell_pressure_fixed]);
        }
      }

      for (auto it = local.mf_cond.begin();
          it != local.mf_cond.end(); ++it) {
        auto idxface = it->GetIdx();
        auto cond = it->GetValue().get();
        if (auto cond_wall =
            dynamic_cast<ConditionFaceReferenceToGlobal*>(cond)) {


          (*cond_wall) = ConditionFaceReferenceToGlobal(
              GetInterpolatedInner(fc_velocity_.iter_curr,
                                   local.ff_global[idxface],
                                   global_mesh));
        }
      }
    }
  }
  FluidSimpleParallel(
        const Mesh& mesh
      , const geom::FieldCell<Vect>& fc_velocity_initial
      , const geom::MapFace<std::shared_ptr<ConditionFaceFluid>>&
        mf_cond
      , const geom::MapCell<std::shared_ptr<ConditionCellFluid>>&
        mc_cond
      , Scal velocity_relaxation_factor
      , Scal pressure_relaxation_factor
      , geom::FieldCell<Scal>* p_fc_density
      , geom::FieldCell<Scal>* p_fc_viscosity
      , geom::FieldCell<Vect>* p_fc_force
      , geom::FieldCell<Scal>* p_fc_volume_source
      , geom::FieldCell<Scal>* p_fc_mass_source
      , double time, double time_step
      , const LinearSolverFactory& linear_factory_velocity
      , const LinearSolverFactory& linear_factory_pressure
      , double convergence_tolerance
      , size_t num_iterations_limit
      , MultiTimer<std::string>* timer
      , size_t num_local_mesh
      , size_t overlap_width
      )
      : FluidSolver<Mesh>(time, time_step, p_fc_density, p_fc_viscosity,
                    p_fc_force, p_fc_volume_source, p_fc_mass_source,
                    convergence_tolerance, num_iterations_limit)
      , global_mesh(mesh)
      , mf_cond_(mf_cond)
      , mc_cond_(mc_cond)
      , velocity_relaxation_factor_(velocity_relaxation_factor)
      , pressure_relaxation_factor_(pressure_relaxation_factor)
      , timer_(timer)
  {
    using namespace fluid_condition;

    fc_velocity_.time_curr = fc_velocity_initial;
    fc_pressure_.time_curr.Reinit(global_mesh);
    ff_vol_flux_.time_curr.Reinit(global_mesh);
    CopyFields(Layers::time_curr, Layers::iter_prev);
    CopyFields(Layers::time_curr, Layers::iter_curr);
    CopyFields(Layers::time_curr, Layers::time_prev);


    // Init local meshes
    InitLocalMeshes(num_local_mesh, overlap_width);

    for (auto& local : local_) {
      FieldCell<Vect> local_velocity_initial;
      CopyToLocal(fc_velocity_initial, local_velocity_initial, local);

      for (IdxFace idxface : local.mesh.Faces()) {
        if (auto ptr = mf_cond_.find(local.ff_global[idxface])) {
          local.mf_cond[idxface] = *ptr;
        } else if (!local.mesh.IsInner(idxface)) {
          local.mf_cond[idxface] =
              std::make_shared<ConditionFaceReferenceToGlobal>(Vect(0., 0.));
//          auto idxcell = local.mesh.GetNeighbourCell(
//              idxface, local.mesh.GetValidNeighbourCellId(idxface));
//          local.mc_cond[idxcell] =
//              std::make_shared<ConditionCellReferenceToGlobal>(
//                  Vect(0., 0.), 0.);
        }

        // Add fixed pressure condition
        local.mc_cond[IdxCell(0)] =
            std::make_shared<solver::fluid_condition::
            GivenPressureFixed<Mesh>>(0.);
      }

      for (IdxCell idxcell : local.mesh.Cells()) {
        if (auto ptr = mc_cond_.find(local.fc_global[idxcell])) {
          local.mc_cond[idxcell] = *ptr;
        }
      }

      std::cout << "set " << local.mf_cond.size() << " conditions"
          << std::endl;

      local.fluid_solver =
          std::make_shared<solver::FluidSimple<Mesh>>(
          local.mesh, local_velocity_initial,
          local.mf_cond, local.mc_cond,
          velocity_relaxation_factor,
          pressure_relaxation_factor,
          &local.fc_density, &local.fc_viscosity, &local.fc_force,
          &local.fc_volume_source, &local.fc_mass_source,
          this->GetTime(), this->GetTimeStep(),
          linear_factory_velocity, linear_factory_pressure,
          convergence_tolerance, num_iterations_limit,
          timer_);
    }


  }
  void StartStep() override {
    this->ClearIterationCount();

    for (auto& local : local_) {
      local.fluid_solver->StartStep();
    }

    CopyFieldsToGlobal(Layers::iter_curr);
  }
  void MakeIteration() override {
    CopyFluidPropertiesToLocal();
    CopyFields(Layers::iter_curr, Layers::iter_prev);

    for (auto& local : local_) {
      local.fluid_solver->MakeIteration();
    }

    CopyFieldsToGlobal(Layers::iter_curr);
    UpdateReferencesToGlobal();
    this->IncIterationCount();
  }
  void FinishStep() override {
    for (auto& local : local_) {
      local.fluid_solver->FinishStep();
    }
    CopyFields(Layers::time_curr, Layers::time_prev);
    CopyFields(Layers::iter_curr, Layers::time_curr);

    if (IsNan(fc_pressure_.time_curr)) {
      throw std::string("Invalid solution");
    }
    this->IncTime();
  }
  double GetConvergenceIndicator() const {
    double res = 0.;
    for (auto& local : local_) {
      res = std::max(res, local.fluid_solver->GetConvergenceIndicator());
    }
    return res;
  }
  const geom::FieldCell<Vect>& GetVelocity() override {
    return fc_velocity_.time_curr;
  }
  const geom::FieldCell<Vect>& GetVelocity(Layers layer) override {
    return fc_velocity_.Get(layer);
  }
  const geom::FieldCell<Scal>& GetPressure() override {
    return fc_pressure_.time_curr;
  }
  const geom::FieldCell<Scal>& GetPressure(Layers layer) override {
    return fc_pressure_.Get(layer);
  }
  const geom::FieldFace<Scal>& GetVolumeFlux() override {
    return ff_vol_flux_.time_curr;
  }
  const geom::FieldFace<Scal>& GetVolumeFlux(Layers layer) override {
    return ff_vol_flux_.Get(layer);
  }
};


} // namespace solver
