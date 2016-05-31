/*
 * heat.hpp
 *
 *  Created on: Jan 31, 2016
 *      Author: Petr Karnakov (petr.karnakov@gmail.com)
 */

#pragma once

#include "mesh.hpp"
#include "linear.hpp"
#include <exception>
#include "solver.hpp"
#include "conv_diff.hpp"

namespace solver {

template <class Mesh>
class HeatSolver : public UnsteadyIterativeSolver {
  const Mesh& mesh;
  static constexpr size_t dim = Mesh::dim;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using Solver = ConvectionDiffusionScalarImplicit<Mesh>;
  geom::MapFace<std::shared_ptr<ConditionFace>> mf_conductivity_cond_shared_;
  geom::MapFace<ConditionFace*> mf_conductivity_cond_;
  std::shared_ptr<Solver> solver_;
  geom::FieldFace<Scal> ff_conductivity_;
  geom::FieldCell<Scal> fc_uniform_one_;
  const geom::FieldCell<Scal>* p_fc_conductivity_;

 public:
  HeatSolver(
      const Mesh& mesh,
      const geom::FieldCell<Scal>& fc_temperature_initial,
      const geom::MapFace<std::shared_ptr<ConditionFace>>&
      mf_temperature_cond_shared,
      Scal relaxation_factor,
      const geom::FieldCell<Scal>* p_fc_conductivity,
      const geom::FieldCell<Scal>* p_fc_source,
      const geom::FieldFace<Scal>* p_ff_vol_flux,
      double time, double time_step,
      const LinearSolverFactory& linear_factory,
      double convergence_tolerance,
      size_t num_iterations_limit,
      bool time_second_order = true)
      : UnsteadyIterativeSolver(time, time_step, convergence_tolerance,
                                num_iterations_limit),
        mesh(mesh),
        fc_uniform_one_(mesh, 1.),
        p_fc_conductivity_(p_fc_conductivity)
  {
    for (auto it = mf_temperature_cond_shared.cbegin() ;
        it != mf_temperature_cond_shared.cend();
        ++it) {
      mf_conductivity_cond_shared_[it->GetIdx()] =
          std::make_shared<solver::ConditionFaceDerivativeFixed<Scal>>(0.);
    }

    mf_conductivity_cond_ = GetPointers(mf_conductivity_cond_shared_);

    // Initialize solver
    solver_ = std::make_shared<Solver>(
        mesh, fc_temperature_initial,
        mf_temperature_cond_shared,
        geom::MapCell<std::shared_ptr<ConditionCell>>() /*empty*/,
        relaxation_factor, &fc_uniform_one_, &ff_conductivity_,
        p_fc_source, p_ff_vol_flux, time, time_step,
        linear_factory, convergence_tolerance,
        num_iterations_limit, time_second_order);
  }
  void StartStep() override {
    this->ClearIterationCount();
    solver_->StartStep();
  }
  void MakeIteration() override {
    ff_conductivity_ =
        Interpolate(*p_fc_conductivity_, mf_conductivity_cond_, mesh);

    solver_->MakeIteration();

    this->IncIterationCount();
  }
  void FinishStep() override {
    this->IncTime();
    solver_->FinishStep();
  }
  double GetConvergenceIndicator() const {
    return solver_->GetConvergenceIndicator();
  }
  const geom::FieldCell<Scal>& GetTemperature() const {
    return solver_->GetField();
  }
  const geom::FieldCell<Scal>& GetTemperature(Layers layer) const {
    return solver_->GetField(layer);
  }
};

} // namespace solver
