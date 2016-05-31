/*
 * chemistry.hpp
 *
 *  Created on: Apr 5, 2016
 *      Author: Petr Karnakov (petr.karnakov@gmail.com)
 */

#pragma once

#include <vector>
#include <cmath>
#include <exception>
#include "solver.hpp"
#include "advection.hpp"

namespace solver {

template <class Mesh>
class Kinetics {
  using IdxCell = geom::IdxCell;
  using Scal = typename Mesh::Scal;
 public:
  virtual std::vector<Scal> GetReactionRate(
      IdxCell idxcell,
      const std::vector<Scal>& molar_concentration) = 0;
  virtual ~Kinetics() {}
};

template <class Mesh>
class KineticsSteady : public Kinetics<Mesh> {
  using IdxCell = geom::IdxCell;
  using Scal = typename Mesh::Scal;
  const Mesh& mesh;
  Scal rate_;
  std::vector<Scal> v_molar_mass_;
  const size_t num_phases_;
 public:
  KineticsSteady(Scal rate, const std::vector<Scal>& v_molar_mass,
                 const Mesh& mesh)
      : mesh(mesh)
      , rate_(rate)
      , v_molar_mass_(v_molar_mass)
      , num_phases_(v_molar_mass_.size())
  {}
  std::vector<Scal> GetReactionRate(
      IdxCell,
      const std::vector<Scal>& v_molar_concentration) override {
    std::vector<Scal> res(num_phases_);

    // First phase converts to the others
    res[0] = -rate_ * v_molar_concentration[0];

    Scal raw_mass_rate = res[0] * v_molar_mass_[0];
    Scal product_mass_rate = -raw_mass_rate / (num_phases_ - 1);

    for (size_t i = 1; i < num_phases_; ++i) {
      res[i] = product_mass_rate / v_molar_mass_[i];
    }

    return res;
  }
};


template <class Mesh>
class KineticsSteadyRadiation : public Kinetics<Mesh> {
  using IdxCell = geom::IdxCell;
  using Scal = typename Mesh::Scal;
  const Mesh& mesh;
  Scal rate_;
  std::vector<Scal> v_molar_mass_;
  const size_t num_phases_;
  const geom::FieldCell<Scal>& fc_radiation_;
 public:
  KineticsSteadyRadiation(
      Scal rate, const std::vector<Scal>& v_molar_mass,
      const geom::FieldCell<Scal>& fc_radiation,
      const Mesh& mesh)
      : mesh(mesh)
      , rate_(rate)
      , v_molar_mass_(v_molar_mass)
      , num_phases_(v_molar_mass_.size())
      , fc_radiation_(fc_radiation)
  {}
  std::vector<Scal> GetReactionRate(
      IdxCell idxcell,
      const std::vector<Scal>& v_molar_concentration) override {
    std::vector<Scal> res(num_phases_);

    // First phase converts to the others
    res[0] = -rate_ * v_molar_concentration[0] * fc_radiation_[idxcell];

    Scal raw_mass_rate = res[0] * v_molar_mass_[0];
    Scal product_mass_rate = -raw_mass_rate / (num_phases_ - 1);

    for (size_t i = 1; i < num_phases_; ++i) {
      res[i] = product_mass_rate / v_molar_mass_[i];
    }

    return res;
  }
};



template <class Mesh>
class KineticsNadirov : public Kinetics<Mesh> {
  using IdxCell = geom::IdxCell;
  using Scal = typename Mesh::Scal;
  const Mesh& mesh;
  std::vector<Scal> v_molar_mass_;
  const size_t num_phases_;
  Scal k_P, k, I_0, E_T, k_B, T, G, D;
 public:
  KineticsNadirov(const std::vector<Scal>& v_molar_mass,
                  Scal k_P, Scal k, Scal I_0, Scal E_T,
                  Scal k_B, Scal T, Scal G, Scal D,
                  const Mesh& mesh)
      : mesh(mesh)
      , v_molar_mass_(v_molar_mass)
      , num_phases_(v_molar_mass_.size())
      , k_P(k_P), k(k), I_0(I_0), E_T(E_T)
      , k_B(k_B), T(T), G(G), D(D)
  {
    if (num_phases_ != 3) {
      std::runtime_error("KineticsNadirov: require 3 phases");
    }
  }
  std::vector<Scal> GetReactionRate(
      IdxCell,
      const std::vector<Scal>& v_molar_concentration) override {
    std::vector<Scal> res(num_phases_);

    Scal I_T = I_0 * std::exp(-E_T / (k_B * T));
    Scal I_P = G * D * 0.01;
    Scal K = k_P * std::sqrt((I_T + I_P) / k);
    const auto& N = v_molar_concentration;
    const auto& M = v_molar_mass_;
    enum { Heavy = 0, Light = 1, Gas = 2 };

    res[Heavy] = -K * (N[Heavy] - (N[Light] + N[Gas]));
    res[Light] = -K * (N[Light] - N[Gas]) +
        (M[Light] - M[Gas]) / M[Heavy] * res[Heavy];
    res[Gas] = -(res[Heavy] * M[Heavy] / M[Gas] +
        res[Light] * M[Light] / M[Gas]);

    return res;
  }
};

template <class Mesh, class Scal = typename Mesh::Scal,
    class Vect = typename Mesh::Vect>
geom::FieldCell<Scal> CalcRadiationField(
    const Mesh& mesh,
    const geom::MapFace<std::shared_ptr<solver::ConditionFace>>&
    mf_cond_radiation_shared,
    const geom::FieldCell<Scal>& fc_absorption_rate,
    Vect radiation_direction,
    geom::FieldCell<Scal> initial_guess) {

  if (radiation_direction.norm() == 0.) {
    return initial_guess;
  }

  radiation_direction /= radiation_direction.norm();

  geom::FieldFace<Scal> ff_flux(mesh);
  for (auto idxface : mesh.Faces()) {
    ff_flux[idxface] =
        radiation_direction.dot(mesh.GetSurface(idxface));
  }

  Scal time_step = 1e16;
  for (auto idxcell : mesh.Cells()) {
    for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
      auto idxface = mesh.GetNeighbourFace(idxcell, i);
      Vect vect = mesh.GetCenter(idxcell) - mesh.GetCenter(idxface);
      Scal vel_dot_vect = std::abs(radiation_direction.dot(vect));
      if (vel_dot_vect != 0.) {
        time_step = std::min(time_step, vect.sqrnorm() / vel_dot_vect);
      }
    }
  }

  geom::FieldCell<Scal> fc_source(mesh, 0.);

  solver::AdvectionSolverExplicit<Mesh, geom::FieldFace<Scal>>
  advection(
      mesh, initial_guess, mf_cond_radiation_shared, &ff_flux,
      &fc_source, 0., time_step);

  Scal error = 1e16;
  size_t count = 0;
  for (; count < 10 && error > 1e-3; ++count) {
    for (auto idxcell : mesh.Cells()) {
      fc_source[idxcell] =
          -fc_absorption_rate[idxcell] * advection.GetField()[idxcell];
    }
    advection.StartStep();
    advection.MakeIteration();
    error = advection.GetConvergenceIndicator();
    advection.FinishStep();
  }
  return advection.GetField();
}

} // namespace solver
