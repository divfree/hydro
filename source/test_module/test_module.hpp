/*******************************************************************
******************       CFD SOLVER 2D        **********************
********************************************************************/

#pragma once

#include "../control/experiment.hpp"
#include "../common/vect.hpp"
#include "../hydro2dmpi/mesh.hpp"
#include "../hydro2dmpi/mesh2d.hpp"
#include "../hydro2dmpi/output.hpp"
#include "../hydro2dmpi/linear.hpp"
#include "../hydro2dmpi/solver.hpp"
#include <memory>

// TODO: Change parameters on events (e.g. certain time moments)

namespace test_module
{

template <class Vect>
Vect GetVect(const column<double>& v);

template <>
geom::Vect<double, 2> GetVect<geom::Vect<double, 2>>(const column<double>& v) {
  return geom::Vect<double, 2>(v[0], v[1]);
}


template <>
geom::Vect<float, 2> GetVect<geom::Vect<float, 2>>(const column<double>& v) {
  return geom::Vect<float, 2>(float(v[0]), float(v[1]));
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

  using Expr = solver::Expression<Scal, IdxCell, 1 + dim * 2>;
  std::shared_ptr<solver::LinearSolver<Scal, IdxCell, Expr>> linear_;
  FieldCell<Expr> fc_system_, fc_interpolation_;
  FieldCell<solver::Expression<Scal, IdxCell, 9>>
  fc_restriction_;
  FieldCell<solver::Expression<Scal, IdxCell, 180>>
  fc_system_coarse_;
 public:
  hydro(TExperiment* _ex);
  ~hydro() {}
  void step();
  void write_results(bool force=false);
  Mesh mesh, mesh_output;
  output::Content content, content_scalar;
  std::shared_ptr<output::Session> session, session_scalar;
  FieldCell<IdxCell> output_to_mesh_;
  FieldCell<IdxCell> mesh_to_coarse_;

  FieldCell<Scal> fc_u_, fc_u_coarse_;
  FieldCell<Scal> fc_f_;

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

  // Prepare output mesh
  size_t output_factor = 4;
  geom::InitUniformMesh(mesh_output, domain, mesh_size * output_factor);
  output_to_mesh_.Reinit(mesh_output);
  for (auto idxcell : mesh_output.Cells()) {
    MIdx midxoutput = mesh_output.GetBlockCells().GetMIdx(idxcell);
    output_to_mesh_[idxcell] =
        mesh.GetBlockCells().GetIdx(midxoutput / output_factor);
  }

  P_int.set("cells_number", static_cast<int>(mesh.GetNumCells()));

/*
  switch (linear_id) {
    case solver::LinearSolverId::lu: {
      linear_ = std::make_shared<
          solver::LuDecomposition<Scal, IdxCell, Expr>>();
      break;
    }
    case solver::LinearSolverId::gauss_seidel: {
      linear_ = std::make_shared<
          solver::GaussSeidel<Scal, IdxCell, Expr>>();
      break;
    }
    case solver::LinearSolverId::pardiso: {
      linear_ = std::make_shared<
          solver::Pardiso<Scal, IdxCell, Expr>>();
      break;
    }
  }*/

/*
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
*/

  // Initial conditions and RHS
  fc_u_.Reinit(mesh, 0.);
  fc_f_.Reinit(mesh);
  const Scal freqhigh = static_cast<Scal>(P_double["freqhigh"]);
  const Scal freqlow = static_cast<Scal>(P_double["freqlow"]);
  const Scal ahigh = static_cast<Scal>(P_double["ahigh"]);
  const Scal alow = static_cast<Scal>(P_double["alow"]);
  const Scal pi = static_cast<Scal>(M_PI);
  for (auto idxcell : mesh.Cells()) {
    Vect x = mesh.GetCenter(idxcell);
    Vect mhigh = x * 2 * pi / domain_size * freqhigh;
    Vect mlow = x * 2 * pi / domain_size * freqlow;
    fc_f_[idxcell] =
        ahigh * std::sin(mhigh[0]) * std::sin(mhigh[1]) +
        alow * std::sin(mlow[0]) * std::sin(mlow[1]);
  }

  // Init system
  fc_system_.Reinit(mesh);
  solver::DerivativeInnerFacePlain<Mesh, Expr> derivative(mesh);
  for (auto idxcell : mesh.Cells()) {
    Expr& eqn = fc_system_[idxcell];
    if (mesh.IsInner(idxcell)) {
      for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
        auto idxface = mesh.GetNeighbourFace(idxcell, i);
        eqn += derivative.GetExpression(idxface) * mesh.GetArea(idxface) *
            mesh.GetOutwardFactor(idxcell, i);
      }
      eqn.SetConstant(-fc_f_[idxcell] * mesh.GetVolume(idxcell));
    } else {
      eqn.InsertTerm(1., idxcell);
    }
  }

  // Coarse grid cell subset
  fc_u_coarse_.Reinit(mesh, 0.);
  const auto& b_cells = mesh.GetBlockCells();
  std::vector<IdxCell> v_coarse_cells;
  for (auto midx : b_cells) {
    if (midx[0] % 2 == 0 && midx[1] % 2 == 0) {
      v_coarse_cells.push_back(b_cells.GetIdx(midx));
    }
  }

  // Init nearset-neighbour refs
  mesh_to_coarse_.Reinit(mesh, v_coarse_cells[0]);
  for (auto idxcell : mesh.Cells()) {
    Vect x = mesh.GetCenter(idxcell);
    for (auto idxcoarse : v_coarse_cells) {
      Vect dold = mesh.GetCenter(mesh_to_coarse_[idxcell]) - x;
      Vect dnew = mesh.GetCenter(idxcoarse) - x;
      if (dnew[0] >=0 && dnew[1] >= 0 &&
          (!(dold[0] >=0 && dold[1] >= 0) ||
          std::max(std::abs(dnew[0]), std::abs(dnew[1])) <
          std::max(std::abs(dold[0]), std::abs(dold[1])))) {
        mesh_to_coarse_[idxcell] = idxcoarse;
      }
    }
  }

  // Init interpolation operator
  fc_interpolation_.Reinit(mesh);
  fc_restriction_.Reinit(mesh);
  for (auto idxcenter : v_coarse_cells) {
    MIdx center = b_cells.GetMIdx(idxcenter);
    geom::BlockGeneric<size_t, dim> b_area(MIdx(-1), MIdx(3));
    for (auto offset : b_area) {
      MIdx side = center + offset;
      if (b_cells.IsInside(side)) {
        IdxCell idxside = b_cells.GetIdx(side);
        if (idxside == idxcenter) {
          fc_interpolation_[idxside].InsertTerm(
              Scal(b_area.size() - 1), idxcenter);
        } else {
          fc_interpolation_[idxside].InsertTerm(
              1., idxcenter);
        }
      }
    }
  }
  Normalize(fc_interpolation_);

  // Init restriction operator
  solver::Transpose(fc_interpolation_, fc_restriction_);
  Normalize(fc_restriction_);

  for (auto idxcell : mesh.Cells()) {
    fc_u_coarse_[idxcell] = fc_restriction_[idxcell].Evaluate(fc_u_);
  }

  // Init coarse system
  fc_system_coarse_.Reinit(mesh);
  auto tmp = fc_system_coarse_;
  decltype(fc_system_coarse_) fc_interpolation_ext(fc_interpolation_);
  for (auto idxcell : mesh.Cells()) {
    tmp[idxcell] = fc_system_[idxcell].Evaluate(fc_interpolation_ext);
  }
  for (auto idxcell : mesh.Cells()) {
    fc_system_coarse_[idxcell] = fc_restriction_[idxcell].Evaluate(tmp);
  }
  std::cout << "fc_system_coarse_" << std::endl;
  for (auto idxcell : mesh.Cells()) {
    std::cout << fc_system_coarse_[idxcell] << std::endl;
  }

  // TODO: EntryFunctionField
  content = {
      std::make_shared<output::EntryFunction<Scal, IdxNode, Mesh>>(
          "x", mesh_output,
          [this](IdxNode idx) { return mesh_output.GetNode(idx)[0]; }),
      std::make_shared<output::EntryFunction<Scal, IdxNode, Mesh>>(
          "y", mesh_output,
          [this](IdxNode idx) { return mesh_output.GetNode(idx)[1]; }),
      std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
          "u", mesh_output,
          [this](IdxCell idx) { return fc_u_[output_to_mesh_[idx]]; }),
      std::make_shared<output::EntryFunction<Scal, IdxCell, Mesh>>(
          "u_coarse", mesh_output,
          [this](IdxCell idx) {
              return fc_u_coarse_[mesh_to_coarse_[output_to_mesh_[idx]]]; })
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

  session = std::make_shared<output::SessionTecplotBinaryStructured<Mesh>>(
      content, P_string[_plt_title], mesh_output, P_string["filename_field"]);

  session_scalar = std::make_shared<output::SessionTecplotAsciiScalar<Scal>>(
      content_scalar, P_string[_plt_title],
      "zone_title", P_string["filename_scalar"]);

  last_frame_time_ = 0;
  last_frame_scalar_time_ = 0;
  session->Write(0., "initial");
  session_scalar->Write();
}

template <class Mesh>
void hydro<Mesh>::StepMultigrid() {
  Scal relaxation_fine = P_double["relaxation_fine"];
  Scal relaxation_coarse = P_double["relaxation_coarse"];

  // Calc residual
  FieldCell<Scal> residual(mesh);
  for (auto idxcell : mesh.Cells()) {
    residual[idxcell] = fc_system_[idxcell].Evaluate(fc_u_);
  }

  // Interpolate residual to coarse
  FieldCell<Scal> residual_coarse(mesh);
  for (auto idxcell : mesh.Cells()) {
    residual_coarse[idxcell] = fc_restriction_[idxcell].Evaluate(residual);
  }

  // Calc delta for coarse
  FieldCell<Scal> delta_coarse(mesh);
  for (auto idxcell : mesh.Cells()) {
    auto& eqn = fc_system_coarse_[idxcell];
    Scal sum = 0.;
    Scal diag_coeff = 0.;
    for (size_t i = 0; i < eqn.size(); ++i) {
      auto& term = eqn[i];
      if (term.idx.GetRaw() < idxcell.GetRaw()) {
        sum += term.coeff * delta_coarse[term.idx];
      } else if (term.idx == idxcell) {
        diag_coeff = term.coeff;
      }
    }
    delta_coarse[idxcell] = -(residual_coarse[idxcell] + sum) / diag_coeff;
    delta_coarse[idxcell] *= relaxation_coarse;
  }

  // First correction
  for (auto idxcell : mesh.Cells()) {
    fc_u_[idxcell] += fc_interpolation_[idxcell].Evaluate(delta_coarse);
  }

  // Second residual
  for (auto idxcell : mesh.Cells()) {
    residual[idxcell] = fc_system_[idxcell].Evaluate(fc_u_);
  }

  // Second correction
  FieldCell<Scal> delta(mesh);
  for (auto idxcell : mesh.Cells()) {
    auto& eqn = fc_system_[idxcell];
    Scal sum = 0.;
    Scal diag_coeff = 0.;
    for (size_t i = 0; i < eqn.size(); ++i) {
      auto& term = eqn[i];
      if (term.idx.GetRaw() < idxcell.GetRaw()) {
        sum += term.coeff * delta[term.idx];
      } else if (term.idx == idxcell) {
        diag_coeff = term.coeff;
      }
    }
    delta[idxcell] = -(residual[idxcell] + sum) / diag_coeff;
    delta[idxcell] *= relaxation_fine;
  }

  // Second correction
  for (auto idxcell : mesh.Cells()) {
    fc_u_[idxcell] += delta[idxcell];
  }

  for (auto idxcell : mesh.Cells()) {
    fc_u_coarse_[idxcell] = fc_restriction_[idxcell].Evaluate(fc_u_);
  }

}

template <class Scal, class Idx, class Expr>
class LuDecompositionLocal : public solver::LinearSolver<Scal, Idx, Expr> {
  template <class T>
  using Field = geom::FieldGeneric<T, Idx>;

 public:
  Field<Scal> Solve(const Field<Expr>& system) override {
    Field<Scal> res(system.GetRange());

    // forward step
    for(size_t i = 0; i < system.size(); ++i) {
      Scal sum = 0;
      const auto &expr = system[Idx(i)];
      size_t k = 0;
      while (k < expr.size() && expr[k].idx.GetRaw() < i) {
        sum += expr[k].coeff * res[expr[k].idx];
        ++k;
      }
      assert(k < expr.size() && expr[k].idx.GetRaw() == i);
      Scal coeff_diag = expr[k].coeff;
      res[Idx(i)] = (-expr.GetConstant() - sum) / coeff_diag;
    }

    // backward step
    for (size_t i = system.size(); i > 0; ) {
      --i;
      Scal sum = 0;
      const auto &expr = system[Idx(i)];
      size_t k = expr.size() - 1;
      while (k >= 0 && expr[k].idx.GetRaw() > i) {
        sum += expr[k].coeff * res[expr[k].idx];
        --k;
      }
      assert(k >= 0 && expr[k].idx.GetRaw() == i);
      Scal coeff_diag = expr[k].coeff;
      res[Idx(i)] -= sum / coeff_diag;
    }

    return res;
  }
};

template <class Mesh>
void hydro<Mesh>::step() {
  ex->timer_.Push("step");

  // StepMultigrid();

  linear_ = std::make_shared<
      LuDecompositionLocal<Scal, IdxCell, Expr>>();

  FieldCell<Scal> tmp(mesh);

  for (auto idxcell : mesh.Cells()) {
    auto& expr = fc_system_[idxcell];
    tmp[idxcell] = expr.GetConstant();
    expr.SetConstant(expr.Evaluate(fc_u_));
  }
  auto corr = linear_->Solve(fc_system_);

  // Restore constant terms
  for (auto idxcell : mesh.Cells()) {
    auto& expr = fc_system_[idxcell];
    expr.SetConstant(tmp[idxcell]);
  }

  for (auto idxcell : mesh.Cells()) {
    fc_u_[idxcell] += corr[idxcell];
  }

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
    session->Write(time, "step");
    flog<<"Frame "<<(P_int["current_frame"]++)<<": t="<<time<<endl;
  }

  const size_t max_frame_scalar_index = P_int["max_frame_scalar_index"];
  const double frame_scalar_duration = total_time / max_frame_scalar_index;

  if (force || time >= last_frame_scalar_time_ + frame_scalar_duration) {
    last_frame_scalar_time_ = time;
    session_scalar->Write();
    flog<<"Frame_scalar "<<(P_int["current_frame_scalar"]++)<<": t="<<time<<endl;
  }
}

}
