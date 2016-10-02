/*
 * mesh1d.hpp
 *
 *  Created on: Oct 02, 2016
 *      Author: Petr Karnakov (petr.karnakov@gmail.com)
 */

#pragma once

#include "mesh.hpp"

namespace geom {

namespace geom1d {


template <class Scal>
class MeshStructured : public MeshGeneric<Scal, 1> {
 public:
  static constexpr size_t dim = 1;
  using Vect = geom::Vect<Scal, dim>;
  using Direction = geom::Direction<dim>;
  using MIdx = geom::Vect<IntIdx, dim>;
  using BlockNodes = geom::BlockNodes<dim>;
  using BlockCells = geom::BlockCells<dim>;
  using BlockFaces = geom::BlockFaces<dim>;
  static constexpr size_t kCellNumNeighbourFaces = 2;
  static constexpr size_t kCellNumNeighbourNodes = 2;
  static constexpr size_t kFaceNumNeighbourNodes = 1;
  static constexpr size_t kFaceNumNeighbourCells = 2;

 private:
  // b:Block, fc:FieldCell, ff:FieldFace, fn:FieldNode
  BlockNodes b_nodes_;
  FieldNode<Vect> fn_node_;
  BlockCells b_cells_;
  BlockFaces b_faces_;
  FieldCell<Vect> fc_center_;
  FieldCell<Scal> fc_volume_;
  FieldFace<Vect> ff_center_;
  FieldFace<Scal> ff_area_;
  FieldFace<Vect> ff_surface_;
  std::array<IntIdx, kCellNumNeighbourFaces> cell_neighbour_cell_offset_;
  FieldCell<std::array<IdxFace, kCellNumNeighbourFaces>> fc_neighbour_face_;
  FieldCell<IdxNode> fc_neighbour_node_base_;
  std::array<IntIdx, kCellNumNeighbourNodes> cell_neighbour_node_offset_;
  FieldFace<Direction> ff_direction_;
  FieldFace<std::array<IdxCell, kFaceNumNeighbourCells>> ff_neighbour_cell_;
  FieldFace<std::array<IdxNode, kFaceNumNeighbourNodes>> ff_neighbour_node_;
  FieldCell<bool> fc_is_inner_;
  FieldFace<bool> ff_is_inner_;
  FieldCell<bool> fc_is_excluded_;
  FieldFace<bool> ff_is_excluded_;

 public:
  MeshStructured() {}
  MeshStructured(const BlockNodes& b_nodes, const FieldNode<Vect>& fn_node);
  const BlockCells& GetBlockCells() const {
    return b_cells_;
  }
  const BlockFaces& GetBlockFaces() const {
    return b_faces_;
  }
  const BlockNodes& GetBlockNodes() const {
    return b_nodes_;
  }
  Vect GetCenter(IdxCell idx) const override {
    return fc_center_[idx];
  }
  Vect GetCenter(IdxFace idx) const override {
    return ff_center_[idx];
  }
  Vect GetSurface(IdxFace idx) const override {
    return ff_surface_[idx];
  }
  Vect GetNode(IdxNode idx) const override {
    return fn_node_[idx];
  }
  Scal GetVolume(IdxCell idx) const override {
    return fc_volume_[idx];
  }
  Scal GetArea(IdxFace idx) const override {
    return ff_area_[idx];
  }
  size_t GetNumCells() const override {
    return b_cells_.size();
  }
  size_t GetNumFaces() const override {
    return b_faces_.size();
  }
  size_t GetNumNodes() const override {
    return b_nodes_.size();
  }
  // n = 0 .. 1
  IdxCell GetNeighbourCell(IdxCell idx, size_t n) const override {
    if (fc_is_inner_[idx]) {
      idx.AddRaw(cell_neighbour_cell_offset_[n]);
    } else {
      MIdx base = b_cells_.GetMIdx(idx);
      MIdx offset = (std::vector<MIdx>{MIdx{-1}, MIdx{1}})[n];
      if (b_cells_.IsInside(base + offset)) {
        idx.AddRaw(cell_neighbour_cell_offset_[n]);
      } else {
        idx = IdxCell::None();
      }
    }
    return idx;
  }
  // n same as for GetNeighbourCell()
  IdxFace GetNeighbourFace(IdxCell idxcell, size_t n) const override {
    return fc_neighbour_face_[idxcell][n];
  }
  Scal GetOutwardFactor(IdxCell, size_t n) const override {
    return (n == 0) ? -1. : 1.;
  }
  // n same as for GetNeighbourCell()
  Vect GetOutwardSurface(IdxCell idxcell, size_t n) const override {
    return GetSurface(GetNeighbourFace(idxcell, n)) *
        GetOutwardFactor(idxcell, n);
  }
  // n = 0 .. 1
  IdxNode GetNeighbourNode(IdxCell idxcell, size_t n) const override {
    IdxNode idxnode = fc_neighbour_node_base_[idxcell];
    idxnode.AddRaw(cell_neighbour_node_offset_[n]);
    return idxnode;
  }
  Direction GetDirection(IdxFace idx) const {
    return ff_direction_[idx];
  }
  IdxCell GetNeighbourCell(IdxFace idx, size_t n) const override {
    return ff_neighbour_cell_[idx][n];
  }
  IdxNode GetNeighbourNode(IdxFace idx, size_t n) const override {
    return ff_neighbour_node_[idx][n];
  }
  size_t GetNumNeighbourFaces(IdxCell) const override {
    return kCellNumNeighbourFaces;
  }
  size_t GetNumNeighbourNodes(IdxCell) const override {
    return kCellNumNeighbourNodes;
  }
  size_t GetNumNeighbourCells(IdxFace) const override {
    return kFaceNumNeighbourCells;
  }
  size_t GetNumNeighbourNodes(IdxFace) const override {
    return kFaceNumNeighbourNodes;
  }
  bool IsInner(IdxCell idxcell) const override {
    return fc_is_inner_[idxcell];
  }
  bool IsInner(IdxFace idxface) const override {
    return ff_is_inner_[idxface];
  }
  bool IsInside(IdxCell idxcell, Vect vect) const override {
    for (size_t i = 0; i < GetNumNeighbourFaces(idxcell); ++i) {
      IdxFace idxface = GetNeighbourFace(idxcell, i);
      if (GetOutwardSurface(idxcell, i).dot(vect - GetCenter(idxface)) > 0) {
        return false;
      }
    }
    return true;
  }
  void ExcludeCells(const std::vector<IdxCell>& list) {
    for (auto idxcell : list) {
      fc_is_excluded_[idxcell] = true;
      fc_is_inner_[idxcell] = false;
    }
    for (auto idxface : this->Faces()) {
      for (size_t i = 0; i < GetNumNeighbourCells(idxface); ++i) {
        IdxCell idxcell = GetNeighbourCell(idxface, i);
        if (!idxcell.IsNone() && fc_is_excluded_[idxcell]) {
          ff_is_inner_[idxface] = false;
          ff_neighbour_cell_[idxface][i] = IdxCell::None();
        }
      }
    }
    for (auto idxface : this->Faces()) {
      bool found_valid = false;
      for (size_t i = 0; i < GetNumNeighbourCells(idxface); ++i) {
        IdxCell idxcell = GetNeighbourCell(idxface, i);
        if (!idxcell.IsNone()) {
          found_valid = true;
        }
      }
      if (!found_valid) {
        ff_is_excluded_[idxface] = true;
      }
    }
  }
  bool IsExcluded(IdxCell idxcell) const {
    return fc_is_excluded_[idxcell];
  }
  bool IsExcluded(IdxFace idxface) const {
    return ff_is_excluded_[idxface];
  }

 private:
   Vect CalcCenter(IdxCell idxcell) const {
     Vect res = Vect::kZero;
     for (size_t i = 0; i < GetNumNeighbourNodes(idxcell); ++i) {
       res += GetNode(GetNeighbourNode(idxcell, i));
     }
     return res / Scal(GetNumNeighbourNodes(idxcell));
   }
   Vect CalcCenter(IdxFace idxface) const {
     Vect res = Vect::kZero;
     for (size_t i = 0; i < GetNumNeighbourNodes(idxface); ++i) {
       res += GetNode(GetNeighbourNode(idxface, i));
     }
     return res / Scal(GetNumNeighbourNodes(idxface));
   }
   Scal CalcArea(IdxFace idxface) const {
     return GetNode(GetNeighbourNode(idxface, 1)).dist(
         GetNode(GetNeighbourNode(idxface, 0)));
   }
   Vect CalcSurface(IdxFace idxface) const {
     Vect v = GetNode(GetNeighbourNode(idxface, 1)) -
         GetNode(GetNeighbourNode(idxface, 0));
     return Vect(v[0]);
   }
   Scal CalcVolume(IdxCell idxcell) const {
     Scal res = 0.;
     for (size_t i = 0; i < GetNumNeighbourFaces(idxcell); ++i) {
       IdxFace idxface = GetNeighbourFace(idxcell, i);
       res += GetCenter(idxface)[0] * GetOutwardSurface(idxcell, i)[0];
     }
     return res;
   }
};

template <class Scal>
MeshStructured<Scal>::MeshStructured(const BlockNodes& b_nodes,
                                          const FieldNode<Vect>& fn_node)
    : b_nodes_(b_nodes)
    , fn_node_(fn_node)
    , b_cells_(b_nodes_.GetBegin(), b_nodes_.GetDimensions() - MIdx(1))
    , b_faces_(b_cells_.GetDimensions())
    , fc_center_(b_cells_)
    , fc_volume_(b_cells_)
    , ff_center_(b_faces_)
    , ff_area_(b_faces_)
    , ff_surface_(b_faces_)
    , fc_neighbour_face_(b_cells_)
    , fc_neighbour_node_base_(b_cells_)
    , ff_direction_(b_faces_)
    , ff_neighbour_cell_(b_faces_)
    , ff_neighbour_node_(b_faces_)
    , fc_is_inner_(b_cells_)
    , ff_is_inner_(b_faces_)
    , fc_is_excluded_(b_cells_, false)
    , ff_is_excluded_(b_faces_, false)
{
  MIdx mb = b_cells_.GetBegin(), me = b_cells_.GetEnd();

  // Base for cell neighbours
  for (auto midx : b_cells_) {
    IdxCell idxcell(b_cells_.GetIdx(midx));
    fc_neighbour_node_base_[idxcell] = b_nodes_.GetIdx(midx);

    auto nface = [this, midx](IntIdx i, Direction dir) {
      return b_faces_.GetIdx(midx + MIdx(i), dir);
    };
    fc_neighbour_face_[idxcell] = {{
        nface(0, Direction::i),
        nface(1, Direction::i)}};

    fc_is_inner_[idxcell] = (mb < midx && midx + MIdx(1) < me);
  }

  // Offset for cell neighbour cells
  {
    auto offset = [this](IntIdx i) {
      return static_cast<IntIdx>(b_cells_.GetIdx(MIdx(i)).GetRaw() -
          b_cells_.GetIdx(MIdx(0)).GetRaw());
    };

    cell_neighbour_cell_offset_ = {{offset(-1), offset(1)}};
  }

  // Offset for cell neighbour nodes
  {
    auto offset = [this](IntIdx i) {
      return static_cast<IntIdx>(b_nodes_.GetIdx(MIdx(i)).GetRaw() -
          b_nodes_.GetIdx(MIdx(0)).GetRaw());
    };

    cell_neighbour_node_offset_ = {{offset(0), offset(1)}};
  }

  // Base for i-face neighbours
  for (Direction dir : {Direction::i}) {
    for (auto midx : BlockCells(mb, me - mb + dir)) {
      IdxFace idxface = b_faces_.GetIdx(midx, dir);
      ff_direction_[idxface] = dir;
      ff_neighbour_cell_[idxface] = {{
          b_cells_.GetIdx(midx - dir),
          b_cells_.GetIdx(midx)}};
      if (dir == Direction::i) {
        ff_neighbour_node_[idxface] = {{b_nodes_.GetIdx(midx)}};
      }
      if (midx[dir] == mb[dir]) {
        ff_neighbour_cell_[idxface][0] = IdxCell::None();
      }
      if (midx[dir] == me[dir]) {
        ff_neighbour_cell_[idxface][1] = IdxCell::None();
      }
    }
  }

  // Face centers, area
  for (auto idx : this->Faces()) {
    ff_center_[idx] = CalcCenter(idx);
    ff_area_[idx] = CalcArea(idx);
    ff_surface_[idx] = CalcSurface(idx);
  }

  // Cell centers, volume
  for (auto idx : this->Cells()) {
    fc_center_[idx] = CalcCenter(idx);
    fc_volume_[idx] = CalcVolume(idx);
  }

  for (auto idxface : this->Faces()) {
    ff_is_inner_[idxface] = (!GetNeighbourCell(idxface, 0).IsNone() &&
        !GetNeighbourCell(idxface, 1).IsNone());
  }
}

} // namespace geom1d

} // namespace geom
