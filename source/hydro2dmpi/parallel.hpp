#pragma once

#include "mesh.hpp"
#include "solver.hpp"
#include <mpi.h>
#include <exception>

namespace solver {

template <class Mesh>
class ParallelTools {
 public:
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using IdxNode = geom::IdxNode;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  static constexpr auto MPI_Scal =
      (sizeof(Scal) == sizeof(double) ? MPI_DOUBLE : MPI_FLOAT);
  struct Processor {
    int rank;
    Mesh mesh;
    geom::FieldCell<IdxCell> fc_origin;
    geom::FieldCell<bool> fc_is_active;
    geom::FieldFace<IdxFace> ff_origin;
    geom::FieldFace<bool> ff_is_active;
  };

 private:
  const Mesh& global_mesh_;
  int num_processors_;
  int current_rank_;
  std::vector<Processor> v_processor_;
  void InitLocalMesh(const Mesh& global_mesh_, int rank,
                       int overlap_width = 2) {
    using MIdx = typename Mesh::MIdx;
    using Rect = geom::Rect<MIdx>;
    using BlockNodes = typename Mesh::BlockNodes;
    using BlockFaces = typename Mesh::BlockFaces;
    using BlockCells = typename Mesh::BlockCells;
    int num_pieces = num_processors_ - 1;
    int piece_idx = rank - 1;  // local mesh idx (mesh piece number)

    auto& proc = v_processor_[rank];

    MIdx domain_size = global_mesh_.GetBlockCells().GetDimensions();

    // Rect containing MIdx of active cells
    Rect cells_active(MIdx(0), domain_size);
    cells_active.lb[0] = piece_idx * domain_size[0] / num_pieces;
    cells_active.rt[0] = (piece_idx + 1) * domain_size[0] / num_pieces;

    std::cout << "\npart = " << piece_idx << std::endl;

    std::cout << "begin = " << cells_active.lb
        << " end = " << cells_active.rt << std::endl;

    // Rect containing MIdx of all cells
    Rect cells_ext = cells_active;
    cells_ext.lb[0] = std::max<IntIdx>(0,
                                       cells_active.lb[0] - overlap_width);
    cells_ext.rt[0] = std::min<IntIdx>(domain_size[0],
                                       cells_active.rt[0] + overlap_width);

    // Initialize local nodes
    BlockNodes bnodes(cells_ext.GetDimensions() + MIdx(1));
    geom::FieldNode<Vect> fn_node(bnodes);
    for (auto midx : bnodes) {
      IdxNode idxnode = bnodes.GetIdx(midx);
      MIdx global_midx = cells_ext.lb + midx;
      IdxNode global_idxnode = global_mesh_.GetBlockNodes().GetIdx(global_midx);
      fn_node[idxnode] = global_mesh_.GetNode(global_idxnode);
    }

    // Create mesh
    proc.mesh = Mesh(bnodes, fn_node);

    // Fill references to global cells (origin)
    BlockCells bcells = proc.mesh.GetBlockCells();
    BlockCells global_bcells = global_mesh_.GetBlockCells();
    proc.fc_origin.Reinit(proc.mesh);
    proc.fc_is_active.Reinit(proc.mesh);
    for (IdxCell idxcell : proc.mesh.Cells()) {
      MIdx midx = bcells.GetMIdx(idxcell);
      MIdx global_midx = cells_ext.lb + midx;
      proc.fc_origin[idxcell] = global_bcells.GetIdx(global_midx);
      proc.fc_is_active[idxcell] =
          (cells_active.lb <= global_midx && global_midx < cells_active.rt);
    }

    // Fill references to global faces
    BlockFaces bfaces = proc.mesh.GetBlockFaces();
    BlockFaces global_bfaces = global_mesh_.GetBlockFaces();
    proc.ff_origin.Reinit(proc.mesh);
    for (IdxFace idxface : proc.mesh.Faces()) {
      MIdx midx = bfaces.GetMIdx(idxface);
      auto dir = bfaces.GetDirection(idxface);
      MIdx global_midx = cells_ext.lb + midx;
      proc.ff_origin[idxface] = global_bfaces.GetIdx(global_midx, dir);
    }

    // Traverse all active cells and mark their neighbour faces active
    proc.ff_is_active.Reinit(proc.mesh, false);
    for (IdxCell idxcell : proc.mesh.Cells()) {
      if (proc.fc_is_active[idxcell]) {
        for (size_t k = 0;
             k < proc.mesh.GetNumNeighbourFaces(idxcell); ++k) {
          proc.ff_is_active[proc.mesh.GetNeighbourFace(idxcell, k)] = true;
        }
      }
    }
  }

 public:
  ParallelTools(const Mesh& mesh)
      : global_mesh_(mesh) {
    MPI_Comm_size(MPI_COMM_WORLD, &num_processors_);
    MPI_Comm_rank(MPI_COMM_WORLD, &current_rank_);

    v_processor_.resize(num_processors_);
    for (int rank = 1; rank < num_processors_; ++rank) {
      InitLocalMesh(mesh, rank);
    }
  }
  int GetNumProcessors() const {
    return num_processors_;
  }
  int GetRank() const {
    return current_rank_;
  }
  const Mesh& GetLocalMesh() const {
    return v_processor_[current_rank_].mesh;
  }
  const Mesh& GetLocalMesh(int rank) const {
    return v_processor_[rank].mesh;
  }
  const Processor& GetProcessor(int rank) const {
    return v_processor_[rank];
  }
  void SendLocal(const geom::FieldCell<Scal>& fc_local, int dest) {
    if (dest == current_rank_) {
      throw std::runtime_error("Can't send a message to myself");
    }
    int tag = 0;
    MPI_Send(fc_local.data(), fc_local.size(),
             MPI_Scal, dest, tag, MPI_COMM_WORLD);
  }
  void RecvGlobal(geom::FieldCell<Scal>& fc_global) {
    std::vector<geom::FieldCell<Scal>> v_fc_buf(num_processors_);
    for (int rank = 0; rank < num_processors_; ++rank) {
      if (rank != current_rank_) {
        auto& proc = v_processor_[rank];
        auto& fc_buf = v_fc_buf[rank];
        v_fc_buf[rank].Reinit(proc.mesh);
        int tag = 0;
        MPI_Recv(fc_buf.data(), fc_buf.size(),
                 MPI_Scal, rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (auto idxcell : proc.mesh.Cells()) {
          if (proc.fc_is_active[idxcell]) {
            fc_global[proc.fc_origin[idxcell]] = fc_buf[idxcell];
          }
        }
      }
    }
  }
  void SendGlobal(const geom::FieldCell<Scal>& fc_global) {
    std::vector<geom::FieldCell<Scal>> v_fc_buf(num_processors_);
    for (int rank = 0; rank < num_processors_; ++rank) {
      if (rank != current_rank_) {
        auto& proc = v_processor_[rank];
        auto& fc_buf = v_fc_buf[rank];
        v_fc_buf[rank].Reinit(proc.mesh);

        for (auto idxcell : proc.mesh.Cells()) {
          fc_buf[idxcell] = fc_global[proc.fc_origin[idxcell]];
        }

        int tag = 0;
        MPI_Send(fc_buf.data(), fc_buf.size(),
                 MPI_Scal, rank, tag, MPI_COMM_WORLD);
      }
    }
  }
  void RecvLocal(geom::FieldCell<Scal>& fc_local, int source) {
    if (source == current_rank_) {
      throw std::runtime_error("Can't receive a message from myself");
    }
    int tag = 0;
    MPI_Recv(fc_local.data(), fc_local.size(),
             MPI_Scal, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
};

} // namespace solver
