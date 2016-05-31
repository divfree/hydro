/*
 * test.cpp
 *
 *  Created on: Jan 30, 2016
 *      Author: Petr Karnakov (petr.karnakov@gmail.com)
 */

#include "test.hpp"
#include <iostream>
#include <string>
#include <memory>
#include <fstream>
#include "time.h"

/*
void TestMesh3d() {
  using Scal = double;
  using Mesh = geom::geom3d::MeshStructured<Scal>;
  using Vect = typename Mesh::Vect;
  using MIdx = geom::geom3d::MIdx;
  const MIdx mesh_size(1);
  geom::geom3d::BlockNodes b_nodes(mesh_size + MIdx(1));
  geom::FieldNode<Vect> fn_node(b_nodes);
  for (auto midx : b_nodes) {
    geom::IdxNode idxnode = b_nodes.GetIdx(midx);
    Vect unit = Vect(midx) / Vect(mesh_size);
    fn_node[idxnode] = unit;
  }
  Mesh mesh = Mesh(b_nodes, fn_node);

  std::cout << "Mesh3d" << std::endl;
  std::cout << "Nodes" << std::endl;
  for (auto idxnode : mesh.Nodes()) {
    std::cout << idxnode.GetRaw()
        << " " << mesh.GetNode(idxnode)
        << std::endl;
  }
  std::cout << "Cells" << std::endl;
  for (auto idxcell : mesh.Cells()) {
    std::cout << idxcell.GetRaw()
        << " " << mesh.GetCenter(idxcell)
        << " " << mesh.GetVolume(idxcell)
        << std::endl;
    std::cout << "neighbour faces:" << std::endl;
    for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
      auto idxface = mesh.GetNeighbourFace(idxcell, i);
      std::cout << "\t" << i << " " << idxface.GetRaw() << std::endl;
    }
  }
  auto& b_faces = mesh.GetBlockFaces();
  std::cout << "Faces" << std::endl;
  for (auto idxface : mesh.Faces()) {
    std::cout << idxface.GetRaw()
        << " " << GetLetter(mesh.GetDirection(idxface))
        << " " << mesh.GetCenter(idxface)
        << " " << mesh.GetSurface(idxface)
        << " " << b_faces.GetMIdx(idxface)
        << " " << b_faces.GetIdx(b_faces.GetMIdx(idxface), mesh.GetDirection(idxface)).GetRaw()
        << std::endl;
    std::cout << "neighbour nodes:" << std::endl;
    for (size_t i = 0; i < mesh.GetNumNeighbourNodes(idxface); ++i) {
      std::cout << "\t" << i << " " << mesh.GetNeighbourNode(idxface, i).GetRaw() << std::endl;
    }
    std::cout << "neighbour cells:" << std::endl;
    for (size_t i = 0; i < mesh.GetNumNeighbourCells(idxface); ++i) {
      auto idxcell = mesh.GetNeighbourCell(idxface, i);
      std::cout << "\t" << i << " " << idxcell.GetRaw() << std::endl;
    }
  }

  auto& b_cells = mesh.GetBlockCells();
  std::cout << "Faces block" << std::endl;
  MIdx mb = b_cells.GetBegin(), me = b_cells.GetEnd();
  using Direction = geom::geom3d::Direction;
  for (Direction dir : {Direction::i, Direction::j, Direction::k}) {
    for (auto midx : geom::geom3d::BlockGeneric<size_t>(mb,
                                          me - mb + GetMIdx(dir))) {
      IdxFace idxface = b_faces.GetIdx(midx, dir);
      std::cout << midx << " " << GetLetter(dir)
          << " " << idxface.GetRaw()
          << " " << GetLetter(b_faces.GetDirection(idxface))
          << std::endl;
    }
  }
}
*/

/*
void Stat(const Mesh& mesh, std::ostream& out) {
  out << "Mesh statistics\n\n";

  out << "Nodes, number = " << mesh.GetNumNodes() <<":\n";

  for (auto idx : mesh.Nodes()) {
    out << idx.GetRaw() << ": " << mesh.GetNode(idx) << "\n";
  }
  std::cout << std::endl;

  out << "Cells, number = " << mesh.GetNumCells() <<":\n";

  for (auto idx : mesh.Cells()) {
    out << idx.GetRaw() << ": c=" << mesh.GetCenter(idx)
        << ", vol=" << mesh.GetVolume(idx)
        << ", c0=" << mesh.GetNeighbourCell(idx, 0).GetRaw()
        << ", c1=" << mesh.GetNeighbourCell(idx, 1).GetRaw()
        << ", c2=" << mesh.GetNeighbourCell(idx, 2).GetRaw()
        << ", c3=" << mesh.GetNeighbourCell(idx, 3).GetRaw()
        << ", f0=" << mesh.GetNeighbourFace(idx, 0).GetRaw()
        << ", f1=" << mesh.GetNeighbourFace(idx, 1).GetRaw()
        << ", f2=" << mesh.GetNeighbourFace(idx, 2).GetRaw()
        << ", f3=" << mesh.GetNeighbourFace(idx, 3).GetRaw()
        << "\n";
  }
  std::cout << std::endl;

  out << "Faces, number = " << mesh.GetNumFaces() <<":\n";

  for (auto idx : mesh.Faces()) {
    out << idx.GetRaw() << ": c=" << mesh.GetCenter(idx) << ", area="
        << mesh.GetArea(idx) << ", surf=" << mesh.GetSurface(idx)
        << ", n0=" << mesh.GetNeighbourNode(idx, 0).GetRaw()
        << ", n1=" << mesh.GetNeighbourNode(idx, 1).GetRaw()
        << "\n";
  }
  std::cout << std::endl;
}
*/
void Test() {
  long int t_start = clock();

  long int t_stop = clock();
  std::cout << "Execution time (ms): "
      << ((t_stop - t_start) * 1000) / CLOCKS_PER_SEC;


}




