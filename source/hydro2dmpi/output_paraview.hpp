#pragma once

#include <string>
#include <iostream>
#include <memory>
#include "mesh.hpp"
#include <cassert>
#include <fstream>
#include <fstream>
#include "output.hpp"

namespace output {

class SessionParaview : public Session {
 protected:
  Content content_;
  std::string title_;
  std::string filename_;
  std::ofstream out_;
 public:
  SessionParaview(const Content& content,
                  std::string title,
                  std::string filename)
      : content_(content),
        title_(title),
        filename_(filename),
        out_(filename_) {
    out_.sync_with_stdio(false);
  }
};

template <class Mesh>
class SessionParaviewStructured : public SessionParaview {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using MIdx = typename Mesh::MIdx;
  using FieldCell = geom::FieldCell<Scal>;
  using FieldNode = geom::FieldNode<Scal>;

  const Mesh& mesh;
  void WriteDataArrayHeader(std::string name, size_t num_components) {
    out_ << "        <DataArray "
        << "Name=\"" << name << "\" "
        << "NumberOfComponents=\"" << num_components << "\" "
        << "type=\"Float32\" format=\"ascii\">\n";
  }
  void WriteDataArrayFooter() {
    out_ << "        </DataArray>\n";
  }
  void WriteField(EntryField<FieldCell>* entry) {
    auto& field = entry->GetField();

    WriteDataArrayHeader(entry->GetName(), 1);
    for (auto idxcell : mesh.Cells()) {
      out_ << field[idxcell] << " ";
    }
    out_ << "\n";
    WriteDataArrayFooter();
  }
  void WriteField(EntryField<FieldNode>* entry) {
    auto& field = entry->GetField();

    WriteDataArrayHeader(entry->GetName(), 1);
    for (auto idxnode : mesh.Nodes()) {
      out_ << field[idxnode] << " ";
    }
    out_ << "\n";
    WriteDataArrayFooter();
  }
 public:
  SessionParaviewStructured(const Content& content,
                            std::string title,
                            std::string filename,
                            const Mesh& mesh)
      : SessionParaview(content, title, filename),
        mesh(mesh) {
    out_ << "<VTKFile type=\"StructuredGrid\" "
        << "version=\"0.1\" byte_order=\"LittleEndian\">\n";

    MIdx size = mesh.GetBlockCells().GetDimensions();
    out_ << "  <StructuredGrid timestep=\"5\" WholeExtent=\""
        << "0 " << size[0] << " 0 " << size[1] << " 0 " << size[2] << "\">\n";

    out_ << "    <Piece Extent=\""
        << "0 " << size[0] << " 0 " << size[1] << " 0 " << size[2] << "\">\n";
  }
  ~SessionParaviewStructured() {
    out_ << "    </Piece>\n";
    out_ << "  </StructuredGrid>\n";
    out_ << "</VTKFile>\n";
  }
  void Write(double time, std::string title) override {
    out_ << "      <PointData>\n";
    for (auto& entry_generic : content_) {
      if (auto entry = dynamic_cast<EntryField<FieldNode>*>(
          entry_generic.get())) {
        entry->Prepare();
        WriteField(entry);
      }
    }
    out_ << "      </PointData>\n";

    out_ << "      <CellData>\n";
    for (auto& entry_generic : content_) {
      if (auto entry = dynamic_cast<EntryField<FieldCell>*>(
          entry_generic.get())) {
        entry->Prepare();
        WriteField(entry);
      }
    }
    out_ << "      </CellData>\n";

    out_ << "      <Points>\n";
    WriteDataArrayHeader("mesh", Mesh::dim);
    for (auto idxnode : mesh.Nodes()) {
      Vect p = mesh.GetNode(idxnode);
      for (size_t i = 0; i < Mesh::dim; ++i) {
        out_ << p[i] << " ";
      }
    }
    WriteDataArrayFooter();
    out_ << "      </Points>\n";

  }
};

} // namespace output
