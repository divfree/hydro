/*
 * output.hpp
 *
 *  Created on: Jan 30, 2016
 *      Author: Petr Karnakov (petr.karnakov@gmail.com)
 */

#pragma once

#include <string>
#include <iostream>
#include <memory>
#include "mesh.hpp"
#include <cassert>
#include <fstream>

namespace tecplotio {
#include "TECIO.h"
}

namespace output {

class EntryGeneric {
  std::string name_;
 public:
  EntryGeneric(std::string name)
      : name_(name)
  {}
  virtual std::string GetName() {
    return name_;
  }
  virtual void Prepare() = 0;
  virtual ~EntryGeneric() {}
};

template <class FieldType>
class EntryField : public EntryGeneric {
 public:
  EntryField(std::string name)
      : EntryGeneric(name)
  {}
  virtual const FieldType& GetField() = 0;
};

template <class FieldType>
class EntryFieldCopy : public EntryField<FieldType> {
  const FieldType& field_;
 public:
  EntryFieldCopy(std::string name, const FieldType& field)
      : EntryField<FieldType>(name)
      , field_(field)
  {}
  void Prepare() override {}
  const FieldType& GetField() override {
    return field_;
  }
};

template <class Vect, class IdxType>
class EntryExtractScalar :
    public EntryField<geom::FieldGeneric<typename Vect::value_type, IdxType>> {
  using VectField = geom::FieldGeneric<Vect, IdxType>;
  using ScalarField = geom::FieldGeneric<typename Vect::value_type, IdxType>;
  const VectField& vect_field_;
  ScalarField scalar_field_;
  size_t component_number_;
 public:
  EntryExtractScalar(std::string name, const VectField& vect_field,
                     size_t component_number)
      : EntryField<ScalarField>(name)
      , vect_field_(vect_field)
      , scalar_field_(vect_field_.size())
      , component_number_(component_number)
  {}
  void Prepare() override {
    for (size_t i = 0; i < vect_field_.size(); ++i) {
      IdxType idx(i);
      scalar_field_[idx] = vect_field_[idx][component_number_];
    }
  }
  const ScalarField& GetField() override {
    return scalar_field_;
  }
};

template <class Value, class Idx, class Mesh>
class EntryFunction :
    public EntryField<geom::FieldGeneric<Value, Idx>> {
  using ScalarField = geom::FieldGeneric<Value, Idx>;
  using Function = std::function<Value (Idx)>;
  const Mesh& mesh_;
  ScalarField scalar_field_;
  Function function_;
 public:
  EntryFunction(std::string name, const Mesh& mesh, const Function& function)
      : EntryField<ScalarField>(name)
      , mesh_(mesh)
      , scalar_field_(mesh_)
      , function_(function)
  {}
  void Prepare() override {
    for (auto idx : geom::Range<Idx>(mesh_)) {
      scalar_field_[idx] = function_(idx);
    }
  }
  const ScalarField& GetField() override {
    return scalar_field_;
  }
};

template <class Value>
class EntryScalar : public EntryGeneric {
 public:
  EntryScalar(std::string name)
      : EntryGeneric(name)
  {}
  virtual Value GetValue() = 0;
};

template <class Value>
class EntryScalarFunction :
    public EntryScalar<Value> {
  using Function = std::function<Value ()>;
  Function function_;
 public:
  EntryScalarFunction(std::string name, const Function& function)
      : EntryScalar<Value>(name)
      , function_(function)
  {}
  void Prepare() override {}
  Value GetValue() override {
    return function_();
  }
};

using Content = std::vector<std::shared_ptr<EntryGeneric>>;

class Session {
 public:
  virtual void Write(double time = 0., std::string title = "") = 0;
  virtual ~Session() {}
};

namespace plain {

template <class Scal>
class SessionPlain : public Session {
  std::ostream& out_;
  Content content_;
  void WriteHeader() {
    out_ << "Output session, plain format" << std::endl;
    out_ << "All fields: ";
    for (auto& entry_generic : content_) {
      out_ << entry_generic->GetName() << " ";
    }
    out_ << "\n";
    out_ << "Cell fields: ";
    for (auto& entry_generic : content_)
    if (auto entry = dynamic_cast<EntryField<geom::FieldCell<Scal>>*>(
        entry_generic.get())) {
      out_ << entry->GetName() << " ";
    }
    out_ << "\n";
    out_ << "Face fields: ";
    for (auto& entry_generic : content_)
    if (auto entry = dynamic_cast<EntryField<geom::FieldFace<Scal>>*>(
        entry_generic.get())) {
      out_ << entry->GetName() << " ";
    }
    out_ << "\n";
    out_ << "Node fields: ";
    for (auto& entry_generic : content_)
    if (auto entry = dynamic_cast<EntryField<geom::FieldNode<Scal>>*>(
        entry_generic.get())) {
      out_ << entry->GetName() << " ";
    }
    out_ << "\n";
    out_ << std::endl;
  }
  void WriteFooter() {
    out_ << "Output session finished" << std::endl;
  }
  template <class FieldType>
  bool TryWriteField(EntryGeneric* entry_generic) {
    if (auto entry = dynamic_cast<EntryField<FieldType>*>(
        entry_generic)) {
      auto& field = entry->GetField();
      for (size_t i = 0; i < field.size(); ++i) {
        out_ << field[typename FieldType::IdxType(i)] << " ";
      }
      out_ << "\n";
      return true;
    }
    return false;
  }
 public:
  SessionPlain(const Content& content, std::ostream& out)
      : out_(out)
      , content_(content)
  {
    WriteHeader();
  }
  void Write(double time = 0., std::string title = "") override {
    out_ << "Time = " << time
        << ", Title = " << title << std::endl;
    for (auto& entry_generic : content_) {
      entry_generic->Prepare();
      TryWriteField<geom::FieldCell<Scal>>(entry_generic.get());
      TryWriteField<geom::FieldFace<Scal>>(entry_generic.get());
      TryWriteField<geom::FieldNode<Scal>>(entry_generic.get());
    }
    out_ << std::endl;
  }
  ~SessionPlain() {
    WriteFooter();
  }
};

} // namespace plain

namespace tecplot {

class SessionTecplotAscii : public Session {
 protected:
  Content content_;
  std::string title_;
  std::ostream& out_;
  std::ofstream output_file_;
  void WriteFileHeader() {
    out_ << "Title=\"" << title_ << "\"\n";
    out_ << "Variables=";
    for (auto& entry_generic : content_) {
      out_ << "\"" << entry_generic->GetName() << "\" ";
    }
    out_ << "\n" << std::endl;
    out_.flush();
  }
  SessionTecplotAscii(const Content& content, std::string title,
                      std::ostream& out)
      : content_(content)
      , title_(title)
      , out_(out)
  {
    WriteFileHeader();
  }
  SessionTecplotAscii(const Content& content, std::string title,
                      std::string filename)
      : content_(content)
      , title_(title)
      , out_(output_file_)
  {
    output_file_.open(filename);
    WriteFileHeader();
  }
};

template <class Mesh>
class SessionTecplotAsciiGeneric : public SessionTecplotAscii {
  const Mesh& mesh_;

 public:
  SessionTecplotAsciiGeneric(const Content& content, std::string title,
                             const Mesh& mesh,
                             std::ostream& out)
      : SessionTecplotAscii(content, title, out)
      , mesh_(mesh)
  {}
  void Write(double /*time = 0.*/, std::string /*title = ""*/) override {
    out_ << "tecplot with generic mesh: write" << std::endl;
  }
};

template <class Mesh>
class SessionTecplotAsciiStructured : public SessionTecplotAscii {
  const Mesh& mesh_;
  using Scal = typename Mesh::Scal;
  using MIdx = typename Mesh::MIdx;
  const size_t max_line_count_;
  using FieldCell = geom::FieldCell<Scal>;
  using FieldNode = geom::FieldNode<Scal>;

  void WriteField(EntryField<FieldCell>* entry) {
    size_t line_count = 0;
    auto& b_cells = mesh_.GetBlockCells();

    for (auto midx : b_cells) {
      out_ << entry->GetField()[b_cells.GetIdx(midx)] << " ";
      if (++line_count >= max_line_count_) {
        out_ << "\n";
        line_count = 0;
      }
    }

    if (line_count) {
      out_ << "\n";
    }
    out_ << std::endl;
  }
  void WriteField(EntryField<FieldNode>* entry) {
    size_t line_count = 0;
    auto& b_nodes = mesh_.GetBlockNodes();

    for (auto midx : b_nodes) {
      out_ << entry->GetField()[b_nodes.GetIdx(midx)] << " ";
      if (++line_count >= max_line_count_) {
        out_ << "\n";
        line_count = 0;
      }
    }

    if (line_count) {
      out_ << "\n";
    }
    out_ << std::endl;
  }
  void WriteZoneHeader(double time, std::string title) {
    auto& b_nodes = mesh_.GetBlockNodes();
    out_ << "Zone I=" << b_nodes.GetDimensions()[0]
        << ", J=" << b_nodes.GetDimensions()[1]
        << ", Datapacking=BLOCK";
    out_ << ", VarLocation=([";
    bool need_comma = false;
    for (size_t i = 0; i < content_.size(); ++i) {
      if (dynamic_cast<EntryField<FieldCell>*>(content_[i].get())) {
        if (need_comma) {
          out_ << ",";
        }
        need_comma = true;
        out_ << i + 1;
      }
    }
    out_ << "]=CELLCENTERED)\n";
    out_ << "T=\"" << title << "\"\n";
    out_ << "SolutionTime=" << time << "\n";
  }
 public:
  SessionTecplotAsciiStructured(const Content& content, std::string title,
                                const Mesh& mesh,
                                std::ostream& out)
      : SessionTecplotAscii(content, title, out)
      , mesh_(mesh)
      , max_line_count_(100)
  {}
  void Write(double time = 0., std::string title = "") override {
    WriteZoneHeader(time, title);

    for (auto& entry_generic : content_) {
      entry_generic->Prepare();
      if (auto entry = dynamic_cast<EntryField<FieldCell>*>(
          entry_generic.get())) {
        WriteField(entry);
      } else if (auto entry = dynamic_cast<EntryField<FieldNode>*>(
          entry_generic.get())) {
        WriteField(entry);
      }
    }

    out_ << std::endl;
  }
};

class SessionTecplotBinary : public Session {
 protected:
  Content content_;
  std::string title_;
  std::string filename_;
  void WriteFileHeader() {
    std::string variables;
    for (auto& entry_generic : content_) {
      variables += entry_generic->GetName() + " ";
    }

    using namespace tecplotio;

    INTEGER4 Debug      = 0;
    INTEGER4 OutputIsDouble   = 0;
    INTEGER4 FileType   = 0;
    INTEGER4 fileFormat = 0; // plt

    TECINI142((char*)title_.c_str(),
              (char*)variables.c_str(),
              (char*)filename_.c_str(),
              (char*)".",
              &fileFormat,
              &FileType,
              &Debug,
              &OutputIsDouble);

  }
  SessionTecplotBinary(const Content& content, std::string title,
                      std::string filename)
      : content_(content)
      , title_(title)
      , filename_(filename)
  {
    WriteFileHeader();
  }
  ~SessionTecplotBinary() {
    using namespace tecplotio;
    TECEND142();
  }
};

template <class Mesh>
class SessionTecplotBinaryStructured : public SessionTecplotBinary {
  const Mesh& mesh_;
  using Scal = typename Mesh::Scal;
  using MIdx = typename Mesh::MIdx;
  using FieldCell = geom::FieldCell<Scal>;
  using FieldNode = geom::FieldNode<Scal>;

  void WriteField(EntryField<FieldCell>* entry) {
    using namespace tecplotio;

    auto& field = entry->GetField();
    INTEGER4 size = static_cast<INTEGER4>(field.size());

    INTEGER4 FieldIsDouble;
    if (sizeof(Scal) == 4) {
      FieldIsDouble = 0;
    } else if (sizeof(Scal) == 8) {
      FieldIsDouble = 1;
    } else {
      // sizeof(Scal) should be either 4 or 8
      assert(false);
    }

    std::vector<Scal> field_raw(size);

    auto& b_cells = mesh_.GetBlockCells();
    size_t m = 0;
    for (auto midx : b_cells) {
      field_raw[m++] = field[b_cells.GetIdx(midx)];
    }
    TECDAT142(&size, field_raw.data(), &FieldIsDouble);
  }
  void WriteField(EntryField<FieldNode>* entry) {
    using namespace tecplotio;

    auto& field = entry->GetField();
    INTEGER4 size = static_cast<INTEGER4>(field.size());

    INTEGER4 FieldIsDouble;
    if (sizeof(Scal) == 4) {
      FieldIsDouble = 0;
    } else if (sizeof(Scal) == 8) {
      FieldIsDouble = 1;
    } else {
      // sizeof(Scal) should be either 4 or 8
      assert(false);
    }

    std::vector<Scal> field_raw(size);

    auto& b_nodes = mesh_.GetBlockNodes();
    size_t m = 0;
    for (auto midx : b_nodes) {
      field_raw[m++] = field[b_nodes.GetIdx(midx)];
    }
    TECDAT142(&size, field_raw.data(), &FieldIsDouble);
  }
  void WriteZoneHeader(double time, std::string zone_title) {
    auto& b_nodes = mesh_.GetBlockNodes();

    using namespace tecplotio;
    std::vector<INTEGER4> value_location(content_.size());
    for (size_t i = 0; i < content_.size(); ++i) {
      value_location[i] =
          (dynamic_cast<EntryField<FieldCell>*>(content_[i].get()) ? 0 : 1);
    }

    INTEGER4 ICellMax                 = 0;
    INTEGER4 JCellMax                 = 0;
    INTEGER4 KCellMax                 = 0;
    INTEGER4 StrandID                 = 1;
    INTEGER4 ParentZn                 = 0;
    INTEGER4 IsBlock                  = 1;      /* Block */
    INTEGER4 NFConns                  = 0;
    INTEGER4 FNMode                   = 0;
    INTEGER4 TotalNumFaceNodes        = 1;
    INTEGER4 TotalNumBndryFaces       = 1;
    INTEGER4 TotalNumBndryConnections = 1;
    INTEGER4 ShrConn                  = 0;

    /*Ordered Zone Parameters*/
    INTEGER4 IMax = static_cast<INTEGER4>(b_nodes.GetDimensions()[0]);
    INTEGER4 JMax = static_cast<INTEGER4>(b_nodes.GetDimensions()[1]);
    INTEGER4 KMax = static_cast<INTEGER4>((Mesh::dim == 2) ? 1 : b_nodes.GetDimensions()[2]);

    /*  Ordered Zone */
    INTEGER4 ZoneType = 0;
    TECZNE142((char*)zone_title.c_str(),
              &ZoneType,
              &IMax,
              &JMax,
              &KMax,
              &ICellMax,
              &JCellMax,
              &KCellMax,
              &time,
              &StrandID,
              &ParentZn,
              &IsBlock,
              &NFConns,
              &FNMode,
              &TotalNumFaceNodes,
              &TotalNumBndryFaces,
              &TotalNumBndryConnections,
              NULL,
              value_location.data(),
              NULL,
              &ShrConn);
  }
 public:
  SessionTecplotBinaryStructured(const Content& content, std::string title,
                                const Mesh& mesh,
                                std::string filename)
      : SessionTecplotBinary(content, title, filename)
      , mesh_(mesh)
  {}
  void Write(double time = 0., std::string title = "") override {
    WriteZoneHeader(time, title);

    for (auto& entry_generic : content_) {
      entry_generic->Prepare();
      if (auto entry = dynamic_cast<EntryField<FieldCell>*>(
          entry_generic.get())) {
        WriteField(entry);
      } else if (auto entry = dynamic_cast<EntryField<FieldNode>*>(
          entry_generic.get())) {
        WriteField(entry);
      }
    }
  }
};


template <class Scal>
class SessionTecplotAsciiScalar : public SessionTecplotAscii {
  void WriteZoneHeader(std::string title) {
    out_ << "Zone ZoneType=ORDERED DataPacking=POINT\n";
    out_ << "T=\"" << title << "\"\n";
  }
 public:
  SessionTecplotAsciiScalar(const Content& content, std::string title,
                                std::string zonetitle,
                                std::ostream& out)
      : SessionTecplotAscii(content, title, out)
  {
    WriteZoneHeader(zonetitle);
  }
  SessionTecplotAsciiScalar(const Content& content, std::string title,
                                std::string zonetitle,
                                std::string filename)
      : SessionTecplotAscii(content, title, filename)
  {
    WriteZoneHeader(zonetitle);
  }
  void Write(double /*time = 0.*/, std::string /*title = ""*/) override {
    for (auto& entry_generic : content_) {
      entry_generic->Prepare();
      if (auto entry = dynamic_cast<EntryScalar<Scal>*>(
          entry_generic.get())) {
        out_ << entry->GetValue() << " ";
      } else {
        throw std::runtime_error(
            "SessionTecplotAsciiScalar: Unknown entry type");
      }
    }

    out_ << std::endl;
  }
};


} // namespace tecplot

template <class Mesh>
using SessionPlain = plain::SessionPlain<Mesh>;

using SessionTecplotAscii = tecplot::SessionTecplotAscii;

template <class Mesh>
using SessionTecplotAsciiGeneric = tecplot::SessionTecplotAsciiGeneric<Mesh>;

template <class Mesh>
using SessionTecplotAsciiStructured =
    tecplot::SessionTecplotAsciiStructured<Mesh>;

using SessionTecplotBinary = tecplot::SessionTecplotBinary;

template <class Mesh>
using SessionTecplotBinaryStructured =
    tecplot::SessionTecplotBinaryStructured<Mesh>;

template <class Scal>
using SessionTecplotAsciiScalar =
    tecplot::SessionTecplotAsciiScalar<Scal>;

} // namespace output
