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

template <class Mesh>
using SessionPlain = plain::SessionPlain<Mesh>;

} // namespace output
