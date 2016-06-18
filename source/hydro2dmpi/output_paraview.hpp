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
    out_ << "test" << std::endl;
  }
};

template <class Mesh>
class SessionParaviewStructured : public SessionParaview {
  const Mesh& mesh_;
 public:
  SessionParaviewStructured(const Content& content,
                            std::string title,
                            std::string filename,
                            const Mesh& mesh)
      : SessionParaview(content, title, filename),
        mesh_(mesh) {
    out_ << mesh_.GetNumCells() << std::endl;
  }
  void Write(double time, std::string title) override {
    file_ << time << " " << title << std::endl;
  }
};

} // namespace output
