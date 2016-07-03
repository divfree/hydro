#pragma once

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>

namespace logger {

class Logger {
 protected:
  const std::string default_prefix_;
  std::shared_ptr<std::ostream> p_stream_;
  const bool no_eol_ = false;
 public:
  class LogStream {
    Logger* parent_;
    const std::string prefix_;
    std::shared_ptr<std::stringstream> buf_;
    const bool no_eol_;
   public:
    LogStream() = delete;
    LogStream(LogStream&&) = default;
    LogStream(Logger* parent, bool no_eol)
        : parent_(parent), buf_(std::make_shared<std::stringstream>()),
          no_eol_(no_eol) {}
    LogStream(Logger* parent, std::string prefix, bool no_eol)
        : parent_(parent), prefix_(prefix),
          buf_(std::make_shared<std::stringstream>()),
          no_eol_(no_eol) {}
    template <class T>
    LogStream& operator<<(const T& value) {
      (*buf_) << value;
      return *this;
    }
    ~LogStream() {
      (*parent_->p_stream_) << prefix_ << (*buf_).str();
      if (!no_eol_) {
        (*parent_->p_stream_) << std::endl;
      } else {
        (*parent_->p_stream_).flush();
      }
    }
  };

  Logger(bool no_eol = false)
      : p_stream_(&std::cout, [](std::ostream*){}),
        no_eol_(no_eol) {}
  Logger(std::shared_ptr<std::ostream> p_stream, bool no_eol = false)
      : p_stream_(p_stream),
        no_eol_(no_eol) {}
  Logger(std::string prefix, bool no_eol = false)
      : default_prefix_(prefix),
        p_stream_(&std::cout, [](std::ostream*){}),
        no_eol_(no_eol) {}
  Logger(std::shared_ptr<std::ostream> p_stream, std::string prefix,
         bool no_eol = false)
      : default_prefix_(prefix),
        p_stream_(p_stream),
        no_eol_(no_eol) {}
  LogStream operator()() {
    return LogStream(this, default_prefix_, no_eol_);
  }
  LogStream operator()(const std::string prefix) {
    return LogStream(this, prefix, no_eol_);
  }
  void SetStream(std::shared_ptr<std::ostream> p_stream) {
    p_stream_ = p_stream;
  }
};

} // namespace logger
