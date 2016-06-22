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
 public:
  class LogStream {
    Logger* parent_;
    const std::string prefix_;
    std::stringstream buf_;
   public:
    LogStream() = delete;
    LogStream(LogStream&&) = default;
    LogStream(Logger* parent)
        : parent_(parent) {}
    LogStream(Logger* parent, std::string prefix)
        : parent_(parent), prefix_(prefix) {}
    template <class T>
    LogStream& operator<<(const T& value) {
      buf_ << value;
      return *this;
    }
    ~LogStream() {
      (*parent_->p_stream_) << prefix_ << buf_.rdbuf() << std::endl;
    }
  };

  Logger()
      : p_stream_(&std::cout, [](std::ostream*){}) {}
  Logger(std::shared_ptr<std::ofstream> p_stream)
      : p_stream_(p_stream) {}
  Logger(std::string prefix)
      : default_prefix_(prefix),
        p_stream_(&std::cout, [](std::ostream*){}) {}
  Logger(std::shared_ptr<std::ofstream> p_stream, std::string prefix)
      : default_prefix_(prefix),
        p_stream_(p_stream) {}
  LogStream operator()() {
    return LogStream(this, default_prefix_);
  }
  LogStream operator()(const std::string prefix) {
    return LogStream(this, prefix);
  }
};

} // namespace logger
