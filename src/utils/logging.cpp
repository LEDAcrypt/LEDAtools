#include "logging.hpp"
#include <spdlog/sinks/basic_file_sink.h>
#include <string>

void configure_logger(const std::optional<std::string> filename) {
  // Initialize the logger
  const std::string ff =
      filename.has_value() ? filename.value() : "logs/default.log";
  auto logger = spdlog::basic_logger_mt("default_logger", ff);
  spdlog::set_default_logger(logger);

  // Retrieve the environment variable for log level
  const char *log_level_env = std::getenv("LOG_LEVEL");

  if (log_level_env) {
    std::string log_level_str(log_level_env);

    // Configure the log level based on the environment variable
    if (log_level_str == "trace") {
      spdlog::set_level(spdlog::level::trace);
    } else if (log_level_str == "debug") {
      spdlog::set_level(spdlog::level::debug);
    } else if (log_level_str == "info") {
      spdlog::set_level(spdlog::level::info);
    } else if (log_level_str == "warn") {
      spdlog::set_level(spdlog::level::warn);
    } else if (log_level_str == "err") {
      spdlog::set_level(spdlog::level::err);
    } else if (log_level_str == "critical") {
      spdlog::set_level(spdlog::level::critical);
    } else {
      spdlog::set_level(spdlog::level::info); // Default level
    }
  } else {
    spdlog::set_level(spdlog::level::err); // Default level if environment
                                            // variable is not set
  }
}

std::string optional_to_string(const std::optional<uint32_t> &opt) {
  if (opt) {
    return std::to_string(*opt);
  } else {
    return "Not Initialized";
  }
}
