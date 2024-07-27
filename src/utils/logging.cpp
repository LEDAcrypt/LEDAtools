#include "logging.hpp"
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <string>
#include <vector>
#include <iostream>

Logger::LoggerManager &Logger::LoggerManager::getInstance() {
  static LoggerManager instance;
  return instance;
}

void Logger::LoggerManager::setup_logger(
    const std::string &logger_name, spdlog::level::level_enum console_level,
    spdlog::level::level_enum file_level, const std::string &pattern) {
  // Creating sinks, console
  auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  console_sink->set_level(console_level);
  console_sink->set_pattern(pattern);

  // ... and file
  // TODO change hard-coded dir
  auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
      "logs/" + logger_name + ".log", true);
  file_sink->set_level(file_level);
  file_sink->set_pattern(pattern);

  std::vector<spdlog::sink_ptr> sinks{console_sink, file_sink};

  std::shared_ptr<spdlog::logger> logger;
  auto it = loggers.find(logger_name);
  if (it != loggers.end()) {
    std::cout << "Logger already present " << logger_name << std::endl;
    logger = it->second;
  } else {
    std::cout << "Creating logger " << logger_name << std::endl;
    logger = std::make_shared<spdlog::logger>(logger_name, sinks.begin(),
                                              sinks.end());
    loggers[logger_name] = logger;
    spdlog::register_logger(logger);
  }
  // logger->set_level(spdlog::level::trace); // Set to the most verbose level
  // logger->flush_on(spdlog::level::err);
}

std::shared_ptr<spdlog::logger>
Logger::LoggerManager::get_logger(const std::string &logger_name) {
  auto it = loggers.find(logger_name);
  if (it != loggers.end()) {
    return it->second;
  }
  // Logger not found, so set it up
  setup_logger(logger_name, spdlog::level::info,
               spdlog::level::info); // Default levels

  return loggers[logger_name];
}
