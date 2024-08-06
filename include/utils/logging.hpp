#pragma once

#include <map>
#include <memory>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <string>
#include <vector>
#include <optional>

namespace Logger {
  class LoggerManager {
  public:
    static LoggerManager& getInstance();

    void setup_logger(const std::string& logger_name,
                      spdlog::level::level_enum console_level,
                      spdlog::level::level_enum file_level,
                      const std::string& pattern = "[%Y-%m-%d %H:%M:%S.%e] [%n] [%l] [%s:%#] %v");
    std::shared_ptr<spdlog::logger> get_logger(const std::string& logger_name);

    template <typename T>
    std::string optional_to_string(const std::optional<T>& opt);
    
    template <typename T>
    std::string array_to_string(const std::vector<T>& vec);
    
    template <typename T>
    std::string array_to_string(const T* array, size_t size);

  private:
    LoggerManager() = default;
    std::map<std::string, std::shared_ptr<spdlog::logger>> loggers;
  };

  template <typename T>
  inline std::string LoggerManager::optional_to_string(const std::optional<T>& opt) {
    if (opt) {
      return std::to_string(*opt); // Convert the value to string if present
    } else {
      return "None"; // Represent the absence of value
    }
  }

  template <>
  inline std::string LoggerManager::optional_to_string(const std::optional<std::string>& opt) {
    if (opt) {
      return *opt; // Return the string directly if present
    } else {
      return "None"; // Represent the absence of value
    }
  }

  template <typename T>
  inline std::string LoggerManager::array_to_string(const std::vector<T>& vec) {
    std::string result = "[";
    for (size_t i = 0; i < vec.size(); ++i) {
      result += std::to_string(vec[i]);
      if (i < vec.size() - 1) {
        result += ", ";
      }
    }
    result += "]";
    return result;
  }

  template <typename T>
  inline std::string LoggerManager::array_to_string(const T* array, size_t size) {
    std::string result = "[";
    for (size_t i = 0; i < size; ++i) {
      result += std::to_string(array[i]);
      if (i < size - 1) {
        result += ", ";
      }
    }
    result += "]";
    return result;
  }
}

