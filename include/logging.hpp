#pragma once
#include <optional>
#include <spdlog/spdlog.h>
#include <string>
#include <vector>

void configure_logger(const std::optional<std::string> filename);

std::string optional_to_string(const std::optional<uint32_t> &opt);

template <typename T> std::string array_to_string(const std::vector<T> &vec) {
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

template <typename T> std::string array_to_string(const T *array, size_t size) {
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
