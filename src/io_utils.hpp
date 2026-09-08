#pragma once

#include <fstream>
#include <stdexcept>
#include <string>

inline void closeChecked(std::ofstream &stream, const std::string &error) {
  stream.close();
  if (!stream)
    throw std::runtime_error(error);
}
