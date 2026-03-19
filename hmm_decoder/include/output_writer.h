#ifndef OUTPUT_WRITER_H
#define OUTPUT_WRITER_H

#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include "data_structures.h"

void write_results_to_text(
    const std::pair<std::string, std::vector<HetLocation>>& result,
    std::ofstream& file);

#endif // OUTPUT_WRITER_H
