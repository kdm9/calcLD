#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <sstream>

// Calculate mean of a range of doubles
double calculate_mean(const double* values, size_t size) {
    if (size == 0) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < size; ++i) {
        sum += values[i];
    }
    return sum / size;
}

// Determine the maximum number of columns needed for output
size_t determine_max_columns(const std::vector<std::vector<double>>& values, int n) {
    if (values.empty() || values[0].empty()) {
        return 0;
    }
    
    // Find the maximum length among all value vectors
    size_t max_length = 0;
    for (const auto& row : values) {
        max_length = std::max(max_length, row.size());
    }
    
    // Calculate how many bins we'll need
    return (max_length + n - 1) / n; // Ceiling division
}

// Process a bin of data
void process_bin(const std::string& chrom, 
                 const std::vector<double>& positions, 
                 const std::vector<std::vector<double>>& values, 
                 int n = 100,
                 bool print_header = false,
                 size_t max_columns = 0) {
    
    // If printing header
    if (print_header) {
        std::cout << "Chromosome\tPosition";
        for (size_t i = 0; i < max_columns; ++i) {
            std::cout << "\tBin" << (i+1);
        }
        std::cout << std::endl;
        return;
    }
    
    // Calculate mean position
    double mean_pos = 0.0;
    if (!positions.empty()) {
        mean_pos = std::accumulate(positions.begin(), positions.end(), 0.0) / positions.size();
    }
    
    // Output chromosome and mean position
    std::cout << chrom << "\t" << mean_pos;
    
    // If no values, output empty columns to maintain consistent column count
    if (values.empty() || values[0].empty()) {
        for (size_t i = 0; i < max_columns; ++i) {
            std::cout << "\t";
        }
        std::cout << std::endl;
        return;
    }
    
    size_t value_length = values[0].size();
    size_t num_rows = values.size();
    
    // Pre-allocate a buffer for bin means to avoid repeated allocations
    std::vector<double> bin_means(num_rows);
    
    // Count of bins processed
    size_t bins_processed = 0;
    
    for (size_t i = 0; i < value_length; i += n) {
        size_t end = std::min(i + n, value_length);
        size_t valid_means = 0;
        
        // Calculate mean for each row's segment
        for (size_t j = 0; j < num_rows; ++j) {
            const std::vector<double>& row = values[j];
            if (i < row.size()) {  // Ensure we don't access beyond the row's bounds
                double row_mean = 0.0;
                size_t count = 0;
                
                for (size_t k = i; k < end && k < row.size(); ++k) {
                    row_mean += row[k];
                    count++;
                }
                
                if (count > 0) {
                    bin_means[valid_means++] = row_mean / count;
                }
            }
        }
        
        // Calculate mean of means
        double mean_val = 0.0;
        if (valid_means > 0) {
            mean_val = calculate_mean(bin_means.data(), valid_means);
        }
        std::cout << "\t" << mean_val;
        bins_processed++;
    }
    
    // Add empty columns if needed to maintain consistent column count
    for (size_t i = bins_processed; i < max_columns; ++i) {
        std::cout << "\t";
    }
    
    std::cout << std::endl;
}

// Bin N function
void bin_n(const std::string& flatld_file, int n = 10) {
    std::ifstream file(flatld_file);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << flatld_file << std::endl;
        return;
    }
    
    std::string line;
    std::string current_chrom;
    std::vector<double> positions;
    std::vector<std::vector<double>> bin_values;
    size_t nlines = 0;
    
    // Pre-allocate memory to reduce reallocations
    positions.reserve(n);
    bin_values.reserve(n);
    
    // Determine max columns from header
    size_t max_columns = 0;
    
    // Read header
    if (std::getline(file, line)) {
        // Count columns in header (minus the first two which are chrom and pos)
        std::istringstream header_iss(line);
        std::string dummy;
        
        // Skip first two columns (chrom and pos)
        header_iss >> dummy >> dummy;
        
        // Count remaining columns
        size_t header_columns = 0;
        while (header_iss >> dummy) {
            header_columns++;
        }
        
        // Calculate max columns based on bin size
        max_columns = (header_columns + n - 1) / n; // Ceiling division
    }
    
    // Reset file position to beginning
    //file.clear();
    //file.seekg(0);
    //
    //// Skip header
    //if (std::getline(file, line)) {
    //    // Skip header processing
    //}
    
    // Print header
    process_bin("", {}, {}, n, true, max_columns);
    
    // Reset for data processing
    bin_values.clear();
    positions.clear();
    current_chrom = "";
    nlines = 0;
    
    // Second pass to process data with consistent column count
    while (std::getline(file, line)) {
        if (++nlines % 1000 == 0) {
            std::cerr << "\rProcessed " << nlines << " lines";
            std::cerr.flush();
        }
        
        std::istringstream iss(line);
        std::string chrom;
        double pos;
        
        // Parse chromosome and position
        if (!(iss >> chrom >> pos)) {
            continue; // Skip malformed lines
        }
        
        // Process bin if chromosome changes or we've reached n positions
        if (chrom != current_chrom || positions.size() == static_cast<size_t>(n)) {
            if (!current_chrom.empty()) {
                process_bin(current_chrom, positions, bin_values, n, false, max_columns);
            }
            positions.clear();
            bin_values.clear();
            current_chrom = chrom;
        }
        
        positions.push_back(pos);
        
        // Parse remaining values
        std::vector<double> values;
        values.reserve(100); // Reserve space for a reasonable number of values
        double val;
        while (iss >> val) {
            values.push_back(val);
        }
        
        bin_values.push_back(std::move(values)); // Use move semantics to avoid copying
    }
    
    // Process the last bin
    if (!current_chrom.empty() && !positions.empty()) {
        process_bin(current_chrom, positions, bin_values, n, false, max_columns);
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <flatld_file> [bin_size]" << std::endl;
        return 1;
    }
    
    int bin_size = 10;
    if (argc > 2) {
        bin_size = std::stoi(argv[2]);
    }
    
    bin_n(argv[1], bin_size);
    return 0;
}
