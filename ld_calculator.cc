#include "ld_calculator.hh"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <thread>
#include <filesystem>
#include <iomanip>

LDCalculator::LDCalculator(const std::string& vcf_file, 
                           size_t window_size,
                           const std::string& output_file,
                           const std::vector<std::string>& samples,
                           const int min_mac,
                           bool write_header)
    : vcf_file_(vcf_file),
      sample_list_(samples),
      min_mac_(min_mac),
      window_size_(window_size),
      write_header_(write_header) {
    
    // Open output file
    output_file_.open(output_file);
    if (!output_file_.is_open()) {
        throw std::runtime_error("Failed to open output file: " + output_file);
    }
    
    // Write header with fixed number of columns only if requested
    if (write_header_) {
        output_file_ << "chrom\tposition";
        for (size_t i = 0; i < window_size_ - 1; ++i) {
            output_file_ << "\tr2_" << (i + 1);
        }
        output_file_ << "\n";
    }
}

LDCalculator::~LDCalculator() {
    if (output_file_.is_open()) {
        output_file_.close();
    }
}

void LDCalculator::calculate_ld(const GenomeRegion &region) {
    // Process each variant
    GenomeRegion full_region = region.extend(std::min(static_cast<size_t>(1000000), region.size())); // Double region
    VcfReader vcf_reader(vcf_file_, window_size_, full_region.to_string(), sample_list_, true, min_mac_); // Enable normalization and set MAC threshold
    size_t n = 0;
    size_t extended = 0;
    std::string last_chrom;

    if (region.start > 1) { // Not the first window, so skip the first N SNPs
        for (size_t i = 0; i < window_size_ && vcf_reader.has_more(); ++i) {
            vcf_reader.next();
        }
        // auto pos_it = vcf_reader.get_positions().begin();
        //std::cerr << "prefilled buffer " << region.to_string() << " " << full_region.to_string() << " " << i << " " << window_size_ << " " << *pos_it << std::endl;
    }
    
    // Now process variants one by one
    do {  // as we pre-filled, we want at very least a single run through these in case the whole chrom only has N snps
        const auto& genotypes = vcf_reader.get_genotypes();
        const auto& chromosomes = vcf_reader.get_chromosomes();
        const auto& positions = vcf_reader.get_positions();
        
        if (genotypes.empty()) continue;
        
        // Get iterators for the first elements
        auto geno_it = genotypes.begin();
        auto chrom_it = chromosomes.begin();
        auto pos_it = positions.begin();
        
        // Current variant is the first one
        const std::string& current_chrom = *chrom_it;
        
        int current_pos = *pos_it;
        
        assert(last_chrom.empty() || last_chrom == current_chrom);
        last_chrom = current_chrom;
        
        // If we're past the region end, increment extension
        if (current_pos > region.end) {
            extended++;
        }
        
        // Start output line
        output_file_ << current_chrom << "\t" << current_pos;
        
        // Calculate r² with previous variants in the buffer
        auto geno_it2 = geno_it;
        auto chrom_it2 = chrom_it;
        auto pos_it2 = pos_it;
        ++geno_it2;
        ++chrom_it2;
        ++pos_it2;
        
        // Count how many r² values we've written
        size_t values_written = 0;
        size_t i = 1;
        
        while (geno_it2 != genotypes.end() && i < window_size_) {
            // Only calculate LD for variants on the same chromosome
            if (*chrom_it2 != current_chrom) {
                ++geno_it2;
                ++chrom_it2;
                ++pos_it2;
                ++i;
                continue;
            }
            
            double r2 = calculate_r2(*geno_it, *geno_it2);
            output_file_ << "\t" << r2;
            values_written++;
            
            ++geno_it2;
            ++chrom_it2;
            ++pos_it2;
            ++i;
        }
        
        // Fill remaining columns with tabs to ensure fixed width
        for (size_t j = values_written; j < window_size_ - 1; ++j) {
            output_file_ << "\t";
        }
        
        output_file_ << "\n";
        n++;
        if (n%100 == 0) {
            std::cerr << "Processed " << n << " SNPs in " << region.to_string() << "...\r";
        }
    } while (vcf_reader.next() & extended < window_size_);
    
    if (n > 0) {
        std::cerr << "\nCompleted processing " << n << " SNPs on region "  << region.to_string()
                  << " (extended by " << extended << " SNPs at region end)" << std::endl;
    }
}

double LDCalculator::calculate_r2(const std::vector<float>& genotype1, const std::vector<float>& genotype2) {
    // Ensure vectors are of the same size
    if (genotype1.size() != genotype2.size() || genotype1.empty()) {
        return 0.0;
    }
    
    // With normalization enabled, we shouldn't have missing values, so we can optimize
    // by avoiding the creation of temporary vectors
    
    const size_t n = genotype1.size();
    if (n < 2) return 0.0; // Need at least 2 samples for correlation
    
    // Calculate means, covariance, and variances in a single pass
    double sum1 = 0.0, sum2 = 0.0;
    
    // First pass to calculate means
    for (size_t i = 0; i < n; ++i) {
        sum1 += genotype1[i];
        sum2 += genotype2[i];
    }
    
    double mean1 = sum1 / n;
    double mean2 = sum2 / n;
    
    // Second pass to calculate covariance and variances
    double covariance = 0.0;
    double variance1 = 0.0;
    double variance2 = 0.0;
    
    for (size_t i = 0; i < n; ++i) {
        double diff1 = genotype1[i] - mean1;
        double diff2 = genotype2[i] - mean2;
        
        covariance += diff1 * diff2;
        variance1 += diff1 * diff1;
        variance2 += diff2 * diff2;
    }
    
    // Avoid division by zero
    if (variance1 <= 0.0 || variance2 <= 0.0) return 0.0;
    
    // Calculate correlation coefficient
    double r = covariance / std::sqrt(variance1 * variance2);
    
    // Return r²
    return r * r;
}
std::set<std::string> LDCalculator::get_chromosomes(const std::string& vcf_file) {
    std::set<std::string> chromosomes;
    
    // Create a temporary VCF reader to get chromosomes from header
    try {
        VcfReader reader(vcf_file, 1); // Minimal buffer size
        
        // Get chromosomes directly from the header
        auto chrom_map = reader.get_chromosomes_from_header();
        
        // Extract just the chromosome names
        for (const auto& [chrom, size] : chrom_map) {
            chromosomes.insert(chrom);
        }
        
        if (chromosomes.empty()) {
            std::cerr << "Warning: No chromosomes found in VCF header, falling back to scanning variants" << std::endl;
            
            // Fallback: scan variants if header doesn't have chromosome info
            std::string last_chrom;
            while (reader.next()) {
                const auto& chroms = reader.get_chromosomes();
                if (chroms.empty()) continue;
                
                std::string current_chrom = *chroms.begin();
                if (current_chrom != last_chrom) {
                    chromosomes.insert(current_chrom);
                    last_chrom = current_chrom;
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error scanning chromosomes: " << e.what() << std::endl;
    }
    
    return chromosomes;
}

std::map<std::string, int> LDCalculator::get_chromosome_sizes(const std::string& vcf_file) {
    std::map<std::string, int> chrom_sizes;
    
    try {
        VcfReader reader(vcf_file, 1); // Minimal buffer size
        chrom_sizes = reader.get_chromosomes_from_header();
    } catch (const std::exception& e) {
        std::cerr << "Error getting chromosome sizes: " << e.what() << std::endl;
    }
    
    return chrom_sizes;
}

std::vector<GenomeRegion> LDCalculator::divide_genome_into_regions(
    const std::string& vcf_file, int region_size) {
    
    std::vector<GenomeRegion> regions;
    auto chrom_sizes = get_chromosome_sizes(vcf_file);
    
    // Sort chromosomes to ensure consistent ordering
    std::vector<std::string> sorted_chroms;
    for (const auto& [chrom, size] : chrom_sizes) {
        sorted_chroms.push_back(chrom);
    }
    
    // Sort chromosomes naturally (e.g., chr1, chr2, ..., chr10, ...)
    std::sort(sorted_chroms.begin(), sorted_chroms.end(), 
        [](const std::string& a, const std::string& b) {
            // Extract numeric part if present
            auto extract_num = [](const std::string& s) -> int {
                std::string num_part;
                for (char c : s) {
                    if (std::isdigit(c)) num_part += c;
                }
                return num_part.empty() ? 0 : std::stoi(num_part);
            };
            
            int num_a = extract_num(a);
            int num_b = extract_num(b);
            
            if (num_a != num_b) return num_a < num_b;
            return a < b;
        });
    
    // Divide each chromosome into regions
    for (const auto& chrom : sorted_chroms) {
        int size = chrom_sizes[chrom];
        
        // Skip chromosomes with unknown or zero size
        if (size <= 0) continue;
        
        // Create regions for this chromosome
        for (int start = 1; start <= size; start += region_size) {
            int end = std::min(start + region_size - 1, size);
            regions.emplace_back(chrom, start, end, size);
        }
    }
    
    return regions;
}
