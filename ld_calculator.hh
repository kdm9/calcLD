#ifndef LD_CALCULATOR_H
#define LD_CALCULATOR_H

#include "vcf_reader.hh"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <list>
#include <mutex>
#include <set>
#include <map>

// Structure to represent a genomic region
struct GenomeRegion {
    std::string chrom;
    int chrom_len;
    int start;
    int end;
    
    // Constructor
    GenomeRegion(const std::string& c, int s, int e, int l=-1) 
        : chrom(c), chrom_len(l), start(s), end(e) {}
    //
    // Constructor
    GenomeRegion(const std::string& reg) {
        const auto i=reg.find(":");
        const auto j=reg.find("-", i);
        assert(i!=std::string::npos && j!=std::string::npos);
        chrom = reg.substr(0, i);
        start = std::stoi(reg.substr(i+1, j));
        end = std::stoi(reg.substr(j+1, std::string::npos));
    }
    
    
    // Convert to region string format (e.g., "chr1:1000-2000")
    std::string to_string() const {
        return chrom + ":" + std::to_string(start) + "-" + std::to_string(end);
    }

    std::vector<GenomeRegion> subdivide(int length) const {
        std::vector<GenomeRegion> r;
        for (int i = start; i<end; i+= length) {
            int e = std::min(i + length - 1, end);
            r.emplace_back(chrom, i, e);
        }
        return r;
    }
    GenomeRegion extend(int length) const {
        GenomeRegion r(chrom, start, std::min(end+length, chrom_len), chrom_len);
        return r;
    }

    size_t size() const {
        return end-start+1;
    }
};

class LDCalculator {
public:
    /**
     * Constructor for LDCalculator
     * @param vcf_file Path to the VCF file
     * @param window_size Number of variants to keep in the buffer for LD calculation
     * @param output_file Path to the output file
     * @param samples Optional list of samples to subset
     * @param min_mac Minimum minor allele count threshold
     * @param write_header Whether to write the header to the output file
     */
    LDCalculator(const std::string& vcf_file, 
                 size_t window_size,
                 const std::string& output_file,
                 const std::vector<std::string>& samples = {},
                 const int min_mac = 0,
                 bool write_header = true);
    
    /**
     * Destructor
     */
    ~LDCalculator();
    
    /**
     * Calculate LD and write results to output file
     */
    void calculate_ld(const GenomeRegion &region);
    
    /**
     * Get list of all chromosomes in the VCF file
     * @param vcf_file Path to the VCF file
     * @return Set of chromosome names
     */
    static std::set<std::string> get_chromosomes(const std::string& vcf_file);
    
    /**
     * Get chromosome sizes from the VCF file
     * @param vcf_file Path to the VCF file
     * @return Map of chromosome names to their sizes
     */
    static std::map<std::string, int> get_chromosome_sizes(const std::string& vcf_file);
    
    /**
     * Divide genome into regions of specified size
     * @param vcf_file Path to the VCF file
     * @param region_size Size of each region in base pairs
     * @return Vector of GenomeRegion objects
     */
    static std::vector<GenomeRegion> divide_genome_into_regions(
        const std::string& vcf_file, int region_size);
    
private:
    std::string vcf_file_;
    std::vector<std::string> sample_list_;
    int min_mac_;
    std::ofstream output_file_;
    size_t window_size_;
    bool write_header_;
    
    /**
     * Calculate r² (correlation coefficient squared) between two genotype vectors
     * @param genotype1 First genotype vector
     * @param genotype2 Second genotype vector
     * @return r² value
     */
    double calculate_r2(const std::vector<float>& genotype1, const std::vector<float>& genotype2);
};

#endif // LD_CALCULATOR_H
