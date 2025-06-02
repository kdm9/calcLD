#ifndef VCF_READER_H
#define VCF_READER_H

#include <string>
#include <vector>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <map>

class VcfReader {
public:
    /**
     * Constructor for VcfReader
     * @param filename Path to the VCF file
     * @param buffer_size Number of variants to keep in the buffer
     * @param region Optional region string (e.g., "chr1:1000-2000")
     * @param samples Optional list of samples to subset
     * @param normalise Whether to normalize genotypes by imputing missing values
     * @param min_mac Minimum minor allele count threshold
     */
    VcfReader(const std::string& filename, 
              size_t buffer_size,
              const std::string& region = "",
              const std::vector<std::string>& samples = {},
              const bool normalise = true,
              const int min_mac = 0);
    
    /**
     * Destructor
     */
    ~VcfReader();
    
    /**
     * Move to the next variant
     * @return true if successful, false if no more variants
     */
    bool next();
    
    /**
     * Check if there are more variants to read
     * @return true if there are more variants
     */
    bool has_more() const;
    
    /**
     * Get the current buffer of genotypes
     * @return Reference to the deque of genotype vectors
     */
    const std::list<std::vector<float>>& get_genotypes() const;
    
    /**
     * Get the current buffer of chromosomes
     * @return Reference to the list of chromosomes
     */
    const std::list<std::string>& get_chromosomes() const;
    
    /**
     * Get the current buffer of positions
     * @return Reference to the list of positions
     */
    const std::list<int>& get_positions() const;
    
    /**
     * Get the current buffer of info fields
     * @return Reference to the list of info field maps
     */
    const std::list<std::unordered_map<std::string, std::string>>& get_info_fields() const;
    
    /**
     * Get the number of samples in the VCF
     * @return Number of samples
     */
    int get_num_samples() const;
    
    /**
     * Get the sample names
     * @return Vector of sample names
     */
    const std::vector<std::string>& get_sample_names() const;
    
    /**
     * Get list of all chromosomes and their sizes from the VCF header
     * @return Map of chromosome names to their sizes
     */
    std::map<std::string, int> get_chromosomes_from_header() const;

private:
    // File handling
    std::string filename_;
    bcf_srs_t* sr_ = nullptr;
    bcf_hdr_t* header_ = nullptr;
    bool normalise_variants_ = false;
    
    // Buffer parameters
    size_t buffer_size_;
    bool has_more_ = true;
    
    // Sample handling
    std::vector<std::string> sample_names_;
    std::unordered_set<int> sample_indices_;
    int num_samples_ = 0;
    int min_mac_ = 0;  // Minimum minor allele count threshold
    
    // Data buffers - using list for better insertion/deletion performance
    std::list<std::vector<float>> genotypes_;
    std::list<std::string> chromosomes_;
    std::list<int> positions_;
    std::list<std::unordered_map<std::string, std::string>> info_fields_;
    
    // Helper methods
    void initialize_reader(const std::string& region, const std::vector<std::string>& samples);
    bool read_variant();
    std::vector<float> extract_genotypes(bcf1_t* record);
    std::unordered_map<std::string, std::string> extract_info_fields(bcf1_t* record);
};

#endif // VCF_READER_H
