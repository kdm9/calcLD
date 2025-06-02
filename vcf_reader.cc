#include "vcf_reader.hh"
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <boost/math/statistics/univariate_statistics.hpp>

VcfReader::VcfReader(const std::string& filename, 
                     size_t buffer_size,
                     const std::string& region,
                     const std::vector<std::string>& samples,
                     const bool normalise,
                     const int min_mac) 
    : filename_(filename), buffer_size_(buffer_size), normalise_variants_(normalise), min_mac_(min_mac) {
    
    initialize_reader(region, samples);
    
    // Fill the initial buffer
    //for (size_t i = 0; i < buffer_size_ && has_more_; ++i) {
    //    if (!read_variant()) {
    //        break;
    //    }
    //}
}

VcfReader::~VcfReader() {
    if (sr_) {
        bcf_sr_destroy(sr_);
        sr_ = nullptr;
    }
    // header is cleaned up by bcf_sr_destroy
}

void VcfReader::initialize_reader(const std::string& region, const std::vector<std::string>& samples) {
    // Initialize the synced reader
    sr_ = bcf_sr_init();
    
    // Set region if provided
    if (!region.empty()) {
        if (bcf_sr_set_regions(sr_, region.c_str(), 0) != 0) {
            throw std::runtime_error("Failed to set region: " + region);
        }
    }
    
    // Open the VCF file
    if (bcf_sr_add_reader(sr_, filename_.c_str()) != 1) {
        std::string error = "Failed to open VCF file: " + filename_;
        if (sr_->errnum) {
            error += " - " + std::string(bcf_sr_strerror(sr_->errnum));
        }
        bcf_sr_destroy(sr_);
        sr_ = nullptr;
        throw std::runtime_error(error);
    }
    
    // Get the header
    header_ = bcf_sr_get_header(sr_, 0);
    if (!header_) {
        bcf_sr_destroy(sr_);
        sr_ = nullptr;
        throw std::runtime_error("Failed to get VCF header");
    }
    
    // Get sample names from the header
    num_samples_ = bcf_hdr_nsamples(header_);
    for (int i = 0; i < num_samples_; ++i) {
        sample_names_.push_back(header_->samples[i]);
    }
    
    // Process sample subset if provided
    if (!samples.empty()) {
        std::unordered_set<std::string> sample_set(samples.begin(), samples.end());
        std::vector<std::string> filtered_samples;
        
        for (int i = 0; i < num_samples_; ++i) {
            if (sample_set.find(sample_names_[i]) != sample_set.end()) {
                sample_indices_.insert(i);
                filtered_samples.push_back(sample_names_[i]);
            }
        }
        
        if (sample_indices_.empty()) {
            bcf_sr_destroy(sr_);
            sr_ = nullptr;
            throw std::runtime_error("None of the specified samples found in VCF");
        }
        
        sample_names_ = filtered_samples;
        num_samples_ = static_cast<int>(sample_names_.size());
    } else {
        // If no subset, use all samples
        for (int i = 0; i < num_samples_; ++i) {
            sample_indices_.insert(i);
        }
    }
}

bool VcfReader::next() {
    if (!has_more_) {
        return false;
    }
    
    // Read the next variant
    if (!read_variant()) {
        return false;
    }
    
    // Remove the oldest variant if buffer is full
    if (genotypes_.size() > buffer_size_) {
        genotypes_.pop_back();
        chromosomes_.pop_back();
        positions_.pop_back();
        if (!info_fields_.empty()) {
            info_fields_.pop_back();
        }
    }
    
    return true;
}

bool VcfReader::read_variant() {
    if (!sr_ || !has_more_) {
        return false;
    }
    
    while (true) {
        // Read the next record
        if (bcf_sr_next_line(sr_) <= 0) {
            has_more_ = false;
            return false;
        }
        
        bcf1_t* record = bcf_sr_get_line(sr_, 0);
        if (!record) {
            has_more_ = false;
            return false;
        }
        
        // Extract genotypes first to check MAC
        std::vector<float> genos = extract_genotypes(record);
        
        // Check if variant passes the minor allele count threshold
        if (min_mac_ > 0) {
            int alt_count = 0;
            int total_alleles = 0;
            
            for (const float& val : genos) {
                if (val != -1) {  // Skip missing values
                    alt_count += static_cast<int>(val);
                    total_alleles += 2;  // Assuming diploid
                }
            }
            
            // Calculate minor allele count
            int mac = std::min(alt_count, total_alleles - alt_count);
            
            // Skip variant if it doesn't meet the threshold
            if (mac < min_mac_) {
                continue;  // Try the next variant
            }
        }
        
        // Extract chromosome
        const char* chrom = bcf_seqname(header_, record);
        chromosomes_.push_front(chrom);
        
        // Extract position (0-based to 1-based)
        positions_.push_front(record->pos + 1);
        
        // Store the genotypes
        genotypes_.push_front(genos);
        
        // Extract INFO fields - only if needed
        if (!info_fields_.empty() || info_fields_.size() < buffer_size_) {
            info_fields_.push_front(extract_info_fields(record));
        }
        
        return true;
    }
}

std::vector<float> VcfReader::extract_genotypes(bcf1_t* record) {
    // Prepare result vector (initialize with missing values)
    std::vector<float> result(num_samples_, -1);
    
    // Get the genotypes
    int32_t* gt_arr = nullptr;
    int n_gt = 0;
    int n_alleles = record->n_allele;
    
    // Extract the genotype field
    int ret = bcf_get_genotypes(header_, record, &gt_arr, &n_gt);
    
    if (ret <= 0 || !gt_arr) {
        // No genotypes found, return vector of missing values
        free(gt_arr);
        return result;
    }
    
    // Calculate ploidy (usually 2 for diploid)
    int ploidy = n_gt / bcf_hdr_nsamples(header_);
    
    // Process each sample
    int result_idx = 0;
    for (int i = 0; i < bcf_hdr_nsamples(header_); ++i) {
        // Skip samples not in our subset
        if (!sample_indices_.empty() && sample_indices_.find(i) == sample_indices_.end()) {
            continue;
        }
        
        // Calculate dosage of non-reference allele
        int dosage = 0;
        bool has_missing = false;
        
        for (int j = 0; j < ploidy; ++j) {
            int idx = i * ploidy + j;
            if (idx >= n_gt) break;
            
            // Check for missing genotype
            if (gt_arr[idx] == bcf_int32_vector_end || bcf_gt_is_missing(gt_arr[idx])) {
                has_missing = true;
                break;
            }
            
            // Get allele index (0 = ref, 1+ = alt)
            int allele_idx = bcf_gt_allele(gt_arr[idx]);
            
            // Count non-reference alleles
            if (allele_idx > 0 && allele_idx < n_alleles) {
                dosage++;
            }
        }
        
        // Set the result
        result[result_idx++] = has_missing ? -1 : dosage;
    }

    if (normalise_variants_) {
        // Normalize genotypes by imputing missing values with the mean - more efficient implementation
        float sum = 0.0f;
        int count = 0;
    
        // First pass: calculate mean in one loop
        for (const float& val : result) {
            if (val != -1) {
                sum += val;
                count++;
            }
        }
    
        // Calculate mean and impute missing values
        if (count > 0) {
            float mean = sum / count;
            int imputed_value = std::round(mean);
        
            // Replace missing values with the mean
            for (float& val : result) {
                if (val == -1) {
                    val = imputed_value;
                }
            }
        }
    }
    
    free(gt_arr);
    return result;
}

std::unordered_map<std::string, std::string> VcfReader::extract_info_fields(bcf1_t* record) {
    std::unordered_map<std::string, std::string> info;
    
    // Iterate through all INFO fields in the header
    for (int i = 0; i < header_->n[BCF_DT_ID]; ++i) {
        bcf_idpair_t* idpair = &header_->id[BCF_DT_ID][i];
        if (!idpair->key || !idpair->val) continue;
        
        bcf_info_t* info_field = bcf_get_info(header_, record, idpair->key);
        if (!info_field) continue;
        
        // Convert the info field to string based on its type
        std::string value;
        switch (info_field->type) {
            case BCF_BT_NULL:
                value = "true"; // Flag type
                break;
            case BCF_BT_INT8:
            case BCF_BT_INT16:
            case BCF_BT_INT32:
                value = std::to_string(info_field->v1.i);
                break;
            case BCF_BT_FLOAT:
                value = std::to_string(info_field->v1.f);
                break;
            case BCF_BT_CHAR:
                value = std::string((char*)info_field->vptr, info_field->len);
                break;
            default:
                // For complex types, just indicate presence
                value = "present";
                break;
        }
        
        info[idpair->key] = value;
    }
    
    return info;
}

bool VcfReader::has_more() const {
    return has_more_;
}

const std::list<std::vector<float>>& VcfReader::get_genotypes() const {
    return genotypes_;
}

const std::list<std::string>& VcfReader::get_chromosomes() const {
    return chromosomes_;
}

const std::list<int>& VcfReader::get_positions() const {
    return positions_;
}

const std::list<std::unordered_map<std::string, std::string>>& VcfReader::get_info_fields() const {
    return info_fields_;
}

int VcfReader::get_num_samples() const {
    return num_samples_;
}

const std::vector<std::string>& VcfReader::get_sample_names() const {
    return sample_names_;
}

std::map<std::string, int> VcfReader::get_chromosomes_from_header() const {
    std::map<std::string, int> chromosomes;
    
    if (!header_) {
        return chromosomes;
    }

    int32_t nctg = header_->n[BCF_DT_CTG];
    for (int32_t i = 0; i < nctg; i++) {
      bcf_idpair_t *ctg = header_->id[BCF_DT_CTG];
      chromosomes[(ctg[i].key)]=  (ctg[i].val->info[0]);
    }
    return chromosomes;
}
