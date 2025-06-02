#include "ld_calculator.hh"
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include <filesystem>
#include <mutex>
#include <future>
#include <getopt.h>

namespace fs = std::filesystem;

// Function to display usage information
void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " [OPTIONS] <vcf_file> <output_dir>\n"
              << "Options:\n"
              << "  -w, --window-size=SIZE    Set window size for LD calculation (default: 1000)\n"
              << "  -r, --region=REGION       Restrict analysis to a specific region (format: chr:start-end)\n"
              << "  -t, --threads=NUM         Number of threads to use (default: all available cores)\n"
              << "  -m, --min-mac=COUNT       Minimum minor allele count threshold (default: 1)\n"
              << "  -s, --region-size=SIZE    Size of each genomic region in bp (default: 5000000)\n"
              << "  -h, --help                Display this help message and exit\n"
              << std::endl;
}

// Function to process a single region
void process_region(const std::string& vcf_file, 
                    size_t window_size, 
                    const std::string& output_dir,
                    const GenomeRegion& region,
                    const int chunk_id,
                    const std::vector<std::string>& samples,
                    int min_mac,
                    bool write_header) {
    try {
        // Create output filename based on chunk ID with leading zeros
        std::stringstream ss;
        ss << std::setw(6) << std::setfill('0') << chunk_id;
        std::string output_file = (fs::path(output_dir) / (ss.str() + ".tsv")).string();
        
        // Create calculator for this region
        LDCalculator calculator(vcf_file, window_size, output_file, samples, min_mac, write_header);
        
        // Calculate LD for this region
        calculator.calculate_ld(region);
    } catch (const std::exception& e) {
        std::cerr << "Error processing region " << region.to_string() << ": " << e.what() << std::endl;
    }
}

int main(int argc, char *argv[])
{
    // Default parameter values
    size_t window_size = 1000;
    std::string region = "";
    unsigned int thread_count = std::thread::hardware_concurrency();
    int min_mac = 0;
    int region_size = 5000000; // Default region size: 5Mb
    
    // Define long options
    static struct option long_options[] = {
        {"window-size", required_argument, 0, 'w'},
        {"region",      required_argument, 0, 'r'},
        {"threads",     required_argument, 0, 't'},
        {"min-mac",     required_argument, 0, 'm'},
        {"region-size", required_argument, 0, 's'},
        {"help",        no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };
    
    // Parse command line options
    int opt;
    int option_index = 0;
    
    while ((opt = getopt_long(argc, argv, "w:r:t:m:s:h", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'w':
                window_size = std::stoul(optarg);
                break;
            case 'r':
                region = optarg;
                break;
            case 't':
                thread_count = std::stoul(optarg);
                if (thread_count == 0) thread_count = 1;
                break;
            case 'm':
                min_mac = std::stoi(optarg);
                if (min_mac < 0) {
                    std::cerr << "Error: Minimum minor allele count must be non-negative" << std::endl;
                    return 1;
                }
                break;
            case 's':
                region_size = std::stoi(optarg);
                if (region_size <= 0) {
                    std::cerr << "Error: Region size must be positive" << std::endl;
                    return 1;
                }
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            case '?':
                // getopt_long already printed an error message
                print_usage(argv[0]);
                return 1;
            default:
                abort();
        }
    }
    
    // Check for required positional arguments
    if (optind + 1 >= argc) {
        std::cerr << "Error: Missing required arguments" << std::endl;
        print_usage(argv[0]);
        return 1;
    }
    
    std::string vcf_file = argv[optind];
    std::string output_dir = argv[optind + 1];
    
    try {
        // Create output directory if it doesn't exist
        fs::create_directories(output_dir);
        
        std::cerr << "Calculating LD for " << vcf_file << " with window size " << window_size << std::endl;
        std::cerr << "Using " << thread_count << " threads" << std::endl;
        std::cerr << "Using region size of " << region_size << " bp" << std::endl;
        if (!region.empty()) {
            std::cerr << "Using region: " << region << std::endl;
        }
        if (min_mac > 0) {
            std::cerr << "Using minimum minor allele count threshold: " << min_mac << std::endl;
        }
        
        std::vector<GenomeRegion> regions;
        if (!region.empty()) {
            // If region is specified, split it according to region size
            GenomeRegion reg(region);
            regions = reg.subdivide(region_size);
        } else {
            // Divide genome into regions
            regions = LDCalculator::divide_genome_into_regions(vcf_file, region_size);
        }
        
        if (regions.empty()) {
            std::cerr << "No regions found in VCF file" << std::endl;
            return 1;
        }
        
        std::cerr << "Calculating in parallel across " << regions.size() << " regions" << std::endl;
        
        // Create thread pool
        std::vector<std::future<void>> futures;
        
        // Process each region in a separate thread
        for (size_t i = 0; i < regions.size(); ++i) {
            // Only the first chunk should write the header
            bool write_header = (i == 0);
            
            // If we have too many active threads, wait for some to complete
            while (futures.size() >= thread_count) {
                // Wait for at least one thread to complete
                for (auto it = futures.begin(); it != futures.end(); ) {
                    if (it->wait_for(std::chrono::milliseconds(1)) == std::future_status::ready) {
                        it = futures.erase(it);
                        break;
                    } else {
                        ++it;
                    }
                }
                
                // If we still have too many threads, sleep briefly
                if (futures.size() >= thread_count) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(100));
                    continue;
                }
            }
            
            // Launch a new thread for this region
            futures.push_back(std::async(std::launch::async, 
                                        process_region, 
                                        vcf_file, 
                                        window_size, 
                                        output_dir, 
                                        regions[i], 
                                        i,
                                        std::vector<std::string>(),
                                        min_mac,
                                        write_header));
        }
        
        // Wait for all threads to complete
        for (auto& future : futures) {
            future.wait();
        }
        
        std::cerr << "LD calculation complete. Results written to " << output_dir << std::endl;
        std::cerr << "To combine results: cat " << output_dir << "/*.tsv > combined.ld" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
