#include <exception>
#include <string>
#include <vector>
#include <iostream>
#include <numeric>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

namespace KDMBCF {

using namespace std;

class BCFReader
{
public:
  template <typename T> using Mat2d = vector<vector<T>>;
  BCFReader (const string &filename_) {
    filename = filename_;
    if((bcf = bcf_open(filename.c_str(), "r")) == NULL) {
      throw runtime_error("Unable to open file.");
    }
    if (hts_set_threads(bcf, 1) != 0) {
      throw runtime_error("Unable to open file (set threads).");
    }
    if((header = bcf_hdr_read(bcf)) == NULL) {
      throw runtime_error("Unable to read header.");
    }
    nsamp = bcf_hdr_nsamples(header);
  }
  ~BCFReader ()
  {
    bcf_hdr_destroy(header);
    bcf_close(bcf);
  }

  vector<string> get_sample_names()
  {
      vector<string> samples;
      for (int32_t i = 0; i < nsamp; i++) {
          samples.push_back(header->samples[i]);
      }
      return samples;
  }

  bool set_samples(const vector<string> &samples)
  {
      samplecsv = "";
      for (size_t i = 0; i < samples.size(); i++) {
          samplecsv += samples[i];
          if (i < samples.size() - 1)
              samplecsv += ",";
      }
      if (bcf_hdr_set_samples(header, samplecsv.c_str(), 0) != 0) {
          return false;
      }
      nsamp = bcf_hdr_nsamples(header);
      if (nsamp != samples.size()) return false;
      return true;
  }

  bool get_contig_names_lengths(vector<string> &names, vector<int32_t> &lengths)
  {
    int32_t nctg = header->n[BCF_DT_CTG];
    names.clear();
    lengths.clear();
    for (int32_t i = 0; i < nctg; i++) {
      bcf_idpair_t *ctg = header->id[BCF_DT_CTG];
      names.emplace_back(ctg[i].key);
      lengths.emplace_back(ctg[i].val->info[0]);
    }
    return nctg > 0;
  }

  bool read_chunk(const string &region)
  {
    bool ret = true;
    CHROM.clear();
    POS.clear();
    GT.clear();
#ifdef KDM_BCF_READER_USE_AD
    AD_ref.clear();
    AD_alt.clear();
#endif

    bcf_srs_t *sr = bcf_sr_init();
    //bcf_sr_set_threads(sr, 1);
    bcf_sr_set_opt(sr, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_BOTH_REF);
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    if (samplecsv.size() > 0) {
      if (bcf_sr_set_samples(sr, samplecsv.c_str(), 0) != 1) {
        cerr << "Failed to set samples for " << filename << ": " << samplecsv << endl;
        ret = false;
      }
    }
    if (region.size() > 0) {
      if (bcf_sr_set_regions(sr, region.c_str(), 0) != 0) {
        cerr << "Failed to set region for " << filename << endl;
        ret = false;
      }
    }
    if (bcf_sr_add_reader(sr, filename.c_str()) != 1) {
      cerr << "Failed to add reader for " << filename << endl;
      ret = false;
    }
    hts_set_threads(sr->readers[0].file, 1);

    while (ret && bcf_sr_next_line(sr)) {
        bcf1_t *record = bcf_sr_get_line(sr,0);
        process_record(record);
    }
    if (sr->errnum) {
      cerr << "Error: " << bcf_sr_strerror(sr->errnum) << endl;
      ret = false;
    }
    bcf_sr_destroy(sr);
    return ret;
  }

  vector< string > CHROM;
  vector< uint64_t > POS;
  Mat2d<int32_t> GT;
#ifdef KDM_BCF_READER_USE_AD
  Mat2d<int32_t> AD_ref;
  Mat2d<int32_t> AD_alt;
#endif

protected:
  htsFile *bcf = NULL;
  bcf_hdr_t *header = NULL;
  int32_t *buffer = NULL;
  int32_t buffersz = 0;
  int32_t nsamp = 0;
  string filename;
  string samplecsv;

  void process_record(bcf1_t * record)
  {
      vector<int32_t> GT1;
#ifdef KDM_BCF_READER_USE_AD
      vector<int32_t> AD_ref1, AD_alt1;
#endif
      if (!bcf_is_snp(record) ||   record->n_allele > 2) return;
      if (!get_GT(GT1, record)) return;
      GT.push_back(GT1);
#ifdef KDM_BCF_READER_USE_AD
      if (!get_AD(AD_ref1, AD_alt1)) return;
      AD_ref.push_back(AD_ref1);
      AD_alt.push_back(AD_alt1);
#endif
      CHROM.emplace_back(bcf_hdr_id2name(header, record->rid));
      POS.push_back(record->pos);
  }

  int32_t fill_buffer(const char *tag, bcf1_t *record) {
    int32_t nentry = bcf_get_format_int32(header, record, tag, &buffer, &buffersz);
    int32_t ploid = nentry/nsamp;
    if (nentry < 0) return nentry;
    if (ploid != 2) return 0;
    return nentry;
  }

  bool get_GT(vector<int32_t> &GT1, bcf1_t *record)
  {
    int32_t n = fill_buffer("GT", record);
    if (n < 1) return false;
    GT1.resize(nsamp);
    for (int32_t i=0; i<nsamp; i++) {
      int32_t a=bcf_gt_allele(buffer[i]), b=bcf_gt_allele(buffer[i+1]);
      int32_t altcnt = (a<0 || b< 0) ? -1 : a + b;
      GT1[i] = altcnt;
    }
    return true;
  }

  bool get_AD(vector<int32_t> &AD_ref1, vector<int32_t> &AD_alt1, bcf1_t *record)
  {
    int32_t n = fill_buffer("AD", record);
    int32_t ploid = n/nsamp;
    if (n < 1) return false;

    AD_ref1.resize(nsamp);
    AD_alt1.resize(nsamp);
    for (int32_t i=0; i<nsamp; i++) {
      int32_t *ADs = buffer + i*ploid;
      AD_ref1[i] = ADs[0];
      AD_alt1[i] = ADs[1];
    }
    return true;
  }

};


} /* namespace gusld */
