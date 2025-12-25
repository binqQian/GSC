
#pragma once

#include <iostream>
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "defs.h"
#include "gsc_params.h"
#include "bit_memory.h"
#include "queues.h"
#include "variant.h"
#include <string>
#include <vector>
#include "utils.h"
#include "file_handle.h"
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <memory>
#include <cstdint>
#include "parallel_vcf_reader.h"
#include "adaptive_format_known_fields.h"

class CompressionReader {

    
	htsFile* in_file = nullptr;
    vector<htsFile*> merge_files;
    bcf_hdr_t * vcf_hdr= nullptr;
    // vector<bcf_hdr_t*> vcf_hdrs;
	bcf1_t * vcf_record;
    file_type in_type;
    std::string in_file_name;
    bool in_open;
    bool vcf_hdr_read;
    compress_mode_t compress_mode;
    compression_backend_t backend_;
    bool merge_flag;
    bool merge_failure_flag;

    uint32_t no_samples;
    uint32_t ploidy;
    uint64_t vec_len;
    uint64_t no_vec;

    // GT column tiling metadata
    uint32_t max_block_cols;            // Max haplotypes per column block (from params)
    uint32_t n_col_blocks;              // Number of column blocks
    uint32_t total_haplotypes;          // Total number of haplotypes
    vector<uint32_t> col_block_sizes;   // Haplotypes per column block
    vector<uint64_t> col_block_vec_lens; // vec_len for each column block

    // Tiled mode: per-column-block buffers and counters
    vector<CBitMemory> col_bv_buffers;   // One buffer per column block
    vector<uint32_t> col_vec_read_in_block;  // Vectors read in current block for each column

    vector<string> samples_list;
    vector<variant_desc_t> v_vcf_data_compress;
    int32_t *cur_g_data = nullptr;
    int ncur_g_data = 0;
    int *gt_data = nullptr;
    int temp;
    int64_t cur_pos;
    int32_t tmpi;
    size_t chunk_size;  

    CBitMemory bv;
    fixed_field_block fixed_field_block_buf;
    int64_t block_max_size;
    uint32_t no_vec_in_block, vec_read_in_block, block_id;
    uint32_t no_fixed_fields,vec_read_fixed_fields,fixed_fields_id;
    GtBlockQueue * Gt_queue = nullptr;
    VarBlockQueue<fixed_field_block>* Var_queue = nullptr;


    void *dst_int = nullptr;
    void *dst_real = nullptr;
    void *dst_str = nullptr;
    void *dst_flag = nullptr;
    int  ndst_int = 0;
    int  ndst_real = 0;
    int  ndst_str = 0;
    int  ndst_flag = 0;
    std::vector<int> FilterIdToFieldId;
    std::vector<int> InfoIdToFieldId;
    std::vector<int> FormatIdToFieldId;
    int no_flt_keys, no_info_keys, no_fmt_keys;
    uint32_t no_keys;
    vector<key_desc> keys;
    int key_gt_id;
    uint32_t no_actual_variants;
    vector<uint32_t> actual_variants;
    vector<uint32_t> v_size;
	vector<uint8_t> v_data;	
    vector<CBuffer> v_o_buf;
    vector<int> v_buf_ids_size;
	vector<int> v_buf_ids_data;
    File_Handle_2 * file_handle2 = nullptr;
    PartQueue<SPackage> * part_queue = nullptr;

    bool use_adaptive_format_ = false;  // Enable adaptive FORMAT compression
    bool adaptive_format_primary_ = false;  // Omit legacy known FORMAT fields when true
    adaptive_format_part_backend_t adaptive_format_part_backend_ = adaptive_format_part_backend_t::auto_select;
    gsc::AdaptiveKnownFieldsIndices adaptive_known_indices_;
    gsc::AdaptiveKnownFieldsDicts adaptive_known_dicts_;

    // Stream IDs for adaptive FORMAT compression
    int adaptive_format_stream_id_ = -1;
    std::vector<uint8_t> adaptive_format_buffer_;  // Buffer for adaptive FORMAT data
    int adaptive_format_ad_dict_stream_id_ = -1;
    int adaptive_format_pl_dict_stream_id_ = -1;
    int adaptive_format_pid_dict_stream_id_ = -1;

    unordered_map<int, unordered_set<int>> field_order_graph;
    bool field_order_flag;
    // unordered_map<int, unordered_set<int>> graph;
    unordered_map<int, int> inDegree;
    vector<int> order;
    size_t no_chrom_num;
    string no_chrom;
    vector<pair<std::string,uint32_t>> where_chrom;
    string cur_chrom;
    vector<int64_t> chunks_min_pos;
    bool start_flag;

	    // Parallel VCF reading
	    gsc::ParallelVCFReader* parallel_reader_;
	    int num_parse_threads_;
	    bool use_parallel_reading_;
        uint64_t max_memory_mb_;

    // int temp_count = 0;
    #ifdef LOG_INFO
	unordered_map<int, unordered_set<int>> distinct_values;
    #endif
    bool ReadFile();
    bool setBitVector();
	void addVariant(int * gt_data, int ngt_data,variant_desc_t &desc);
    bool GetVariantFromRec(bcf1_t* rec, vector<field_desc>& fields);
    bool GetFilterInfoFormatKeys(int &no_flt_keys, int &no_info_keys, int &no_fmt_keys, vector<key_desc> &keys);
    void ProcessFixedVariants(bcf1_t *vcf_record, variant_desc_t &desc);
	bool SetVariantOtherFields(bcf1_t *vcf_rec, vector<field_desc> &fields);
    vector<int> topoSort(unordered_map<int, unordered_set<int>>& graph);
    vector<int> topo_sort(unordered_map<int, unordered_set<int>> &graph,unordered_map<int, int> inDegree);


    
public:

	    CompressionReader() {
        in_open = false;
        vcf_hdr_read = false;
        no_samples = 0;
        ploidy = 0;
        max_block_cols = 0;
        n_col_blocks = 1;
        no_chrom_num = 0;
        no_vec = 0;
        start_flag = true;
        field_order_flag = false;
	        parallel_reader_ = nullptr;
	        num_parse_threads_ = 1;
	        use_parallel_reading_ = false;
            max_memory_mb_ = 0;
	    }
	    CompressionReader(const GSC_Params & params) {
        in_open = false;
        vcf_hdr_read = false;
        no_samples = 0;
        ploidy = params.ploidy;
        max_block_cols = params.max_block_cols;
        n_col_blocks = 1;
        no_chrom_num = 0;
        in_file_name = params.in_file_name;
        in_type = params.in_type;
        compress_mode = params.compress_mode;
        backend_ = params.backend;
        merge_flag = params.merge_file_flag;
        merge_failure_flag = false;
        v_vcf_data_compress.reserve(no_variants_in_buf);
        actual_variants.reserve(no_variants_in_buf);
        no_vec = 0;
        start_flag = true;
        field_order_flag = false;
	        parallel_reader_ = nullptr;
	        num_parse_threads_ = (params.no_parse_threads == 0) ? 1 : static_cast<int>(params.no_parse_threads);
	        use_parallel_reading_ = (num_parse_threads_ > 1);
            max_memory_mb_ = params.max_memory_mb;

		        // Adaptive FORMAT compression (known-field codec; see adaptive_format_known_fields.*)
		        use_adaptive_format_ = (params.compress_mode == compress_mode_t::lossless_mode) &&
		                              (params.adaptive_format_mode != adaptive_format_mode_t::off);
	        adaptive_format_primary_ = (params.adaptive_format_mode == adaptive_format_mode_t::primary);
	        adaptive_format_part_backend_ = params.adaptive_format_part_backend;
	    }
    
  
    
    ~CompressionReader() {

    if(in_open)
    {
        bcf_hdr_destroy(vcf_hdr);
        vcf_hdr = nullptr;
        hts_close(in_file);
        in_file = nullptr;
        in_open = false;
        // bcf_destroy1(vcf_record);
    }

    if (parallel_reader_) {
        delete parallel_reader_;
        parallel_reader_ = nullptr;
    }
        // test.close();
    //    var_out.close();
    }
          
    void setQueue(GtBlockQueue * _queue)
    {
        Gt_queue = _queue;
    };
    void setPartQueue(PartQueue<SPackage> *_part_queue){
        part_queue = _part_queue;
    }
    uint64_t getNoVec()
    {
        return no_vec;
    }

    bool OpenForReading(string & file_name);
    uint32_t GetSamples(vector<string> &s_list);
    bool GetHeader(string &v_header);
    void InitVarinats(File_Handle_2 *_file_handle2);
	bool ProcessInVCF();
	uint32_t setNoVecBlock(GSC_Params & params);
    void initializeColumnBlocks();
	void GetWhereChrom(vector<pair<std::string,uint32_t>> &_where_chrom,vector<int64_t> &chunks_min_pos);
    uint32_t GetOtherFieldsBlockSum();
    void GetOtherField(vector<key_desc> &_keys,uint32_t &_no_keys,int &_key_gt_id);
    vector<uint32_t> GetActualVariants();
    void UpdateKeys(vector<key_desc> &_keys);
    void CloseFiles();

    // Get column block tiling information
    uint32_t GetNumColumnBlocks() const { return n_col_blocks; }
    const vector<uint32_t>& GetColumnBlockSizes() const { return col_block_sizes; }
    uint32_t GetTotalHaplotypes() const { return total_haplotypes; }
    uint32_t GetMaxBlockCols() const { return max_block_cols; }
    const vector<uint64_t>& GetColumnBlockVecLens() const { return col_block_vec_lens; }
};
