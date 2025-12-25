#pragma once

#include <stdio.h>
// #include <filesystem>
#include <iostream>
#include <string>
#include <memory>
#include <sdsl/bit_vectors.hpp>
#include<future>
#include "bit_memory.h"
#include "defs.h"
#include "variant.h"
#include "file_handle.h"
#include "queues.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "utils.h"
#include "compression_strategy.h"
#include "vint_code.h"
#include "fmt_field_processor.h"
#define MMAP

#ifdef MMAP
#include <cpp-mmf/memory_mapped_file.hpp>
#endif

class DecompressionReader {

    friend class Decompressor;

    struct fixed_fields_rb_dir_entry
    {
        fixed_fields_row_block_meta meta;
        uint32_t chrom_off = 0, chrom_size = 0;
        uint32_t pos_off = 0, pos_size = 0;
        uint32_t id_off = 0, id_size = 0;
        uint32_t ref_off = 0, ref_size = 0;
        uint32_t alt_off = 0, alt_size = 0;
        uint32_t qual_off = 0, qual_size = 0;
        uint32_t gt_off = 0, gt_size = 0; // v2+: per-row_block GT index segment
    };
    
	
    // sdsl vectors and ranks
    sdsl::rrr_vector<> rrr_copy_bit_vector[2];  //for copy vectors (even and odd)
    sdsl::rrr_vector<> rrr_zeros_bit_vector[2];  //for zeros vectors (even and odd)

	uint8_t * buf = nullptr;
	uint64_t buf_pos;
#ifdef MMAP

    memory_mapped_file::read_only_mmf *fm;
#endif

	
    uint64_t vec_len;
    // uint32_t max_no_vec_in_block;
    uint32_t ploidy;
    uint32_t max_block_rows = 0;
    uint32_t max_block_cols = 0;
    uint32_t n_samples;

    // GT column tiling metadata
    uint32_t n_col_blocks = 1;           // Number of column blocks
    uint32_t total_haplotypes = 0;       // Total number of haplotypes
    vector<pair<uint32_t, uint32_t>> col_block_ranges; // (start_haplotype, n_haplotypes) per column block
    vector<uint64_t> col_block_vec_lens; // vec_len for each column block
    bool useLegacyPath = true;           // Backward compatibility flag
    vector<uint32_t> chunk_block_offsets;
    uint64_t no_vec;
    uint64_t no_copy;
    uint32_t used_bits_cp;
    uint32_t bm_comp_cp_size;
    CBitMemory bm_comp_copy_orgl_id; //id of copied vector (where is original)


	vector<pair<std::string,uint32_t>> d_where_chrom;
	vector<int64_t> chunks_min_pos;
	map<int, chunk_stream> chunks_streams;
    // GT column tiling: 2D permutation map (row_block, col_block) -> perm
    map<pair<uint32_t, uint32_t>, vector<uint8_t>> vint_last_perm_2d;
	map<uint32_t,vector<uint8_t>> vint_last_perm;  // Legacy format (deprecated)

	size_t p_chrom;
	size_t p_pos;
	size_t p_id;
	size_t p_ref;
	size_t p_alt;
	size_t p_gt;
	size_t p_qual;
	uint32_t no_variants;

	uint32_t no_threads;

   	vector<uint8_t> all_v_header, comp_v_header;
	vector<uint8_t> all_v_samples, comp_v_samples;

	DecompressPartQueue<uint32_t> *decomp_part_queue; 
	vector<thread> part_decompress_thread;
	File_Handle_2 * file_handle2 = nullptr;
	string temp_file2_fname ;
	string fname;

	uint32_t no_keys;
	vector<key_desc> keys;
	int key_gt_id;
	vector<uint32_t> actual_variants;
	vector<CBuffer> v_i_buf;
	vector<SPackage*> v_packages;
	mutex m_packages;
	condition_variable cv_packages;

	// sort_field_block sort_field_block_compress;
	// sort_field_block fixed_field_block_io;

	fixed_field_block fixed_field_block_compress;
	fixed_field_block fixed_field_block_io;

    // New fixed-fields chunk format (row_block directory)
    bool has_fixed_fields_rb_dir = false;
    uint64_t fixed_fields_chunk_start = 0;
    uint32_t fixed_fields_chunk_version = 0;
    uint32_t fixed_fields_total_variants = 0;
    uint32_t fixed_fields_row_block_count = 0;
    uint32_t fixed_fields_gt_off = 0;
    uint32_t fixed_fields_gt_size = 0;
    std::vector<fixed_fields_rb_dir_entry> fixed_fields_rb_dir;

	vector<uint32_t> v_coder_part_ids;
    vector<std::unique_ptr<CompressionStrategy>> field_size_codecs;
    vector<std::unique_ptr<CompressionStrategy>> field_data_codecs;
    compression_backend_t backend = compression_backend_t::bsc;

	int64_t prev_pos;
	// int id_block = 0;

	
	

	
	void initDecoderParams();
	bool decompress_other_fileds(SPackage* pck);
	void Decoder(vector<uint8_t>& v_tmp, vector<uint8_t>& v_data);

public:

    DecompressionReader()
    {
		
        // no_variants = 0;
        prev_pos = 0;

#ifdef MMAP
        fm = nullptr;

#endif
    }

    DecompressionReader(const GSC_Params &params) : DecompressionReader() { backend = params.backend; }
    
    ~DecompressionReader()
    {
#ifdef MMAP
        if(fm)
            delete fm;

#endif
	
	for (auto p : v_packages)
		if (p)
			delete p;

	if (file_handle2)
		delete file_handle2;

	if(decomp_part_queue)
		delete decomp_part_queue;
    }

	bool OpenReading(const string &in_file_name, const bool &_decompression_mode_type);
	bool OpenReadingPart2(const string &in_file_name);
	void InitDecompressParams();
	void decompress_meta(vector<string> &v_samples, string &header);
	bool readFixedFields();

	bool Decoder(vector<block_t> &v_blocks,vector<vector<vector<uint32_t>>> &s_perm,vector<uint8_t> &gt_index,uint32_t cur_chunk_id);
    bool DecoderByRange(vector<block_t> &v_blocks, vector<vector<vector<uint32_t>>> &s_perm,
                        vector<uint8_t> &gt_index, uint32_t cur_chunk_id,
                        int64_t range_1, int64_t range_2, uint32_t &variants_before);
	bool setStartChunk(uint32_t start_chunk);
	uint32_t getActualPos(uint32_t chunk_id);

	void GetVariants(vector<field_desc> &fields);

	// bool writeOtherFields(FILE* f,CompOtherFields<int,uint8_t,uint8_t> * fields_queue);
	// bool readOtherFields(FILE* f,BlockingQueue<int,uint8_t,uint8_t> *fields_queue);

	

	// bool SetVariantOtherFields(vector<field_desc> &fields);

	// void SetNoKeys(int &_no_keys,vector<key_desc> &_keys,uint32_t &_other_fields_id);
	// void GetNoKeys(int &_no_keys,vector<key_desc> &_keys,size_t &_no_other_fields);
	// bool setBuffer(vector<vector<uint32_t>> &_v_size,vector<vector<uint8_t>> &_v_data);

	void close();

	void SetNoSamples(uint32_t _no_samples);

	void SetNoThreads(uint32_t _no_threads);


	bool GetHeader(string &_v_header);
	bool SetHeader(string &_v_header);

	bool SetSamples(vector<string> &_v_samples);
	bool GetSamples(vector<string> &_v_samples);
	

	void out_perm(vector<uint32_t> &perm,vector<variant_desc_t> &v_vcf_data_io);





	
};



