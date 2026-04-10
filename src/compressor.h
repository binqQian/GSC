#pragma once

#include <iostream>
#include <vector>
#include "gsc_params.h"
#include "defs.h"
#include <sdsl/bit_vectors.hpp>
#include "bit_memory.h"
#include "file_handle.h"
#include "compression_reader.h"
#include "variant.h"
#include "block_processing.h"
#include "queues.h"
#include <tuple>
#include "compression_strategy.h"
#include "fmt_compress/fmt_field_processor.h"
#include "field_stats.h"
// #include <filesystem>
using namespace std;

// 压缩总控：负责串联读取、GT 块处理、fixed fields 打包和最终写盘。
class Compressor
{
    // 当前任务的参数快照。
    GSC_Params params;

    // CompOtherFields<int,uint8_t,uint8_t> comp_other_fields;
    // CompVarBlockQueue<fixed_field_block> comp_var_block_queue;
    // 线程和中间队列：GT 处理与 fixed fields 压缩在这里衔接。
    CompVarBlockQueue<fixed_field_block> comp_sort_block_queue;
    vector<thread> part_compress_thread;
    vector<thread> block_process_thread;

    // 输出文件和临时文件句柄。
    FILE *comp = nullptr;
    bool is_stdout = false;
    FILE *temp_file = nullptr;
    uint64_t sdsl_offset;
    uint64_t other_fields_offset;
    bool mode_type;
    string fname ;
    char *temp_file1_fname = nullptr;
    string temp_file2_fname;
    File_Handle_2 * file_handle2 = nullptr;


    // zero/copy 位图和 rank 结构，用于解压端快速判断向量类型。
    sdsl::bit_vector zeros_only_bit_vector[2];
    sdsl::bit_vector copy_bit_vector[2];
    sdsl::bit_vector unique;
    sdsl::rank_support_v5<> rank_unique;
    sdsl::rank_support_v5<> rank_zeros_only_vector[2];
    sdsl::rank_support_v5<> rank_copy_bit_vector[2];


    uint64_t copy_no = 0;
    uint64_t  unique_no = 0;
    vector<uint32_t> comp_pos_copy;
    vector<bool> all_zeros;
    vector<bool> all_copies;
    uint32_t no_blocks = 0;
    // uint64_t start = 0;
    int64_t prev_pos = 0;
    // uint32_t sort_field_block_id = 0;
    uint32_t fixed_field_block_id = 0;
    // sort_field_block sort_field_block_io,sort_field_block_compress,comp_sort_field_block;
    // 每个 chunk 内先按 row_block 聚合 fixed fields，再统一压缩写盘。
    fixed_field_chunk fixed_chunk_io;
    fixed_field_block fixed_field_block_compress;


    
    bool end_of_processing = false;
    uint32_t no_curr_chrom_block =  0;
    vector<int64_t> chunks_min_pos;
    uint32_t cur_chunk_actual_pos = 0;
    // GT 列分块后的置换表；旧格式只保留单维映射。
    map<pair<uint32_t, uint32_t>, vector<uint8_t>> vint_last_perm_2d;
    map<uint32_t,vector<uint8_t>> vint_last_perm;  // Legacy format (deprecated)

    // GT 列分块元数据。
    uint32_t n_col_blocks = 1;
    vector<uint32_t> col_block_sizes;
    vector<uint64_t> col_block_vec_lens;
    uint32_t total_haplotypes = 0;
    uint32_t max_block_cols = 0;

    map<int, chunk_stream> chunks_streams;   

    CBitMemory bm_comp_copy_orgl_id;
    uint32_t used_bits_cp;
    vector<pair<std::string, uint32_t>> where_chrom;

    // tiled 模式下，先把同一个 row_block 的所有列块 GT 索引拼到一起，
    // 等这个 row_block 完整后再写入 fixed_chunk_io.gt_row_blocks。
    std::vector<uint8_t> pending_gt_row_block;
    uint32_t pending_gt_row_block_id = 0;
    bool has_pending_gt_row_block = false;

    uint64_t no_vec;
    size_t block_size;
    mutex mtx_gt_block;
	condition_variable cv_gt_block;
    int cur_block_id = 0;
    uint32_t cur_col_block_id = 0;
    bool use_legacy_perm = true;



    vector<uint8_t> all_v_header, comp_v_header;
	vector<uint8_t> all_v_samples, comp_v_samples;


    uint64_t Meta_comp_size = 0;
    uint64_t CHORM_comp_size = 0;
    uint64_t POS_comp_size = 0;
    uint64_t ID_comp_size = 0;
    uint64_t REF_comp_size = 0;
    uint64_t ALT_comp_size = 0;
    uint64_t QUAL_comp_size = 0;
    uint64_t GT_comp_size = 0;

    // 字段级统计和 codec 选择，主要服务于 lossless 其他字段压缩。
    CompressionStatistics comp_stats_;
    bool enable_field_stats_ = false;

    vector<std::unique_ptr<CompressionStrategy>> field_size_codecs;
    vector<std::unique_ptr<CompressionStrategy>> field_data_codecs;

    
    bool input_pos;
    size_t toal_all_size = 0;
    uint32_t no_keys;
    vector<key_desc> keys;
    int key_gt_id;

    // FORMAT 特殊字段的字典与处理器。
    std::unique_ptr<fmt_compress::FmtDictionaries> fmt_dictionaries_;
    std::unique_ptr<fmt_compress::FmtFieldProcessor> fmt_processor_;
    bool isFmtSpecialField(const std::string& name) const;

    bool OpenForWriting(const string &out_file_name);
    char bits_used(unsigned int n);
    void compressReplicatedRow();;
    bool writeCompressFlie();

    bool compress_other_fields_to_payloads(SPackage& pck, vector<uint8_t>& data_payload, vector<uint8_t>& size_payload, vector<uint8_t>& scratch);
    void compress_INT_fileds(SPackage& pck, vector<uint8_t>& v_compressed, vector<uint8_t>& v_tmp);
    void lock_gt_block_process(int &_block_id, uint32_t _col_block_id);
    bool check_gt_block_process(int &_block_id);
    void unlock_gt_block_process(uint32_t _col_block_id);
    void Encoder(vector<uint8_t>& v_data, vector<uint8_t>& v_tmp);
    bool compress_meta(vector<string> v_samples,const string& v_header);
	    void InitCompressParams();
	    bool compressFixedFields(fixed_field_block &field_block_io);
	    bool compressFixedFieldsChunk(fixed_field_chunk &chunk_io);
	    bool compressFixedFieldsChunkToPayload(const fixed_field_chunk &chunk_io, std::vector<uint8_t> &payload);
	    bool OpenTempFile(const string &out_file_name);
	    bool writeTempFlie(fixed_field_block &fixed_field_block_io);
	    bool writeTempChunkRB(const fixed_field_chunk &chunk_io,
	                          const std::vector<fixed_field_block> &row_blocks_comp,
	                          const std::vector<std::vector<uint8_t>> &gt_row_blocks_comp);
public:
    ~Compressor()
    {
	    if (file_handle2)
		    delete file_handle2;
        // if(fname)
        //     free(fname);
        if(temp_file1_fname)
            free(temp_file1_fname);
        // if(comp)
        //     fclose(comp);
           
    }

    Compressor(GSC_Params &_params)
    {
        params = _params;
        // curr_vec_id = 0;
        copy_no = 0;
        unique_no = 0;
        no_blocks = 0;
        input_pos = true;
        all_zeros.reserve(no_variants_in_buf);
        all_copies.reserve(no_variants_in_buf);
        comp_pos_copy.reserve(no_variants_in_buf);

        chunks_min_pos.reserve(no_variants_in_buf);
        
    }
    // 压缩主入口：组织读取、并行处理、临时文件写入和最终封包。
    bool CompressProcess();

    // Get compression statistics
    const CompressionStatistics& GetCompressionStats() const { return comp_stats_; }

private:
    // Log compression statistics at debug level
    void LogCompressionStats();
};


