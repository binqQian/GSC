#include "compressor.h"
#include "logger.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>

// **************************************************************************************************
// write the compressed file
bool Compressor::writeCompressFlie()
{
    auto logger = LogManager::Instance().Logger();
    const char *backend_name = "bsc";
    switch (params.backend)
    {
    case compression_backend_t::bsc:
        backend_name = "bsc";
        break;
    case compression_backend_t::zstd:
        backend_name = "zstd";
        break;
    case compression_backend_t::brotli:
        backend_name = "brotli";
        break;
    }
    logger->info("Using {} to compress the genotype part.", backend_name);
    // cout << "Gt_index_original_size:" << toal_all_size<< endl;

    uint32_t curr_no_blocks = 0;
    uint32_t where_chrom_size = static_cast<uint32_t>(where_chrom.size());
    uint32_t chunks_streams_size = static_cast<uint32_t>(chunks_streams.size());
    fwrite(&chunks_streams_size, sizeof(uint32_t), 1, comp);

    for (auto cur_chunk : chunks_streams)
    {
        fwrite(&cur_chunk.second.cur_chunk_actual_pos, sizeof(uint32_t), 1, comp);
        fwrite(&cur_chunk.second.offset, sizeof(size_t), 1, comp);
    }
    fwrite(&params.ploidy, sizeof(uint8_t), 1, comp);
    fwrite(&params.max_block_rows, sizeof(uint32_t), 1, comp);
    fwrite(&params.max_block_cols, sizeof(uint32_t), 1, comp);

    fwrite(&params.vec_len, sizeof(params.vec_len), 1, comp);

    fwrite(&no_vec, sizeof(no_vec), 1, comp);

    fwrite(&copy_no, sizeof(copy_no), 1, comp);
    fwrite(&used_bits_cp, sizeof(char), 1, comp);
    fwrite(&bm_comp_copy_orgl_id.mem_buffer_pos, sizeof(int32_t), 1, comp);
    fwrite(bm_comp_copy_orgl_id.mem_buffer, 1, bm_comp_copy_orgl_id.mem_buffer_pos, comp);
    // fwrite(&no_blocks, sizeof(no_blocks), 1, comp);
    // fwrite(&params.max_no_vec_in_block, sizeof(params.max_no_vec_in_block), 1, comp);
    fwrite(&params.n_samples, sizeof(params.n_samples), 1, comp);

    uint32_t chunks_min_pos_size = static_cast<uint32_t>(chunks_min_pos.size());
    fwrite(&chunks_min_pos_size, sizeof(uint32_t), 1, comp);
    fwrite(&chunks_min_pos[0], sizeof(int64_t), chunks_min_pos_size, comp);

    fwrite(&where_chrom_size, sizeof(uint32_t), 1, comp);
    // fwrite(&all_indexes_pos[0], all_indexes_pos.size() * sizeof(uint64_t), 1, comp);
    for (const auto &elem : where_chrom)
    {
        size_t chrom_size = elem.first.size();
        fwrite(&chrom_size, sizeof(size_t), 1, comp);
        fwrite(elem.first.data(), sizeof(char), chrom_size, comp);
        fwrite(&elem.second, sizeof(int), 1, comp);
    }

    // Write GT column tiling metadata
    fwrite(&n_col_blocks, sizeof(uint32_t), 1, comp);
    for (uint32_t cb = 0; cb < n_col_blocks; ++cb)
    {
        uint32_t start_hap = cb * max_block_cols;
        uint32_t block_size = col_block_sizes[cb];
        fwrite(&start_hap, sizeof(uint32_t), 1, comp);
        fwrite(&block_size, sizeof(uint32_t), 1, comp);
    }

    // Write permutations (legacy: chunk-based last block only, tiled: per row/col block)
    if (!use_legacy_perm)
    {
        // New format: 2D permutations
        uint32_t vint_last_perm_size = static_cast<uint32_t>(vint_last_perm_2d.size());
        fwrite(&vint_last_perm_size, sizeof(uint32_t), 1, comp);
        for (auto it = vint_last_perm_2d.begin(); it != vint_last_perm_2d.end(); ++it)
        {
            fwrite(&it->first.first, sizeof(uint32_t), 1, comp);  // row_block_id
            fwrite(&it->first.second, sizeof(uint32_t), 1, comp); // col_block_id
            uint32_t data_size = static_cast<uint32_t>(it->second.size());
            fwrite(&data_size, sizeof(uint32_t), 1, comp);
            fwrite(it->second.data(), sizeof(uint8_t), data_size, comp);
        }
    }
    else
    {
        // Legacy format: write row_block_id and data size only (no col_block_id)
        uint32_t vint_last_perm_size = static_cast<uint32_t>(vint_last_perm.size());
        fwrite(&vint_last_perm_size, sizeof(uint32_t), 1, comp);
        for (const auto &data : vint_last_perm)
        {
            fwrite(&data.first, sizeof(uint32_t), 1, comp); // row_block_id (chunk_id in legacy)
            uint32_t data_size = static_cast<uint32_t>(data.second.size());
            fwrite(&data_size, sizeof(uint32_t), 1, comp);
            fwrite(data.second.data(), sizeof(uint8_t), data_size, comp);
        }
    }

    bm_comp_copy_orgl_id.Close();
    Meta_comp_size += comp_v_header.size();
    uint32_t comp_size = static_cast<uint32_t>(comp_v_header.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    fwrite(comp_v_header.data(), 1, comp_v_header.size(), comp);

    Meta_comp_size += comp_v_samples.size();
    comp_size = static_cast<uint32_t>(comp_v_samples.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    fwrite(comp_v_samples.data(), 1, comp_v_samples.size(), comp);

    uint64_t size;
    uint8_t *temp_buffer = NULL;
    while (fread(&size, sizeof(size), 1, temp_file) == 1)
    {
        chunks_streams[curr_no_blocks + 1].offset = ftell(comp);

        temp_buffer = (uint8_t *)realloc(temp_buffer, size);
        if (!temp_buffer)
        {
            perror("Memory allocation failed");
            fclose(temp_file);
            return 1;
        }

        if (fread(temp_buffer, 1, size, temp_file) != size)
        {
            perror("Error reading data block");
            free(temp_buffer);
            fclose(temp_file);
            return 1;
        }
        fwrite(temp_buffer, 1, size, comp);
        curr_no_blocks++;
    }

    // while (comp_sort_block_queue.Pop(fixed_field_block_id,fixed_field_block_io))
    // {
    //     // std::cerr<<"fixed_fixed_field_block_id: "<<fixed_fixed_field_block_id<<":"<<sort_fixed_field_block_id<<endl;
    //     assert(fixed_field_block_id == curr_no_blocks);
    //     chunks_streams[curr_no_blocks+1].offset = ftell(comp);
    //     curr_no_blocks++;

    //     fwrite(&fixed_field_block_io.no_variants, sizeof(uint32_t), 1, comp);
    //     CHORM_comp_size += fixed_field_block_io.chrom.size();
    //     comp_size = static_cast<uint32_t>(fixed_field_block_io.chrom.size());
    //     fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    //     fwrite(fixed_field_block_io.chrom.data(), 1, fixed_field_block_io.chrom.size(), comp);

    //     POS_comp_size += fixed_field_block_io.pos.size();
    //     comp_size = static_cast<uint32_t>(fixed_field_block_io.pos.size());
    //     fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    //     fwrite(fixed_field_block_io.pos.data(), 1, fixed_field_block_io.pos.size(), comp);

    //     ID_comp_size += fixed_field_block_io.id.size();
    //     comp_size = static_cast<uint32_t>(fixed_field_block_io.id.size());
    //     fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    //     fwrite(fixed_field_block_io.id.data(), 1, fixed_field_block_io.id.size(), comp);

    //     REF_comp_size += fixed_field_block_io.ref.size();
    //     comp_size = static_cast<uint32_t>(fixed_field_block_io.ref.size());
    //     fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    //     fwrite(fixed_field_block_io.ref.data(), 1, fixed_field_block_io.ref.size(), comp);

    //     ALT_comp_size += fixed_field_block_io.alt.size();
    //     comp_size = static_cast<uint32_t>(fixed_field_block_io.alt.size());
    //     fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    //     fwrite(fixed_field_block_io.alt.data(), 1, fixed_field_block_io.alt.size(), comp);

    //     QUAL_comp_size += fixed_field_block_io.qual.size();
    //     comp_size = static_cast<uint32_t>(fixed_field_block_io.qual.size());
    //     fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    //     fwrite(fixed_field_block_io.qual.data(), 1, fixed_field_block_io.qual.size(), comp);

    //     GT_comp_size += fixed_field_block_io.gt_block.size();
    //     comp_size = static_cast<uint32_t>(fixed_field_block_io.gt_block.size());
    //     fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    //     fwrite(fixed_field_block_io.gt_block.data(), 1, fixed_field_block_io.gt_block.size(), comp);
    // }
    // std::cerr << "Meta_comp_size:    " << Meta_comp_size <<"\tByte"<<endl;
    // std::cerr << "CHORM_comp_size:   " << CHORM_comp_size <<"\tByte"<< endl;
    // std::cerr << "POS_comp_size:     " << POS_comp_size <<"\tByte"<< endl;
    // std::cerr << "ID_comp_size:      " << ID_comp_size <<"\tByte"<< endl;
    // std::cerr << "REF_comp_size:     " << REF_comp_size <<"\tByte"<< endl;
    // std::cerr << "ALT_comp_size:     " << ALT_comp_size <<"\tByte"<< endl;
    // std::cerr << "QUAL_comp_size:    " << QUAL_comp_size <<"\tByte"<< endl;
    // cout << "GT_index_comp_size:      " << GT_comp_size <<"\tByte"<< endl;
    free(temp_buffer);
    fclose(temp_file);
    if (remove(temp_file1_fname) != 0)
    {
        perror("Error deleting temp1 file");
        return EXIT_FAILURE;
    }
    // remove(temp_file1_fname);
    other_fields_offset = ftell(comp);

    FILE *other_f = fopen(temp_file2_fname.c_str(), "rb");

    if (other_f)
    {

        const size_t buffer_size = 1024;
        char buffer[buffer_size];
        size_t bytes_read;
        while ((bytes_read = fread(buffer, 1, buffer_size, other_f)) > 0)
        {
            fwrite(buffer, 1, bytes_read, comp);
        }
        fclose(other_f);
        other_f = nullptr;
        if (remove(temp_file2_fname.c_str()) != 0)
        {
            perror("Error deleting temp2 file");
            return EXIT_FAILURE;
        }

        mode_type = true;
    }

    sdsl_offset = ftell(comp);

    // sdsl_offset = end_offset;
    fseek(comp, 0, SEEK_SET);
    fwrite(&mode_type, sizeof(mode_type), 1, comp);
    fwrite(&other_fields_offset, sizeof(other_fields_offset), 1, comp);
    fwrite(&sdsl_offset, sizeof(sdsl_offset), 1, comp);
    logger->debug("mode_type: {} other_fields_offset: {} sdsl_offset: {}", static_cast<int>(mode_type), other_fields_offset, sdsl_offset);
    fseek(comp, 21, SEEK_SET);
    for (auto cur_chunk : chunks_streams)
    {
        fwrite(&cur_chunk.second.cur_chunk_actual_pos, sizeof(uint32_t), 1, comp);
        fwrite(&cur_chunk.second.offset, sizeof(size_t), 1, comp);
    }
    fseek(comp, 0, SEEK_END);

    if (comp && comp != stdout)
    {
        fclose(comp);
        comp = nullptr;
    }
    if (is_stdout)
    {

        for (int v = 0; v < 2; v++)
        {
            sdsl::rrr_vector<> rrr_bit_vector(zeros_only_bit_vector[v]);
            rrr_bit_vector.serialize(std::cout);
            sdsl::util::clear(rrr_bit_vector);
        }
        for (int v = 0; v < 2; v++)
        {
            sdsl::rrr_vector<> rrr_bit_vector(copy_bit_vector[v]);
            rrr_bit_vector.serialize(std::cout);
            sdsl::util::clear(rrr_bit_vector);
        }
    }
    else
    {

        sdsl::osfstream out(fname, std::ios::binary | std::ios::app);
        if (!out)
        {
            if (sdsl::util::verbose)
            {
                logger->error("ERROR: store_to_file not successful for: `{}`", fname);
            }
            exit(1);
        }
        sdsl::rrr_vector<> rrr_bit_vector[5];
        for (int v = 0; v < 2; v++)
        {
            sdsl::rrr_vector<> rrr_bit_vector(zeros_only_bit_vector[v]);
            rrr_bit_vector.serialize(out);
            sdsl::util::clear(rrr_bit_vector);
        }
        for (int v = 0; v < 2; v++)
        {
            rrr_bit_vector[v] = sdsl::rrr_vector<>(copy_bit_vector[v]);
            rrr_bit_vector[v].serialize(out);
            sdsl::util::clear(rrr_bit_vector[v]);
        }

        out.close();
    }

    if (sdsl::util::verbose)
    {
        logger->debug("store_to_file: `{}`", fname);
    }
    // sdsl::osfstream out(fname, std::ios::binary | std::ios::app);
    // if (!out)
    // {
    //     if (sdsl::util::verbose)
    //     {
    //     std::cerr << "ERROR: store_to_file not successful for: `" << fname << "`" << std::endl;
    //     }
    //     exit(1);
    // }
    // sdsl::rrr_vector<> rrr_bit_vector[5];

    // for (int v = 0; v < 2; v++)
    // {
    //     rrr_bit_vector[v] = sdsl::rrr_vector<>(zeros_only_bit_vector[v]);
    //     rrr_bit_vector[v].serialize(out);
    //     sdsl::util::clear(rrr_bit_vector[v]);
    // }
    // for (int v = 0; v < 2; v++)
    // {
    //     rrr_bit_vector[v + 2] = sdsl::rrr_vector<>(copy_bit_vector[v]);
    //     rrr_bit_vector[v + 2].serialize(out);
    //     sdsl::util::clear(rrr_bit_vector[v + 2]);
    // }

    // out.close();

    // if (sdsl::util::verbose)
    // {
    // std::cerr << "INFO: store_to_file: `" << fname << "`" << std::endl;
    // }

    // std::cerr << "genotype compress file (" << params.out_file_name + ".gsc"<< ") created." << std::endl;

    return 0;
}
// open file for writing
//  *******************************************************************************************************************
bool Compressor::OpenForWriting(const string &out_file_name)
{
    auto logger = LogManager::Instance().Logger();
    if (out_file_name != "-")
    {
        fname = out_file_name;
        // fname = (char *)malloc(strlen(out_file_name.c_str()) + 5);
        // snprintf(fname, strlen(out_file_name.c_str()) + 5, "%s.gsc", out_file_name.c_str());
        comp = fopen(fname.c_str(), "wb");
        if (!comp)
        {

            logger->error("storing archive not successful for: `{}`", fname);
            exit(1);
        }
    }
    else
    {

        comp = stdout;
        is_stdout = true;
    }
    if (setvbuf(comp, nullptr, _IOFBF, 64 << 20) != 0)
    {
        logger->error("Buffer setup failed for: `{}`", fname);
        // if (fname != nullptr) {
        //     free(fname);
        // }
        exit(1);
    }

    mode_type = false;

    fwrite(&mode_type, sizeof(mode_type), 1, comp);

    other_fields_offset = ftell(comp) + sizeof(uint64_t);
    fwrite(&other_fields_offset, sizeof(other_fields_offset), 1, comp);

    sdsl_offset = ftell(comp) + sizeof(uint64_t);

    if (fwrite(&sdsl_offset, sizeof(sdsl_offset), 1, comp) != 1)
    {
        logger->error("Write operation failed for: `{}`", fname);
    }

    // setvbuf(comp, nullptr, _IOFBF, 64 << 20);
    // sdsl_offset = ftell(comp) + sizeof(uint64_t);

    // fwrite(&sdsl_offset, sizeof(sdsl_offset), 1, comp);

    int id = (int)chunks_streams.size();

    chunks_streams[id] = chunk_stream(0, 0);

    return true;
}
//*******************************************************************************************************************
bool Compressor::OpenTempFile(const string &out_file_name)
{
    auto logger = LogManager::Instance().Logger();
    temp_file1_fname = (char *)malloc(strlen(out_file_name.c_str()) + 5);

    snprintf(temp_file1_fname, strlen(out_file_name.c_str()) + 5, "%s.temp", out_file_name.c_str());

    temp_file = fopen(temp_file1_fname, "wb");
    if (!temp_file)
    {

        logger->error("storing archive not successful for: `{}`", temp_file1_fname);
        exit(1);
    }

    setvbuf(temp_file, nullptr, _IOFBF, 64 << 20);

    if (params.compress_mode == compress_mode_t::lossless_mode)
    {
        if (file_handle2)
            delete file_handle2;
        file_handle2 = new File_Handle_2(false);

        temp_file2_fname = out_file_name + ".com_tmp_gsc";

        if (!file_handle2->Open(temp_file2_fname))
        {
            logger->error("Cannot open {}", temp_file2_fname);
            return false;
        }
    }
    return true;
}
// *******************************************************************************************************************
// compess pragma entry
bool Compressor::CompressProcess()
{
    // it is important to initialize the library before using it
    auto logger = LogManager::Instance().Logger();
    InitializeCompressionBackend(params.backend);

    enable_field_stats_ = false;
    if (const char* env = std::getenv("GSC_LOG_LEVEL"))
    {
        std::string val(env);
        std::transform(val.begin(), val.end(), val.begin(), ::tolower);
        enable_field_stats_ = (val == "debug");
    }
    if (enable_field_stats_)
        comp_stats_.Reset();

    MyBarrier my_barrier(3);

    unique_ptr<CompressionReader> compression_reader(new CompressionReader(params));

    if (!compression_reader->OpenForReading(params.in_file_name))
    {
        logger->error("Cannot open: {}", params.in_file_name);
        return false;
    }
    if (!OpenForWriting(params.out_file_name))
        return false;

    if (!OpenTempFile(params.out_file_name))
        return false;

    string header;

    vector<string> v_samples;

    compression_reader->GetHeader(header);

    uint32_t no_samples = compression_reader->GetSamples(v_samples);
    // for(int i =0;i<v_samples.size();i++)
    //     std::cerr<<v_samples[i]<<endl;
    logger->info("Number of samples: {}", no_samples);

    if (!no_samples)
    {
        logger->error("The number of genotype samples is zero and cannot be compressed!");
        return false;
    }

    // compress meta data
    compress_meta(v_samples, header);

    compression_reader->setNoVecBlock(params);

    // Get column block tiling information BEFORE creating threads
    n_col_blocks = compression_reader->GetNumColumnBlocks();
    col_block_sizes = compression_reader->GetColumnBlockSizes();
    col_block_vec_lens = compression_reader->GetColumnBlockVecLens();
    total_haplotypes = compression_reader->GetTotalHaplotypes();
    max_block_cols = compression_reader->GetMaxBlockCols();
    use_legacy_perm = (n_col_blocks == 1) &&
                      (params.max_block_cols == 0 || params.max_block_cols >= total_haplotypes) &&
                      (params.max_block_rows == 0 || params.max_block_rows == total_haplotypes);

    logger->info("Number of GT threads: {}", params.no_gt_threads);
    GtBlockQueue inGtBlockQueue(max((int)(params.no_blocks * params.no_gt_threads), 8));

    // VarBlockQueue<fixed_fixed_field_block> inVarBlockQueue(max((int)params.no_threads * 2, 8));
    VarBlockQueue<fixed_field_chunk> sortVarBlockQueue(max((int)params.no_threads * 2, 8));
    compression_reader->setQueue(&inGtBlockQueue);
    if (enable_field_stats_)
    {
        comp_stats_.total_samples = no_samples;
        comp_stats_.ploidy = params.ploidy;
    }

    PartQueue<SPackage> part_queue(max((int)params.no_threads * 2, 8));

    if (params.compress_mode == compress_mode_t::lossless_mode)
    {

        compression_reader->setPartQueue(&part_queue);

        compression_reader->InitVarinats(file_handle2);

        compression_reader->GetOtherField(keys, no_keys, key_gt_id);

        InitCompressParams();

        part_compress_thread.reserve(params.no_threads);

        for (uint32_t i = 0; i < params.no_threads; ++i)
        {
            part_compress_thread.emplace_back(thread([&]()
                                                     {
                                                         SPackage pck;
                                                         vector<uint8_t> v_compressed;
                                                         vector<uint8_t> v_tmp;

                                                         auto fo = [this](SPackage &pck) -> bool
                                                         { return check_coder_compressor(pck); };

                                                         while (true)
                                                         {

                                                             if (!part_queue.Pop<SPackage>(pck, fo))
                                                                 break;
                                                             compress_other_fileds(pck, v_compressed, v_tmp);
                                                         }
                                                     }));
        }
    }
    block_size = params.var_in_block * 2;

    unique_ptr<thread> compress_thread(new thread([&]
                                                  {
                                                      fixed_field_chunk fixed_field_chunk_process;
                                                      while (true)
                                                      {

                                                          if (!sortVarBlockQueue.Pop(fixed_field_block_id, fixed_field_chunk_process))
                                                          {
                                                              break;
                                                          }

                                                          compressFixedFieldsChunk(fixed_field_chunk_process);
                                                      }
                                                  }));

    // create multiple threads to handle individual blocks
    block_process_thread.reserve(params.no_gt_threads);
    string prev_chrom = "";
    int chunk_id = 0;
    for (uint32_t i = 0; i < params.no_gt_threads; ++i)
        block_process_thread.emplace_back(thread([&]()
                                                 {
                                                     auto logger = LogManager::Instance().Logger();

                                                     int block_id = 0;
                                                     uint32_t col_block_id = 0; // NEW: column block ID
                                                     unsigned long num_rows;
                                                     unsigned char *data = nullptr;
                                                     vector<variant_desc_t> v_vcf_data_io;
                                                     vector<uint32_t> origin_of_copy;
                                                     origin_of_copy.reserve(no_variants_in_buf);
                                                     vector<uint8_t> samples_indexes; // Index of the location where the block storing 1 is stored.
                                                     vector<uint32_t> perm;
                                                     BlockProcess block_process(params);

                                                     while (true)
                                                     {
                                                         if (!inGtBlockQueue.Pop(block_id, col_block_id, data, num_rows, v_vcf_data_io))
                                                         {
                                                             break;
                                                         }

                                                         vector<bool> zeros_only(num_rows, false);
                                                         vector<bool> copies(num_rows, false);
                                                         block_process.SetCurBlock(num_rows, data);
                                                         // Gets the sparse encoding for each block

                                                        if (num_rows)
                                                        {
                                                            // Calculate column block parameters from precomputed tiling
                                                            if (col_block_id >= col_block_sizes.size() || col_block_id >= col_block_vec_lens.size())
                                                            {
                                                                logger->error("Invalid column block id: {} (n_col_blocks={}, sizes={}, vec_lens={})",
                                                                              col_block_id, n_col_blocks, col_block_sizes.size(), col_block_vec_lens.size());
                                                                continue;
                                                            }
                                                            uint32_t col_block_size = col_block_sizes[col_block_id];
                                                            uint32_t col_vec_len = static_cast<uint32_t>(col_block_vec_lens[col_block_id]);

                                                            if (col_block_size == 0 || col_vec_len == 0)
                                                            {
                                                                logger->error("Empty column block detected: block_id={}, col_block_id={}, size={}, vec_len={}",
                                                                              block_id, col_block_id, col_block_size, col_vec_len);
                                                                continue;
                                                            }

                                                            logger->debug("Processing GT block: block_id={}, col_block_id={}, num_rows={}, col_block_size={}, col_vec_len={}",
                                                                         block_id, col_block_id, num_rows, col_block_size, col_vec_len);

                                                            // Prepare permutation buffer for this column block
                                                             perm.resize(col_block_size);
                                                             for (size_t i_p = 0; i_p < perm.size(); i_p++)
                                                                 perm[i_p] = static_cast<uint32_t>(i_p);

                                                             // if(num_rows == block_size)
                                                             block_process.ProcessSquareBlock(col_block_size, col_vec_len, perm, zeros_only, copies, origin_of_copy, samples_indexes, true);
                                                             // else
                                                             //     block_process.ProcessLastBlock(zeros_only, copies, origin_of_copy,samples_indexes);

                                                            if (use_legacy_perm && num_rows == block_size)
                                                                block_process.ProcessVariant(perm, v_vcf_data_io);

                                                             logger->debug("Finished GT block: block_id={}, col_block_id={}, perm_size={}", block_id, col_block_id, perm.size());
                                                         }

                                                         if (data != nullptr)
                                                             delete[] data;
                                                         // Gets the sparse encoding for each block
                                                         lock_gt_block_process(block_id, col_block_id);
                                                         {
                                                             if (num_rows)
                                                             {
                                                                 if (use_legacy_perm)
                                                                 {
                                                                     if (num_rows % block_size)
                                                                         vint_last_perm.emplace(chunk_id, vint_code::EncodeArray(perm));
                                                                 }
                                                                 else
                                                                 {
                                                                     vint_last_perm_2d.emplace(make_pair(block_id, col_block_id), vint_code::EncodeArray(perm));
                                                                 }
                                                                 if (v_vcf_data_io.empty())
                                                                 {
                                                                     if (!has_pending_gt_row_block)
                                                                     {
                                                                         pending_gt_row_block.clear();
                                                                         pending_gt_row_block_id = static_cast<uint32_t>(block_id);
                                                                         has_pending_gt_row_block = true;
                                                                     }
                                                                     block_process.AddGtIndexBlock(pending_gt_row_block, all_zeros, all_copies, comp_pos_copy,
                                                                                                  zeros_only, copies, origin_of_copy, samples_indexes);
                                                                 }
                                                             }
                                                             // if(!cur_block_id)
                                                             //     prev_chrom = v_vcf_data_io[0].chrom;
                                                             // if(num_rows != block_size)
                                                             //     block_process.ProcessLastPerm(perm,vint_last_perm);

                                                             // In Tiled mode, only the last column block has variant descriptors
                                                            // Skip variant descriptor processing for column blocks with empty v_vcf_data_io
                                                            if (!v_vcf_data_io.empty() && prev_chrom != v_vcf_data_io[0].chrom)
                                                             {

                                                                 prev_chrom = v_vcf_data_io[0].chrom;
                                                                 if (no_curr_chrom_block)
                                                                 {

                                                                     size_t gt_sz = 0;
                                                                     for (const auto &b : fixed_chunk_io.gt_row_blocks)
                                                                         gt_sz += b.size();
                                                                     toal_all_size += gt_sz;
                                                                     sortVarBlockQueue.Push(chunk_id, std::move(fixed_chunk_io));
                                                                     int id = (int)chunks_streams.size();
                                                                     chunks_streams[id] = chunk_stream(cur_chunk_actual_pos, 0);

                                                                     no_curr_chrom_block = 0;
                                                                     chunk_id++;
                                                                     fixed_chunk_io.Clear();
                                                                 }
                                                                 if (num_rows)
                                                                 {

                                                                     cur_chunk_actual_pos += (uint32_t)v_vcf_data_io.size();
                                                                     if (!has_pending_gt_row_block)
                                                                     {
                                                                         pending_gt_row_block.clear();
                                                                         pending_gt_row_block_id = static_cast<uint32_t>(block_id);
                                                                         has_pending_gt_row_block = true;
                                                                     }
                                                                     block_process.AddGtIndexBlock(pending_gt_row_block, all_zeros, all_copies, comp_pos_copy,
                                                                                                  zeros_only, copies, origin_of_copy, samples_indexes);
                                                                     fixed_field_block row_block_fixed;
                                                                     row_block_fixed.Clear();
                                                                     int64_t prev_pos_rb = 0;
                                                                     block_process.addFixedFieldsBlock(row_block_fixed, v_vcf_data_io, prev_pos_rb);
                                                                     fixed_fields_row_block_meta meta;
                                                                     meta.variant_count = row_block_fixed.no_variants;
                                                                     meta.first_pos = (int64_t)v_vcf_data_io.front().pos;
                                                                     meta.last_pos = (int64_t)v_vcf_data_io.back().pos;
                                                                     fixed_chunk_io.no_variants += meta.variant_count;
                                                                     fixed_chunk_io.row_blocks.emplace_back(std::move(row_block_fixed));
                                                                     fixed_chunk_io.row_meta.emplace_back(meta);
                                                                     fixed_chunk_io.gt_row_blocks.emplace_back(std::move(pending_gt_row_block));
                                                                     pending_gt_row_block.clear();
                                                                     has_pending_gt_row_block = false;
                                                                     no_curr_chrom_block++;

                                                                     if (no_curr_chrom_block == params.no_blocks)
                                                                     {
                                                                         size_t gt_sz = 0;
                                                                         for (const auto &b : fixed_chunk_io.gt_row_blocks)
                                                                             gt_sz += b.size();
                                                                         toal_all_size += gt_sz;
                                                                         sortVarBlockQueue.Push(chunk_id, std::move(fixed_chunk_io));
                                                                         int id = (int)chunks_streams.size();
                                                                         chunks_streams[id] = chunk_stream(cur_chunk_actual_pos, 0);
                                                                         no_curr_chrom_block = 0;
                                                                         chunk_id++;
                                                                         fixed_chunk_io.Clear();
                                                                     }
                                                                 }
                                                                 // Note: num_rows == 0 is a termination marker
                                                                 // sortVarBlockQueue.Complete() is now called by the main thread after all workers finish
                                                             }
                                                             else if (!v_vcf_data_io.empty())
                                                             {

                                                                 cur_chunk_actual_pos += (uint32_t)v_vcf_data_io.size();
                                                                 if (!has_pending_gt_row_block)
                                                                 {
                                                                     pending_gt_row_block.clear();
                                                                     pending_gt_row_block_id = static_cast<uint32_t>(block_id);
                                                                     has_pending_gt_row_block = true;
                                                                 }
                                                                 block_process.AddGtIndexBlock(pending_gt_row_block, all_zeros, all_copies, comp_pos_copy,
                                                                                              zeros_only, copies, origin_of_copy, samples_indexes);
                                                                 fixed_field_block row_block_fixed;
                                                                 row_block_fixed.Clear();
                                                                 int64_t prev_pos_rb = 0;
                                                                 block_process.addFixedFieldsBlock(row_block_fixed, v_vcf_data_io, prev_pos_rb);
                                                                 fixed_fields_row_block_meta meta;
                                                                 meta.variant_count = row_block_fixed.no_variants;
                                                                 meta.first_pos = (int64_t)v_vcf_data_io.front().pos;
                                                                 meta.last_pos = (int64_t)v_vcf_data_io.back().pos;
                                                                 fixed_chunk_io.no_variants += meta.variant_count;
                                                                 fixed_chunk_io.row_blocks.emplace_back(std::move(row_block_fixed));
                                                                 fixed_chunk_io.row_meta.emplace_back(meta);
                                                                 fixed_chunk_io.gt_row_blocks.emplace_back(std::move(pending_gt_row_block));
                                                                 pending_gt_row_block.clear();
                                                                 has_pending_gt_row_block = false;
                                                                 no_curr_chrom_block++;

                                                                 if (no_curr_chrom_block == params.no_blocks)
                                                                 {
                                                                     size_t gt_sz = 0;
                                                                     for (const auto &b : fixed_chunk_io.gt_row_blocks)
                                                                         gt_sz += b.size();
                                                                     toal_all_size += gt_sz;
                                                                     sortVarBlockQueue.Push(chunk_id, std::move(fixed_chunk_io));
                                                                     int id = (int)chunks_streams.size();
                                                                     chunks_streams[id] = chunk_stream(cur_chunk_actual_pos, 0);
                                                                     no_curr_chrom_block = 0;
                                                                     chunk_id++;
                                                                     fixed_chunk_io.Clear();
                                                                 }
                                                             }
                                                         }
                                                         unlock_gt_block_process(col_block_id);
                                                     }
                                                 }));

    if (!compression_reader->ProcessInVCF())
        return false;

    no_vec = compression_reader->getNoVec();

    // Wait for all GT processing threads to finish
    for (uint32_t i = 0; i < params.no_gt_threads; ++i)
        block_process_thread[i].join();

    // After all GT blocks are processed, push any remaining data and complete the queue
    {
        lock_guard<mutex> lock(mtx_gt_block);
        if (no_curr_chrom_block > 0)
        {
            // Push the last chunk if there's remaining data
            size_t gt_sz = 0;
            for (const auto &b : fixed_chunk_io.gt_row_blocks)
                gt_sz += b.size();
            toal_all_size += gt_sz;
            sortVarBlockQueue.Push(chunk_id, std::move(fixed_chunk_io));
            int id = (int)chunks_streams.size();
            chunks_streams[id] = chunk_stream(cur_chunk_actual_pos, 0);
            fixed_chunk_io.Clear();
        }
        // Always complete the queue to unblock the compress thread
        sortVarBlockQueue.Complete();
    }

    // obtain the variation description

    compression_reader->GetWhereChrom(where_chrom, chunks_min_pos);

    if (params.compress_mode == compress_mode_t::lossless_mode)
    {

        for (uint32_t i = 0; i < params.no_threads; ++i)
            part_compress_thread[i].join();

        auto stream_id = file_handle2->RegisterStream("part2_params");
        vector<uint8_t> v_desc;
        vector<uint32_t> actual_variants = compression_reader->GetActualVariants();
        compression_reader->UpdateKeys(keys);
        append(v_desc, static_cast<uint32_t>(actual_variants.size()));

        for (uint32_t i = 0; i < actual_variants.size(); ++i)
        {
            append(v_desc, actual_variants[i]);
        }
        append(v_desc, no_keys);
        append(v_desc, key_gt_id);
        for (uint32_t i = 0; i < no_keys; ++i)
        {
            // std::cerr<<"key_id:"<<keys[i].key_id<<":"<<keys[i].actual_field_id<<endl;
            append(v_desc, keys[i].key_id);
            append(v_desc, keys[i].actual_field_id);
            append(v_desc, static_cast<uint64_t>(keys[i].keys_type));
            append(v_desc, keys[i].type);
            append_str(v_desc, keys[i].name);  // Save field name for FMT compression
        }
        append(v_desc, static_cast<uint32_t>(params.backend));

	        // Serialize FORMAT dictionaries (AD/PL/PID) for lossless special codecs.
	        // Backward compatible: older files won't have this tail.
	        {
	            std::vector<uint8_t> dict_blob;
	            if (fmt_dictionaries_)
	                fmt_dictionaries_->serialize(dict_blob);
	            if (fmt_dictionaries_)
	            {
	                auto logger = LogManager::Instance().Logger();
	                logger->info("FMT dict: AD={} PL={} PID={} (blob {} bytes)",
	                             fmt_dictionaries_->getADCount(),
	                             fmt_dictionaries_->getPLCount(),
	                             fmt_dictionaries_->getPIDCount(),
	                             dict_blob.size());
	            }
	            append(v_desc, static_cast<uint32_t>(1)); // fmt_dict_version
	            append(v_desc, static_cast<uint32_t>(dict_blob.size()));
	            v_desc.insert(v_desc.end(), dict_blob.begin(), dict_blob.end());
	        }

        vector<uint8_t> v_desc_compressed;
        auto codec = MakeCompressionStrategy(params.backend, p_bsc_fixed_fields);
        codec->Compress(v_desc, v_desc_compressed);
        // prepend backend id for safe detection during decompression
        std::vector<uint8_t> v_desc_payload;
        v_desc_payload.reserve(v_desc_compressed.size() + 1);
        v_desc_payload.push_back(static_cast<uint8_t>(params.backend));
        v_desc_payload.insert(v_desc_payload.end(), v_desc_compressed.begin(), v_desc_compressed.end());
        file_handle2->AddParamsPart(stream_id, v_desc_payload);
        file_handle2->Close();
    }

    // GT threads already joined earlier, just wait for compress thread
    compress_thread->join();

    if (temp_file)
        fclose(temp_file);
    temp_file = fopen(temp_file1_fname, "rb");
    logger->info("Complete and process the chunking.");

    compressReplicatedRow();

    writeCompressFlie();

    // Set total variants for statistics and log compression statistics
    if (enable_field_stats_)
    {
        uint64_t total_actual_variants = 0;
        for (uint32_t v : compression_reader->GetActualVariants())
            total_actual_variants += v;
        comp_stats_.total_variants = total_actual_variants;
        LogCompressionStats();
    }

    return true;
}
// process the all_zeros and all_copies
//  **************************************************************************************************
void Compressor::compressReplicatedRow()
{
    zeros_only_bit_vector[0] = sdsl::bit_vector(no_vec / 2 + no_vec % 2, 0);
    zeros_only_bit_vector[1] = sdsl::bit_vector(no_vec / 2 + no_vec % 2, 0);
    copy_bit_vector[0] = sdsl::bit_vector(no_vec / 2 + no_vec % 2, 0);
    copy_bit_vector[1] = sdsl::bit_vector(no_vec / 2 + no_vec % 2, 0);

    // comp_pos_copy = new uint32_t[all_copy_num]();

    unique = sdsl::bit_vector(no_vec, 0);

    for (uint64_t i = 0; i < all_zeros.size(); i++)
    {
        zeros_only_bit_vector[i % 2][i / 2] = all_zeros[i];
        copy_bit_vector[i % 2][i / 2] = all_copies[i];
    }
    comp_pos_copy.shrink_to_fit();
    copy_no = comp_pos_copy.size();

    // std::cerr << "no_vec:" << no_vec << endl;
    // rank_copy_bit_vector[0] = sdsl::rank_support_v5<>(&copy_bit_vector[0]);
    // rank_copy_bit_vector[1] = sdsl::rank_support_v5<>(&copy_bit_vector[1]);
    // rank_zeros_only_vector[0] = sdsl::rank_support_v5<>(&zeros_only_bit_vector[0]);
    // rank_zeros_only_vector[1] = sdsl::rank_support_v5<>(&zeros_only_bit_vector[1]);

    uint32_t n_copies = 0;

    uint64_t i = 0, zeros_no = 0, uqq = 0;

    for (i = 0; i < no_vec; i++)
    {
        if (zeros_only_bit_vector[i % 2][i / 2])
        {
            zeros_no++;
            continue;
        }
        else if (!copy_bit_vector[i % 2][i / 2])
        {
            unique[i] = 1;
            uqq++;
        }
        else
            n_copies++;
    }
    rank_unique = sdsl::rank_support_v5<>(&unique);
    uint64_t curr_non_copy_vec_id = 0;
    uint64_t copy_id = 0, origin_unique_id;
    uint32_t max_diff_copy = 0;

    for (uint64_t i = 0; i < no_vec; i++)
    {
        if (zeros_only_bit_vector[i % 2][i / 2])
            continue;
        if (!copy_bit_vector[i % 2][i / 2])
        {
            curr_non_copy_vec_id++;

            continue;
        }
        // Here curr_non_copy_vec_id is a fake curr_non_copy_vec_id (it is ud of the next non_copy_vec_id)

        origin_unique_id = rank_unique(comp_pos_copy[copy_id]);

        // Store difference -1, to not waste one value
        comp_pos_copy[copy_id] = curr_non_copy_vec_id - origin_unique_id - 1;

        if (comp_pos_copy[copy_id] > max_diff_copy)
        {
            max_diff_copy = comp_pos_copy[copy_id];
        }
        copy_id++;
    }
    used_bits_cp = bits_used(max_diff_copy);
    bm_comp_copy_orgl_id.Create(copy_no * 4);

    for (i = 0; i < copy_no; i++)
    {
        // std::cerr<<comp_pos_copy[i]<<endl;
        bm_comp_copy_orgl_id.PutBits(comp_pos_copy[i], (int32_t)used_bits_cp);
    }
    bm_comp_copy_orgl_id.FlushPartialWordBuffer();

    sdsl::util::clear(rank_unique);
}

// get the bits of the number
//  **************************************************************************************************
char Compressor::bits_used(unsigned int n)
{
    char bits = 0;
    while (n)
    {
        n = n >> 1;
        bits++;
    }
    return bits;
}
// ************************************************************************************
void Compressor::lock_coder_compressor(SPackage &pck)
{
    unique_lock<mutex> lck(mtx_v_coder);
    cv_v_coder.wait(lck, [&, this]
                    {
		int sid = pck.key_id;

		return (int) v_coder_part_ids[sid] == pck.part_id; });
}

// ************************************************************************************
bool Compressor::check_coder_compressor(SPackage &pck)
{
    unique_lock<mutex> lck(mtx_v_coder);
    int sid = pck.key_id;

    return (int)v_coder_part_ids[sid] == pck.part_id;
}

// ************************************************************************************
void Compressor::unlock_coder_compressor(SPackage &pck)
{
    lock_guard<mutex> lck(mtx_v_coder);
    int sid = pck.key_id;

    ++v_coder_part_ids[sid];
    cv_v_coder.notify_all();
}
void Compressor::compress_other_fileds(SPackage &pck, vector<uint8_t> &v_compressed, vector<uint8_t> &v_tmp)
{

    lock_coder_compressor(pck);
    CompressionStrategy *cbsc_size = field_size_codecs[pck.key_id].get();

    // Special lossless compression for selected FORMAT fields (AD/PL/PGT/PID).
    // Note: SPackage stores per-variant payloads concatenated in v_data; v_size stores per-variant lengths.
    const std::string &field_name = keys[pck.key_id].name;
    const bool is_fmt_field = (keys[pck.key_id].keys_type == key_type_t::fmt);
    const bool stats_enabled = enable_field_stats_;
    const bool is_gt_field = is_fmt_field && (field_name == "GT");

    auto addFieldStats = [&](uint64_t raw_bytes, uint64_t compressed_bytes) {
        if (!stats_enabled || is_gt_field || (raw_bytes == 0 && compressed_bytes == 0))
            return;
        if (keys[pck.key_id].keys_type == key_type_t::fmt)
        {
            comp_stats_.AddFormatRaw(field_name, raw_bytes);
            comp_stats_.AddFormatCompressed(field_name, compressed_bytes);
        }
        else if (keys[pck.key_id].keys_type == key_type_t::info)
        {
            comp_stats_.info.AddRaw(raw_bytes);
            comp_stats_.info.AddCompressed(compressed_bytes);
        }
        else if (keys[pck.key_id].keys_type == key_type_t::flt)
        {
            comp_stats_.filter.AddRaw(raw_bytes);
            comp_stats_.filter.AddCompressed(compressed_bytes);
        }
    };
	    const bool is_sparse_field = (field_name == "PGT" || field_name == "PID");
	    const bool is_ad_field = (field_name == "AD");
	    const bool is_dp_field = (field_name == "DP");
	    const bool is_min_dp_field = (field_name == "MIN_DP");
	    const bool is_gq_field = (field_name == "GQ");
	    const bool is_pl_field = (field_name == "PL");
	    const bool is_pid_field = (field_name == "PID");
	    const uint32_t n_samples = params.n_samples;

    if (pck.v_data.size())
    {

        v_compressed.clear();
        // lzma2::lzma2_compress(pck.v_data, v_compressed, 10, 3);

        CompressionStrategy *cbsc = field_data_codecs[pck.key_id].get();
        uint64_t raw_data_bytes = 0;

	        if (is_fmt_field && (is_dp_field || is_min_dp_field || is_gq_field) && keys[pck.key_id].type == BCF_HT_INT && n_samples > 0)
	        {
	            // DP/MIN_DP/GQ use record-wise byte streams (each record starts with its own marker);
	            // keep bytes as-is (no transpose), and add a stream-level codec marker.
	            vector<uint8_t> to_compress;
	            to_compress.reserve(1 + pck.v_data.size());
	            to_compress.push_back(1); // codec_id: special
            to_compress.insert(to_compress.end(), pck.v_data.begin(), pck.v_data.end());
            raw_data_bytes = to_compress.size();
            cbsc->Compress(to_compress, v_compressed);
        }
	        else if (is_fmt_field && is_sparse_field && keys[pck.key_id].type == BCF_HT_STR && n_samples > 0)
	        {
	            bool ok = true;
	            vector<uint8_t> sparse_data;
	            sparse_data.reserve(pck.v_data.size() / 2);

	            uint8_t vint_buf[16];
	            size_t in_off = 0;

	            if (is_pid_field && fmt_dictionaries_)
	            {
		                // PID: sparse + dictionary; store only non-missing samples.
		                // Record encoding: [count][delta_pos][tag]...
		                // tag: 0 => raw cell bytes [len][bytes...], else dict_id = tag-1 (dict stores bytes without bcf_str_vector_end).
		                for (size_t rec_idx = 0; rec_idx < pck.v_size.size(); ++rec_idx)
		                {
	                    const uint32_t rec_bytes = pck.v_size[rec_idx];
	                    if (rec_bytes == 0)
	                        continue;
	                    if (rec_bytes % n_samples != 0)
	                    {
	                        ok = false;
	                        break;
	                    }
	                    const uint32_t stride = rec_bytes / n_samples;
	                    const uint8_t *rec_ptr = pck.v_data.data() + in_off;

	                    vector<uint32_t> positions;
	                    positions.reserve(n_samples / 4);
	                    vector<uint32_t> tags;
	                    tags.reserve(n_samples / 4);
	                    vector<pair<uint32_t, vector<uint8_t>>> raws;
	                    raws.reserve(n_samples / 32);

	                    for (uint32_t s = 0; s < n_samples; ++s)
	                    {
	                        const uint8_t *cell = rec_ptr + (size_t)s * stride;
	                        const bool is_missing =
	                            (stride > 0 && (cell[0] == bcf_str_missing ||
	                                            (cell[0] == '.' && (stride == 1 || cell[1] == bcf_str_vector_end))));
	                        if (is_missing)
	                            continue;

		                        uint32_t term_pos = 0;
		                        while (term_pos < stride && cell[term_pos] != bcf_str_vector_end)
		                            ++term_pos;
		                        const bool has_term = (term_pos < stride);
		                        const uint32_t raw_len = has_term ? (term_pos + 1) : stride;
		                        if (!has_term || term_pos > 510)
		                        {
		                            // No terminator in stride, or too long for dictionary item: store raw.
		                            vector<uint8_t> raw(cell, cell + raw_len);
		                            positions.push_back(s);
		                            tags.push_back(0);
		                            raws.emplace_back(s, std::move(raw));
		                            continue;
		                        }

		                        fmt_compress::PIDItem item;
		                        memset(item.data_, 0, sizeof(item.data_));
		                        const uint16_t ulen = static_cast<uint16_t>(term_pos);
		                        memcpy(item.data_, &ulen, 2);
		                        if (term_pos)
		                            memcpy(item.data_ + 2, cell, term_pos);
		                        uint32_t id = fmt_dictionaries_->getPIDItemId(item);
		                        positions.push_back(s);
		                        tags.push_back(id + 1);
		                    }

	                    uint8_t vlen = fmt_compress::VintCodec::encode((uint64_t)positions.size(), vint_buf);
	                    sparse_data.insert(sparse_data.end(), vint_buf, vint_buf + vlen);

	                    uint32_t prev = 0;
	                    size_t raw_idx = 0;
	                    for (size_t i = 0; i < positions.size(); ++i)
	                    {
	                        const uint32_t pos = positions[i];
	                        const uint32_t delta = pos - prev;
	                        prev = pos;
	                        vlen = fmt_compress::VintCodec::encode(delta, vint_buf);
	                        sparse_data.insert(sparse_data.end(), vint_buf, vint_buf + vlen);
	                        const uint32_t tag = tags[i];
	                        vlen = fmt_compress::VintCodec::encode(tag, vint_buf);
	                        sparse_data.insert(sparse_data.end(), vint_buf, vint_buf + vlen);

	                        if (tag == 0)
	                        {
	                            // Find raw payload for this pos.
	                            while (raw_idx < raws.size() && raws[raw_idx].first != pos)
	                                ++raw_idx;
	                            if (raw_idx >= raws.size())
	                            {
	                                ok = false;
	                                break;
	                            }
	                            const auto &raw = raws[raw_idx].second;
	                            vlen = fmt_compress::VintCodec::encode((uint64_t)raw.size(), vint_buf);
	                            sparse_data.insert(sparse_data.end(), vint_buf, vint_buf + vlen);
	                            sparse_data.insert(sparse_data.end(), raw.begin(), raw.end());
	                            ++raw_idx;
	                        }
	                    }
	                    if (!ok)
	                        break;

	                    in_off += rec_bytes;
	                }
	            }
	            else
	            {
	                // PGT (and fallback PID): sparse storage per variant, store full stride bytes.
	                // Each record encoding: [count][pos][stride-bytes]...
	                for (size_t rec_idx = 0; rec_idx < pck.v_size.size(); ++rec_idx)
	                {
	                    const uint32_t rec_bytes = pck.v_size[rec_idx];
	                    if (rec_bytes == 0)
	                        continue;
	                    if (rec_bytes % n_samples != 0)
	                    {
	                        ok = false;
	                        break;
	                    }
	                    const uint32_t stride = rec_bytes / n_samples;
	                    const uint8_t *rec_ptr = pck.v_data.data() + in_off;

	                    vector<uint32_t> positions;
	                    positions.reserve(n_samples / 4);
	                    for (uint32_t s = 0; s < n_samples; ++s)
	                    {
	                        const uint8_t *cell = rec_ptr + (size_t)s * stride;
	                        const bool is_missing =
	                            (stride > 0 && (cell[0] == bcf_str_missing ||
	                                            (cell[0] == '.' && (stride == 1 || cell[1] == bcf_str_vector_end))));
	                        if (!is_missing)
	                            positions.push_back(s);
	                    }

	                    uint8_t vlen = fmt_compress::VintCodec::encode((uint64_t)positions.size(), vint_buf);
	                    sparse_data.insert(sparse_data.end(), vint_buf, vint_buf + vlen);
	                    for (uint32_t pos : positions)
	                    {
	                        vlen = fmt_compress::VintCodec::encode(pos, vint_buf);
	                        sparse_data.insert(sparse_data.end(), vint_buf, vint_buf + vlen);
	                        const uint8_t *cell = rec_ptr + (size_t)pos * stride;
	                        sparse_data.insert(sparse_data.end(), cell, cell + stride);
	                    }
	                    in_off += rec_bytes;
	                }
	            }

	            vector<uint8_t> to_compress;
	            if (ok && in_off == pck.v_data.size())
	            {
	                to_compress.reserve(1 + sparse_data.size());
	                to_compress.push_back(is_pid_field ? 2 : 1); // codec_id: 2=PID dict, 1=legacy sparse
	                to_compress.insert(to_compress.end(), sparse_data.begin(), sparse_data.end());
	            }
	            else
	            {
	                to_compress.reserve(1 + pck.v_data.size());
	                to_compress.push_back(0); // codec_id: legacy raw bytes
	                to_compress.insert(to_compress.end(), pck.v_data.begin(), pck.v_data.end());
	            }
	            raw_data_bytes = to_compress.size();
	            cbsc->Compress(to_compress, v_compressed);
	        }
		        else if (is_fmt_field && is_ad_field && keys[pck.key_id].type == BCF_HT_INT && n_samples > 0)
		        {
		            // AD compression per variant: 2-bit tip encoding per sample.
		            // Record encoding:
		            //   codec_id=1: [tip-bytes][values...]
		            //   codec_id=2: [tip-bytes][values...], and for Tip=11 each sample starts with tag (0=raw, else dict_id+1)
		            //   codec_id=3: baseline+delta for AD[0], plus optional dict tags for Tip=11.
		            //   codec_id=4: baseline+delta for sum(AD) (variant-count blocks), plus optional dict tags for Tip=11.
		            //
		            // codec_id=4 stream layout (after codec_id):
		            //   [block_size_vint]
		            //   for each variant record:
		            //     if (rec_idx % block_size == 0): [baseline_sum[sample]=sum(AD) vint] * n_samples
		            //     [tip-bytes][payload_len_vint][payload...]
		            //
		            // Tips (2-bit) for codec_id=4:
		            //   00: all zeros
		            //   01: ref-only (AD[0]>0, AD[1..]=0) => store zigzag(delta_sum) vs baseline_sum, then update baseline_sum
		            //   10: alt-only (biallelic, AD[0]=0, AD[1]>0) => store zigzag(delta_sum), then update baseline_sum
		            //   11: general case => store tag then:
		            //         tag=0: raw values (u32 bits via vint)
		            //         tag=1: sum+alt encoding (biallelic, both non-zero): [delta_sum][alt]
		            //         tag>=2: dict_id = tag-2
		            vector<uint8_t> ad_compressed;
		            uint8_t vint_buf[16];
		            const uint8_t codec_id = 4;
		            uint32_t block_size = 16384;
		            if (const char *bs = getenv("GSC_AD_BASELINE_BLOCK"))
		            {
	                const long v = strtol(bs, nullptr, 10);
	                if (v >= 64 && v <= 1'000'000)
	                    block_size = (uint32_t)v;
	            }

	            bool ok = true;
	            ad_compressed.reserve(pck.v_data.size() / 2);

	            auto append_vint_u64 = [&](uint64_t v, std::vector<uint8_t>& out) {
	                uint8_t vlen = fmt_compress::VintCodec::encode(v, vint_buf);
	                out.insert(out.end(), vint_buf, vint_buf + vlen);
	            };
	            auto append_vint_u32 = [&](uint32_t u, std::vector<uint8_t>& out) {
	                append_vint_u64((uint64_t)u, out);
	            };
	            auto append_zigzag_i64 = [&](int64_t diff, std::vector<uint8_t>& out) {
	                const uint64_t ud = (uint64_t)diff;
	                const uint64_t zz = (ud << 1) ^ (uint64_t)(diff >> 63);
	                append_vint_u64(zz, out);
	            };
	            auto append_vint_u32_bits = [&](int32_t v, std::vector<uint8_t>& out) {
	                uint32_t u = 0;
	                memcpy(&u, &v, sizeof(u));
	                append_vint_u64((uint64_t)u, out);
	            };
		            auto try_build_ad_item = [&](const int32_t* sv, uint32_t cnt, fmt_compress::ADItem& item) -> bool {
		                if (cnt == 0)
		                    return false;
		                // ADItem payload capacity: 46 bytes (bytes 2..47).
	                uint32_t max_v = 0;
	                for (uint32_t j = 0; j < cnt; ++j)
	                {
	                    uint32_t u = 0;
	                    memcpy(&u, &sv[j], sizeof(u));
	                    if (u > max_v) max_v = u;
	                }
	                uint8_t type = 0;
	                if (max_v < 0xFFu) type = 0;
	                else if (max_v < 0xFFFFu) type = 1;
	                else type = 2;

	                const size_t bytes_per = (type == 0) ? 1u : (type == 1) ? 2u : 4u;
	                if (2u + (size_t)cnt * bytes_per > sizeof(item.data_))
	                    return false;

	                memset(item.data_, 0, sizeof(item.data_));
	                item.data_[1] = type;
	                uint8_t* ptr = item.data_ + 2;
	                for (uint32_t j = 0; j < cnt; ++j)
	                {
	                    uint32_t u = 0;
	                    memcpy(&u, &sv[j], sizeof(u));
	                    if (type == 0)
	                    {
	                        *ptr++ = (uint8_t)u;
	                    }
	                    else if (type == 1)
	                    {
	                        uint16_t w = (uint16_t)u;
	                        memcpy(ptr, &w, sizeof(w));
	                        ptr += 2;
	                    }
	                    else
	                    {
	                        memcpy(ptr, &u, sizeof(u));
	                        ptr += 4;
	                    }
	                }
	                item.data_[0] = (uint8_t)(ptr - item.data_ - 1);
	                return true;
		            };

		            // Stream header: block size.
		            append_vint_u32(block_size, ad_compressed);
		            std::vector<uint32_t> baseline_sum(n_samples, 0);

		            size_t in_off = 0;
		            for (size_t rec_idx = 0; rec_idx < pck.v_size.size(); ++rec_idx)
		            {
	                const uint32_t total_elems = pck.v_size[rec_idx];
	                if (total_elems == 0)
	                    continue;
	                if (total_elems % n_samples != 0)
	                {
	                    ok = false;
	                    break;
	                }
		                const uint32_t per_sample = total_elems / n_samples;
		                const int32_t *vals = reinterpret_cast<const int32_t *>(pck.v_data.data() + in_off);

		                // Start of block: write baselines (sum(AD)) for all samples.
		                if (block_size > 0 && (rec_idx % block_size) == 0)
		                {
		                    for (uint32_t s = 0; s < n_samples; ++s)
		                    {
		                        const int32_t *sv = vals + (size_t)s * per_sample;
		                        uint64_t sum = 0;
		                        bool ok_sum = true;
		                        for (uint32_t j = 0; j < per_sample; ++j)
		                        {
		                            const int32_t v = sv[j];
		                            if (v == bcf_int32_missing || v == bcf_int32_vector_end || v < 0)
		                            {
		                                ok_sum = false;
		                                break;
		                            }
		                            sum += (uint64_t)(uint32_t)v;
		                        }
		                        const uint32_t b = ok_sum ? (uint32_t)std::min<uint64_t>(sum, (uint64_t)UINT32_MAX) : 0u;
		                        baseline_sum[s] = b;
		                        append_vint_u32(b, ad_compressed);
		                    }
		                }

		                vector<uint8_t> tips;
	                tips.reserve((size_t)n_samples * 2);
	                vector<uint8_t> payload;
	                payload.reserve((size_t)total_elems);

		                for (uint32_t s = 0; s < n_samples; ++s)
		                {
		                    const int32_t *sv = vals + (size_t)s * per_sample;
		                    bool has_special = false;
	                    for (uint32_t j = 0; j < per_sample; ++j)
	                    {
	                        const int32_t v = sv[j];
	                        if (v == bcf_int32_missing || v == bcf_int32_vector_end || v < 0)
	                        {
	                            has_special = true;
	                            break;
	                        }
	                    }

		                    if (has_special || per_sample == 0)
		                    {
		                        tips.push_back(1);
		                        tips.push_back(1);
		                        // tag=0 => raw values follow
		                        append_vint_u64(0, payload);
		                        for (uint32_t j = 0; j < per_sample; ++j)
		                            append_vint_u32_bits(sv[j], payload);
		                        // Keep baseline_sum[s] unchanged (missing/special values should not pollute the baseline).
		                    }
		                    else
		                    {
		                        // Sum(AD) baseline+delta with biallelic fast paths.
		                        uint64_t sum64 = 0;
		                        for (uint32_t j = 0; j < per_sample; ++j)
		                            sum64 += (uint64_t)(uint32_t)sv[j];
		                        const uint32_t sum = (uint32_t)std::min<uint64_t>(sum64, (uint64_t)UINT32_MAX);

		                        bool all_alts_zero = true;
		                        for (uint32_t j = 1; j < per_sample; ++j)
		                        {
		                            if (sv[j] != 0)
		                            {
		                                all_alts_zero = false;
		                                break;
		                            }
		                        }
		                        const bool all_zero = (sum == 0);
		                        const bool ref_only = (!all_zero && all_alts_zero);
		                        const bool alt_only = (per_sample == 2 && sv[0] == 0 && sv[1] != 0);
		                        const bool biallelic_both = (per_sample == 2 && sv[0] != 0 && sv[1] != 0);

		                        if (all_zero)
		                        {
		                            tips.push_back(0);
		                            tips.push_back(0);
		                            baseline_sum[s] = 0;
		                        }
		                        else if (ref_only)
		                        {
		                            tips.push_back(0);
		                            tips.push_back(1);
		                            const int64_t diff = (int64_t)sum - (int64_t)baseline_sum[s];
		                            append_zigzag_i64(diff, payload);
		                            baseline_sum[s] = sum;
		                        }
		                        else if (alt_only)
		                        {
		                            tips.push_back(1);
		                            tips.push_back(0);
		                            const int64_t diff = (int64_t)sum - (int64_t)baseline_sum[s];
		                            append_zigzag_i64(diff, payload);
		                            baseline_sum[s] = sum;
		                        }
		                        else if (biallelic_both)
		                        {
		                            tips.push_back(1);
		                            tips.push_back(1);
		                            // tag=1 => sum+alt encoding
		                            append_vint_u64(1, payload);
		                            const int64_t diff = (int64_t)sum - (int64_t)baseline_sum[s];
		                            append_zigzag_i64(diff, payload);
		                            append_vint_u32((uint32_t)sv[1], payload);
		                            baseline_sum[s] = sum;
		                        }
		                        else
		                        {
		                            tips.push_back(1);
		                            tips.push_back(1);
		                            bool stored_dict = false;
		                            uint32_t dict_id = 0;
		                            if (fmt_dictionaries_)
		                            {
		                                fmt_compress::ADItem item;
		                                if (try_build_ad_item(sv, per_sample, item))
		                                {
		                                    dict_id = fmt_dictionaries_->getADItemId(item);
		                                    stored_dict = true;
		                                }
		                            }
		                            // tag=0 raw, tag=1 reserved for sum+alt, tag>=2 => dict_id=tag-2
		                            append_vint_u64(stored_dict ? ((uint64_t)dict_id + 2u) : 0u, payload);
		                            if (!stored_dict)
		                            {
		                                for (uint32_t j = 0; j < per_sample; ++j)
		                                    append_vint_u32_bits(sv[j], payload);
		                            }
		                            baseline_sum[s] = sum;
		                        }
		                    }
		                }

	                const size_t tip_bytes = ((size_t)n_samples * 2 + 7) / 8;
	                const size_t tip_start = ad_compressed.size();
	                ad_compressed.resize(tip_start + tip_bytes, 0);
	                for (size_t i = 0; i < tips.size(); ++i)
	                {
	                    if (!tips[i])
	                        continue;
	                    ad_compressed[tip_start + (i / 8)] |= (uint8_t)(1u << (i % 8));
	                }

	                // codec_id=3: write payload length for robust decoding.
	                append_vint_u32((uint32_t)payload.size(), ad_compressed);
	                ad_compressed.insert(ad_compressed.end(), payload.begin(), payload.end());

	                in_off += (size_t)total_elems * 4;
	            }

	            vector<uint8_t> to_compress;
	            if (ok && in_off == pck.v_data.size())
	            {
	                to_compress.reserve(1 + ad_compressed.size());
	                to_compress.push_back(codec_id); // codec_id: AD special
	                to_compress.insert(to_compress.end(), ad_compressed.begin(), ad_compressed.end());
	            }
	            else
	            {
	                // Fallback to legacy payload (transposed bytes) with an explicit codec marker.
	                Encoder(pck.v_data, v_tmp);
                to_compress.reserve(1 + v_tmp.size());
                to_compress.push_back(0); // codec_id: legacy
                to_compress.insert(to_compress.end(), v_tmp.begin(), v_tmp.end());
            }
            raw_data_bytes = to_compress.size();
            cbsc->Compress(to_compress, v_compressed);
        }
		        else if (is_fmt_field && is_pl_field && keys[pck.key_id].type == BCF_HT_INT && n_samples > 0)
		        {
		            // PL compression per variant: 2-bit tip encoding + per-sample payload.
		            // Record encoding:
		            //   codec_id=1: [tip-bytes][payload...]
		            //   codec_id=2: [tip-bytes][payload...], and for Tip=11 each sample starts with tag (0=raw, else dict_id+1)
		            //   codec_id=3: like codec_id=2, plus per-record optional permutation data (for PL[0]!=0 normalization on len==3)
		            //   codec_id=4: block baseline for a + delta(a) and residual(b-11*a) for biallelic len==3, plus perm sidecar.
		            // Type 00: all zeros
		            // Type 01: pattern 1 (store a,b)
		            // Type 10: pattern 2 (store a; b=15*a)
		            // Type 11: store full vector
		            vector<uint8_t> pl_compressed;
		            uint8_t vint_buf[16];
		            const uint8_t codec_id = 4;
		            uint32_t block_size = 16384;
		            if (const char *bs = getenv("GSC_PL_BLOCK_SIZE"))
		            {
		                const long v = strtol(bs, nullptr, 10);
		                if (v >= 64 && v <= 1'000'000)
		                    block_size = (uint32_t)v;
		            }

		            bool ok = true;
		            pl_compressed.reserve(pck.v_data.size() / 2);

		            auto append_vint_u64 = [&](uint64_t v, std::vector<uint8_t>& out) {
		                uint8_t vlen = fmt_compress::VintCodec::encode(v, vint_buf);
		                out.insert(out.end(), vint_buf, vint_buf + vlen);
		            };
		            auto append_vint_u32 = [&](uint32_t u, std::vector<uint8_t>& out) {
		                append_vint_u64((uint64_t)u, out);
		            };
		            auto append_vint_u32_bits = [&](int32_t v, std::vector<uint8_t>& out) {
		                uint32_t u = 0;
		                memcpy(&u, &v, sizeof(u));
		                append_vint_u64((uint64_t)u, out);
		            };
		            auto append_zigzag_i64 = [&](int64_t diff, std::vector<uint8_t>& out) {
		                const uint64_t ud = (uint64_t)diff;
		                const uint64_t zz = (ud << 1) ^ (uint64_t)(diff >> 63);
		                append_vint_u64(zz, out);
		            };
		            auto abs_i64 = [](int64_t v) -> uint64_t {
		                return (uint64_t)(v < 0 ? -v : v);
		            };
		            auto set_bit = [](std::vector<uint8_t>& buf, size_t bit) {
		                buf[bit / 8] |= (uint8_t)(1u << (bit % 8));
		            };
	            auto pack_2bit_values = [&](const std::vector<uint8_t>& vals, std::vector<uint8_t>& out) {
	                // vals contain 2-bit values (0..3). Pack in order, little-endian within a byte.
	                const size_t bits = vals.size() * 2;
	                const size_t bytes = (bits + 7) / 8;
	                out.assign(bytes, 0);
	                for (size_t i = 0; i < vals.size(); ++i)
	                {
	                    const uint8_t v = vals[i] & 3u;
	                    const size_t b0 = i * 2;
	                    if (v & 1u) out[b0 / 8] |= (uint8_t)(1u << (b0 % 8));
	                    if (v & 2u) out[(b0 + 1) / 8] |= (uint8_t)(1u << ((b0 + 1) % 8));
	                }
	            };
	            auto try_build_pl_item = [&](const std::vector<uint32_t>& pl_vec, fmt_compress::PLItem& item) -> bool {
	                if (pl_vec.empty())
	                    return false;
	                uint32_t max_v = 0;
	                for (uint32_t u : pl_vec)
	                    if (u != fmt_compress::kMissingValue && u > max_v) max_v = u;
	                uint8_t type = 0;
	                if (max_v < 0xFFu) type = 0;
	                else if (max_v < 0xFFFFu) type = 1;
	                else type = 2;
	                const size_t bytes_per = (type == 0) ? 1u : (type == 1) ? 2u : 4u;
	                if (2u + pl_vec.size() * bytes_per > sizeof(item.data_))
	                    return false;
	                memset(item.data_, 0, sizeof(item.data_));
	                item.data_[1] = type;
	                uint8_t* ptr = item.data_ + 2;
	                for (uint32_t u : pl_vec)
	                {
	                    if (type == 0)
	                    {
	                        *ptr++ = (uint8_t)u;
	                    }
	                    else if (type == 1)
	                    {
	                        uint16_t w = (uint16_t)u;
	                        memcpy(ptr, &w, sizeof(w));
	                        ptr += 2;
	                    }
	                    else
	                    {
	                        memcpy(ptr, &u, sizeof(u));
	                        ptr += 4;
	                    }
	                }
		                item.data_[0] = (uint8_t)(ptr - item.data_ - 1);
		                return true;
		            };

		            // Stream header: block size.
		            append_vint_u32(block_size, pl_compressed);
		            std::vector<uint32_t> baseline_a(n_samples, 0);

            size_t in_off = 0;
            for (size_t rec_idx = 0; rec_idx < pck.v_size.size(); ++rec_idx)
            {
                const uint32_t total_elems = pck.v_size[rec_idx];
                if (total_elems == 0)
                    continue;
                if (total_elems % n_samples != 0)
                {
                    ok = false;
                    break;
                }
                const uint32_t per_sample = total_elems / n_samples;
                const int32_t *vals = reinterpret_cast<const int32_t *>(pck.v_data.data() + in_off);

	                // Start of block: write baselines (a) for all samples.
	                if (block_size > 0 && (rec_idx % block_size) == 0)
	                {
	                    for (uint32_t s = 0; s < n_samples; ++s)
	                    {
	                        uint32_t a = 0;
	                        if (per_sample == 3)
	                        {
	                            const int32_t *sv = vals + (size_t)s * per_sample;
	                            bool has_special = false;
	                            for (uint32_t j = 0; j < 3; ++j)
	                            {
	                                const int32_t v = sv[j];
	                                if (v == bcf_int32_missing || v == bcf_int32_vector_end || v < 0)
	                                {
	                                    has_special = true;
	                                    break;
	                                }
	                            }
	                            if (!has_special)
	                            {
	                                const uint32_t p0 = (uint32_t)sv[0];
	                                const uint32_t p1 = (uint32_t)sv[1];
	                                const uint32_t p2 = (uint32_t)sv[2];
	                                uint32_t min_v = p0;
	                                uint32_t min_pos = 0;
	                                if (p1 < min_v) { min_v = p1; min_pos = 1; }
	                                if (p2 < min_v) { min_v = p2; min_pos = 2; }
	                                if (min_v == 0 && min_pos == 1)
	                                    a = p0; // [a,0,b] -> [0,a,b]
	                                else
	                                    a = p1; // [0,a,b] or [b,a,0]
	                            }
	                        }
	                        baseline_a[s] = a;
	                        append_vint_u32(a, pl_compressed);
	                    }
	                }

	                vector<uint8_t> tips;
	                tips.reserve((size_t)n_samples * 2);
	                std::vector<uint8_t> perm_present;
	                std::vector<uint8_t> perm_vals;
	                const size_t perm_bytes = ((size_t)n_samples + 7) / 8;
	                perm_present.assign(perm_bytes, 0);
	                perm_vals.reserve(n_samples / 16);
	                vector<uint8_t> payload;
	                payload.reserve((size_t)total_elems);

	                for (uint32_t s = 0; s < n_samples; ++s)
	                {
	                    const int32_t *sv = vals + (size_t)s * per_sample;
	                    bool has_special = false;
	                    vector<uint32_t> pl_vec;
	                    pl_vec.resize(per_sample);
	                    for (uint32_t j = 0; j < per_sample; ++j)
	                    {
	                        const int32_t v = sv[j];
	                        if (v == bcf_int32_missing || v == bcf_int32_vector_end || v < 0)
	                        {
	                            has_special = true;
	                            break;
	                        }
	                        pl_vec[j] = static_cast<uint32_t>(v);
	                    }

	                    if (has_special || per_sample == 0)
	                    {
	                        tips.push_back(1);
	                        tips.push_back(1);
	                        append_vint_u64(0, payload); // tag=0 => raw values (u32 bits)
	                        for (uint32_t j = 0; j < per_sample; ++j)
	                            append_vint_u32_bits(sv[j], payload);
	                        continue;
	                    }

	                    // Optional PL[0]!=0 normalization for biallelic diploid (len==3): if min is 0 but at pos!=0,
	                    // permute to [0,a,b] so it can be encoded by the Type1/Type2 pattern model.
	                    const uint32_t* pl_used = pl_vec.data();
	                    uint32_t tmp_pl[3];
	                    bool used_perm = false;
	                    uint8_t perm_val = 0;
	                    if (per_sample == 3 && pl_vec[0] != 0)
	                    {
	                        uint32_t min_v = pl_vec[0];
	                        uint32_t min_pos = 0;
	                        for (uint32_t j = 1; j < 3; ++j)
	                        {
	                            if (pl_vec[j] < min_v)
	                            {
	                                min_v = pl_vec[j];
	                                min_pos = j;
	                            }
	                        }
	                        if (min_v == 0 && min_pos > 0)
	                        {
	                            tmp_pl[0] = 0;
	                            if (min_pos == 1)
	                            {
	                                tmp_pl[1] = pl_vec[0];
	                                tmp_pl[2] = pl_vec[2];
	                            }
	                            else
	                            {
	                                tmp_pl[1] = pl_vec[1];
	                                tmp_pl[2] = pl_vec[0];
	                            }
	                            uint32_t a_tmp = 0, b_tmp = 0;
	                            const uint8_t t_tmp = fmt_compress::checkPlPattern(tmp_pl, 3, a_tmp, b_tmp);
	                            if (t_tmp != 3)
	                            {
	                                pl_used = tmp_pl;
	                                used_perm = true;
	                                perm_val = (uint8_t)min_pos;
	                            }
	                        }
	                    }

	                    uint32_t a_val = 0, b_val = 0;
	                    const uint8_t type = fmt_compress::checkPlPattern(pl_used, per_sample, a_val, b_val);
	                    if (type == 0)
	                    {
	                        tips.push_back(0);
	                        tips.push_back(0);
	                        if (per_sample == 3)
	                            baseline_a[s] = 0;
	                    }
	                    else if (type == 1)
	                    {
	                        if (used_perm)
	                        {
	                            set_bit(perm_present, s);
	                            perm_vals.push_back(perm_val);
	                        }
	                        if (per_sample == 3)
	                        {
	                            const int64_t delta_a = (int64_t)a_val - (int64_t)baseline_a[s];
	                            const int64_t residual = (int64_t)b_val - 11ll * (int64_t)a_val;
	                            if (abs_i64(delta_a) <= 63 && abs_i64(residual) <= 63)
	                            {
	                                tips.push_back(1);
	                                tips.push_back(0);
	                                append_zigzag_i64(delta_a, payload);
	                                append_zigzag_i64(residual, payload);
	                            }
	                            else
	                            {
	                                tips.push_back(1);
	                                tips.push_back(1);
	                                append_vint_u64(1, payload); // tag=1 => pattern1 (a,b)
	                                append_vint_u32(a_val, payload);
	                                append_vint_u32(b_val, payload);
	                            }
	                            baseline_a[s] = a_val;
	                        }
	                        else
	                        {
	                            tips.push_back(1);
	                            tips.push_back(1);
	                            append_vint_u64(1, payload); // tag=1 => pattern1 (a,b)
	                            append_vint_u32(a_val, payload);
	                            append_vint_u32(b_val, payload);
	                            baseline_a[s] = a_val;
	                        }
	                    }
	                    else if (type == 2)
	                    {
	                        if (used_perm)
	                        {
	                            set_bit(perm_present, s);
	                            perm_vals.push_back(perm_val);
	                        }
	                        tips.push_back(0);
	                        tips.push_back(1);
	                        const int64_t delta_a = (int64_t)a_val - (int64_t)baseline_a[s];
	                        append_zigzag_i64(delta_a, payload);
	                        baseline_a[s] = a_val;
	                    }
	                    else
	                    {
	                        tips.push_back(1);
	                        tips.push_back(1);
	                        bool stored_dict = false;
	                        uint32_t dict_id = 0;
	                        if (fmt_dictionaries_)
	                        {
	                            fmt_compress::PLItem item;
	                            if (try_build_pl_item(pl_vec, item))
	                            {
	                                dict_id = fmt_dictionaries_->getPLItemId(item);
	                                stored_dict = true;
	                            }
	                        }
	                        // tag=0 raw (vint values), tag=1 reserved for pattern1 (a,b), tag>=2 => dict_id=tag-2
	                        append_vint_u64(stored_dict ? ((uint64_t)dict_id + 2u) : 0u, payload);
	                        if (!stored_dict)
	                        {
	                            for (uint32_t j = 0; j < per_sample; ++j)
	                                append_vint_u64((uint64_t)pl_vec[j], payload);
	                        }
	                        if (per_sample == 3)
	                        {
	                            uint32_t a = pl_vec[1];
	                            uint32_t min_v = pl_vec[0];
	                            uint32_t min_pos = 0;
	                            for (uint32_t j = 1; j < 3; ++j)
	                                if (pl_vec[j] < min_v) { min_v = pl_vec[j]; min_pos = j; }
	                            if (min_v == 0 && min_pos == 1) a = pl_vec[0];
	                            baseline_a[s] = a;
	                        }
	                    }
	                }

                const size_t tip_bytes = ((size_t)n_samples * 2 + 7) / 8;
                const size_t tip_start = pl_compressed.size();
                pl_compressed.resize(tip_start + tip_bytes, 0);
                for (size_t i = 0; i < tips.size(); ++i)
                {
                    if (!tips[i])
                        continue;
                    pl_compressed[tip_start + (i / 8)] |= (uint8_t)(1u << (i % 8));
                }

	                // codec_id=4: optional permutation section for this record.
	                bool has_perm = false;
	                for (size_t i = 0; i < perm_present.size(); ++i)
	                {
	                    if (perm_present[i] != 0)
	                    {
	                        has_perm = true;
	                        break;
	                    }
	                }
	                pl_compressed.push_back(has_perm ? 1u : 0u);
	                if (has_perm)
	                {
	                    pl_compressed.insert(pl_compressed.end(), perm_present.begin(), perm_present.end());
	                    std::vector<uint8_t> packed_perm;
	                    pack_2bit_values(perm_vals, packed_perm);
	                    pl_compressed.insert(pl_compressed.end(), packed_perm.begin(), packed_perm.end());
	                }

	                // codec_id=4: write payload length for robust decoding.
	                append_vint_u32((uint32_t)payload.size(), pl_compressed);
	                pl_compressed.insert(pl_compressed.end(), payload.begin(), payload.end());

	                in_off += (size_t)total_elems * 4;
	            }

	            vector<uint8_t> to_compress;
	            if (ok && in_off == pck.v_data.size())
	            {
	                to_compress.reserve(1 + pl_compressed.size());
	                to_compress.push_back(codec_id); // codec_id: PL special
	                to_compress.insert(to_compress.end(), pl_compressed.begin(), pl_compressed.end());
	            }
	            else
	            {
                // Fallback to legacy payload (transposed bytes) with an explicit codec marker.
                Encoder(pck.v_data, v_tmp);
                to_compress.reserve(1 + v_tmp.size());
                to_compress.push_back(0); // codec_id: legacy
                to_compress.insert(to_compress.end(), v_tmp.begin(), v_tmp.end());
            }
            raw_data_bytes = to_compress.size();
            cbsc->Compress(to_compress, v_compressed);
        }
        else if (keys[pck.key_id].type != BCF_HT_INT)
        {
            raw_data_bytes = pck.v_data.size();
            cbsc->Compress(pck.v_data, v_compressed);
            // zstd::zstd_compress(pck.v_data, v_compressed);
        }
        else
        {
            Encoder(pck.v_data, v_tmp);
            raw_data_bytes = v_tmp.size();
            cbsc->Compress(v_tmp, v_compressed);
            // zstd::zstd_compress(v_tmp, v_compressed);
        }

        // cbsc->Compress(pck.v_data, v_compressed);

        file_handle2->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed);
        addFieldStats(raw_data_bytes, v_compressed.size());
    }
    else
    {
        v_compressed.clear();

        file_handle2->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed);
    }

    v_compressed.clear();

    // lzma2::lzma2_compress(pck.v_size, v_compressed, 10, 1);
    v_tmp.resize(pck.v_size.size() * 4);

    copy_n((uint8_t *)pck.v_size.data(), v_tmp.size(), v_tmp.data());

    cbsc_size->Compress(v_tmp, v_compressed);
    // zstd::zstd_compress(v_tmp, v_compressed);

    file_handle2->AddPartComplete(pck.stream_id_size, pck.part_id, v_compressed);
    addFieldStats(v_tmp.size(), v_compressed.size());

    unlock_coder_compressor(pck);
}
// void Compressor::compress_INT_fileds(SPackage& pck, vector<uint8_t>& v_compressed, vector<uint8_t>& v_tmp)
// {

//     vector<uint32_t> vec_temp1;
//     vector<uint32_t> vec_temp2;
//     lock_coder_compressor(pck);
// 	if (pck.v_data.size())
// 	{
// 		vec_temp1.resize(pck.v_data.size()/4);
//         v_compressed.clear();
//         for (size_t i = 0; i < vec_temp1.size(); i++) {
//             vec_temp1[i] = ((uint32_t)pck.v_data[i*4] ) |
//                         ((uint32_t)pck.v_data[i*4+1] << 8) |
//                         ((uint32_t)pck.v_data[i*4+2] << 16) |
//                         ((uint32_t)pck.v_data[i*4+3]<< 24);

//         }

//         // lzma2::lzma2_compress(pck.v_data, v_compressed, 10, 3);

//         FastPForCompress::FastPFor_Compress(vec_temp1, vec_temp2);

//         v_compressed.resize(vec_temp2.size()*4);

//         copy_n((uint8_t*)vec_temp2.data(),vec_temp2.size()*4,v_compressed.data());
// 		file_handle2->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed);
// 	}
// 	else
// 	{
// 		v_compressed.clear();

// 		file_handle2->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed);
// 	}

//     v_compressed.clear();
//     lzma2::lzma2_compress(pck.v_size, v_compressed, 10, 1);
//     // FastPForCompress::FastPFor_Compress(pck.v_size, vec_temp2);
//     // v_compressed.resize(vec_temp2.size()*4);
// 	// copy_n((uint8_t*)vec_temp2.data(),vec_temp2.size()*4,v_compressed.data());

// 	file_handle2->AddPartComplete(pck.stream_id_size, pck.part_id, v_compressed);

// 	unlock_coder_compressor(pck);
// }

void Compressor::Encoder(vector<uint8_t> &v_data, vector<uint8_t> &v_tmp)
{
    v_tmp.resize(v_data.size());
    size_t size = v_data.size() / 4;
    for (size_t i = 0; i < v_data.size() / 4; i++)
    {
        v_tmp[i] = v_data[i * 4];
        v_tmp[i + size] = v_data[i * 4 + 1];
        v_tmp[i + size * 2] = v_data[i * 4 + 2];
        v_tmp[i + size * 3] = v_data[i * 4 + 3];
    }
    // v_data = tmp;
    // uint8_t count = 1;
    // for (size_t i = 1; i < v_data.size(); ++i) {
    //     if(v_data[i] == v_data[i-1] && count < 255)
    //     {
    //         count++;
    //     }
    //     else
    //     {

    //         v_tmp.emplace_back(count);
    //         v_tmp.emplace_back(v_data[i-1]);
    //         count = 1;
    //     }

    // }
    // v_tmp.emplace_back(count);
    // v_tmp.emplace_back(v_data[v_data.size() - 1]);

    // uint8_t count = 0;
    // for (size_t i = 0; i < v_data.size(); ++i) {
    //     if (v_data[i] == 0) {
    //         count++;
    //         if(count == 255)
    //         {
    //             v_tmp.emplace_back(0);
    //             v_tmp.emplace_back(count);
    //             count = 0;
    //         }
    //     }
    //     else {
    //         if(count)
    //         {
    //             v_tmp.emplace_back(0);
    //             v_tmp.emplace_back(count);
    //             count = 0;
    //         }
    //         v_tmp.emplace_back(v_data[i]);
    //     }
    // }
    // if(count)
    // {
    //     v_tmp.emplace_back(0);
    //     v_tmp.emplace_back(count);
    // }

    v_tmp.shrink_to_fit();
    // for (size_t i = 0; i < v_data.size(); i++)
    // {
    //     std::cerr<<(int)v_data[i]<<" ";
    // }
    // std::cerr<<endl;
    // std::cerr<<endl;
    // for(size_t i = 0; i < v_tmp.size(); i++)
    // {
    //     std::cerr<<(int)v_tmp[i]<<" ";
    // }
    // std::cerr<<endl;
    // std::cerr<<endl;
}
// compressor meta data
//  *******************************************************************************************************************************************
bool Compressor::compress_meta(vector<string> v_samples, const string &v_header)
{

    append_str(all_v_header, v_header);

    for (auto &x : v_samples)
        append_str(all_v_samples, x);

    for (auto data : {
             make_tuple(ref(all_v_header), ref(comp_v_header), "header"),
             make_tuple(ref(all_v_samples), ref(comp_v_samples), "samples"),
         })
    {
        auto codec = MakeCompressionStrategy(params.backend, p_bsc_meta);
        codec->Compress(get<0>(data), get<1>(data));
        if (enable_field_stats_)
        {
            comp_stats_.meta.AddRaw(get<0>(data).size());
            comp_stats_.meta.AddCompressed(get<1>(data).size());
        }

        // zstd::zstd_compress(get<0>(data), get<1>(data));

        // LZMACompress::Compress(get<0>(data), get<1>(data), 9);

        // lz4:: lz4_compress(get<0>(data), get<1>(data), 12);
        // BrotliUtils::compressData(get<0>(data), get<1>(data));
        // std::cerr<<get<0>(data).size()<<":"<<get<1>(data).size()<<endl;
        // lzma2::lzma2_compress(get<0>(data), get<1>(data), get<2>(data), 10);
        // std::cerr << get<2>(data) << " size: " << get<1>(data).size() << endl;

        // fh.WriteUInt(get<1>(data).size(), 4, f);
        // fh.Write(get<1>(data).data(), get<1>(data).size(), f);
    }

    return true;
}
// init bsc compress params
//******************************************************************************************************************************************
void Compressor::InitCompressParams()
{

    v_coder_part_ids.resize(no_keys, 0);
    field_data_codecs.resize(no_keys);
    field_size_codecs.resize(no_keys);

    for (uint32_t i = 0; i < no_keys; ++i)
    {
        field_size_codecs[i] = MakeCompressionStrategy(params.backend, p_bsc_size);
        bsc_params_t param = p_bsc_text;

        switch (keys[i].type)
        {
        case BCF_HT_FLAG:
            param = p_bsc_flag;
            break;
        case BCF_HT_INT:
            param = p_bsc_int;
            break;
        case BCF_HT_REAL:
            param = p_bsc_real;
            break;
        case BCF_HT_STR:
            param = p_bsc_text;
            break;
        }
        field_data_codecs[i] = MakeCompressionStrategy(params.backend, param);
    }

    // Initialize FMT field compression
    fmt_dictionaries_ = std::make_unique<fmt_compress::FmtDictionaries>();
    fmt_processor_ = std::make_unique<fmt_compress::FmtFieldProcessor>(fmt_dictionaries_.get());
}

bool Compressor::isFmtSpecialField(const std::string& name) const {
    return name == "AD" || name == "DP" || name == "GQ" ||
           name == "PL" || name == "PGT" || name == "PID";
}

void Compressor::lock_gt_block_process(int &_block_id, uint32_t _col_block_id)
{
    unique_lock<mutex> lck(mtx_gt_block);
    cv_gt_block.wait(lck, [&, this]
                     {
        if (n_col_blocks <= 1) {
            return cur_block_id == _block_id;
        }
        return cur_block_id == _block_id && cur_col_block_id == _col_block_id; });
}

// ************************************************************************************
bool Compressor::check_gt_block_process(int &_block_id)
{
    unique_lock<mutex> lck(mtx_gt_block);

    return (int)cur_block_id == _block_id;
}

// ************************************************************************************
void Compressor::unlock_gt_block_process(uint32_t _col_block_id)
{
    lock_guard<mutex> lck(mtx_gt_block);
    (void)_col_block_id; // suppress unused warning; ordering handled via cur_col_block_id
    if (n_col_blocks <= 1)
    {
        ++cur_block_id;
    }
    else
    {
        // Advance column block cursor; roll over to next row block when last column is processed
        ++cur_col_block_id;
        if (cur_col_block_id >= n_col_blocks)
        {
            cur_col_block_id = 0;
            ++cur_block_id;
        }
    }
    cv_gt_block.notify_all();
}
bool Compressor::compressFixedFields(fixed_field_block &fixed_field_block_io)
{

    fixed_field_block_compress.no_variants = fixed_field_block_io.no_variants;
    auto codec = MakeCompressionStrategy(params.backend, p_bsc_fixed_fields);
    for (auto data : {
             make_tuple(ref(fixed_field_block_io.chrom), ref(fixed_field_block_compress.chrom), "chrom"),
             make_tuple(ref(fixed_field_block_io.id), ref(fixed_field_block_compress.id), "id"),
             make_tuple(ref(fixed_field_block_io.alt), ref(fixed_field_block_compress.alt), "alt"),
             make_tuple(ref(fixed_field_block_io.qual), ref(fixed_field_block_compress.qual), "qual"),
             make_tuple(ref(fixed_field_block_io.pos), ref(fixed_field_block_compress.pos), "pos"),
             make_tuple(ref(fixed_field_block_io.ref), ref(fixed_field_block_compress.ref), "ref"),
             // make_tuple(ref(fixed_field_block_io.gt_block), ref(fixed_field_block_compress.gt_block), "GT"),
         })
    {

        // BSC

        codec->Compress(get<0>(data), get<1>(data));

        // ZSTD

        // zstd::zstd_compress(get<0>(data), get<1>(data));

        // LZMA

        // LZMACompress::Compress(get<0>(data), get<1>(data), 9);

        // lz4:: lz4_compress(get<0>(data), get<1>(data), 12);
        // BrotliUtils::compressData(get<0>(data), get<1>(data));

        // std::cerr<<get<0>(data).size()<<":"<<get<1>(data).size()<<endl;
    }
    auto gt_codec = MakeCompressionStrategy(params.backend, p_bsc_fixed_fields);
    gt_codec->Compress(fixed_field_block_io.gt_block, fixed_field_block_compress.gt_block);
    uint8_t marker = 0;
    switch (params.backend)
    {
    case compression_backend_t::bsc:
        marker = 0;
        break;
    case compression_backend_t::zstd:
        marker = 1;
        break;
    case compression_backend_t::brotli:
        marker = 2;
        break;
    default:
        marker = 0;
        break;
    }
    fixed_field_block_compress.gt_block.emplace_back(marker);
    writeTempFlie(fixed_field_block_compress);
    // comp_sort_block_queue.Push(fixed_field_block_id,fixed_field_block_compress);
    fixed_field_block_compress.Clear();
    return true;
}

bool Compressor::compressFixedFieldsChunk(fixed_field_chunk &chunk_io)
{
    if (chunk_io.row_blocks.size() != chunk_io.row_meta.size())
        return false;
    if (chunk_io.row_blocks.size() != chunk_io.gt_row_blocks.size())
        return false;

    const uint32_t row_block_count = static_cast<uint32_t>(chunk_io.row_blocks.size());
    std::vector<fixed_field_block> row_blocks_comp(row_block_count);

    auto codec = MakeCompressionStrategy(params.backend, p_bsc_fixed_fields);
    for (uint32_t i = 0; i < row_block_count; ++i)
    {
        const fixed_field_block &rb = chunk_io.row_blocks[i];
        fixed_field_block &out = row_blocks_comp[i];
        out.no_variants = rb.no_variants;

        if (enable_field_stats_)
        {
            comp_stats_.chrom.AddRaw(rb.chrom.size());
            comp_stats_.pos.AddRaw(rb.pos.size());
            comp_stats_.id.AddRaw(rb.id.size());
            comp_stats_.ref.AddRaw(rb.ref.size());
            comp_stats_.alt.AddRaw(rb.alt.size());
            comp_stats_.qual.AddRaw(rb.qual.size());
        }

        codec->Compress(rb.chrom, out.chrom);
        codec->Compress(rb.pos, out.pos);
        codec->Compress(rb.id, out.id);
        codec->Compress(rb.ref, out.ref);
        codec->Compress(rb.alt, out.alt);
        codec->Compress(rb.qual, out.qual);

        // Collect compressed field sizes for statistics
        if (enable_field_stats_)
        {
            comp_stats_.chrom.AddCompressed(out.chrom.size());
            comp_stats_.pos.AddCompressed(out.pos.size());
            comp_stats_.id.AddCompressed(out.id.size());
            comp_stats_.ref.AddCompressed(out.ref.size());
            comp_stats_.alt.AddCompressed(out.alt.size());
            comp_stats_.qual.AddCompressed(out.qual.size());
        }
    }

    uint8_t marker = 0;
    {
        switch (params.backend)
        {
        case compression_backend_t::bsc:
            marker = 0;
            break;
        case compression_backend_t::zstd:
            marker = 1;
            break;
        case compression_backend_t::brotli:
            marker = 2;
            break;
        default:
            marker = 0;
            break;
        }
    }

    std::vector<std::vector<uint8_t>> gt_row_blocks_comp(row_block_count);
    auto gt_codec = MakeCompressionStrategy(params.backend, p_bsc_fixed_fields);
    for (uint32_t i = 0; i < row_block_count; ++i)
    {
        if (enable_field_stats_)
            comp_stats_.gt.AddRaw(chunk_io.gt_row_blocks[i].size());

        gt_codec->Compress(chunk_io.gt_row_blocks[i], gt_row_blocks_comp[i]);
        gt_row_blocks_comp[i].emplace_back(marker);

        // Collect GT compressed size for statistics
        if (enable_field_stats_)
            comp_stats_.gt.AddCompressed(gt_row_blocks_comp[i].size());
    }

    return writeTempChunkRB(chunk_io, row_blocks_comp, gt_row_blocks_comp);
}

bool Compressor::writeTempFlie(fixed_field_block &fixed_field_block_io)
{

    uint64_t offset = ftell(temp_file) + sizeof(uint64_t);

    uint64_t start_offset = ftell(temp_file);
    size_t comp_size = 0;
    fwrite(&offset, sizeof(offset), 1, temp_file);

    fwrite(&fixed_field_block_io.no_variants, sizeof(uint32_t), 1, temp_file);
    comp_size = static_cast<uint32_t>(fixed_field_block_io.chrom.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, temp_file);
    fwrite(fixed_field_block_io.chrom.data(), 1, fixed_field_block_io.chrom.size(), temp_file);

    comp_size = static_cast<uint32_t>(fixed_field_block_io.pos.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, temp_file);
    fwrite(fixed_field_block_io.pos.data(), 1, fixed_field_block_io.pos.size(), temp_file);

    comp_size = static_cast<uint32_t>(fixed_field_block_io.id.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, temp_file);
    fwrite(fixed_field_block_io.id.data(), 1, fixed_field_block_io.id.size(), temp_file);

    comp_size = static_cast<uint32_t>(fixed_field_block_io.ref.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, temp_file);
    fwrite(fixed_field_block_io.ref.data(), 1, fixed_field_block_io.ref.size(), temp_file);

    comp_size = static_cast<uint32_t>(fixed_field_block_io.alt.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, temp_file);
    fwrite(fixed_field_block_io.alt.data(), 1, fixed_field_block_io.alt.size(), temp_file);

    comp_size = static_cast<uint32_t>(fixed_field_block_io.qual.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, temp_file);
    fwrite(fixed_field_block_io.qual.data(), 1, fixed_field_block_io.qual.size(), temp_file);

    comp_size = static_cast<uint32_t>(fixed_field_block_io.gt_block.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, temp_file);
    fwrite(fixed_field_block_io.gt_block.data(), 1, fixed_field_block_io.gt_block.size(), temp_file);

    offset = ftell(temp_file) - offset;

    fseek(temp_file, start_offset, SEEK_SET);

    fwrite(&offset, sizeof(offset), 1, temp_file);
    fseek(temp_file, 0, SEEK_END);

    return true;
}

bool Compressor::writeTempChunkRB(const fixed_field_chunk &chunk_io,
                                  const std::vector<fixed_field_block> &row_blocks_comp,
                                  const std::vector<std::vector<uint8_t>> &gt_row_blocks_comp)
{
    const uint32_t row_block_count = static_cast<uint32_t>(row_blocks_comp.size());
    if (row_block_count != chunk_io.row_meta.size())
        return false;
    if (row_block_count != gt_row_blocks_comp.size())
        return false;

    const uint32_t header_bytes = 5u * sizeof(uint32_t); // magic, version, total_variants, row_block_count, flags
    const uint32_t entry_bytes =
        sizeof(uint32_t) + 2u * sizeof(int64_t) + 6u * (2u * sizeof(uint32_t)) + 2u * sizeof(uint32_t); // + (gt_off,gt_size)
    const uint32_t dir_bytes = row_block_count * entry_bytes;

    const uint64_t start_offset = ftell(temp_file);
    uint64_t payload_size_placeholder = 0;
    fwrite(&payload_size_placeholder, sizeof(payload_size_placeholder), 1, temp_file);
    const uint64_t payload_start = ftell(temp_file);

    uint32_t magic = GSC_FIXED_FIELDS_RB_MAGIC;
    uint32_t version = GSC_FIXED_FIELDS_RB_VERSION_LATEST;
    uint32_t total_variants = chunk_io.no_variants;
    uint32_t flags = 0;

    uint32_t cur_off = header_bytes + dir_bytes;

    struct DirEntryOffsets
    {
        uint32_t chrom_off, chrom_size;
        uint32_t pos_off, pos_size;
        uint32_t id_off, id_size;
        uint32_t ref_off, ref_size;
        uint32_t alt_off, alt_size;
        uint32_t qual_off, qual_size;
        uint32_t gt_off, gt_size;
    };
    std::vector<DirEntryOffsets> offsets(row_block_count);

    for (uint32_t i = 0; i < row_block_count; ++i)
    {
        const fixed_field_block &rb = row_blocks_comp[i];
        auto &o = offsets[i];

        o.chrom_off = cur_off;
        o.chrom_size = static_cast<uint32_t>(rb.chrom.size());
        cur_off += o.chrom_size;

        o.pos_off = cur_off;
        o.pos_size = static_cast<uint32_t>(rb.pos.size());
        cur_off += o.pos_size;

        o.id_off = cur_off;
        o.id_size = static_cast<uint32_t>(rb.id.size());
        cur_off += o.id_size;

        o.ref_off = cur_off;
        o.ref_size = static_cast<uint32_t>(rb.ref.size());
        cur_off += o.ref_size;

        o.alt_off = cur_off;
        o.alt_size = static_cast<uint32_t>(rb.alt.size());
        cur_off += o.alt_size;

        o.qual_off = cur_off;
        o.qual_size = static_cast<uint32_t>(rb.qual.size());
        cur_off += o.qual_size;

        o.gt_off = cur_off;
        o.gt_size = static_cast<uint32_t>(gt_row_blocks_comp[i].size());
        cur_off += o.gt_size;
    }

    // Header
    fwrite(&magic, sizeof(uint32_t), 1, temp_file);
    fwrite(&version, sizeof(uint32_t), 1, temp_file);
    fwrite(&total_variants, sizeof(uint32_t), 1, temp_file);
    fwrite(&row_block_count, sizeof(uint32_t), 1, temp_file);
    fwrite(&flags, sizeof(uint32_t), 1, temp_file);

    // Directory
    for (uint32_t i = 0; i < row_block_count; ++i)
    {
        const fixed_fields_row_block_meta &meta = chunk_io.row_meta[i];
        const auto &o = offsets[i];

        fwrite(&meta.variant_count, sizeof(uint32_t), 1, temp_file);
        fwrite(&meta.first_pos, sizeof(int64_t), 1, temp_file);
        fwrite(&meta.last_pos, sizeof(int64_t), 1, temp_file);

        fwrite(&o.chrom_off, sizeof(uint32_t), 1, temp_file);
        fwrite(&o.chrom_size, sizeof(uint32_t), 1, temp_file);
        fwrite(&o.pos_off, sizeof(uint32_t), 1, temp_file);
        fwrite(&o.pos_size, sizeof(uint32_t), 1, temp_file);
        fwrite(&o.id_off, sizeof(uint32_t), 1, temp_file);
        fwrite(&o.id_size, sizeof(uint32_t), 1, temp_file);
        fwrite(&o.ref_off, sizeof(uint32_t), 1, temp_file);
        fwrite(&o.ref_size, sizeof(uint32_t), 1, temp_file);
        fwrite(&o.alt_off, sizeof(uint32_t), 1, temp_file);
        fwrite(&o.alt_size, sizeof(uint32_t), 1, temp_file);
        fwrite(&o.qual_off, sizeof(uint32_t), 1, temp_file);
        fwrite(&o.qual_size, sizeof(uint32_t), 1, temp_file);

        fwrite(&o.gt_off, sizeof(uint32_t), 1, temp_file);
        fwrite(&o.gt_size, sizeof(uint32_t), 1, temp_file);
    }

    // Data region in the same order used for offsets above.
    for (uint32_t i = 0; i < row_block_count; ++i)
    {
        const fixed_field_block &rb = row_blocks_comp[i];
        fwrite(rb.chrom.data(), 1, rb.chrom.size(), temp_file);
        fwrite(rb.pos.data(), 1, rb.pos.size(), temp_file);
        fwrite(rb.id.data(), 1, rb.id.size(), temp_file);
        fwrite(rb.ref.data(), 1, rb.ref.size(), temp_file);
        fwrite(rb.alt.data(), 1, rb.alt.size(), temp_file);
        fwrite(rb.qual.data(), 1, rb.qual.size(), temp_file);
        fwrite(gt_row_blocks_comp[i].data(), 1, gt_row_blocks_comp[i].size(), temp_file);
    }

    const uint64_t payload_end = ftell(temp_file);
    const uint64_t payload_size = payload_end - payload_start;

    fseek(temp_file, start_offset, SEEK_SET);
    fwrite(&payload_size, sizeof(payload_size), 1, temp_file);
    fseek(temp_file, 0, SEEK_END);

    return true;
}

// Log field compression statistics at debug level
void Compressor::LogCompressionStats()
{
    auto logger = LogManager::Instance().Logger();
    if (!enable_field_stats_)
        return;

    // Helper function to format size with appropriate unit
    auto formatSize = [](uint64_t bytes) -> std::string {
        if (bytes < 1024) {
            return std::to_string(bytes) + " B";
        } else if (bytes < 1024 * 1024) {
            return std::to_string(bytes / 1024) + " KB";
        } else if (bytes < 1024ULL * 1024 * 1024) {
            char buf[32];
            snprintf(buf, sizeof(buf), "%.2f MB", bytes / (1024.0 * 1024));
            return std::string(buf);
        } else {
            char buf[32];
            snprintf(buf, sizeof(buf), "%.2f GB", bytes / (1024.0 * 1024 * 1024));
            return std::string(buf);
        }
    };

    // Helper function to format ratio as percentage
    auto formatRatio = [](double ratio) -> std::string {
        char buf[32];
        snprintf(buf, sizeof(buf), "%.2f%%", ratio * 100);
        return std::string(buf);
    };

    // Helper to calculate percentage of total
    auto calcPercent = [](uint64_t part, uint64_t total) -> double {
        if (total == 0) return 0.0;
        return static_cast<double>(part) / total * 100.0;
    };

    logger->debug("===============================================================================");
    logger->debug("                    Field Compression Statistics");
    logger->debug("===============================================================================");

    uint64_t total_raw = comp_stats_.GetTotalRawSize();
    uint64_t total_compressed = comp_stats_.GetTotalCompressedSize();

    logger->debug("Total variants: {}", comp_stats_.total_variants);
    logger->debug("Total samples: {}", comp_stats_.total_samples);
    logger->debug("Ploidy: {}", comp_stats_.ploidy);
    logger->debug("");

    logger->debug("{:<12} {:>15} {:>15} {:>12} {:>10}",
                  "Field", "Raw Size", "Compressed", "Ratio", "% of Total");
    logger->debug("-------------------------------------------------------------------------------");

    // Fixed fields
    auto logField = [&](const FieldCompressionStats& field) {
        uint64_t raw = field.raw_size.load();
        uint64_t comp = field.compressed_size.load();
        if (raw > 0 || comp > 0) {
            logger->debug("{:<12} {:>15} {:>15} {:>12} {:>9.2f}%",
                          field.name,
                          formatSize(raw),
                          formatSize(comp),
                          formatRatio(field.GetCompressionRatio()),
                          calcPercent(comp, total_compressed));
        }
    };

    logField(comp_stats_.chrom);
    logField(comp_stats_.pos);
    logField(comp_stats_.id);
    logField(comp_stats_.ref);
    logField(comp_stats_.alt);
    logField(comp_stats_.qual);
    logField(comp_stats_.filter);
    logField(comp_stats_.info);
    logField(comp_stats_.gt);
    logField(comp_stats_.meta);

    // FORMAT fields
    {
        std::lock_guard<std::mutex> lock(comp_stats_.format_mutex);
        for (const auto& kv : comp_stats_.format_fields) {
            logField(kv.second);
        }
    }

    logger->debug("-------------------------------------------------------------------------------");
    logger->debug("{:<12} {:>15} {:>15} {:>12}",
                  "TOTAL",
                  formatSize(total_raw),
                  formatSize(total_compressed),
                  formatRatio(total_compressed > 0 && total_raw > 0 ?
                             static_cast<double>(total_compressed) / total_raw : 0.0));
    logger->debug("===============================================================================");

    // GT field detailed statistics
    uint64_t gt_raw = comp_stats_.gt.raw_size.load();
    uint64_t gt_comp = comp_stats_.gt.compressed_size.load();
    if (gt_raw > 0 && comp_stats_.total_variants > 0 && comp_stats_.total_samples > 0) {
        uint64_t total_genotypes = static_cast<uint64_t>(comp_stats_.total_variants) *
                                   comp_stats_.total_samples * comp_stats_.ploidy;
        double bits_per_gt_raw = (gt_raw * 8.0) / total_genotypes;
        double bits_per_gt_comp = (gt_comp * 8.0) / total_genotypes;

        logger->debug("");
        logger->debug("GT Field Details:");
        logger->debug("  Total genotypes: {}", total_genotypes);
        logger->debug("  Bits/genotype (raw): {:.2f}", bits_per_gt_raw);
        logger->debug("  Bits/genotype (compressed): {:.4f}", bits_per_gt_comp);
    }
}
