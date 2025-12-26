#include "decompression_reader.h"
#include "logger.h"
#include "fmt_utils.h"
#include <algorithm>
#include <chrono>
#include <limits>
#include <unistd.h>

using namespace std::chrono;
using namespace std;

// *******************************************************************************************************************************
void DecompressionReader::out_perm(vector<uint32_t> &perm, vector<variant_desc_t> &v_vcf_data_io)
{

	vector<sblock> temp;
	for (size_t i = 0; i < v_vcf_data_io.size(); i++)
	{
		sblock sb(v_vcf_data_io[i].pos, i, v_vcf_data_io[i].ref);
		temp.emplace_back(sb);
	}
	my_merge_sort(temp, temp.size());
	// my_merge_sort(temp);
	for (size_t i = 0; i < v_vcf_data_io.size(); i++)
	{
		perm[i] = static_cast<uint32_t>(temp[i].p);
	}
	for (size_t i = 0; i < v_vcf_data_io.size(); i++)
	{

		if (atoi(temp[i].s_ref.c_str()))
		{

			temp[i].s_ref = temp[i].s_ref.substr(to_string(atoi(temp[i].s_ref.c_str())).length());
		}
		v_vcf_data_io[i].pos = temp[i].val;
		v_vcf_data_io[i].ref = temp[i].s_ref;
	}
}

// *******************************************************************************************************************************
bool DecompressionReader::OpenReading(const string &in_file_name, const bool &_decompression_mode_type)
{

	auto logger = LogManager::Instance().Logger();
	InitializeCompressionBackend(backend);
	fname = in_file_name;
	if (in_file_name == "-")
	{
		std::istream *in_stream = &std::cin;

		fname = in_file_name + ".decom_tmp_gsc";
		{
			std::ofstream temp_file(fname, std::ios::binary | std::ios::trunc);
			if (!temp_file)
			{
				logger->error("Failed to create temporary file: {}", fname);
				return 1;
			}
			const std::streamsize bufferSize = 1024 * 1024;
			char buffer[bufferSize];
			while (in_stream->read(buffer, bufferSize) || in_stream->gcount() > 0)
			{
				temp_file.write(buffer, in_stream->gcount());
			}

			temp_file.close();
		}
	}
	// sdsl vectors
	sdsl::isfstream in(fname, std::ios::binary | std::ios::in);
	if (!in)
	{
		// if (sdsl::util::verbose)
		{
			logger->error("Could not load file `{}`", fname);
		}
		exit(1);
	}

	uint64_t other_fields_offset;
	uint64_t sdsl_offset;
	bool mode_type;
	in.read((char *)&mode_type, sizeof(bool));
	in.read((char *)&other_fields_offset, sizeof(uint64_t));
	in.read((char *)&sdsl_offset, sizeof(uint64_t));
	file_mode_type = mode_type;
	if (mode_type)
	{
		std::string base = fname;
		size_t slash = base.find_last_of("/\\");
		if (slash != std::string::npos)
			base = base.substr(slash + 1);

		temp_file2_fname = "tmp/gsc_part2_" + std::to_string(getpid()) + "_" + base + ".tmp";
		in.seekg(other_fields_offset, std::ios::beg);
		std::ofstream out(temp_file2_fname, std::ios::binary | std::ios::trunc);
		if (!out)
		{
			// Fallback to the legacy behavior (next to the input file).
			temp_file2_fname = fname + ".tmp";
			out.open(temp_file2_fname, std::ios::binary | std::ios::trunc);
			if (!out)
			{
				logger->error("Cannot open destination file: {}", temp_file2_fname);
				return 1;
			}
		}
		const std::streamsize bufferSize = 1024 * 1024;
		char buffer[bufferSize];
		std::streamsize bytesToRead = sdsl_offset - other_fields_offset;

		while (bytesToRead > 0 && in)
		{
			std::streamsize readSize = std::min(bytesToRead, bufferSize);
			in.read(buffer, readSize);
			std::streamsize bytesRead = in.gcount();
			out.write(buffer, bytesRead);
			bytesToRead -= bytesRead;
		}
		if (in.bad())
		{
			logger->error("Error reading source file.");
			return 1;
		}

		out.close();
	}
	logger->info("Mode: {} other_fields_offset: {} sdsl_offset: {}", mode_type, other_fields_offset, sdsl_offset);
	in.seekg(sdsl_offset, std::ios::beg);

	rrr_zeros_bit_vector[0].load(in);
	rrr_zeros_bit_vector[1].load(in);
	rrr_copy_bit_vector[0].load(in);
	rrr_copy_bit_vector[1].load(in);
	in.seekg(sizeof(mode_type) + sizeof(other_fields_offset) + sizeof(sdsl_offset), std::ios::beg); // 17 = sizeof()
	uint64_t FileStartPosition = in.tellg();

	in.close();
	if (sdsl::util::verbose)
	{
		logger->info("Load file `{}`", fname);
	}

	// rest of archive
	buf_pos = 0;

#ifndef MMAP
	{
		uint64_t arch_size = 0;
		FILE *comp = fopen(fname.c_str(), "rb");

		if (!comp)
		{
			logger->error("Input file ({}) error", fname);
			exit(1);
		}
		arch_size = other_fields_offset - FileStartPosition;
		fseek(comp, FileStartPosition, SEEK_SET);

		buf = new uchar[arch_size];
		fread(buf, sizeof(uchar) * arch_size, 1, comp);

		fclose(comp);

		logger->info("File start position: {}", FileStartPosition);
		logger->info("Archive size: {}", arch_size);
	}
#else // if BOOST
	/*    const boost::interprocess::mode_t mode = boost::interprocess::read_only;
		fm = new boost::interprocess::file_mapping(fname.c_str(), mode);
		region = new boost::interprocess::mapped_region(*fm, mode, 0, 0);
		buf = reinterpret_cast<unsigned char*>(region->get_address()) + FileStartPosition;
		arch_size = region->get_size() - FileStartPosition;*/
	fm = new memory_mapped_file::read_only_mmf(fname.c_str());
	if (!fm->is_open())
	{
		logger->error("No file: {}", fname);
		exit(1);
	}
	buf = (uint8_t *)fm->data() + FileStartPosition;
	// arch_size = other_fields_offset - FileStartPosition;

#endif
	uint32_t chunks_streams_size;
	memcpy(&chunks_streams_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	// std::cerr<<"chunks_streams_size: "<<chunks_streams_size<<endl;
	for (uint32_t i = 0; i < chunks_streams_size; i++)
	{
		size_t offset;
		uint32_t cur_chunk_actual_pos;

		memcpy(&cur_chunk_actual_pos, buf + buf_pos, sizeof(uint32_t));
		buf_pos = buf_pos + sizeof(uint32_t);

		memcpy(&offset, buf + buf_pos, sizeof(size_t));
		buf_pos = buf_pos + sizeof(size_t);

		chunks_streams[i].cur_chunk_actual_pos = cur_chunk_actual_pos;
		chunks_streams[i].offset = offset;
	}
	memcpy(&ploidy, buf + buf_pos, sizeof(uint8_t));
	buf_pos = buf_pos + sizeof(uint8_t);

	memcpy(&max_block_rows, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);

	memcpy(&max_block_cols, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);

	memcpy(&vec_len, buf + buf_pos, sizeof(uint64_t));
	buf_pos = buf_pos + sizeof(uint64_t);

	memcpy(&no_vec, buf + buf_pos, sizeof(uint64_t));
	buf_pos = buf_pos + sizeof(uint64_t);

	memcpy(&no_copy, buf + buf_pos, sizeof(uint64_t));
	buf_pos = buf_pos + sizeof(uint64_t);

	memcpy(&used_bits_cp, buf + buf_pos, sizeof(char));
	buf_pos = buf_pos + sizeof(char);

	memcpy(&bm_comp_cp_size, buf + buf_pos, sizeof(int));
	buf_pos = buf_pos + sizeof(int);

	bm_comp_copy_orgl_id.Open(buf + buf_pos, bm_comp_cp_size);
	buf_pos = buf_pos + sizeof(uint8_t) * bm_comp_cp_size;

	// memcpy(&max_no_vec_in_block, buf + buf_pos, sizeof(uint32_t));
	// buf_pos = buf_pos + sizeof(uint32_t);

	memcpy(&n_samples, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);

	uint32_t chunks_min_pos_size;
	memcpy(&chunks_min_pos_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	// std::cerr<<"chunks_min_pos_size: "<<chunks_min_pos_size<<endl;
	chunks_min_pos.resize(chunks_min_pos_size);
	memcpy(&chunks_min_pos[0], buf + buf_pos, chunks_min_pos_size * sizeof(int64_t));
	buf_pos = buf_pos + chunks_min_pos_size * sizeof(int64_t);
	// for(uint32_t i = 0; i < chunks_min_pos_size; i++){
	// 	std::cerr<<chunks_min_pos[i]<<endl;
	// }
	uint32_t where_chrom_size;
	memcpy(&where_chrom_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);

	d_where_chrom.resize(where_chrom_size);
	// std::cerr<<where_chrom_size<<endl;
	for (size_t i = 0; i < where_chrom_size; ++i)
	{
		size_t chrom_size;
		memcpy(&chrom_size, buf + buf_pos, sizeof(size_t));
		buf_pos = buf_pos + sizeof(size_t);
		// std::cerr<<chrom_size<<endl;
		d_where_chrom[i].first.resize(chrom_size);
		memcpy(&d_where_chrom[i].first[0], buf + buf_pos, chrom_size * sizeof(char));
		buf_pos = buf_pos + chrom_size * sizeof(char);

		memcpy(&d_where_chrom[i].second, buf + buf_pos, sizeof(uint32_t));
		buf_pos = buf_pos + sizeof(uint32_t);
	}

	// Read GT column tiling metadata (with backward compatibility for old files)
	// Try to read column block metadata (may not exist in old files)
	if (buf_pos + sizeof(uint32_t) <= sdsl_offset)
	{
		uint32_t saved_pos = buf_pos;
		memcpy(&n_col_blocks, buf + buf_pos, sizeof(uint32_t));
		buf_pos += sizeof(uint32_t);

		// Validate that n_col_blocks is reasonable (between 1 and 10000)
		if (n_col_blocks > 0 && n_col_blocks <= 10000)
		{
			col_block_ranges.reserve(n_col_blocks);
			col_block_vec_lens.reserve(n_col_blocks);
			total_haplotypes = n_samples * ploidy;

			for (uint32_t cb = 0; cb < n_col_blocks; ++cb)
			{
				uint32_t start, size;
				memcpy(&start, buf + buf_pos, sizeof(uint32_t));
				buf_pos += sizeof(uint32_t);
				memcpy(&size, buf + buf_pos, sizeof(uint32_t));
				buf_pos += sizeof(uint32_t);
				col_block_ranges.emplace_back(start, size);
				col_block_vec_lens.push_back((size + 7) / 8);
			}
		}
		else
		{
			// Invalid n_col_blocks, restore position (old file format)
			buf_pos = saved_pos;
			n_col_blocks = 1;
			total_haplotypes = n_samples * ploidy;
			col_block_ranges.emplace_back(0, total_haplotypes);
			col_block_vec_lens.push_back(vec_len);
		}
	}
	else
	{
		// Old format: no column tiling metadata
		n_col_blocks = 1;
		total_haplotypes = n_samples * ploidy;
		col_block_ranges.emplace_back(0, total_haplotypes);
		col_block_vec_lens.push_back(vec_len);
	}

	// Set useLegacyPath flag before reading permutations
	if (max_block_rows == 0)
		max_block_rows = total_haplotypes;
	if (max_block_cols == 0)
		max_block_cols = total_haplotypes;
	bool legacy_cols = (max_block_cols >= total_haplotypes);
	bool legacy_rows = (max_block_rows == total_haplotypes);
	useLegacyPath = legacy_cols && legacy_rows && (n_col_blocks == 1);

	// Read vint_last_perm (supports both old and new formats)
	uint32_t vint_last_perm_size;
	memcpy(&vint_last_perm_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	// std::cerr<<"vint_last_perm_size: "<<vint_last_perm_size<<endl;
	for (uint32_t i = 0; i < vint_last_perm_size; ++i)
	{
		uint32_t data_size;
		uint32_t first, second = 0;

		memcpy(&first, buf + buf_pos, sizeof(uint32_t));
		buf_pos = buf_pos + sizeof(uint32_t);

		// New format: read col_block_id (second parameter)
		// Old format (legacy path) stores only row_block_id
		if (!useLegacyPath)
		{
			memcpy(&second, buf + buf_pos, sizeof(uint32_t));
			buf_pos = buf_pos + sizeof(uint32_t);
		}

		memcpy(&data_size, buf + buf_pos, sizeof(uint32_t));
		buf_pos = buf_pos + sizeof(uint32_t);
		// std::cerr<<first<<":"<<data_size<<endl;
		vector<uint8_t> data(data_size);

		memcpy(&data[0], buf + buf_pos, data_size * sizeof(uint8_t));
		buf_pos = buf_pos + data_size * sizeof(uint8_t);

		// Store in 2D map
		vint_last_perm_2d.emplace(make_pair(first, second), data);

		// For backward compatibility: also store in legacy map
		if (useLegacyPath)
		{
			vint_last_perm.emplace(first, data);
		}
	}

	// std::cerr<<"vint_last_perm_size: "<<vint_last_perm_size<<endl;
	uint32_t comp_size;
	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);

	comp_v_header.resize(comp_size);
	memcpy(&comp_v_header[0], buf + buf_pos, comp_size * sizeof(uint8_t));
	buf_pos = buf_pos + comp_size * sizeof(uint8_t);

	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);

	comp_v_samples.resize(comp_size);
	memcpy(&comp_v_samples[0], buf + buf_pos, comp_size * sizeof(uint8_t));
	buf_pos = buf_pos + comp_size * sizeof(uint8_t);

	// Precompute block offsets per chunk for tiled permutations
	uint32_t row_block_variants = max_block_rows ? max_block_rows : total_haplotypes;
	chunk_block_offsets.clear();
	if (row_block_variants)
	{
		chunk_block_offsets.resize(chunks_streams_size);
		uint32_t blocks_before = 0;
		uint32_t prev_variants = 0;
		for (uint32_t chunk_id = 0; chunk_id < chunks_streams_size; ++chunk_id)
		{
			auto it = chunks_streams.find(chunk_id);
			uint32_t cur_variants = 0;
			if (it != chunks_streams.end())
				cur_variants = it->second.cur_chunk_actual_pos;
			uint32_t chunk_variants = cur_variants - prev_variants;
			uint32_t blocks_in_chunk = (chunk_variants + row_block_variants - 1) / row_block_variants;
			chunk_block_offsets[chunk_id] = blocks_before;
			blocks_before += blocks_in_chunk;
			prev_variants = cur_variants;
		}
	}

	buf = buf - FileStartPosition;

	return true;
}
void DecompressionReader::SetNoThreads(uint32_t _no_threads)
{
	no_threads = _no_threads;
}
bool DecompressionReader::OpenReadingPart2(const string &in_file_name, bool start_threads)
{

	auto logger = LogManager::Instance().Logger();
	vector<uint8_t> part2_params_data;
	vector<uint8_t> v_desc;
	size_t pos_part2_params = 0;
	uint64_t temp = 0;
	if (file_handle2)
		delete file_handle2;
	file_handle2 = new File_Handle_2(true);

	if (!file_handle2->Open(temp_file2_fname))
	{
		logger->error("Can't open the temporary file: {}", temp_file2_fname);
		logger->warn("The .gsc file is generated through lossy compression mode and can only be decompressed using a lossy mode.");
		return false;
	}
	int stream_id = file_handle2->GetStreamId("part2_params");

	if (stream_id < 0)
	{
		logger->error("Corrupted part2!");
		exit(0);
	}
	if (!file_handle2->GetPart(stream_id, part2_params_data))
	{
		logger->error("Corrupted part2!");
		exit(0);
	}
	if (part2_params_data.empty())
	{
		logger->error("Corrupted part2: empty params payload.");
		return false;
	}
	const uint8_t backend_id = part2_params_data[0];
	if (backend_id <= static_cast<uint8_t>(compression_backend_t::brotli))
	{
		backend = static_cast<compression_backend_t>(backend_id);
	}
	InitializeCompressionBackend(backend);
	auto codec = MakeCompressionStrategy(backend, p_bsc_fixed_fields);
	std::vector<uint8_t> compressed_payload(part2_params_data.begin() + 1, part2_params_data.end());
	codec->Decompress(compressed_payload, v_desc);
	uint32_t actual_variants_size = 0;
	read(v_desc, pos_part2_params, actual_variants_size);
	actual_variants.resize(actual_variants_size);

	for (uint32_t i = 0; i < actual_variants_size; i++)
	{
		read(v_desc, pos_part2_params, actual_variants[i]);
	}
	read(v_desc, pos_part2_params, no_keys);
	read(v_desc, pos_part2_params, key_gt_id);
	keys.resize(no_keys);

	for (uint32_t i = 0; i < no_keys; i++)
	{
		read(v_desc, pos_part2_params, keys[i].key_id);
		read(v_desc, pos_part2_params, keys[i].actual_field_id);
		read(v_desc, pos_part2_params, temp);
		keys[i].keys_type = (key_type_t)temp;
		read(v_desc, pos_part2_params, keys[i].type);
		read_str(v_desc, pos_part2_params, keys[i].name);  // Read field name for FMT decompression
		// std::cerr<<keys[i].key_id<<" "<<keys[i].actual_field_id<<endl;
	}
		if (pos_part2_params < v_desc.size())
		{
			uint32_t backend_id = 0;
			if (pos_part2_params + sizeof(uint32_t) <= v_desc.size())
			{
				read(v_desc, pos_part2_params, backend_id);
				if (backend_id <= static_cast<uint32_t>(compression_backend_t::brotli))
				{
					backend = static_cast<compression_backend_t>(backend_id);
				}
			}
		}

		// FORMAT dictionaries tail (AD/PL/PID). Backward compatible.
		fmt_dictionaries_ = std::make_unique<fmt_compress::FmtDictionaries>();
	if (pos_part2_params + 2 * sizeof(uint32_t) <= v_desc.size())
	{
			uint32_t dict_ver = 0;
			uint32_t dict_size = 0;
			read(v_desc, pos_part2_params, dict_ver);
			read(v_desc, pos_part2_params, dict_size);
			if (dict_ver == 1 && dict_size > 0)
			{
				if (pos_part2_params + dict_size <= v_desc.size())
				{
					fmt_dictionaries_->unserialize(v_desc.data() + pos_part2_params, dict_size);
					pos_part2_params += dict_size;
				}
				else
				{
					logger->error("Corrupted part2: fmt dict size out of bounds.");
					return false;
				}
			}
	}

	if (!start_threads)
	{
		return true;
	}

	InitDecompressParams();
	for (uint32_t i = 0; i < no_keys; ++i)
	{
		decomp_part_queue->PushQueue(i);
	}
	part_decompress_thread.reserve(no_threads);

	for (uint32_t i = 0; i < no_threads; ++i)
	{
		part_decompress_thread.emplace_back(thread([&]()
													   {
														   while (!decomp_part_queue->IsComplete())
														   {
															   uint32_t p_id;
															   SPackage *pck = new SPackage;
															   if (!decomp_part_queue->PopQueue(p_id))
															   {
																   delete pck;
																   break;
															   }

															   pck->key_id = p_id;
															   pck->stream_id_size = file_handle2->GetStreamId("key_" + to_string(p_id) + "_size");
															   pck->stream_id_data = file_handle2->GetStreamId("key_" + to_string(p_id) + "_data");

															   if (decompress_other_fileds(pck))
															   {

																   lock_guard<mutex> lck(m_packages);

																   v_packages[pck->key_id] = pck;
															   }
															   else
															   {

																   pck->v_size.clear();
																   pck->v_data.clear();

																   lock_guard<mutex> lck(m_packages);
																   v_packages[pck->key_id] = pck;
															   }
															   cv_packages.notify_all();
														   } }));
	}

	return true;
}
//**********************************************************************************************************************
void DecompressionReader::close()
{
	if (decomp_part_queue)
	{
		decomp_part_queue->Complete();
		for (auto &t : part_decompress_thread)
			t.join();
	}

	if (!temp_file2_fname.empty())
	{
		if (remove(temp_file2_fname.c_str()) != 0)
		{
			perror("Error deleting temp2 file");
		}
	}

	// file_handle2->Close();
}
//**********************************************************************************************************************
bool DecompressionReader::decompress_other_fileds(SPackage *pck)
{
	auto logger = LogManager::Instance().Logger();
	vector<uint8_t> v_compressed;
	vector<uint8_t> v_tmp;
	CompressionStrategy *cbsc_size = field_size_codecs[pck->key_id].get();
	CompressionStrategy *cbsc = field_data_codecs[pck->key_id].get();

	if (!file_handle2->GetPart(pck->stream_id_size, v_compressed))
		return false;

	cbsc_size->Decompress(v_compressed, v_tmp);

	pck->v_size.resize(v_tmp.size() / 4);

	copy_n(v_tmp.data(), v_tmp.size(), (uint8_t *)pck->v_size.data());

	v_compressed.clear();

	v_tmp.clear();

	file_handle2->GetPart(pck->stream_id_data, v_compressed);

	// Check for special FMT fields
	const std::string& field_name = keys[pck->key_id].name;
		bool is_sparse_field = (field_name == "PGT" || field_name == "PID");
		bool is_ad_field = (field_name == "AD");
		bool is_dp_field = (field_name == "DP");
		bool is_min_dp_field = (field_name == "MIN_DP");
		bool is_gq_field = (field_name == "GQ");
		bool is_pl_field = (field_name == "PL");
		bool is_fmt_field = (keys[pck->key_id].keys_type == key_type_t::fmt);
		uint32_t n_samples_u32 = n_samples;


	// New lossless codecs for selected FORMAT fields (AD/PL/PGT/PID).
	// After decompression, the first byte is a codec marker:
		//  - 0: legacy payload (INT=transposed bytes, STR=raw bytes)
		//  - 1: special payload (record-wise stream)
		if (v_compressed.size() && is_fmt_field && n_samples_u32 > 0)
		{
				if ((is_dp_field || is_min_dp_field || is_gq_field) && keys[pck->key_id].type == BCF_HT_INT)
				{
					cbsc->Decompress(v_compressed, v_tmp);
					if (v_tmp.empty())
					{
						pck->v_data.clear();
					return true;
				}

				const uint8_t codec_id = v_tmp[0];
				if (codec_id == 0)
				{
					vector<uint8_t> transposed(v_tmp.begin() + 1, v_tmp.end());
					pck->v_data.resize(transposed.size());
					Decoder(transposed, pck->v_data);
					return true;
				}
					if (codec_id != 1)
					{
						logger->error("Unsupported codec_id={} for FMT {} (INT)", (int)codec_id, field_name);
						return false;
					}

					pck->v_data.assign(v_tmp.begin() + 1, v_tmp.end());
					return true;
				}

			if (is_sparse_field && keys[pck->key_id].type == BCF_HT_STR)
			{
				cbsc->Decompress(v_compressed, v_tmp);
				if (v_tmp.empty())
				{
				pck->v_data.clear();
				return true;
			}

				const uint8_t codec_id = v_tmp[0];
				if (codec_id == 0)
				{
					pck->v_data.assign(v_tmp.begin() + 1, v_tmp.end());
					return true;
				}
				if (codec_id != 1 && codec_id != 2)
				{
					logger->error("Unsupported codec_id={} for FMT {} (STR)", (int)codec_id, field_name);
					return false;
				}
				if (codec_id == 2 && field_name != "PID")
				{
					logger->error("codec_id=2 is only supported for PID, got {}", field_name);
					return false;
				}

			auto decode_vint = [&](size_t &off, uint64_t &val) -> bool {
				if (off >= v_tmp.size())
					return false;
				uint8_t len = fmt_compress::VintCodec::decode(v_tmp.data() + off, val);
				if (len == 0 || off + len > v_tmp.size())
					return false;
				off += len;
				return true;
			};

			size_t total_bytes = 0;
			for (size_t i = 0; i < pck->v_size.size(); ++i)
				total_bytes += pck->v_size[i];
			pck->v_data.assign(total_bytes, 0);

			size_t in_off = 1;
			size_t out_off = 0;
			for (size_t rec_idx = 0; rec_idx < pck->v_size.size(); ++rec_idx)
			{
				const uint32_t rec_bytes = pck->v_size[rec_idx];
				if (rec_bytes == 0)
					continue;
				if (rec_bytes % n_samples_u32 != 0)
				{
					logger->error("FMT {} STR record_bytes={} not divisible by n_samples={}", field_name, rec_bytes, n_samples_u32);
					return false;
				}
				const uint32_t stride = rec_bytes / n_samples_u32;

				uint8_t *rec_out = pck->v_data.data() + out_off;
				memset(rec_out, 0, rec_bytes);
				if (stride > 0)
				{
					for (uint32_t s = 0; s < n_samples_u32; ++s)
						rec_out[(size_t)s * stride] = bcf_str_missing;
				}

					uint64_t count = 0;
					if (!decode_vint(in_off, count))
						return false;
					uint32_t prev_pos = 0;
					for (uint64_t i = 0; i < count; ++i)
					{
						uint64_t pos = 0;
						if (codec_id == 2)
						{
							uint64_t delta = 0;
							if (!decode_vint(in_off, delta))
								return false;
							pos = (uint64_t)prev_pos + delta;
							prev_pos = (uint32_t)pos;
						}
						else
						{
							if (!decode_vint(in_off, pos))
								return false;
						}
						if (pos >= n_samples_u32)
						{
							logger->error("FMT {} STR decoded pos={} out of range n_samples={}", field_name, pos, n_samples_u32);
							return false;
						}

						if (codec_id == 1)
						{
							if (in_off + stride > v_tmp.size())
							{
								logger->error("FMT {} STR payload truncated (need {}, have {})", field_name, stride, v_tmp.size() - in_off);
								return false;
							}
							memcpy(rec_out + (size_t)pos * stride, v_tmp.data() + in_off, stride);
							in_off += stride;
						}
						else
						{
							uint64_t tag = 0;
							if (!decode_vint(in_off, tag))
								return false;

							uint8_t *cell = rec_out + (size_t)pos * stride;
							if (tag == 0)
							{
								uint64_t raw_len = 0;
								if (!decode_vint(in_off, raw_len))
									return false;
								if (raw_len > stride)
								{
									logger->error("FMT PID raw_len={} exceeds stride={}", raw_len, stride);
									return false;
								}
								if (in_off + raw_len > v_tmp.size())
								{
									logger->error("FMT PID raw payload truncated (need {}, have {})", raw_len, v_tmp.size() - in_off);
									return false;
								}
								memcpy(cell, v_tmp.data() + in_off, (size_t)raw_len);
								in_off += (size_t)raw_len;
							}
							else
							{
								if (!fmt_dictionaries_)
								{
									logger->error("FMT PID codec_id=2 requires dictionaries");
									return false;
								}
								const uint32_t id = (uint32_t)(tag - 1);
								const uint8_t *ptr = fmt_dictionaries_->getPIDItemPtr(id);
								if (!ptr)
								{
									logger->error("FMT PID dict id={} not found", id);
									return false;
								}
								uint16_t len = 0;
								memcpy(&len, ptr, 2);
								if ((uint32_t)len + 1u > stride)
								{
									logger->error("FMT PID dict len={} exceeds stride={}", (uint32_t)len, stride);
									return false;
								}
								if (len)
									memcpy(cell, ptr + 2, len);
								cell[len] = bcf_str_vector_end;
							}
						}
					}

				out_off += rec_bytes;
			}

			if (out_off != total_bytes || in_off != v_tmp.size())
			{
				logger->error("FMT {} STR decode mismatch (out_off={}, total_bytes={}, in_off={}, tmp_size={})",
					field_name, out_off, total_bytes, in_off, v_tmp.size());
				return false;
			}
			return true;
		}

		if ((is_ad_field || is_pl_field) && keys[pck->key_id].type == BCF_HT_INT)
		{
			cbsc->Decompress(v_compressed, v_tmp);
			if (v_tmp.empty())
			{
				pck->v_data.clear();
				return true;
			}

				const uint8_t codec_id = v_tmp[0];
				if (codec_id == 0)
				{
					vector<uint8_t> transposed(v_tmp.begin() + 1, v_tmp.end());
					pck->v_data.resize(transposed.size());
					Decoder(transposed, pck->v_data);
					return true;
				}
				if (codec_id != 1 && codec_id != 2)
				{
					logger->error("Unsupported codec_id={} for FMT {} (INT)", (int)codec_id, field_name);
					return false;
				}

			auto decode_vint = [&](size_t &off, uint64_t &val) -> bool {
				if (off >= v_tmp.size())
					return false;
				uint8_t len = fmt_compress::VintCodec::decode(v_tmp.data() + off, val);
				if (len == 0 || off + len > v_tmp.size())
					return false;
				off += len;
				return true;
			};

			auto get_bit = [](const uint8_t *buf, size_t bit) -> uint8_t {
				return (buf[bit / 8] >> (bit % 8)) & 1u;
			};

			size_t total_elems = 0;
			for (size_t i = 0; i < pck->v_size.size(); ++i)
				total_elems += pck->v_size[i];
			pck->v_data.assign(total_elems * 4, 0);

			const size_t tip_bytes = ((size_t)n_samples_u32 * 2 + 7) / 8;
			size_t in_off = 1;
			size_t out_off = 0;

			for (size_t rec_idx = 0; rec_idx < pck->v_size.size(); ++rec_idx)
			{
				const uint32_t total = pck->v_size[rec_idx];
				if (total == 0)
					continue;
				if (total % n_samples_u32 != 0)
				{
					logger->error("FMT {} INT total_elems={} not divisible by n_samples={}", field_name, total, n_samples_u32);
					return false;
				}
				const uint32_t per_sample = total / n_samples_u32;
				if (in_off + tip_bytes > v_tmp.size())
				{
					logger->error("FMT {} INT payload truncated in tips", field_name);
					return false;
				}
				const uint8_t *tips = v_tmp.data() + in_off;
				in_off += tip_bytes;

				int32_t *out = reinterpret_cast<int32_t *>(pck->v_data.data() + out_off);
				for (uint32_t s = 0; s < n_samples_u32; ++s)
				{
					const uint8_t tip0 = get_bit(tips, (size_t)s * 2);
					const uint8_t tip1 = get_bit(tips, (size_t)s * 2 + 1);
					int32_t *sv = out + (size_t)s * per_sample;

					if (is_ad_field)
					{
						if (tip0 == 0 && tip1 == 0)
						{
							for (uint32_t j = 0; j < per_sample; ++j) sv[j] = 0;
						}
						else if (tip0 == 0 && tip1 == 1)
						{
							if (per_sample > 0) sv[0] = 2;
							for (uint32_t j = 1; j < per_sample; ++j) sv[j] = 0;
						}
						else if (tip0 == 1 && tip1 == 0)
						{
							uint64_t v = 0;
							if (!decode_vint(in_off, v)) return false;
							if (per_sample > 0) sv[0] = (int32_t)(uint32_t)v;
							for (uint32_t j = 1; j < per_sample; ++j) sv[j] = 0;
							}
							else
							{
								if (codec_id == 1)
								{
									for (uint32_t j = 0; j < per_sample; ++j)
									{
										uint64_t v = 0;
										if (!decode_vint(in_off, v)) return false;
										sv[j] = (int32_t)(uint32_t)v;
									}
								}
								else
								{
									uint64_t tag = 0;
									if (!decode_vint(in_off, tag)) return false;
									if (tag == 0)
									{
										for (uint32_t j = 0; j < per_sample; ++j)
										{
											uint64_t v = 0;
											if (!decode_vint(in_off, v)) return false;
											sv[j] = (int32_t)(uint32_t)v;
										}
									}
									else
									{
										if (!fmt_dictionaries_)
										{
											logger->error("FMT AD codec_id=2 requires dictionaries");
											return false;
										}
										const uint32_t id = (uint32_t)(tag - 1);
										const uint8_t *ptr = fmt_dictionaries_->getADItemPtr(id);
										if (!ptr)
										{
											logger->error("FMT AD dict id={} not found", id);
											return false;
										}
										const uint8_t item_total = (uint8_t)(ptr[0] + 1);
										const uint8_t type = ptr[1];
										const uint8_t *data = ptr + 2;
										const uint32_t bytes_per = (type == 0) ? 1u : (type == 1) ? 2u : 4u;
										if (2u + per_sample * bytes_per > item_total)
										{
											logger->error("FMT AD dict item too short (need {}, have {})", 2u + per_sample * bytes_per, (uint32_t)item_total);
											return false;
										}
										for (uint32_t j = 0; j < per_sample; ++j)
										{
											uint32_t u = 0;
											if (type == 0)
												u = data[j];
											else if (type == 1)
												u = reinterpret_cast<const uint16_t *>(data)[j];
											else
												u = reinterpret_cast<const uint32_t *>(data)[j];
											sv[j] = (int32_t)u;
										}
									}
								}
							}
						}
						else
						{
						if (tip0 == 0 && tip1 == 0)
						{
							for (uint32_t j = 0; j < per_sample; ++j) sv[j] = 0;
						}
						else if (tip0 == 0 && tip1 == 1)
						{
							uint64_t a = 0, b = 0;
							if (!decode_vint(in_off, a) || !decode_vint(in_off, b)) return false;
							const int32_t ai = (int32_t)(uint32_t)a;
							const int32_t bi = (int32_t)(uint32_t)b;
							if (per_sample > 0) sv[0] = 0;
							if (per_sample > 1) sv[1] = ai;
							if (per_sample > 2) sv[2] = bi;
							uint32_t pos = 3;
							uint32_t b_count = 2;
							while (pos < per_sample)
							{
								sv[pos++] = ai;
								for (uint32_t k = 0; k < b_count && pos < per_sample; ++k) sv[pos++] = bi;
								++b_count;
							}
						}
						else if (tip0 == 1 && tip1 == 0)
						{
							uint64_t a = 0;
							if (!decode_vint(in_off, a)) return false;
							const int32_t ai = (int32_t)(uint32_t)a;
							const int32_t bi = (int32_t)(uint32_t)(a * 15);
							if (per_sample > 0) sv[0] = 0;
							if (per_sample > 1) sv[1] = ai;
							if (per_sample > 2) sv[2] = bi;
							uint32_t pos = 3;
							uint32_t b_count = 2;
							while (pos < per_sample)
							{
								sv[pos++] = ai;
								for (uint32_t k = 0; k < b_count && pos < per_sample; ++k) sv[pos++] = bi;
								++b_count;
							}
							}
							else
							{
								if (codec_id == 1)
								{
									for (uint32_t j = 0; j < per_sample; ++j)
									{
										uint64_t v = 0;
										if (!decode_vint(in_off, v)) return false;
										sv[j] = (int32_t)(uint32_t)v;
									}
								}
								else
								{
									uint64_t tag = 0;
									if (!decode_vint(in_off, tag)) return false;
									if (tag == 0)
									{
										for (uint32_t j = 0; j < per_sample; ++j)
										{
											uint64_t v = 0;
											if (!decode_vint(in_off, v)) return false;
											sv[j] = (int32_t)(uint32_t)v;
										}
									}
									else
									{
										if (!fmt_dictionaries_)
										{
											logger->error("FMT PL codec_id=2 requires dictionaries");
											return false;
										}
										const uint32_t id = (uint32_t)(tag - 1);
										const uint8_t *ptr = fmt_dictionaries_->getPLItemPtr(id);
										if (!ptr)
										{
											logger->error("FMT PL dict id={} not found", id);
											return false;
										}
										const uint8_t item_total = (uint8_t)(ptr[0] + 1);
										const uint8_t type = ptr[1];
										const uint8_t *data = ptr + 2;
										const uint32_t bytes_per = (type == 0) ? 1u : (type == 1) ? 2u : 4u;
										if (2u + per_sample * bytes_per > item_total)
										{
											logger->error("FMT PL dict item too short (need {}, have {})", 2u + per_sample * bytes_per, (uint32_t)item_total);
											return false;
										}
										for (uint32_t j = 0; j < per_sample; ++j)
										{
											uint32_t u = 0;
											if (type == 0)
												u = data[j];
											else if (type == 1)
												u = reinterpret_cast<const uint16_t *>(data)[j];
											else
												u = reinterpret_cast<const uint32_t *>(data)[j];
											sv[j] = (int32_t)u;
										}
									}
								}
							}
						}
					}

				out_off += (size_t)total * 4;
			}

			if (out_off != total_elems * 4 || in_off != v_tmp.size())
			{
				logger->error("FMT {} INT decode mismatch (out_off={}, exp={}, in_off={}, tmp_size={})",
					field_name, out_off, total_elems * 4, in_off, v_tmp.size());
				return false;
			}
			return true;
		}
	}

	if (v_compressed.size())
	{
		if (false && is_sparse_field && keys[pck->key_id].type == BCF_HT_STR)  // DISABLED for testing
		{
			// Decompress sparse PGT/PID
			cbsc->Decompress(v_compressed, v_tmp);

			// Parse sparse format: [count][pos1][len1][data1]...
			size_t offset = 0;
			uint64_t count;
			offset += fmt_compress::VintCodec::decode(v_tmp.data() + offset, count);

			// Build position map
			std::unordered_map<uint32_t, std::pair<uint32_t, const uint8_t*>> entries;
			for (uint64_t i = 0; i < count; i++)
			{
				uint64_t pos, len;
				offset += fmt_compress::VintCodec::decode(v_tmp.data() + offset, pos);
				offset += fmt_compress::VintCodec::decode(v_tmp.data() + offset, len);
				entries[static_cast<uint32_t>(pos)] = {static_cast<uint32_t>(len), v_tmp.data() + offset};
				offset += len;
			}

			// Reconstruct full data
			size_t num_samples = pck->v_size.size();
			pck->v_data.clear();
			for (size_t i = 0; i < num_samples; i++)
			{
				auto it = entries.find(static_cast<uint32_t>(i));
				if (it != entries.end())
				{
					pck->v_data.insert(pck->v_data.end(), it->second.second, it->second.second + it->second.first);
					pck->v_size[i] = it->second.first;
				}
				else
				{
					pck->v_data.push_back('.');
					pck->v_size[i] = 1;
				}
			}
		}
		else if (false && is_ad_field && keys[pck->key_id].type == BCF_HT_INT)  // DISABLED for testing
		{
			// Decompress AD with 2-bit tip encoding
			cbsc->Decompress(v_compressed, v_tmp);

			size_t offset = 0;
			uint64_t num_samples;
			offset += fmt_compress::VintCodec::decode(v_tmp.data() + offset, num_samples);

			// Read tips
			size_t tip_bytes = (num_samples * 2 + 7) / 8;
			std::vector<uint8_t> tips;
			tips.reserve(num_samples * 2);
			for (size_t i = 0; i < num_samples * 2; i++)
			{
				size_t byte_idx = i / 8;
				size_t bit_idx = i % 8;
				tips.push_back((v_tmp[offset + byte_idx] >> bit_idx) & 1);
			}
			offset += tip_bytes;

			// Read values
			uint64_t value_count;
			offset += fmt_compress::VintCodec::decode(v_tmp.data() + offset, value_count);
			std::vector<uint32_t> values(value_count);
			for (uint64_t i = 0; i < value_count; i++)
			{
				uint64_t val;
				offset += fmt_compress::VintCodec::decode(v_tmp.data() + offset, val);
				values[i] = static_cast<uint32_t>(val);
			}

			// Reconstruct AD data
			pck->v_data.clear();
			size_t val_idx = 0;
			for (size_t i = 0; i < num_samples; i++)
			{
				uint8_t tip0 = tips[i * 2];
				uint8_t tip1 = tips[i * 2 + 1];
				// v_size stores element count, not byte count (see WriteInt in bit_memory.h)
				uint32_t elem_count = pck->v_size[i];

				if (tip0 == 0 && tip1 == 0)
				{
					// All zeros
					for (uint32_t j = 0; j < elem_count; j++)
					{
						int32_t val = 0;
						pck->v_data.insert(pck->v_data.end(), reinterpret_cast<uint8_t*>(&val), reinterpret_cast<uint8_t*>(&val) + 4);
					}
				}
				else if (tip0 == 0 && tip1 == 1)
				{
					// First = 2, rest = 0
					int32_t val = 2;
					pck->v_data.insert(pck->v_data.end(), reinterpret_cast<uint8_t*>(&val), reinterpret_cast<uint8_t*>(&val) + 4);
					val = 0;
					for (uint32_t j = 1; j < elem_count; j++)
					{
						pck->v_data.insert(pck->v_data.end(), reinterpret_cast<uint8_t*>(&val), reinterpret_cast<uint8_t*>(&val) + 4);
					}
				}
				else if (tip0 == 1 && tip1 == 0)
				{
					// First = stored value, rest = 0
					int32_t val = static_cast<int32_t>(values[val_idx++]);
					pck->v_data.insert(pck->v_data.end(), reinterpret_cast<uint8_t*>(&val), reinterpret_cast<uint8_t*>(&val) + 4);
					val = 0;
					for (uint32_t j = 1; j < elem_count; j++)
					{
						pck->v_data.insert(pck->v_data.end(), reinterpret_cast<uint8_t*>(&val), reinterpret_cast<uint8_t*>(&val) + 4);
					}
				}
				else
				{
					// Full array stored
					for (uint32_t j = 0; j < elem_count; j++)
					{
						int32_t val = static_cast<int32_t>(values[val_idx++]);
						pck->v_data.insert(pck->v_data.end(), reinterpret_cast<uint8_t*>(&val), reinterpret_cast<uint8_t*>(&val) + 4);
					}
				}
			}

			// Debug: verify AD data size (v_size stores element count, multiply by 4 for bytes)
			size_t expected_size = 0;
			for (size_t k = 0; k < pck->v_size.size(); k++) expected_size += pck->v_size[k] * 4;
			logger->info("AD decompress: num_samples={}, v_size.size()={}, expected_size={}, v_data.size()={}, val_idx={}, values.size()={}",
				num_samples, pck->v_size.size(), expected_size, pck->v_data.size(), val_idx, values.size());
			if (pck->v_data.size() != expected_size) {
				logger->error("AD data size mismatch! expected={}, got={}", expected_size, pck->v_data.size());
			}
		}
		else if (false && is_pl_field && keys[pck->key_id].type == BCF_HT_INT)  // DISABLED for testing
		{
			// Decompress PL with pattern recognition
			cbsc->Decompress(v_compressed, v_tmp);

			size_t offset = 0;
			uint64_t num_samples;
			offset += fmt_compress::VintCodec::decode(v_tmp.data() + offset, num_samples);

			// Read tips
			size_t tip_bytes = (num_samples * 2 + 7) / 8;
			std::vector<uint8_t> tips;
			tips.reserve(num_samples * 2);
			for (size_t i = 0; i < num_samples * 2; i++)
			{
				size_t byte_idx = i / 8;
				size_t bit_idx = i % 8;
				tips.push_back((v_tmp[offset + byte_idx] >> bit_idx) & 1);
			}
			offset += tip_bytes;

			// Read a_values
			uint64_t a_count;
			offset += fmt_compress::VintCodec::decode(v_tmp.data() + offset, a_count);
			std::vector<uint32_t> a_values(a_count);
			for (uint64_t i = 0; i < a_count; i++)
			{
				uint64_t val;
				offset += fmt_compress::VintCodec::decode(v_tmp.data() + offset, val);
				a_values[i] = static_cast<uint32_t>(val);
			}

			// Read b_values
			uint64_t b_count;
			offset += fmt_compress::VintCodec::decode(v_tmp.data() + offset, b_count);
			std::vector<uint32_t> b_values(b_count);
			for (uint64_t i = 0; i < b_count; i++)
			{
				uint64_t val;
				offset += fmt_compress::VintCodec::decode(v_tmp.data() + offset, val);
				b_values[i] = static_cast<uint32_t>(val);
			}

			// Read all_values
			uint64_t all_count;
			offset += fmt_compress::VintCodec::decode(v_tmp.data() + offset, all_count);
			std::vector<uint32_t> all_values(all_count);
			for (uint64_t i = 0; i < all_count; i++)
			{
				uint64_t val;
				offset += fmt_compress::VintCodec::decode(v_tmp.data() + offset, val);
				all_values[i] = static_cast<uint32_t>(val);
			}

			// Reconstruct PL data
			pck->v_data.clear();
			size_t a_idx = 0, b_idx = 0, all_idx = 0;
			for (size_t i = 0; i < num_samples; i++)
			{
				uint8_t tip0 = tips[i * 2];
				uint8_t tip1 = tips[i * 2 + 1];
				// v_size stores element count, not byte count (see WriteInt in bit_memory.h)
				uint32_t elem_count = pck->v_size[i];

				if (tip0 == 0 && tip1 == 0)
				{
					// All zeros
					for (uint32_t j = 0; j < elem_count; j++)
					{
						int32_t val = 0;
						pck->v_data.insert(pck->v_data.end(), reinterpret_cast<uint8_t*>(&val), reinterpret_cast<uint8_t*>(&val) + 4);
					}
				}
				else if (tip0 == 0 && tip1 == 1)
				{
					// Pattern [0, a, b, a, b, b, ...]
					uint32_t a = a_values[a_idx++];
					uint32_t b = b_values[b_idx++];
					// Generate pattern based on elem_count
					std::vector<int32_t> pattern;
					pattern.push_back(0);
					uint32_t b_repeat = 1;
					while (pattern.size() < elem_count)
					{
						pattern.push_back(static_cast<int32_t>(a));
						for (uint32_t k = 0; k < b_repeat && pattern.size() < elem_count; k++)
						{
							pattern.push_back(static_cast<int32_t>(b));
						}
						b_repeat++;
					}
					for (size_t j = 0; j < pattern.size(); j++)
					{
						pck->v_data.insert(pck->v_data.end(), reinterpret_cast<uint8_t*>(&pattern[j]), reinterpret_cast<uint8_t*>(&pattern[j]) + 4);
					}
				}
				else if (tip0 == 1 && tip1 == 0)
				{
					// Pattern with b = 15*a
					uint32_t a = a_values[a_idx++];
					uint32_t b = a * 15;
					std::vector<int32_t> pattern;
					pattern.push_back(0);
					uint32_t b_repeat = 1;
					while (pattern.size() < elem_count)
					{
						pattern.push_back(static_cast<int32_t>(a));
						for (uint32_t k = 0; k < b_repeat && pattern.size() < elem_count; k++)
						{
							pattern.push_back(static_cast<int32_t>(b));
						}
						b_repeat++;
					}
					for (size_t j = 0; j < pattern.size(); j++)
					{
						pck->v_data.insert(pck->v_data.end(), reinterpret_cast<uint8_t*>(&pattern[j]), reinterpret_cast<uint8_t*>(&pattern[j]) + 4);
					}
				}
				else
				{
					// Full array stored
					for (uint32_t j = 0; j < elem_count; j++)
					{
						int32_t val = static_cast<int32_t>(all_values[all_idx++]);
						pck->v_data.insert(pck->v_data.end(), reinterpret_cast<uint8_t*>(&val), reinterpret_cast<uint8_t*>(&val) + 4);
					}
				}
			}

			// Debug: verify PL data size (v_size stores element count, multiply by 4 for bytes)
			size_t expected_size = 0;
			for (size_t k = 0; k < pck->v_size.size(); k++) expected_size += pck->v_size[k] * 4;
			logger->info("PL decompress: num_samples={}, v_size.size()={}, expected_size={}, v_data.size()={}, a_idx={}, b_idx={}, all_idx={}",
				num_samples, pck->v_size.size(), expected_size, pck->v_data.size(), a_idx, b_idx, all_idx);
			if (pck->v_data.size() != expected_size) {
				logger->error("PL data size mismatch! expected={}, got={}", expected_size, pck->v_data.size());
			}
		}
		else if (keys[pck->key_id].type == BCF_HT_INT)
		{

			cbsc->Decompress(v_compressed, v_tmp);

			pck->v_data.resize(v_tmp.size());
			Decoder(v_tmp, pck->v_data);
			// std::cerr<<pck->key_id<<":"<<v_compressed.size()<<":"<<v_tmp.size()<<":"<<pck->v_data.size()<<endl;
			// // std::cerr<<pck->v_data.size()<<endl;
		}
		else
		{
			cbsc->Decompress(v_compressed, pck->v_data);
			// std::cerr<<pck->key_id<<":"<<v_compressed.size()<<":"<<pck->v_data.size()<<endl;
		}
	}
	else
	{
		pck->v_data.clear();
		// std::cerr<<pck->key_id<<":"<<v_compressed.size()<<":"<<pck->v_data.size()<<endl;
	}

	return true;
}
void DecompressionReader::Decoder(vector<uint8_t> &v_tmp, vector<uint8_t> &v_data)
{

	size_t size = v_tmp.size() / 4;
	for (size_t i = 0; i < v_data.size() / 4; i++)
	{
		v_data[i * 4] = v_tmp[i];
		v_data[i * 4 + 1] = v_tmp[i + size];
		v_data[i * 4 + 2] = v_tmp[i + size * 2];
		v_data[i * 4 + 3] = v_tmp[i + size * 3];
	}
}
//**********************************************************************************************************************
void DecompressionReader::InitDecompressParams()
{

	v_packages.resize(no_keys, nullptr);
	decomp_part_queue = new DecompressPartQueue<uint32_t>(1);
	v_coder_part_ids.resize(no_keys, 0);
	field_data_codecs.resize(no_keys);
	field_size_codecs.resize(no_keys);
	v_i_buf.resize(no_keys);
	for (uint32_t i = 0; i < no_keys; ++i)
	{
		field_size_codecs[i] = MakeCompressionStrategy(backend, p_bsc_size);
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
		field_data_codecs[i] = MakeCompressionStrategy(backend, param);
	}
}
//**********************************************************************************************************************
	void DecompressionReader::GetVariants(vector<field_desc> &fields)
	{

		// Load and set fields
		for (uint32_t i = 0; i < no_keys; i++)
		{
		if (v_i_buf[i].IsEmpty())
		{

			unique_lock<mutex> lck(m_packages);

			cv_packages.wait(lck, [&, this]
							 { return v_packages[i] != nullptr; });

			v_i_buf[i].SetBuffer(v_packages[i]->v_size, v_packages[i]->v_data);

			delete v_packages[i];

			v_packages[i] = nullptr;

			decomp_part_queue->PushQueue(i);
		}

		switch (keys[i].type)
		{
		case BCF_HT_INT:
			v_i_buf[i].ReadInt(fields[i].data, fields[i].data_size);
			fields[i].present = fields[i].data != nullptr;
			break;
		case BCF_HT_REAL:
			v_i_buf[i].ReadReal(fields[i].data, fields[i].data_size);
			fields[i].present = fields[i].data != nullptr;
			break;
		case BCF_HT_STR:
			v_i_buf[i].ReadText(fields[i].data, fields[i].data_size);
			fields[i].present = fields[i].data != nullptr;
			break;
		case BCF_HT_FLAG:
			uint8_t temp;
			v_i_buf[i].ReadFlag(temp);
			fields[i].present = (bool)temp;
			fields[i].data_size = 0;
			break;
			}
		}

			// Postprocess record-wise FORMAT streams: DP (needs AD), MIN_DP (needs DP), GQ (needs PL).
			int ad_idx = -1;
			int dp_idx = -1;
			int min_dp_idx = -1;
			int pl_idx = -1;
			int gq_idx = -1;
			for (uint32_t i = 0; i < no_keys; ++i)
			{
				if (keys[i].keys_type != key_type_t::fmt)
					continue;
				if (keys[i].name == "AD")
					ad_idx = (int)i;
				else if (keys[i].name == "DP")
					dp_idx = (int)i;
				else if (keys[i].name == "MIN_DP")
					min_dp_idx = (int)i;
				else if (keys[i].name == "PL")
					pl_idx = (int)i;
				else if (keys[i].name == "GQ")
					gq_idx = (int)i;
			}

			// DP
			if (dp_idx >= 0 && fields[dp_idx].present && fields[dp_idx].data && fields[dp_idx].data_size >= 4 && n_samples > 0 &&
				fields[dp_idx].data_size != (size_t)n_samples * 4)
			{
				auto logger = LogManager::Instance().Logger();
				const uint8_t *src = reinterpret_cast<const uint8_t *>(fields[dp_idx].data);
				const size_t src_size = fields[dp_idx].data_size;
				const uint8_t rec_codec = src[0];

				auto decode_vint = [&](size_t &off, uint64_t &val) -> bool {
					if (off >= src_size)
						return false;
					uint8_t len = fmt_compress::VintCodec::decode(src + off, val);
					if (len == 0 || off + len > src_size)
						return false;
					off += len;
					return true;
				};

				if (rec_codec == 0)
				{
					const size_t need = 1 + (size_t)n_samples * 4;
					if (src_size < need)
					{
						logger->error("DP record_codec=0 payload too short (need {}, have {})", need, src_size);
						abort();
					}
					char *out = new char[(size_t)n_samples * 4];
					memcpy(out, src + 1, (size_t)n_samples * 4);
					delete[] fields[dp_idx].data;
					fields[dp_idx].data = out;
					fields[dp_idx].data_size = (uint32_t)((size_t)n_samples * 4);
				}
				else if (rec_codec == 1)
				{
					if (ad_idx < 0 || !fields[ad_idx].present || !fields[ad_idx].data || (fields[ad_idx].data_size % 4 != 0))
					{
						logger->error("DP record_codec=1 requires AD to reconstruct");
						abort();
					}
					const uint32_t total_ad = fields[ad_idx].data_size / 4;
					if (total_ad == 0 || total_ad % n_samples != 0)
					{
						logger->error("AD layout invalid for DP reconstruction (total_ad={}, n_samples={})", total_ad, n_samples);
						abort();
					}
					const uint32_t ad_stride = total_ad / n_samples;
					const int32_t *ad_vals = reinterpret_cast<const int32_t *>(fields[ad_idx].data);

					vector<int32_t> dp_out(n_samples, bcf_int32_missing);
					for (uint32_t s = 0; s < n_samples; ++s)
					{
						const int32_t *sv = ad_vals + (size_t)s * ad_stride;
						bool ad_valid = true;
						int64_t sum = 0;
						for (uint32_t j = 0; j < ad_stride; ++j)
						{
							const int32_t v = sv[j];
							if (v == bcf_int32_missing || v == bcf_int32_vector_end || v < 0)
							{
								ad_valid = false;
								break;
							}
							sum += v;
							if (sum > INT32_MAX)
							{
								ad_valid = false;
								break;
							}
						}
						if (ad_valid)
							dp_out[s] = (int32_t)sum;
					}

					size_t off = 1;
					uint64_t raw_cnt = 0;
					if (!decode_vint(off, raw_cnt))
					{
						logger->error("DP decode: failed to read raw_cnt");
						abort();
					}
					uint32_t pos = 0;
					for (uint64_t i = 0; i < raw_cnt; ++i)
					{
						uint64_t delta = 0, val = 0;
						if (!decode_vint(off, delta) || !decode_vint(off, val))
						{
							logger->error("DP decode: failed in raw entry {}", i);
							abort();
						}
						pos += (uint32_t)delta;
						if (pos >= n_samples)
						{
							logger->error("DP decode: raw pos out of range {}", pos);
							abort();
						}
						dp_out[pos] = (val == 0) ? bcf_int32_missing : (int32_t)((uint32_t)val - 1);
					}

					uint64_t exc_cnt = 0;
					if (!decode_vint(off, exc_cnt))
					{
						logger->error("DP decode: failed to read exc_cnt");
						abort();
					}
					pos = 0;
					for (uint64_t i = 0; i < exc_cnt; ++i)
					{
						uint64_t delta = 0, val = 0;
						if (!decode_vint(off, delta) || !decode_vint(off, val))
						{
							logger->error("DP decode: failed in exc entry {}", i);
							abort();
						}
						pos += (uint32_t)delta;
						if (pos >= n_samples)
						{
							logger->error("DP decode: exc pos out of range {}", pos);
							abort();
						}
						dp_out[pos] = (val == 0) ? bcf_int32_missing : (int32_t)((uint32_t)val - 1);
					}

					for (size_t i = off; i < src_size; ++i)
					{
						if (src[i] != 0)
						{
							logger->error("DP decode: unexpected trailing bytes (off={}, size={})", off, src_size);
							abort();
						}
					}

					char *out = new char[(size_t)n_samples * 4];
					memcpy(out, dp_out.data(), (size_t)n_samples * 4);
					delete[] fields[dp_idx].data;
					fields[dp_idx].data = out;
					fields[dp_idx].data_size = (uint32_t)((size_t)n_samples * 4);
				}
				else
				{
					logger->error("DP unknown record_codec={}", (int)rec_codec);
					abort();
				}
			}

			// MIN_DP
			if (min_dp_idx >= 0 && fields[min_dp_idx].present && fields[min_dp_idx].data && fields[min_dp_idx].data_size >= 4 && n_samples > 0 &&
				fields[min_dp_idx].data_size != (size_t)n_samples * 4)
			{
				auto logger = LogManager::Instance().Logger();
				const uint8_t *src = reinterpret_cast<const uint8_t *>(fields[min_dp_idx].data);
				const size_t src_size = fields[min_dp_idx].data_size;
				const uint8_t rec_codec = src[0];

				auto decode_vint = [&](size_t &off, uint64_t &val) -> bool {
					if (off >= src_size)
						return false;
					uint8_t len = fmt_compress::VintCodec::decode(src + off, val);
					if (len == 0 || off + len > src_size)
						return false;
					off += len;
					return true;
				};

				if (rec_codec == 0)
				{
					const size_t need = 1 + (size_t)n_samples * 4;
					if (src_size < need)
					{
						logger->error("MIN_DP record_codec=0 payload too short (need {}, have {})", need, src_size);
						abort();
					}
					char *out = new char[(size_t)n_samples * 4];
					memcpy(out, src + 1, (size_t)n_samples * 4);
					delete[] fields[min_dp_idx].data;
					fields[min_dp_idx].data = out;
					fields[min_dp_idx].data_size = (uint32_t)((size_t)n_samples * 4);
				}
				else if (rec_codec == 1)
				{
					if (dp_idx < 0 || !fields[dp_idx].present || !fields[dp_idx].data || fields[dp_idx].data_size != (size_t)n_samples * 4)
					{
						logger->error("MIN_DP record_codec=1 requires DP to reconstruct");
						abort();
					}
					vector<int32_t> out(n_samples, bcf_int32_missing);
					memcpy(out.data(), fields[dp_idx].data, (size_t)n_samples * 4);

					size_t off = 1;
					uint64_t exc_cnt = 0;
					if (!decode_vint(off, exc_cnt))
					{
						logger->error("MIN_DP decode: failed to read exc_cnt");
						abort();
					}
					uint32_t pos = 0;
					for (uint64_t i = 0; i < exc_cnt; ++i)
					{
						uint64_t delta = 0, val = 0;
						if (!decode_vint(off, delta) || !decode_vint(off, val))
						{
							logger->error("MIN_DP decode: failed in exc entry {}", i);
							abort();
						}
						pos += (uint32_t)delta;
						if (pos >= n_samples)
						{
							logger->error("MIN_DP decode: pos out of range {}", pos);
							abort();
						}
						out[pos] = (val == 0) ? bcf_int32_missing : (int32_t)((uint32_t)val - 1);
					}

					for (size_t i = off; i < src_size; ++i)
					{
						if (src[i] != 0)
						{
							logger->error("MIN_DP decode: unexpected trailing bytes (off={}, size={})", off, src_size);
							abort();
						}
					}

					char *out_bytes = new char[(size_t)n_samples * 4];
					memcpy(out_bytes, out.data(), (size_t)n_samples * 4);
					delete[] fields[min_dp_idx].data;
					fields[min_dp_idx].data = out_bytes;
					fields[min_dp_idx].data_size = (uint32_t)((size_t)n_samples * 4);
				}
				else
				{
					logger->error("MIN_DP unknown record_codec={}", (int)rec_codec);
					abort();
				}
			}

			// GQ
			if (gq_idx >= 0 && fields[gq_idx].present && fields[gq_idx].data && fields[gq_idx].data_size >= 4 && n_samples > 0 &&
				fields[gq_idx].data_size != (size_t)n_samples * 4)
			{
				auto logger = LogManager::Instance().Logger();
				const uint8_t *src = reinterpret_cast<const uint8_t *>(fields[gq_idx].data);
				const size_t src_size = fields[gq_idx].data_size;
				const uint8_t rec_codec = src[0];

				auto decode_vint = [&](size_t &off, uint64_t &val) -> bool {
					if (off >= src_size)
						return false;
					uint8_t len = fmt_compress::VintCodec::decode(src + off, val);
					if (len == 0 || off + len > src_size)
						return false;
					off += len;
					return true;
				};

				if (rec_codec == 0)
				{
					const size_t need = 1 + (size_t)n_samples * 4;
					if (src_size < need)
					{
						logger->error("GQ record_codec=0 payload too short (need {}, have {})", need, src_size);
						abort();
					}
					char *out = new char[(size_t)n_samples * 4];
					memcpy(out, src + 1, (size_t)n_samples * 4);
					delete[] fields[gq_idx].data;
					fields[gq_idx].data = out;
					fields[gq_idx].data_size = (uint32_t)((size_t)n_samples * 4);
				}
				else if (rec_codec == 1)
				{
					vector<int32_t> out(n_samples, bcf_int32_missing);
					vector<uint8_t> is_raw(n_samples, 0);

					size_t off = 1;
					uint64_t raw_cnt = 0;
					if (!decode_vint(off, raw_cnt))
					{
						logger->error("GQ decode: failed to read raw_cnt");
						abort();
					}
					uint32_t pos = 0;
					for (uint64_t i = 0; i < raw_cnt; ++i)
					{
						uint64_t delta = 0, val = 0;
						if (!decode_vint(off, delta) || !decode_vint(off, val))
						{
							logger->error("GQ decode: failed in raw entry {}", i);
							abort();
						}
						pos += (uint32_t)delta;
						if (pos >= n_samples)
						{
							logger->error("GQ decode: raw pos out of range {}", pos);
							abort();
						}
						is_raw[pos] = 1;
						out[pos] = (val == 0) ? bcf_int32_missing : (int32_t)((uint32_t)val - 1);
					}

					bool has_pl_layout = false;
					uint32_t pl_stride = 0;
					const int32_t *pl_vals = nullptr;
					if (pl_idx >= 0 && fields[pl_idx].present && fields[pl_idx].data && (fields[pl_idx].data_size % 4 == 0))
					{
						const uint32_t total_pl = fields[pl_idx].data_size / 4;
						if (total_pl > 0 && total_pl % n_samples == 0)
						{
							pl_stride = total_pl / n_samples;
							if (pl_stride > 0)
							{
								has_pl_layout = true;
								pl_vals = reinterpret_cast<const int32_t *>(fields[pl_idx].data);
							}
						}
					}

					vector<uint32_t> pl_vec;
					if (has_pl_layout)
						pl_vec.resize(pl_stride);
					for (uint32_t s = 0; s < n_samples; ++s)
					{
						if (is_raw[s] || !has_pl_layout)
							continue;

						const int32_t *sv = pl_vals + (size_t)s * pl_stride;
						bool pl_valid = true;
						for (uint32_t j = 0; j < pl_stride; ++j)
						{
							const int32_t v = sv[j];
							if (v == bcf_int32_missing || v == bcf_int32_vector_end || v < 0)
							{
								pl_valid = false;
								break;
							}
							pl_vec[j] = (uint32_t)v;
						}
						if (!pl_valid)
							continue;

						uint32_t a_val = 0, b_val = 0;
						const uint8_t type = fmt_compress::checkPlPattern(pl_vec.data(), pl_stride, a_val, b_val);
						if (type == 1 || type == 2)
							out[s] = (int32_t)a_val;
						else if (type == 3)
							out[s] = (int32_t)fmt_compress::getSecondMin(pl_vec.data(), pl_stride);
						else
							out[s] = bcf_int32_missing;
					}

					uint64_t exc_cnt = 0;
					if (!decode_vint(off, exc_cnt))
					{
						logger->error("GQ decode: failed to read exc_cnt");
						abort();
					}
					pos = 0;
					for (uint64_t i = 0; i < exc_cnt; ++i)
					{
						uint64_t delta = 0, val = 0;
						if (!decode_vint(off, delta) || !decode_vint(off, val))
						{
							logger->error("GQ decode: failed in exc entry {}", i);
							abort();
						}
						pos += (uint32_t)delta;
						if (pos >= n_samples)
						{
							logger->error("GQ decode: exc pos out of range {}", pos);
							abort();
						}
						out[pos] = (val == 0) ? bcf_int32_missing : (int32_t)((uint32_t)val - 1);
					}

					for (size_t i = off; i < src_size; ++i)
					{
						if (src[i] != 0)
						{
							logger->error("GQ decode: unexpected trailing bytes (off={}, size={})", off, src_size);
							abort();
						}
					}

					char *out_bytes = new char[(size_t)n_samples * 4];
					memcpy(out_bytes, out.data(), (size_t)n_samples * 4);
					delete[] fields[gq_idx].data;
					fields[gq_idx].data = out_bytes;
					fields[gq_idx].data_size = (uint32_t)((size_t)n_samples * 4);
				}
				else
				{
					logger->error("GQ unknown record_codec={}", (int)rec_codec);
					abort();
				}
			}
		}
//**********************************************************************************************************************
void DecompressionReader::decompress_meta(vector<string> &v_samples, string &header)
{
	auto logger = LogManager::Instance().Logger();
	vector<uint8_t> all_v_header;
	vector<uint8_t> all_v_samples;
	size_t p_header;
	size_t p_samples;

	for (auto data : {
			 make_tuple(ref(all_v_header), ref(comp_v_header), ref(p_header), "header"),
			 make_tuple(ref(all_v_samples), ref(comp_v_samples), ref(p_samples), "samples"),
		 })
	{
		auto codec = MakeCompressionStrategy(backend, p_bsc_meta);
		codec->Decompress(get<1>(data), get<0>(data));
		get<2>(data) = 0;
	}
	header.clear();
	read_str(all_v_header, p_header, header);
	v_samples.clear();
	string sample;
	for (uint32_t i = 0; i < n_samples; ++i)
	{
		read_str(all_v_samples, p_samples, sample);
		v_samples.emplace_back(sample);
	}
	logger->info("Decompress header and samples done.");
}
bool DecompressionReader::setStartChunk(uint32_t start_chunk_id)
{

	auto it = chunks_streams.find(start_chunk_id + 1);

	if (it == chunks_streams.end())
		return false;
	buf_pos = it->second.offset;

	return true;
}
// get current chunk actual position
uint32_t DecompressionReader::getActualPos(uint32_t chunk_id)
{
	auto it = chunks_streams.find(chunk_id);
	if (it == chunks_streams.end())
		return 0;

	return it->second.cur_chunk_actual_pos;
}
//**********************************************************************************************************************
// read current chunk fixed fields compressed data
bool DecompressionReader::readFixedFields()
{
	// Reset per-chunk state
	has_fixed_fields_rb_dir = false;
	fixed_fields_rb_dir.clear();
	fixed_fields_chunk_start = 0;
	fixed_fields_chunk_version = 0;
	fixed_fields_total_variants = 0;
	fixed_fields_row_block_count = 0;
	fixed_fields_gt_off = 0;
	fixed_fields_gt_size = 0;

	const uint64_t chunk_start = buf_pos;
	uint32_t magic_or_no_variants = 0;
	memcpy(&magic_or_no_variants, buf + buf_pos, sizeof(uint32_t));

	// New per-row_block fixed fields directory format
	if (magic_or_no_variants == GSC_FIXED_FIELDS_RB_MAGIC)
	{
		buf_pos += sizeof(uint32_t);
		uint32_t version = 0;
		uint32_t flags = 0;
		memcpy(&version, buf + buf_pos, sizeof(uint32_t));
		buf_pos += sizeof(uint32_t);
		memcpy(&fixed_fields_total_variants, buf + buf_pos, sizeof(uint32_t));
		buf_pos += sizeof(uint32_t);
		memcpy(&fixed_fields_row_block_count, buf + buf_pos, sizeof(uint32_t));
		buf_pos += sizeof(uint32_t);
		memcpy(&flags, buf + buf_pos, sizeof(uint32_t));
		buf_pos += sizeof(uint32_t);

		fixed_fields_chunk_version = version;
		if (version == GSC_FIXED_FIELDS_RB_VERSION_V1)
		{
			memcpy(&fixed_fields_gt_off, buf + buf_pos, sizeof(uint32_t));
			buf_pos += sizeof(uint32_t);
			memcpy(&fixed_fields_gt_size, buf + buf_pos, sizeof(uint32_t));
			buf_pos += sizeof(uint32_t);
		}
		else if (version == GSC_FIXED_FIELDS_RB_VERSION_V2)
		{
			fixed_fields_gt_off = 0;
			fixed_fields_gt_size = 0;
		}
		else
		{
			auto logger = LogManager::Instance().Logger();
			logger->error("Unsupported fixed-fields chunk version: {}", version);
			return false;
		}

		fixed_fields_chunk_start = chunk_start;
		has_fixed_fields_rb_dir = true;
		fixed_fields_rb_dir.resize(fixed_fields_row_block_count);

		uint32_t max_end = 0;
		uint32_t header_bytes = 0;
		uint32_t entry_bytes = 0;
		if (version == GSC_FIXED_FIELDS_RB_VERSION_V1)
		{
			header_bytes = 7u * sizeof(uint32_t);
			entry_bytes = sizeof(uint32_t) + 2u * sizeof(int64_t) + 6u * (2u * sizeof(uint32_t));
		}
		else
		{
			header_bytes = 5u * sizeof(uint32_t);
			entry_bytes = sizeof(uint32_t) + 2u * sizeof(int64_t) + 6u * (2u * sizeof(uint32_t)) + 2u * sizeof(uint32_t);
		}
		const uint32_t base_end = header_bytes + fixed_fields_row_block_count * entry_bytes;
		max_end = base_end;
		if (version == GSC_FIXED_FIELDS_RB_VERSION_V1)
			max_end = std::max<uint32_t>(max_end, fixed_fields_gt_off + fixed_fields_gt_size);

		for (uint32_t i = 0; i < fixed_fields_row_block_count; ++i)
		{
			fixed_fields_rb_dir_entry e;
			memcpy(&e.meta.variant_count, buf + buf_pos, sizeof(uint32_t));
			buf_pos += sizeof(uint32_t);
			memcpy(&e.meta.first_pos, buf + buf_pos, sizeof(int64_t));
			buf_pos += sizeof(int64_t);
			memcpy(&e.meta.last_pos, buf + buf_pos, sizeof(int64_t));
			buf_pos += sizeof(int64_t);

			memcpy(&e.chrom_off, buf + buf_pos, sizeof(uint32_t));
			buf_pos += sizeof(uint32_t);
			memcpy(&e.chrom_size, buf + buf_pos, sizeof(uint32_t));
			buf_pos += sizeof(uint32_t);
			memcpy(&e.pos_off, buf + buf_pos, sizeof(uint32_t));
			buf_pos += sizeof(uint32_t);
			memcpy(&e.pos_size, buf + buf_pos, sizeof(uint32_t));
			buf_pos += sizeof(uint32_t);
			memcpy(&e.id_off, buf + buf_pos, sizeof(uint32_t));
			buf_pos += sizeof(uint32_t);
			memcpy(&e.id_size, buf + buf_pos, sizeof(uint32_t));
			buf_pos += sizeof(uint32_t);
			memcpy(&e.ref_off, buf + buf_pos, sizeof(uint32_t));
			buf_pos += sizeof(uint32_t);
			memcpy(&e.ref_size, buf + buf_pos, sizeof(uint32_t));
			buf_pos += sizeof(uint32_t);
			memcpy(&e.alt_off, buf + buf_pos, sizeof(uint32_t));
			buf_pos += sizeof(uint32_t);
			memcpy(&e.alt_size, buf + buf_pos, sizeof(uint32_t));
			buf_pos += sizeof(uint32_t);
			memcpy(&e.qual_off, buf + buf_pos, sizeof(uint32_t));
			buf_pos += sizeof(uint32_t);
			memcpy(&e.qual_size, buf + buf_pos, sizeof(uint32_t));
			buf_pos += sizeof(uint32_t);

			if (version == GSC_FIXED_FIELDS_RB_VERSION_V2)
			{
				memcpy(&e.gt_off, buf + buf_pos, sizeof(uint32_t));
				buf_pos += sizeof(uint32_t);
				memcpy(&e.gt_size, buf + buf_pos, sizeof(uint32_t));
				buf_pos += sizeof(uint32_t);
			}

			max_end = std::max<uint32_t>(max_end, e.chrom_off + e.chrom_size);
			max_end = std::max<uint32_t>(max_end, e.pos_off + e.pos_size);
			max_end = std::max<uint32_t>(max_end, e.id_off + e.id_size);
			max_end = std::max<uint32_t>(max_end, e.ref_off + e.ref_size);
			max_end = std::max<uint32_t>(max_end, e.alt_off + e.alt_size);
			max_end = std::max<uint32_t>(max_end, e.qual_off + e.qual_size);
			if (version == GSC_FIXED_FIELDS_RB_VERSION_V2)
				max_end = std::max<uint32_t>(max_end, e.gt_off + e.gt_size);

			fixed_fields_rb_dir[i] = std::move(e);
		}

		buf_pos = chunk_start + max_end;
		return true;
	}

	// Legacy chunk-level fixed fields format
	memcpy(&fixed_field_block_compress.no_variants, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	uint32_t comp_size;

	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	fixed_field_block_compress.chrom.resize(comp_size);
	memcpy(&fixed_field_block_compress.chrom[0], buf + buf_pos, comp_size * sizeof(uint8_t));
	buf_pos = buf_pos + comp_size * sizeof(uint8_t);

	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	fixed_field_block_compress.pos.resize(comp_size);
	memcpy(&fixed_field_block_compress.pos[0], buf + buf_pos, comp_size * sizeof(uint8_t));
	buf_pos = buf_pos + comp_size * sizeof(uint8_t);

	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	fixed_field_block_compress.id.resize(comp_size);
	memcpy(&fixed_field_block_compress.id[0], buf + buf_pos, comp_size * sizeof(uint8_t));
	buf_pos = buf_pos + comp_size * sizeof(uint8_t);

	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	fixed_field_block_compress.ref.resize(comp_size);
	memcpy(&fixed_field_block_compress.ref[0], buf + buf_pos, comp_size * sizeof(uint8_t));
	buf_pos = buf_pos + comp_size * sizeof(uint8_t);

	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	fixed_field_block_compress.alt.resize(comp_size);
	memcpy(&fixed_field_block_compress.alt[0], buf + buf_pos, comp_size * sizeof(uint8_t));
	buf_pos = buf_pos + comp_size * sizeof(uint8_t);

	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	fixed_field_block_compress.qual.resize(comp_size);
	memcpy(&fixed_field_block_compress.qual[0], buf + buf_pos, comp_size * sizeof(uint8_t));
	buf_pos = buf_pos + comp_size * sizeof(uint8_t);

	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	fixed_field_block_compress.gt_block.resize(comp_size);
	memcpy(&fixed_field_block_compress.gt_block[0], buf + buf_pos, comp_size * sizeof(uint8_t));
	buf_pos = buf_pos + comp_size * sizeof(uint8_t);

	return true;
}
void DecompressionReader::initDecoderParams()
{
	p_chrom = 0;
	p_pos = 0;
	p_id = 0;
	p_ref = 0;
	p_alt = 0;
	p_gt = 0;
	p_qual = 0;
	prev_pos = 0;
	fixed_field_block_io.Clear();
	fixed_field_block_io.Initalize();
	fixed_field_block_io.Clear();
	fixed_field_block_io.Initalize();
}
bool DecompressionReader::Decoder(std::vector<block_t> &v_blocks, std::vector<std::vector<std::vector<uint32_t>>> &s_perm, std::vector<uint8_t> &gt_index, uint32_t cur_chunk_id)
{
	if (has_fixed_fields_rb_dir)
	{
		uint32_t variants_before = 0;
		return DecoderByRange(v_blocks, s_perm, gt_index, cur_chunk_id,
							  std::numeric_limits<int64_t>::min(),
							  std::numeric_limits<int64_t>::max(),
							  variants_before);
	}
	// Initialize decoder parameters
	initDecoderParams();
	no_variants = fixed_field_block_compress.no_variants;

	// Other initializations here...
	std::vector<variant_desc_t> v_vcf_fixed_data_io;
	std::vector<variant_desc_t> v_vcf_sort_data_io;
	uint32_t row_block_variants = max_block_rows ? max_block_rows : total_haplotypes;
	uint32_t block_size = useLegacyPath ? (n_samples * static_cast<uint32_t>(ploidy)) : row_block_variants;
	uint32_t global_block_start = 0;
	if (!useLegacyPath && cur_chunk_id < chunk_block_offsets.size())
		global_block_start = chunk_block_offsets[cur_chunk_id];
	v_vcf_fixed_data_io.reserve(no_variants_in_buf);
	v_vcf_sort_data_io.reserve(no_variants_in_buf);
	std::vector<uint32_t> perm(block_size, 0);

	// Start threads using std::async for better exception handling and cleaner code
	auto fixedFieldHandle = std::async(std::launch::async, [&]()
									   {
        // Code from the original fixed field thread
        // Decompression and processing...
		auto codec = MakeCompressionStrategy(backend, p_bsc_fixed_fields);
		for (auto data : {
        	make_tuple(ref(fixed_field_block_io.chrom), ref(fixed_field_block_compress.chrom), p_chrom, "chrom"),
        	make_tuple(ref(fixed_field_block_io.id), ref(fixed_field_block_compress.id), p_id, "id"),
        	make_tuple(ref(fixed_field_block_io.alt), ref(fixed_field_block_compress.alt), p_alt, "alt"),
    		make_tuple(ref(fixed_field_block_io.qual), ref(fixed_field_block_compress.qual), p_qual, "qual"),


    	})
		{
			codec->Decompress(get<1>(data), get<0>(data));
			// std::cerr<<get<3>(data)<<":"<<get<1>(data).size()<<":"<<get<0>(data).size()<<endl;

		}
		uint32_t i_variant;
		variant_desc_t desc;
		for(i_variant = 0; i_variant < no_variants; i_variant++) 	
		{
			// Load variant description
			read_str(fixed_field_block_io.chrom, p_chrom, desc.chrom);
			read_str(fixed_field_block_io.id, p_id, desc.id);
			read_str(fixed_field_block_io.alt, p_alt, desc.alt);
			read_str(fixed_field_block_io.qual, p_qual, desc.qual);
			v_vcf_fixed_data_io.emplace_back(desc);
		} });

	auto sortFieldHandle = std::async(std::launch::async, [&]()
									  {
        // Code from the original sorted field thread
        // Decompression and processing...
		auto codec = MakeCompressionStrategy(backend, p_bsc_fixed_fields);
		for (auto data : {
			
        	make_tuple(ref(fixed_field_block_io.pos), ref(fixed_field_block_compress.pos), p_pos, "pos"),
        	make_tuple(ref(fixed_field_block_io.ref), ref(fixed_field_block_compress.ref), p_ref, "ref"),
        	// make_tuple(ref(fixed_field_block_io.gt_block), ref(fixed_field_block_compress.gt_block), p_gt, "GT")

    	})
		{
			codec->Decompress(get<1>(data), get<0>(data));
			// std::cerr<<get<3>(data)<<":"<<get<1>(data).size()<<":"<<get<0>(data).size()<<endl;

		}
		uint8_t marker = fixed_field_block_compress.gt_block.back();
		fixed_field_block_compress.gt_block.pop_back();
		switch (marker) {
		case 0: {
			auto gt_codec = MakeCompressionStrategy(compression_backend_t::bsc, p_bsc_fixed_fields);
			gt_codec->Decompress(fixed_field_block_compress.gt_block, fixed_field_block_io.gt_block);
			break;
		}
		case 1:
			zstd::zstd_decompress(fixed_field_block_compress.gt_block, fixed_field_block_io.gt_block);
			break;
		case 2: {
			auto gt_codec = MakeCompressionStrategy(compression_backend_t::brotli, p_bsc_fixed_fields);
			gt_codec->Decompress(fixed_field_block_compress.gt_block, fixed_field_block_io.gt_block);
			break;
		}
		default:
			{
				auto gt_codec = MakeCompressionStrategy(backend, p_bsc_fixed_fields);
				gt_codec->Decompress(fixed_field_block_compress.gt_block, fixed_field_block_io.gt_block);
			}
			break;
		}
		uint32_t i_variant;
		uint32_t block_id_local = 0;
		variant_desc_t desc;
		for (i_variant = 0; i_variant < no_variants; i_variant++)
		{
			int64_t pos = 0;
			// Load variant description
			read(fixed_field_block_io.pos, p_pos, pos);
			pos += prev_pos;
			prev_pos = pos;
			desc.pos = pos;
			read_str(fixed_field_block_io.ref, p_ref, desc.ref);
			desc.filter = ".";
			desc.info = ".";
			v_vcf_sort_data_io.emplace_back(desc);
			if ((i_variant + 1) % block_size == 0)
			{
				if (useLegacyPath)
				{
					out_perm(perm, v_vcf_sort_data_io);
					// For 3D structure: wrap perm in a 2D vector (legacy: single column block)
					vector<vector<uint32_t>> row_perm;
					row_perm.emplace_back(perm);
					s_perm.emplace_back(row_perm);
				}
				else
				{
					vector<vector<uint32_t>> row_perm(n_col_blocks);
					for (uint32_t cb = 0; cb < n_col_blocks; ++cb)
					{
						auto it = vint_last_perm_2d.find(make_pair(global_block_start + block_id_local, cb));
						if (it == vint_last_perm_2d.end())
						{
							row_perm[cb].resize(col_block_ranges[cb].second);
							for (size_t i_p = 0; i_p < row_perm[cb].size(); ++i_p)
								row_perm[cb][i_p] = static_cast<uint32_t>(i_p);
						}
						else
						{
							row_perm[cb] = vint_code::DecodeArray(it->second);
						}
					}
					s_perm.emplace_back(row_perm);
				}
				if (!useLegacyPath)
				{
					for (size_t i_p = 0; i_p < v_vcf_sort_data_io.size(); i_p++)
					{
						if (atoi(v_vcf_sort_data_io[i_p].ref.c_str()))
						{
							v_vcf_sort_data_io[i_p].ref = v_vcf_sort_data_io[i_p].ref.substr(to_string(atoi(v_vcf_sort_data_io[i_p].ref.c_str())).length());
						}
					}
				}
				v_blocks.emplace_back(v_vcf_sort_data_io);
				v_vcf_sort_data_io.clear();
				block_id_local++;
			}
		}
		if (i_variant % block_size)
		{
			if (useLegacyPath)
			{
				auto it = vint_last_perm.find(cur_chunk_id);
				perm = vint_code::DecodeArray(it->second);
				vector<vector<uint32_t>> row_perm;
				row_perm.emplace_back(perm);
				s_perm.emplace_back(row_perm);
			}
			else
			{
				vector<vector<uint32_t>> row_perm(n_col_blocks);
				for (uint32_t cb = 0; cb < n_col_blocks; ++cb)
				{
					auto it = vint_last_perm_2d.find(make_pair(global_block_start + block_id_local, cb));
					if (it == vint_last_perm_2d.end())
					{
						row_perm[cb].resize(col_block_ranges[cb].second);
						for (size_t i_p = 0; i_p < row_perm[cb].size(); ++i_p)
							row_perm[cb][i_p] = static_cast<uint32_t>(i_p);
					}
					else
					{
						row_perm[cb] = vint_code::DecodeArray(it->second);
					}
				}
				s_perm.emplace_back(row_perm);
			}
			for (size_t i_p = 0; i_p < v_vcf_sort_data_io.size(); i_p++)
			{
				if (atoi(v_vcf_sort_data_io[i_p].ref.c_str()))
				{
					v_vcf_sort_data_io[i_p].ref = v_vcf_sort_data_io[i_p].ref.substr(to_string(atoi(v_vcf_sort_data_io[i_p].ref.c_str())).length());
				}
			}
			v_blocks.emplace_back(v_vcf_sort_data_io);
			v_vcf_sort_data_io.clear();
		} });

	// Wait for both tasks to complete
	fixedFieldHandle.wait();
	sortFieldHandle.wait();

	// Merge results
	uint32_t i_variant = 0;

	for (size_t i = 0; i < v_blocks.size(); i++)
	{
		for (size_t j = 0; j < v_blocks[i].data_compress.size(); j++)
		{
			if (i_variant >= v_vcf_fixed_data_io.size()) {
				break;
			}
			v_blocks[i].data_compress[j].chrom = std::move(v_vcf_fixed_data_io[i_variant].chrom);
			v_blocks[i].data_compress[j].id = std::move(v_vcf_fixed_data_io[i_variant].id);
			v_blocks[i].data_compress[j].alt = std::move(v_vcf_fixed_data_io[i_variant].alt);
			v_blocks[i].data_compress[j].qual = std::move(v_vcf_fixed_data_io[i_variant].qual);

			i_variant++;
		}
	}
	// v_vcf_fixed_data_io.clear();
	gt_index = move(fixed_field_block_io.gt_block);

	// Return true or appropriate status
	return true;
}

bool DecompressionReader::DecoderByRange(vector<block_t> &v_blocks, vector<vector<vector<uint32_t>>> &s_perm,
										vector<uint8_t> &gt_index, uint32_t cur_chunk_id,
										int64_t range_1, int64_t range_2, uint32_t &variants_before)
{
	if (!has_fixed_fields_rb_dir)
		return false;

	initDecoderParams();
	no_variants = fixed_fields_total_variants;
	variants_before = 0;

	// Select row_blocks overlapping [range_1, range_2]
	bool any = false;
	uint32_t first_rb = 0;
	uint32_t last_rb = 0;
	uint32_t prefix_variants = 0;
	for (uint32_t rb = 0; rb < fixed_fields_row_block_count; ++rb)
	{
		const auto &meta = fixed_fields_rb_dir[rb].meta;
		if (meta.variant_count == 0)
			continue;

		if (meta.last_pos < range_1)
		{
			prefix_variants += meta.variant_count;
			continue;
		}
		if (meta.first_pos > range_2)
			break;
		if (!any)
		{
			first_rb = rb;
			any = true;
		}
		last_rb = rb;
	}
	variants_before = prefix_variants;

	v_blocks.clear();
	s_perm.clear();
	gt_index.clear();
	if (!any)
		return true;

	// Pre-create codec objects for reuse (avoid repeated allocation)
	auto gt_codec_bsc = MakeCompressionStrategy(compression_backend_t::bsc, p_bsc_fixed_fields);
	auto gt_codec_brotli = MakeCompressionStrategy(compression_backend_t::brotli, p_bsc_fixed_fields);
	auto gt_codec_backend = MakeCompressionStrategy(backend, p_bsc_fixed_fields);

	// Lambda: decompress GT segment directly from pointer (no intermediate memcpy)
	auto decompress_gt_segment = [&](const uint8_t *comp_ptr, uint32_t comp_size, std::vector<uint8_t> &out)
	{
		out.clear();
		if (!comp_size)
			return;
		uint8_t marker = comp_ptr[comp_size - 1];
		const uint8_t *data_ptr = comp_ptr;
		size_t data_size = comp_size - 1;
		switch (marker)
		{
		case 0:
			gt_codec_bsc->DecompressFromPtr(data_ptr, data_size, out);
			break;
		case 1:
			zstd::zstd_decompress_ptr(data_ptr, data_size, out);
			break;
		case 2:
			gt_codec_brotli->DecompressFromPtr(data_ptr, data_size, out);
			break;
		default:
			gt_codec_backend->DecompressFromPtr(data_ptr, data_size, out);
			break;
		}
	};

	// Decode GT index (v1: chunk-level; v2+: per-row_block, only for selected row_blocks)
	if (fixed_fields_chunk_version == GSC_FIXED_FIELDS_RB_VERSION_V1)
	{
		decompress_gt_segment(buf + fixed_fields_chunk_start + fixed_fields_gt_off, fixed_fields_gt_size, gt_index);
	}
	else if (fixed_fields_chunk_version == GSC_FIXED_FIELDS_RB_VERSION_V2)
	{
		// Pre-calculate total GT size for reservation (avoid repeated reallocation)
		size_t estimated_total_size = 0;
		for (uint32_t rb = first_rb; rb <= last_rb; ++rb)
		{
			// Estimate decompressed size as ~4x compressed size (conservative)
			estimated_total_size += fixed_fields_rb_dir[rb].gt_size * 4;
		}
		gt_index.reserve(estimated_total_size);

		std::vector<uint8_t> seg_out;
		for (uint32_t rb = first_rb; rb <= last_rb; ++rb)
		{
			const auto &dir = fixed_fields_rb_dir[rb];
			if (!dir.gt_size)
				continue;
			decompress_gt_segment(buf + fixed_fields_chunk_start + dir.gt_off, dir.gt_size, seg_out);
			if (!seg_out.empty())
				gt_index.insert(gt_index.end(), seg_out.begin(), seg_out.end());
		}
	}
	else
	{
		auto logger = LogManager::Instance().Logger();
		logger->error("DecoderByRange: unsupported fixed-fields chunk version {}", fixed_fields_chunk_version);
		return false;
	}

	// Reuse gt_codec_backend for fixed fields (same backend)
	auto &codec = gt_codec_backend;
	const uint32_t row_block_variants = max_block_rows ? max_block_rows : total_haplotypes;
	const uint32_t block_size = useLegacyPath ? (n_samples * static_cast<uint32_t>(ploidy)) : row_block_variants;
	uint32_t global_block_start = 0;
	if (!useLegacyPath && cur_chunk_id < chunk_block_offsets.size())
		global_block_start = chunk_block_offsets[cur_chunk_id];

	for (uint32_t rb = first_rb; rb <= last_rb; ++rb)
	{
		const auto &dir = fixed_fields_rb_dir[rb];
		const uint32_t vcount = dir.meta.variant_count;
		if (!vcount)
			continue;

		// Decompress fixed fields for this row_block (directly from pointer, no memcpy)
		std::vector<uint8_t> out_chrom, out_id, out_alt, out_qual, out_pos, out_ref;
		const uint8_t *base_ptr = buf + fixed_fields_chunk_start;

		codec->DecompressFromPtr(base_ptr + dir.chrom_off, dir.chrom_size, out_chrom);
		codec->DecompressFromPtr(base_ptr + dir.id_off, dir.id_size, out_id);
		codec->DecompressFromPtr(base_ptr + dir.alt_off, dir.alt_size, out_alt);
		codec->DecompressFromPtr(base_ptr + dir.qual_off, dir.qual_size, out_qual);
		codec->DecompressFromPtr(base_ptr + dir.pos_off, dir.pos_size, out_pos);
		codec->DecompressFromPtr(base_ptr + dir.ref_off, dir.ref_size, out_ref);

		std::vector<variant_desc_t> fixed_desc;
		std::vector<variant_desc_t> sort_desc;
		fixed_desc.resize(vcount);
		sort_desc.resize(vcount);

		size_t p_chrom_l = 0, p_id_l = 0, p_alt_l = 0, p_qual_l = 0;
		for (uint32_t i = 0; i < vcount; ++i)
		{
			read_str(out_chrom, p_chrom_l, fixed_desc[i].chrom);
			read_str(out_id, p_id_l, fixed_desc[i].id);
			read_str(out_alt, p_alt_l, fixed_desc[i].alt);
			read_str(out_qual, p_qual_l, fixed_desc[i].qual);
		}

		size_t p_pos_l = 0, p_ref_l = 0;
		int64_t prev_pos_local = 0;
		for (uint32_t i = 0; i < vcount; ++i)
		{
			int64_t delta = 0;
			read(out_pos, p_pos_l, delta);
			delta += prev_pos_local;
			prev_pos_local = delta;
			sort_desc[i].pos = static_cast<uint64_t>(delta);
			read_str(out_ref, p_ref_l, sort_desc[i].ref);
			sort_desc[i].filter = ".";
			sort_desc[i].info = ".";
		}

		// Permutations
		if (useLegacyPath)
		{
			std::vector<uint32_t> perm(block_size, 0);
			if (vcount == block_size)
			{
				out_perm(perm, sort_desc);
				vector<vector<uint32_t>> row_perm;
				row_perm.emplace_back(perm);
				s_perm.emplace_back(row_perm);
			}
			else
			{
				auto it = vint_last_perm.find(cur_chunk_id);
				if (it == vint_last_perm.end())
				{
					perm.resize(block_size);
					for (size_t i_p = 0; i_p < perm.size(); ++i_p)
						perm[i_p] = static_cast<uint32_t>(i_p);
				}
				else
				{
					perm = vint_code::DecodeArray(it->second);
				}
				vector<vector<uint32_t>> row_perm;
				row_perm.emplace_back(perm);
				s_perm.emplace_back(row_perm);

				for (size_t i_p = 0; i_p < sort_desc.size(); i_p++)
				{
					if (atoi(sort_desc[i_p].ref.c_str()))
					{
						sort_desc[i_p].ref = sort_desc[i_p].ref.substr(to_string(atoi(sort_desc[i_p].ref.c_str())).length());
					}
				}
			}
		}
		else
		{
			vector<vector<uint32_t>> row_perm(n_col_blocks);
			for (uint32_t cb = 0; cb < n_col_blocks; ++cb)
			{
				auto it = vint_last_perm_2d.find(make_pair(global_block_start + rb, cb));
				if (it == vint_last_perm_2d.end())
				{
					row_perm[cb].resize(col_block_ranges[cb].second);
					for (size_t i_p = 0; i_p < row_perm[cb].size(); ++i_p)
						row_perm[cb][i_p] = static_cast<uint32_t>(i_p);
				}
				else
				{
					row_perm[cb] = vint_code::DecodeArray(it->second);
				}
			}
			s_perm.emplace_back(row_perm);

			for (size_t i_p = 0; i_p < sort_desc.size(); i_p++)
			{
				if (atoi(sort_desc[i_p].ref.c_str()))
				{
					sort_desc[i_p].ref = sort_desc[i_p].ref.substr(to_string(atoi(sort_desc[i_p].ref.c_str())).length());
				}
			}
		}

		// Merge fixed fields into sort_desc (file order must already match sort fields)
		for (uint32_t i = 0; i < vcount; ++i)
		{
			sort_desc[i].chrom = std::move(fixed_desc[i].chrom);
			sort_desc[i].id = std::move(fixed_desc[i].id);
			sort_desc[i].alt = std::move(fixed_desc[i].alt);
			sort_desc[i].qual = std::move(fixed_desc[i].qual);
		}

		v_blocks.emplace_back(sort_desc);
	}

	return true;
}
//*******************************************************************************************************************************
// bool DecompressionReader::Decoder(vector<block_t> &v_blocks,vector<vector<uint32_t>> &s_perm,vector<uint8_t> &gt_index,uint32_t cur_chunk_id){

// 	initDecoderParams();
// 	no_variants =  fixed_field_block_compress.no_variants;

// 	// std::cerr<<"no_variants: "<<no_variants<<endl;
// 	vector<variant_desc_t> v_vcf_fixed_data_io;
// 	vector<variant_desc_t> v_vcf_sort_data_io;
// 	uint32_t block_size = n_samples * uint32_t(ploidy);
// 	v_vcf_fixed_data_io.reserve(no_variants_in_buf);
// 	v_vcf_sort_data_io.reserve(no_variants_in_buf);
// 	vector<uint32_t> perm(block_size,0);
// 	unique_ptr<thread> decomp_fixed_field_thread(new thread([&]{
// 		for (auto data : {
//         	make_tuple(ref(fixed_field_block_io.chrom), ref(fixed_field_block_compress.chrom), p_chrom, "chrom"),
//         	make_tuple(ref(fixed_field_block_io.id), ref(fixed_field_block_compress.id), p_id, "id"),
//         	make_tuple(ref(fixed_field_block_io.alt), ref(fixed_field_block_compress.alt), p_alt, "alt"),
//     		make_tuple(ref(fixed_field_block_io.qual), ref(fixed_field_block_compress.qual), p_qual, "qual"),

//     	})
// 		{
// 			CBSCWrapper bsc;
// 			bsc.InitDecompress();
// 			bsc.Decompress(get<1>(data), get<0>(data));
// 			// std::cerr<<get<3>(data)<<":"<<get<1>(data).size()<<":"<<get<0>(data).size()<<endl;

// 		}
// 		uint32_t i_variant;
// 		variant_desc_t desc;
// 		for(i_variant = 0; i_variant < no_variants; i_variant++)
// 		{
// 			// Load variant description
// 			read_str(fixed_field_block_io.chrom, p_chrom, desc.chrom);
// 			read_str(fixed_field_block_io.id, p_id, desc.id);
// 			read_str(fixed_field_block_io.alt, p_alt, desc.alt);
// 			read_str(fixed_field_block_io.qual, p_qual, desc.qual);
// 			v_vcf_fixed_data_io.emplace_back(desc);
// 		}

// 	}));
// 	unique_ptr<thread> decomp_sort_field__thread(new thread([&]{

// 		for (auto data : {

//         	make_tuple(ref(fixed_field_block_io.pos), ref(fixed_field_block_compress.pos), p_pos, "pos"),
//         	make_tuple(ref(fixed_field_block_io.ref), ref(fixed_field_block_compress.ref), p_ref, "ref"),
//         	// make_tuple(ref(fixed_field_block_io.gt_block), ref(fixed_field_block_compress.gt_block), p_gt, "GT")

//     	})
// 		{
// 			CBSCWrapper bsc;
// 			bsc.InitDecompress();
// 			bsc.Decompress(get<1>(data), get<0>(data));
// 			// std::cerr<<get<3>(data)<<":"<<get<1>(data).size()<<":"<<get<0>(data).size()<<endl;

// 		}
// 		if(fixed_field_block_compress.gt_block.back() == 0)
// 		{
// 			fixed_field_block_compress.gt_block.pop_back();
// 			CBSCWrapper bsc;
// 			bsc.InitDecompress();
// 			bsc.Decompress(fixed_field_block_compress.gt_block, fixed_field_block_io.gt_block);

// 		}else{
// 			fixed_field_block_compress.gt_block.pop_back();
// 			zstd::zstd_decompress(fixed_field_block_compress.gt_block, fixed_field_block_io.gt_block);
// 			// std::cerr<<"fixed_field_block_io:"<<fixed_field_block_io.gt_block.size()<<endl;
// 		}
// 		uint32_t i_variant;
// 		variant_desc_t desc;
// 		for(i_variant = 0; i_variant < no_variants; i_variant++)
// 		{
// 			int64_t pos = 0;
// 			// Load variant description

// 			read(fixed_field_block_io.pos, p_pos, pos);

// 			pos += prev_pos;
// 			prev_pos = pos;
// 			desc.pos = pos;
// 			read_str(fixed_field_block_io.ref, p_ref, desc.ref);
// 			desc.filter = ".";
// 			desc.info = ".";
// 			v_vcf_sort_data_io.emplace_back(desc);
// 			if ((i_variant+1) % block_size == 0)
// 			{
// 				out_perm(perm, v_vcf_sort_data_io);
// 				s_perm.emplace_back(perm);
// 				v_blocks.emplace_back(v_vcf_sort_data_io);
// 				v_vcf_sort_data_io.clear();
// 			}
// 		}
// 		if (i_variant % block_size)
// 		{
// 			auto it = vint_last_perm.find(cur_chunk_id);
// 			// std::cerr<<cur_chunk_id<<endl;
// 			// std::cerr<<vint_last_perm.size()<<":"<<it->first<<endl;
// 			// if(it == vint_last_perm.end())
// 			// 	return false;
// 			// for (size_t i_p = 0; i_p < perm.size(); i_p++)
// 			// {
// 			// 	perm[i_p] = i_p;
// 			// }

// 			perm = vint_code::DecodeArray(it->second);
// 			s_perm.emplace_back(perm);

// 			for (size_t i_p = 0; i_p < i_variant % block_size; i_p++)
// 			{
// 				if (atoi(v_vcf_sort_data_io[i_p].ref.c_str()))
// 				{
// 					v_vcf_sort_data_io[i_p].ref = v_vcf_sort_data_io[i_p].ref.substr(to_string(atoi(v_vcf_sort_data_io[i_p].ref.c_str())).length());
// 				}
// 			}

// 			v_blocks.emplace_back(v_vcf_sort_data_io);
// 			v_vcf_sort_data_io.clear();
// 		}
// 	}));
// 	decomp_fixed_field_thread->join();
// 	decomp_sort_field__thread->join();
// 	uint32_t i_variant = 0;

// 	for(size_t i = 0; i < v_blocks.size(); i++)
// 	{
// 		for (size_t j = 0; j < v_blocks[i].data_compress.size(); j++)
// 		{
// 			v_blocks[i].data_compress[j].chrom = std::move(v_vcf_fixed_data_io[i_variant].chrom);
// 			v_blocks[i].data_compress[j].id = std::move(v_vcf_fixed_data_io[i_variant].id);
// 			v_blocks[i].data_compress[j].alt = std::move(v_vcf_fixed_data_io[i_variant].alt);
// 			v_blocks[i].data_compress[j].qual = std::move(v_vcf_fixed_data_io[i_variant].qual);

// 			i_variant++;
// 		}

// 	}
// 	// v_vcf_fixed_data_io.clear();
// 	gt_index = move(fixed_field_block_io.gt_block);
// 	return true;
// }

// // *******************************************************************************************************************************
// void DecompressionReader::SetNoSamples(uint32_t _no_samples)
// {
// 	no_samples = _no_samples;
// }

// // *******************************************************************************************************************************
// bool DecompressionReader::GetHeader(string &_v_header)
// {
// 	_v_header = v_header;
// 	return true;
// }
// // *******************************************************************************************************************************
// bool DecompressionReader::SetHeader(string &_v_header)
// {
// 	v_header = _v_header;
// 	return true;
// }
// // *******************************************************************************************************************************
// bool DecompressionReader::SetSamples(vector<string> &_v_samples)
// {
// 	v_samples = _v_samples;
// 	return true;
// }
// // *******************************************************************************************************************************
// bool DecompressionReader::GetSamples(vector<string> &_v_samples)
// {
// 	_v_samples = v_samples;
// 	return true;
// }
