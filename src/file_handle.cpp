#include "file_handle.h"
#include "defs.h"
#include "logger.h"



// ******************************************************************************
// 构造一个多流句柄；input_mode=false 表示后续会写 footer。
File_Handle_2::File_Handle_2(bool _input_mode)
{
	f = nullptr;
	input_mode = _input_mode;
}

// 析构时兜底关闭文件，确保写模式下 footer 有机会落盘。
// ******************************************************************************
File_Handle_2::~File_Handle_2()
{
    
	if (f){
        
        Close();
    }
		
}

// Open/OpenRange 最终都会走这里。
// 读模式下会先反序列化 footer，建立 stream -> parts 映射。
// ******************************************************************************
bool File_Handle_2::OpenInternal(const string& temp_file2_fname, bool range_enabled, size_t range_begin_in, size_t range_size_in)
{
	auto logger = LogManager::Instance().Logger();
	lock_guard<mutex> lck(mtx);

	if (f)
		fclose(f);

	m_streams.clear();
	use_range = range_enabled;
	range_begin = range_begin_in;
	range_end = range_enabled ? (range_begin_in + range_size_in) : 0;
	// file_name = _file_name+".temp_gsc";

	f = fopen(temp_file2_fname.c_str(), input_mode ? "rb" : "wb");

	if (!f){
        logger->error("Can't open: {} failed", temp_file2_fname);
		return false;
    }
	setvbuf(f, nullptr, _IOFBF, 64 << 20);
	if (input_mode)
	{
		if (!deserialize())
		{
			fclose(f);
			f = nullptr;
			return false;
		}
	}

	f_offset = 0;

	return true;
}

// ******************************************************************************
bool File_Handle_2::Open(const string& temp_file2_fname)
{
	return OpenInternal(temp_file2_fname, false, 0, 0);
}

// ******************************************************************************
bool File_Handle_2::OpenRange(const string& temp_file2_fname, size_t range_begin_in, size_t range_size_in)
{
	if (!input_mode)
		return Open(temp_file2_fname);
	return OpenInternal(temp_file2_fname, true, range_begin_in, range_size_in);
}

// 关闭句柄：
// 写模式会先把 footer 追加到文件尾；读模式直接关闭。
// // ******************************************************************************
bool File_Handle_2::Close()
{
	lock_guard<mutex> lck(mtx);

	if (!f)
		return false;

	if (input_mode)
	{
		
		fclose(f);
		f = nullptr;
	}
	else
	{
        
		serialize();
		fclose(f);
		f = nullptr;
	}

	return true;
}
void File_Handle_2::AddParamsPart(int stream_id,vector<uint8_t> & v_data){
    lock_guard<mutex> lck(mtx);
    // 参数流通常不需要 prepare/complete 两阶段，直接按当前偏移追加。
	m_streams[stream_id].parts.emplace_back(f_offset, v_data.size());;
    if(v_data.size())
		fwrite(v_data.data(), 1, v_data.size(), f);

	f_offset += v_data.size();
}
// 将所有 stream 的 part 目录写到文件尾，并在最后补一个固定 8 字节 footer_size。
// 读取时会先从文件尾倒着找到这段目录。
// // ******************************************************************************
bool File_Handle_2::serialize()
{
    
	size_t footer_size = 0;

	// Store stream part offsets
	footer_size += Write(m_streams.size(), f);
	for (auto& stream : m_streams)
	{
        
		size_t str_size = 0;
		
		footer_size += Write(stream.second.stream_name, f);
		footer_size += Write(stream.second.parts.size(), f);
        for (auto& part : stream.second.parts)
		{
			footer_size += Write(part.offset, f);
			footer_size += Write(part.size, f);
            // std::cerr<<stream.second.stream_name<<":"<<part.offset<<":"<<part.size<<endl;
			str_size += part.size;
		}

#ifdef LOG_INFO
        auto logger = LogManager::Instance().Logger();
		logger->info("{}: {}  raw size: {}   packed size: {}", stream.first, stream.second.stream_name, stream.second.raw_size, str_size);
#endif
	}

	WriteFixed(footer_size, f);
	// std::cerr<<footer_size<<endl;
    // std::std::cerr << "genotype compress file (" << file_name << ") created." << std::endl;
	return true;
}

// 从文件尾读取 footer，恢复每个 stream 的 part 目录。
// range 模式下，footer 位于子区间尾部而不是整个文件尾部。
// // ******************************************************************************
bool File_Handle_2::deserialize()
{
	size_t footer_size;
	if (use_range)
	{
		if (range_end < range_begin + 8)
			return false;
		my_fseek(f, static_cast<long long>(range_end - 8), SEEK_SET);
		read_fixed(footer_size, f);
		my_fseek(f, static_cast<long long>(range_end - 8 - footer_size), SEEK_SET);
	}
	else
	{
		my_fseek(f, -8, SEEK_END);
		read_fixed(footer_size, f);

		my_fseek(f, -(long)(8 + footer_size), SEEK_END);
	}

	// Load stream part offsets
	size_t n_streams;
	read(n_streams, f);
	for (size_t i = 0; i < n_streams; ++i)
	{
		m_streams[(int) i] = stream_t();
		auto& stream_second = m_streams[(int) i];

		read(stream_second.stream_name, f);
		read(stream_second.cur_id, f);

		stream_second.parts.resize(stream_second.cur_id);
		for (size_t j = 0; j < stream_second.cur_id; ++j)
		{
			read(stream_second.parts[j].offset, f);
			read(stream_second.parts[j].size, f);
		}

		stream_second.cur_id = 0;
	}
	
	if (use_range)
		my_fseek(f, static_cast<long long>(range_begin), SEEK_SET);
	else
		my_fseek(f, 0, SEEK_SET);

	return true;
}
// 注册一个新逻辑流，返回内部 stream_id。
// ******************************************************************************
int File_Handle_2::RegisterStream(string stream_name)
{
	lock_guard<mutex> lck(mtx);

	int id = (int) m_streams.size();
	
	m_streams[id] = stream_t();
	m_streams[id].cur_id = 0;
	m_streams[id].stream_name = stream_name;
    
	return id;
}

// 根据 stream_name 反查内部 id。
// ******************************************************************************
int File_Handle_2::GetStreamId(string stream_name)
{
	lock_guard<mutex> lck(mtx);
	for (auto& x : m_streams)
		if (x.second.stream_name == stream_name){
            
            return x.first;
        }
			

	return -1;
}
// bool File_Handle_2::AddPart(int stream_id, vector<uint8_t> &v_data)
// {
	
// 	lock_guard<mutex> lck(mtx);

// 	m_streams[stream_id].parts.emplace_back(f_offset, v_data.size());
	
// 	f_offset += write(metadata, f);
	
// 	if(v_data.size())
// 		fwrite(v_data.data(), 1, v_data.size(), f);

// 	f_offset += v_data.size();

// 	return true;
// }

// 先预留 part 槽位，后续由 AddPartComplete 回填实际 offset/size。
// 这样多个压缩线程可以先并行算 payload，再由单写线程顺序落盘。
// ******************************************************************************
int File_Handle_2::AddPartPrepare(int stream_id)
{
	lock_guard<mutex> lck(mtx);
	m_streams[stream_id].parts.emplace_back(0, 0);
    // std::cerr<<m_streams[stream_id].parts.size()<<endl;

	return (int) m_streams[stream_id].parts.size() - 1;
}


// 把预留 part_id 对应的真实 payload 写入文件，并记录 offset/size。
// ******************************************************************************
bool File_Handle_2::AddPartComplete(int stream_id, int part_id, vector<uint8_t>& v_data)
{
	
	lock_guard<mutex> lck(mtx);
	// std::cerr<<m_streams.size()<<endl;
	m_streams[stream_id].parts[part_id] = part_t(f_offset, v_data.size());

	if (v_data.size())
		fwrite(v_data.data(), 1, v_data.size(), f);
	
	f_offset += v_data.size();

	return true;
}


// ******************************************************************************
size_t File_Handle_2::GetCompressedSize(int stream_id)
{
	lock_guard<mutex> lck(mtx);

	size_t size = 0;

	auto p = m_streams.find(stream_id);

	if (p == m_streams.end())
		return 0;

	for (auto q : p->second.parts)
		size += q.size;

	return size;
}

// 顺序读取某个 stream 的下一个 part。
// 读模式下 cur_id 就是这个 stream 当前的读取游标。
// ******************************************************************************
bool File_Handle_2::GetPart(int stream_id, vector<uint8_t> &v_data)
{
	lock_guard<mutex> lck(mtx);
	
	auto& p = m_streams[stream_id];
	if (p.cur_id >= p.parts.size())
		return false;
	v_data.resize(p.parts[p.cur_id].size);
    
	const size_t offset = use_range ? (range_begin + p.parts[p.cur_id].offset) : p.parts[p.cur_id].offset;
	my_fseek(f, static_cast<long long>(offset), SEEK_SET);

	auto r = fread(v_data.data(), 1, p.parts[p.cur_id].size, f);
    
	p.cur_id++;
    
	if (r != p.parts[p.cur_id-1].size)
		return false;
    
	return r == p.parts[p.cur_id-1].size;
}

// 将某个 stream 的读取游标重置回开头。
// ******************************************************************************
bool File_Handle_2::ResetStreamPartIterator(int stream_id)
{
	lock_guard<mutex> lck(mtx);

	m_streams[stream_id].cur_id = 0;

	return true;
}

// // ******************************************************************************
// size_t File_Handle_2::signature(vector<uint8_t>& v_data)
// {
// 	size_t h = 0;

// 	for (size_t i = 0; i < v_data.size(); i += 8)
// 	{
// 		size_t x = 0;
// 		for (size_t j = 0; j < 8u && i + j < v_data.size(); ++j)
// 			x = (x << 8) + (size_t)v_data[i + j];

// 		x ^= x >> 33;
// 		x *= 0xff51afd7ed558ccdL;
// 		x ^= x >> 33;
// 		x *= 0xc4ceb9fe1a85ec53L;
// 		x ^= x >> 33;

// 		h ^= x;
// 	}

// 	return h;
// }

// 让一个 stream 复用另一个 stream 的 part 目录，仅改逻辑名字。
// ******************************************************************************
bool File_Handle_2::LinkStream(int stream_id, string stream_name, int target_id)
{
	m_streams[stream_id] = m_streams[target_id];
	m_streams[stream_id].stream_name = stream_name;

	return true;
}
size_t File_Handle_2::WriteFixed(size_t x, FILE* file)
{
	// footer_size 固定写 8 字节，方便从文件尾反查。
	fwrite(&x, 1, 8, file);

	return 8;
}

// ******************************************************************************
size_t File_Handle_2::Write(size_t x, FILE* file)
{
	// 非固定宽整数按“字节数 + 大端内容”写入，减少目录体积。
	int no_bytes = 0;

	for (size_t tmp = x; tmp; tmp >>= 8)
		++no_bytes;
	
	putc(no_bytes, file);

	for (int i = no_bytes; i; --i)
		putc((x >> ((i - 1) * 8)) & 0xff, file);

	return no_bytes + 1;
}

size_t File_Handle_2::Write(string s, FILE* file)
{
	fwrite(s.c_str(), 1, s.size(), file);
	putc(0, file);

	return s.size() + 1;
}

// ******************************************************************************
size_t File_Handle_2::read_fixed(size_t& x, FILE* file)
{
	return fread(&x, 8, 1, file) * 8;
}
// ******************************************************************************
size_t File_Handle_2::read(size_t& x, FILE* file)
{
	// 与 Write(size_t) 对应：先读字节数，再按大端拼回整数。
	int no_bytes = getc(file);

	x = 0;

	for (int i = 0; i < no_bytes; ++i)
	{
		x <<= 8;
		x += (size_t)getc(file);
	}

	return no_bytes + 1;
}

// ******************************************************************************
size_t File_Handle_2::read(string& s, FILE* file)
{
	s.clear();

	while (true)
	{
		int c = getc(file);
		if (c == EOF)
			return 0;

		if (c == 0)
			return s.size() + 1;

		s.push_back((char)c);
	}

	return 0;
}
