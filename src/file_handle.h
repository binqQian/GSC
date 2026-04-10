#pragma once

// #include <algorithm>
#include<stdio.h>
#include <vector>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <map> 
#include <mutex>
#include <condition_variable> 
using namespace std;

// lossless part2 的多流容器：
// 每个 stream 由多个 part 组成，文件尾 footer 记录每个 part 的 offset/size。
class File_Handle_2 {
private:
	bool input_mode;
	FILE* f;
	size_t f_offset;
	string file_name;
	size_t range_begin = 0;
	size_t range_end = 0;
	bool use_range = false;

	struct part_t{
		// 一个 part 对应 stream 中的一段连续 payload。
		size_t offset;
		size_t size;

		part_t() : offset(0), size(0)
		{};

		part_t(size_t _offset, size_t _size) : offset(_offset), size(_size)
		{};
	};

	typedef struct {
		// cur_id 在写模式下表示已分配 part 数，在读模式下表示当前读取游标。
		string stream_name;
		size_t cur_id;
		vector<part_t> parts;
	} stream_t;

	map<int, stream_t> m_streams;
	mutex mtx;

	bool serialize();
	bool deserialize();
	// footer 内部使用的轻量写读函数。
	size_t WriteFixed(size_t x, FILE* file);
	size_t Write(size_t x, FILE* file);
	size_t Write(string s, FILE* file);
	size_t read_fixed(size_t& x, FILE* file);
	size_t read(size_t& x, FILE* file);
	size_t read(string& s, FILE* file);
	bool OpenInternal(const string& file_name, bool use_range, size_t range_begin, size_t range_size);

public:
	File_Handle_2(bool _input_mode);
	~File_Handle_2();

	bool Open(const string& temp_file2_fname);
	bool OpenRange(const string& temp_file2_fname, size_t range_begin, size_t range_size);
	bool Close();

	// 写模式下：
	// 1. RegisterStream 注册逻辑流
	// 2. AddPartPrepare 预留 part_id
	// 3. AddPartComplete 按顺序补齐真实数据
	int RegisterStream(string stream_name);
	int GetStreamId(string stream_name);

	bool AddPart(int stream_id, vector<uint8_t> &v_data, size_t metadata = 0);
	int AddPartPrepare(int stream_id);
	bool AddPartComplete(int stream_id, int part_id, vector<uint8_t>& v_data);
	void AddParamsPart(int stream_id,vector<uint8_t> & v_data);

	// 读模式下按 stream 顺序迭代取 part；ResetStreamPartIterator 可回到开头。
	bool GetPart(int stream_id, vector<uint8_t> &v_data);
	void SetRawSize(int stream_id, size_t raw_size);
	size_t GetRawSize(int stream_id);
	size_t GetCompressedSize(int stream_id);
	bool ResetStreamPartIterator(int stream_id);

	bool LinkStream(int stream_id, string stream_name, int target_id);


	size_t GetNoStreams()
	{
		lock_guard<mutex> lck(mtx);

		return m_streams.size();
	}
};

// 简单的带缓冲输出器，支持字节流和位流写出。
class COutFile
{
	const size_t BUFFER_SIZE = 4 << 20;

	FILE *f;
	uint8_t *buffer;
	size_t buffer_pos;
	size_t bit_pos;
	uint8_t byte_buffer;
	bool success;

public:
	COutFile() : f(nullptr), buffer(nullptr)
	{};

	~COutFile()
	{
		if (f)
			Close();
		if (buffer)
			delete[] buffer;
	}

	bool Open(string file_name, const char* write_mode)
	{
		// 这里只负责打开文件和准备大缓冲，不做格式层逻辑。
		if (f)
			return false;
		
		f = fopen(file_name.c_str(), write_mode);
		if (!f)
			return false;
		
		buffer = new uint8_t[BUFFER_SIZE];
		buffer_pos = 0;
		bit_pos = 0;
		byte_buffer = 0;
		success = true;

		return true;
	}

	bool Close()
	{
		// Close 会把残留缓冲刷盘，但不会自动补齐位流尾部。
		if (!f)
			return true;

		if (buffer_pos)
			success &= fwrite(buffer, 1, buffer_pos, f) == buffer_pos;

		if (f)
		{
			fclose(f);
			f = nullptr;
		}
		if (buffer)
		{
			delete[] buffer;
			buffer = nullptr;
		}

		return success;
	}

	void PutByte(uint8_t c)
	{
		
		if (buffer_pos == BUFFER_SIZE)
		{
			success &= fwrite(buffer, 1, BUFFER_SIZE, f) == BUFFER_SIZE;
			buffer_pos = 0;
		}

		buffer[buffer_pos++] = c;
		
	}

	void Put(char c)
	{
		if (buffer_pos == BUFFER_SIZE)
		{
			success &= fwrite(buffer, 1, BUFFER_SIZE, f) == BUFFER_SIZE;
			buffer_pos = 0;
		}

		buffer[buffer_pos++] = c;
	}
	bool PutBit(uint32_t word)
	{
		// 把 bit 先攒到 byte_buffer，凑满 8 位再真正写盘。
		if(bit_pos < 8)
		{
			
			word <<= bit_pos;
			byte_buffer += word;
			++bit_pos;
			
		}
		else
		{
			
			PutByte(byte_buffer);
			bit_pos = 1;
			byte_buffer = word;
			
		}
		
		return true;
	};
	void FlushPartialByteBuffer(){
		// 位流写完后要显式 flush，否则最后不足 1 字节的内容会丢失。
		if(bit_pos){
			
			PutByte(byte_buffer);
			bit_pos = 0;
			byte_buffer = 0;
		}
		return ;
	}
	void Write(const uint8_t *p, size_t n)
	{
		uint8_t *q = (uint8_t *)p;

		while (buffer_pos + n > BUFFER_SIZE)
		{
			size_t small_n = BUFFER_SIZE - buffer_pos;
			memcpy(buffer + buffer_pos, q, small_n);
			success &= fwrite(buffer, 1, BUFFER_SIZE, f) == BUFFER_SIZE;

			buffer_pos = 0;
			n -= small_n;
			q += small_n;
		}

		memcpy(buffer + buffer_pos, q, n);
		buffer_pos += n;
	}
	void Write(const char *p, size_t n)
	{
		char *q = (char *)p;

		while (buffer_pos + n > BUFFER_SIZE)
		{
			size_t small_n = BUFFER_SIZE - buffer_pos;
			memcpy(buffer + buffer_pos, q, small_n);
			success &= fwrite(buffer, 1, BUFFER_SIZE, f) == BUFFER_SIZE;

			buffer_pos = 0;
			n -= small_n;
			q += small_n;
		}

		memcpy(buffer + buffer_pos, q, n);
		buffer_pos += n;
	}
	void WriteUInt(uint64_t x, int no_bytes)
	{
		for (int i = 0; i < no_bytes; ++i)
		{
			PutByte(x & 0xffu);
			x >>= 8;
		}
	}

	void Write(string &s)
	{
		Write((uint8_t*)s.c_str(), s.size());
	}

	void Write(string &s, size_t start_pos, size_t len)
	{
		Write((uint8_t*)s.c_str() + start_pos, len);
	}
};
