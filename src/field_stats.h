#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <mutex>
#include <atomic>

// Structure to hold compression statistics for a single field
struct FieldCompressionStats {
    std::string name;           // Field name (e.g., "CHROM", "POS", "GT")
    std::atomic<uint64_t> raw_size{0};        // Original uncompressed size in bytes
    std::atomic<uint64_t> compressed_size{0}; // Compressed size in bytes

    FieldCompressionStats() = default;
    FieldCompressionStats(const std::string& n) : name(n), raw_size(0), compressed_size(0) {}

    // Copy constructor for atomic members
    FieldCompressionStats(const FieldCompressionStats& other)
        : name(other.name),
          raw_size(other.raw_size.load()),
          compressed_size(other.compressed_size.load()) {}

    // Get compression ratio (compressed/raw), returns 0 if raw_size is 0
    double GetCompressionRatio() const {
        uint64_t raw = raw_size.load();
        if (raw == 0) return 0.0;
        return static_cast<double>(compressed_size.load()) / static_cast<double>(raw);
    }

    // Get bits per element (for GT field, elements = variants * samples * ploidy)
    double GetBitsPerElement(uint64_t element_count) const {
        if (element_count == 0) return 0.0;
        return static_cast<double>(compressed_size.load() * 8) / static_cast<double>(element_count);
    }

    void AddRaw(uint64_t bytes) {
        raw_size.fetch_add(bytes, std::memory_order_relaxed);
    }

    void AddCompressed(uint64_t bytes) {
        compressed_size.fetch_add(bytes, std::memory_order_relaxed);
    }

    void Reset() {
        raw_size.store(0);
        compressed_size.store(0);
    }
};

// Container for all field statistics
struct CompressionStatistics {
    // Fixed fields (VCF standard columns)
    FieldCompressionStats chrom{"CHROM"};
    FieldCompressionStats pos{"POS"};
    FieldCompressionStats id{"ID"};
    FieldCompressionStats ref{"REF"};
    FieldCompressionStats alt{"ALT"};
    FieldCompressionStats qual{"QUAL"};
    FieldCompressionStats filter{"FILTER"};
    FieldCompressionStats info{"INFO"};

    // GT field (genotype) - core compression target
    FieldCompressionStats gt{"GT"};

    // FORMAT fields (sample-level)
    std::unordered_map<std::string, FieldCompressionStats> format_fields;

    // Metadata
    FieldCompressionStats meta{"META"};

    // Total statistics
    uint64_t total_variants{0};
    uint64_t total_samples{0};
    uint32_t ploidy{2};

    mutable std::mutex format_mutex;

    void AddFormatRaw(const std::string& name, uint64_t bytes) {
        std::lock_guard<std::mutex> lock(format_mutex);
        auto it = format_fields.find(name);
        if (it == format_fields.end())
            it = format_fields.emplace(name, FieldCompressionStats{name}).first;
        it->second.AddRaw(bytes);
    }

    void AddFormatCompressed(const std::string& name, uint64_t bytes) {
        std::lock_guard<std::mutex> lock(format_mutex);
        auto it = format_fields.find(name);
        if (it == format_fields.end())
            it = format_fields.emplace(name, FieldCompressionStats{name}).first;
        it->second.AddCompressed(bytes);
    }

    // Calculate total raw and compressed sizes
    uint64_t GetTotalRawSize() const {
        uint64_t total = 0;
        total += chrom.raw_size.load();
        total += pos.raw_size.load();
        total += id.raw_size.load();
        total += ref.raw_size.load();
        total += alt.raw_size.load();
        total += qual.raw_size.load();
        total += filter.raw_size.load();
        total += info.raw_size.load();
        total += gt.raw_size.load();
        total += meta.raw_size.load();
        {
            std::lock_guard<std::mutex> lock(format_mutex);
            for (const auto& kv : format_fields) {
                total += kv.second.raw_size.load();
            }
        }
        return total;
    }

    uint64_t GetTotalCompressedSize() const {
        uint64_t total = 0;
        total += chrom.compressed_size.load();
        total += pos.compressed_size.load();
        total += id.compressed_size.load();
        total += ref.compressed_size.load();
        total += alt.compressed_size.load();
        total += qual.compressed_size.load();
        total += filter.compressed_size.load();
        total += info.compressed_size.load();
        total += gt.compressed_size.load();
        total += meta.compressed_size.load();
        {
            std::lock_guard<std::mutex> lock(format_mutex);
            for (const auto& kv : format_fields) {
                total += kv.second.compressed_size.load();
            }
        }
        return total;
    }

    void Reset() {
        chrom.Reset();
        pos.Reset();
        id.Reset();
        ref.Reset();
        alt.Reset();
        qual.Reset();
        filter.Reset();
        info.Reset();
        gt.Reset();
        meta.Reset();
        {
            std::lock_guard<std::mutex> lock(format_mutex);
            format_fields.clear();
        }
        total_variants = 0;
        total_samples = 0;
        ploidy = 2;
    }
};
