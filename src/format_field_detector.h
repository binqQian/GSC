#pragma once

#include "format_field_codec.h"
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <memory>
#include <vector>
#include <functional>

namespace gsc {

// ============================================================================
// Detector Configuration
// ============================================================================
struct DetectorConfig {
    // Thresholds for codec selection
    float missing_high_threshold = 0.8f;     // High missing ratio threshold
    float predictor_hit_threshold = 0.95f;   // Predictor hit ratio threshold
    float pattern_threshold = 0.9f;          // Pattern matching threshold

    // Observation limits
    uint32_t observe_sample_limit = 256;     // Max samples to observe before freeze
    uint32_t unique_value_limit = 1024;      // Max unique values to track

    // Array detection
    uint32_t expected_array_len = 0;         // Expected array length (from allele count)
};

// ============================================================================
// FormatFieldDetector - Analyzes field values to determine optimal codec
// ============================================================================
class FormatFieldDetector {
public:
    explicit FormatFieldDetector(const std::string& field_name,
                                  const DetectorConfig& config = {});

    // -------------------------------------------------------------------------
    // Observation phase
    // -------------------------------------------------------------------------

    // Feed a single sample's value for analysis
    void feed(const char* value, size_t len);

    // Set expected array length (from allele count)
    void setExpectedArrayLen(uint32_t len);

    // Check if observation limit reached
    bool isComplete() const;

    // -------------------------------------------------------------------------
    // Cross-field predictor support
    // -------------------------------------------------------------------------

    // Register a predictor function that provides expected values
    // predictor(sample_pos) -> expected_value
    void setPredictor(std::function<int64_t(uint32_t)> predictor);

    // Get computed value for this sample (e.g., sum(AD) for DP prediction)
    // Returns INT64_MIN if not applicable
    int64_t getComputedValue(uint32_t sample_pos) const;

    // -------------------------------------------------------------------------
    // Finalization
    // -------------------------------------------------------------------------

    // Compute final features from observations
    FieldFeatures finalize();

    // Create optimal codec based on features
    std::unique_ptr<FormatFieldCodec> createCodec();

    // Get recommended codec type
    CodecType recommendedCodecType() const;

    // -------------------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------------------

    const std::string& fieldName() const { return field_name_; }
    const FieldFeatures& features() const { return features_; }
    bool hasPredictor() const { return predictor_ != nullptr; }

private:
    std::string field_name_;
    DetectorConfig config_;
    FieldFeatures features_;
    bool finalized_ = false;

    // Observation state
    uint32_t sample_count_ = 0;
    uint32_t missing_count_ = 0;

    // Type detection
    bool has_comma_ = false;        // Indicates array type
    bool all_numeric_ = true;       // All values are numeric
    bool all_integer_ = true;       // All values are integers

    // Value range
    int64_t min_val_ = INT64_MAX;
    int64_t max_val_ = INT64_MIN;

    // Array length tracking
    std::vector<uint32_t> array_lengths_;

    // Unique value tracking (with cardinality limit)
    std::unordered_map<std::string, uint32_t> value_freq_;
    uint32_t overflow_count_ = 0;   // Values not tracked due to limit

    // Zero array count (for arrays)
    uint32_t zero_array_count_ = 0;

    // Computed values for cross-field prediction (e.g., sum(AD))
    std::vector<int64_t> computed_values_;

    // Predictor for cross-field correlation
    std::function<int64_t(uint32_t)> predictor_;
    uint32_t predictor_hit_count_ = 0;
    uint32_t predictor_miss_count_ = 0;

    // Pattern detection (for PL-like fields)
    uint32_t pattern_hit_count_ = 0;

    // Helper functions
    void parseValue(const char* value, size_t len);
    void parseArrayValue(const char* value, size_t len);
    void parseScalarValue(const char* value, size_t len);
    bool checkAbPattern(const std::vector<int64_t>& arr);
    int64_t computeArraySum(const std::vector<int64_t>& arr);
    uint8_t computeBitsNeeded(int64_t max_val);
};

// ============================================================================
// FormatFieldManager - Coordinates detection and encoding for all FORMAT fields
// ============================================================================
class FormatFieldManager {
public:
    FormatFieldManager();

    // -------------------------------------------------------------------------
    // Row-level operations
    // -------------------------------------------------------------------------

    // Initialize for a new row with given FORMAT keys
    // allele_count is used to determine expected array lengths
    void initRow(const std::vector<std::string>& format_keys,
                 uint32_t allele_count);

    // Process a single sample's complete FORMAT string
    // sample_str: "0/1:10,5:15:30:..." (GT:AD:DP:GQ:...)
    void processSample(const char* sample_str, size_t len, uint32_t sample_pos);

    // Finalize row: freeze detectors, encode remaining samples, serialize
    void finalizeRow(std::vector<uint8_t>& out);

    // -------------------------------------------------------------------------
    // Decoding
    // -------------------------------------------------------------------------

    // Decode a sample's FORMAT string from serialized data
    std::string decodeSample(const uint8_t* data, size_t len,
                             uint32_t sample_pos,
                             const std::vector<std::string>& format_keys);

    // -------------------------------------------------------------------------
    // Configuration
    // -------------------------------------------------------------------------

    void setConfig(const DetectorConfig& config) { config_ = config; }

    // Enable/disable GT field processing (usually handled separately)
    void setSkipGT(bool skip) { skip_gt_ = skip; }

private:
    DetectorConfig config_;
    bool skip_gt_ = true;  // GT is typically handled separately

    // Current row state
    std::vector<std::string> current_format_keys_;
    uint32_t allele_count_ = 2;
    uint32_t samples_processed_ = 0;
    bool frozen_ = false;

    // Per-field detectors (excluding GT)
    std::unordered_map<std::string, std::unique_ptr<FormatFieldDetector>> detectors_;

    // Per-field codecs (created after freeze)
    std::unordered_map<std::string, std::unique_ptr<FormatFieldCodec>> codecs_;

    // Buffered samples before freeze (kept in FORMAT key order)
    struct BufferedSample {
        std::vector<std::string> values;  // Same size/order as current_format_keys_
    };
    std::vector<BufferedSample> sample_buffer_;

    // Helper functions
    void parseSampleFields(const char* sample_str, size_t len,
                           std::vector<std::string>& values);
    void encodeSampleValues(const std::vector<std::string>& values, uint32_t sample_pos);
    void freezeAndFlush();
    void setupCrossFieldPredictors();
    void setupCrossFieldCodecPredictors();

    // Expected array lengths based on allele count
    uint32_t getExpectedArrayLen(const std::string& field_name) const;
};

} // namespace gsc
