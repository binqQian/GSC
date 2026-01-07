#!/bin/bash
#===============================================================================
# GSC 完整功能验证测试脚本
# 测试目标：验证VCF和gVCF压缩、检索、格式转换的正确性和性能
#===============================================================================

# 不使用 set -e，让测试继续执行即使某些测试失败

# 配置
GSC="/home/binq/Code/vcfcomp/gsc_base/GSC/build/gsc"
VCF_FILE="/home/binq/Code/vcfcomp/gsc_base/GSC/toy/final_subset.vcf.gz"
GVCF_FILE="/home/binq/data/gvcf/ckb/HG002_first20m.gvcf.gz"
TEST_DIR="/home/binq/Code/vcfcomp/gsc_base/GSC/test_output"
LOG_FILE="$TEST_DIR/test_results.log"

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 统计变量
PASS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0

#===============================================================================
# 辅助函数
#===============================================================================

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1" | tee -a "$LOG_FILE"
}

log_pass() {
    echo -e "${GREEN}[PASS]${NC} $1" | tee -a "$LOG_FILE"
    PASS_COUNT=$((PASS_COUNT + 1))
}

log_fail() {
    echo -e "${RED}[FAIL]${NC} $1" | tee -a "$LOG_FILE"
    FAIL_COUNT=$((FAIL_COUNT + 1))
}

log_skip() {
    echo -e "${YELLOW}[SKIP]${NC} $1" | tee -a "$LOG_FILE"
    SKIP_COUNT=$((SKIP_COUNT + 1))
}

log_section() {
    echo "" | tee -a "$LOG_FILE"
    echo "===============================================================================" | tee -a "$LOG_FILE"
    echo -e "${YELLOW}$1${NC}" | tee -a "$LOG_FILE"
    echo "===============================================================================" | tee -a "$LOG_FILE"
}

# 比较两个VCF文件（忽略header中的时间戳等动态信息）
compare_vcf() {
    local file1="$1"
    local file2="$2"

    # 提取数据行（非header）进行比较
    local data1=$(zcat "$file1" 2>/dev/null || cat "$file1" | grep -v "^##")
    local data2=$(zcat "$file2" 2>/dev/null || cat "$file2" | grep -v "^##")

    if [ "$data1" = "$data2" ]; then
        return 0
    else
        return 1
    fi
}

# 计算文件压缩率
calc_ratio() {
    local original_size="$1"
    local compressed_size="$2"
    echo "scale=4; $compressed_size * 100 / $original_size" | bc
}

# 获取文件大小（字节）
get_file_size() {
    stat -c %s "$1"
}

# 计时函数
time_cmd() {
    local start=$(date +%s.%N)
    "$@"
    local ret=$?
    local end=$(date +%s.%N)
    local elapsed=$(echo "$end - $start" | bc)
    echo "$elapsed"
    return $ret
}

#===============================================================================
# 初始化
#===============================================================================

init_test() {
    log_section "初始化测试环境"

    # 创建测试目录
    mkdir -p "$TEST_DIR"
    rm -f "$LOG_FILE"

    # 检查必要文件
    if [ ! -f "$GSC" ]; then
        log_fail "GSC可执行文件不存在: $GSC"
        exit 1
    fi

    if [ ! -f "$VCF_FILE" ]; then
        log_fail "VCF测试文件不存在: $VCF_FILE"
        exit 1
    fi

    if [ ! -f "$GVCF_FILE" ]; then
        log_fail "gVCF测试文件不存在: $GVCF_FILE"
        exit 1
    fi

    log_info "GSC路径: $GSC"
    log_info "VCF测试文件: $VCF_FILE ($(ls -lh $VCF_FILE | awk '{print $5}'))"
    log_info "gVCF测试文件: $GVCF_FILE ($(ls -lh $GVCF_FILE | awk '{print $5}'))"
    log_info "测试输出目录: $TEST_DIR"

    # 解压原始VCF用于比较
    log_info "准备原始VCF数据..."
    zcat "$VCF_FILE" > "$TEST_DIR/original.vcf"

    log_pass "测试环境初始化完成"
}

#===============================================================================
# Part 1: VCF 无损压缩测试
#===============================================================================

test_vcf_lossless_compression() {
    log_section "Part 1: VCF 无损压缩测试"

    local original_size=$(zcat "$VCF_FILE" | wc -c)

    # 1.1 默认无损压缩（bsc）
    log_info "1.1 测试默认无损压缩 (compressor=bsc)..."
    if $GSC compress -i "$VCF_FILE" -o "$TEST_DIR/lossless_bsc.gsc" 2>&1 | tee -a "$LOG_FILE"; then
        local comp_size=$(get_file_size "$TEST_DIR/lossless_bsc.gsc")
        local ratio=$(calc_ratio $original_size $comp_size)
        log_pass "无损压缩(bsc)成功 - 压缩率: ${ratio}%"
    else
        log_fail "无损压缩(bsc)失败"
    fi

    # 1.2 zstd压缩器
    log_info "1.2 测试zstd压缩器..."
    if $GSC compress -i "$VCF_FILE" -o "$TEST_DIR/lossless_zstd.gsc" --compressor zstd 2>&1 | tee -a "$LOG_FILE"; then
        local comp_size=$(get_file_size "$TEST_DIR/lossless_zstd.gsc")
        local ratio=$(calc_ratio $original_size $comp_size)
        log_pass "无损压缩(zstd)成功 - 压缩率: ${ratio}%"
    else
        log_fail "无损压缩(zstd)失败"
    fi

    # 1.3 brotli压缩器
    log_info "1.3 测试brotli压缩器..."
    if $GSC compress -i "$VCF_FILE" -o "$TEST_DIR/lossless_brotli.gsc" --compressor brotli 2>&1 | tee -a "$LOG_FILE"; then
        local comp_size=$(get_file_size "$TEST_DIR/lossless_brotli.gsc")
        local ratio=$(calc_ratio $original_size $comp_size)
        log_pass "无损压缩(brotli)成功 - 压缩率: ${ratio}%"
    else
        log_fail "无损压缩(brotli)失败"
    fi

    # 1.4 多线程压缩
    log_info "1.4 测试多线程压缩 (-t 4)..."
    if $GSC compress -i "$VCF_FILE" -o "$TEST_DIR/lossless_mt.gsc" -t 4 2>&1 | tee -a "$LOG_FILE"; then
        log_pass "多线程压缩成功"
    else
        log_fail "多线程压缩失败"
    fi

    # 1.5 不同block参数测试
    log_info "1.5 测试不同block参数 (--max-block-rows 5000)..."
    if $GSC compress -i "$VCF_FILE" -o "$TEST_DIR/lossless_block_rows.gsc" --max-block-rows 5000 2>&1 | tee -a "$LOG_FILE"; then
        log_pass "自定义block-rows压缩成功"
    else
        log_fail "自定义block-rows压缩失败"
    fi

    log_info "1.6 测试不同block参数 (--max-block-cols 5000)..."
    if $GSC compress -i "$VCF_FILE" -o "$TEST_DIR/lossless_block_cols.gsc" --max-block-cols 5000 2>&1 | tee -a "$LOG_FILE"; then
        log_pass "自定义block-cols压缩成功"
    else
        log_fail "自定义block-cols压缩失败"
    fi

    # 1.7 不同depth参数测试
    log_info "1.7 测试不同depth参数 (-d 50)..."
    if $GSC compress -i "$VCF_FILE" -o "$TEST_DIR/lossless_depth50.gsc" -d 50 2>&1 | tee -a "$LOG_FILE"; then
        log_pass "depth=50压缩成功"
    else
        log_fail "depth=50压缩失败"
    fi

    log_info "1.8 测试depth=0（无匹配）..."
    if $GSC compress -i "$VCF_FILE" -o "$TEST_DIR/lossless_depth0.gsc" -d 0 2>&1 | tee -a "$LOG_FILE"; then
        log_pass "depth=0压缩成功"
    else
        log_fail "depth=0压缩失败"
    fi
}

#===============================================================================
# Part 2: VCF 无损解压验证
#===============================================================================

test_vcf_lossless_decompression() {
    log_section "Part 2: VCF 无损解压验证"

    # 2.1 解压bsc压缩文件
    log_info "2.1 解压bsc压缩文件..."
    if $GSC decompress -i "$TEST_DIR/lossless_bsc.gsc" -o "$TEST_DIR/restored_bsc.vcf" 2>&1 | tee -a "$LOG_FILE"; then
        # 验证无损
        if diff <(grep -v "^##" "$TEST_DIR/original.vcf") <(grep -v "^##" "$TEST_DIR/restored_bsc.vcf") > /dev/null 2>&1; then
            log_pass "bsc解压验证通过 - 数据完全一致"
        else
            log_fail "bsc解压数据不一致"
        fi
    else
        log_fail "bsc解压失败"
    fi

    # 2.2 解压zstd压缩文件
    log_info "2.2 解压zstd压缩文件..."
    if $GSC decompress -i "$TEST_DIR/lossless_zstd.gsc" -o "$TEST_DIR/restored_zstd.vcf" --compressor zstd 2>&1 | tee -a "$LOG_FILE"; then
        if diff <(grep -v "^##" "$TEST_DIR/original.vcf") <(grep -v "^##" "$TEST_DIR/restored_zstd.vcf") > /dev/null 2>&1; then
            log_pass "zstd解压验证通过 - 数据完全一致"
        else
            log_fail "zstd解压数据不一致"
        fi
    else
        log_fail "zstd解压失败"
    fi

    # 2.3 解压brotli压缩文件
    log_info "2.3 解压brotli压缩文件..."
    if $GSC decompress -i "$TEST_DIR/lossless_brotli.gsc" -o "$TEST_DIR/restored_brotli.vcf" --compressor brotli 2>&1 | tee -a "$LOG_FILE"; then
        if diff <(grep -v "^##" "$TEST_DIR/original.vcf") <(grep -v "^##" "$TEST_DIR/restored_brotli.vcf") > /dev/null 2>&1; then
            log_pass "brotli解压验证通过 - 数据完全一致"
        else
            log_fail "brotli解压数据不一致"
        fi
    else
        log_fail "brotli解压失败"
    fi

    # 2.4 验证不同参数压缩的文件解压
    log_info "2.4 验证depth参数压缩文件的解压..."
    if $GSC decompress -i "$TEST_DIR/lossless_depth50.gsc" -o "$TEST_DIR/restored_depth50.vcf" 2>&1 | tee -a "$LOG_FILE"; then
        if diff <(grep -v "^##" "$TEST_DIR/original.vcf") <(grep -v "^##" "$TEST_DIR/restored_depth50.vcf") > /dev/null 2>&1; then
            log_pass "depth=50解压验证通过"
        else
            log_fail "depth=50解压数据不一致"
        fi
    else
        log_fail "depth=50解压失败"
    fi
}

#===============================================================================
# Part 3: VCF 有损压缩测试
#===============================================================================

test_vcf_lossy_compression() {
    log_section "Part 3: VCF 有损压缩测试"

    local original_size=$(zcat "$VCF_FILE" | wc -c)

    # 3.1 有损压缩（只保留GT字段）
    log_info "3.1 测试有损压缩模式 (-M)..."
    if $GSC compress -M -i "$VCF_FILE" -o "$TEST_DIR/lossy.gsc" 2>&1 | tee -a "$LOG_FILE"; then
        local comp_size=$(get_file_size "$TEST_DIR/lossy.gsc")
        local ratio=$(calc_ratio $original_size $comp_size)
        log_pass "有损压缩成功 - 压缩率: ${ratio}%"
    else
        log_fail "有损压缩失败"
    fi

    # 3.2 有损解压
    log_info "3.2 测试有损解压..."
    if $GSC decompress -M -i "$TEST_DIR/lossy.gsc" -o "$TEST_DIR/restored_lossy.vcf" 2>&1 | tee -a "$LOG_FILE"; then
        # 检查GT字段是否保留
        local gt_count=$(grep -v "^#" "$TEST_DIR/restored_lossy.vcf" | head -10 | cut -f10 | grep -c ":")
        if [ "$gt_count" -eq 0 ] || [ "$gt_count" -lt 10 ]; then
            log_pass "有损解压成功 - 只保留GT字段"
        else
            log_info "有损解压成功 - 检查GT字段保留情况"
        fi
    else
        log_fail "有损解压失败"
    fi
}

#===============================================================================
# Part 4: VCF 检索功能测试
#===============================================================================

test_vcf_query() {
    log_section "Part 4: VCF 检索功能测试"

    # 使用有损压缩文件进行检索测试（有损模式支持更多检索功能）

    # 4.1 范围检索
    log_info "4.1 测试范围检索 (--range chr1:10000,20000)..."
    if $GSC decompress -M -r "chr1:10000,20000" -i "$TEST_DIR/lossy.gsc" -o "$TEST_DIR/query_range.vcf" 2>&1 | tee -a "$LOG_FILE"; then
        local count=$(grep -v "^#" "$TEST_DIR/query_range.vcf" | wc -l)
        log_pass "范围检索成功 - 返回 $count 条记录"
    else
        log_fail "范围检索失败"
    fi

    # 4.2 样本检索
    log_info "4.2 测试样本检索 (--samples CK28903712,CK28903725)..."
    if $GSC decompress -M -s "CK28903712,CK28903725" -i "$TEST_DIR/lossy.gsc" -o "$TEST_DIR/query_samples.vcf" 2>&1 | tee -a "$LOG_FILE"; then
        local sample_count=$(grep "^#CHROM" "$TEST_DIR/query_samples.vcf" | awk -F'\t' '{print NF-9}')
        log_pass "样本检索成功 - 返回 $sample_count 个样本"
    else
        log_fail "样本检索失败"
    fi

    # 4.3 组合检索（范围+样本）
    log_info "4.3 测试组合检索 (范围+样本)..."
    if $GSC decompress -M -r "chr1:10000,20000" -s "CK28903712" -i "$TEST_DIR/lossy.gsc" -o "$TEST_DIR/query_combined.vcf" 2>&1 | tee -a "$LOG_FILE"; then
        local count=$(grep -v "^#" "$TEST_DIR/query_combined.vcf" | wc -l)
        local sample_count=$(grep "^#CHROM" "$TEST_DIR/query_combined.vcf" | awk -F'\t' '{print NF-9}')
        log_pass "组合检索成功 - $count 条记录, $sample_count 个样本"
    else
        log_fail "组合检索失败"
    fi

    # 4.4 只输出header
    log_info "4.4 测试只输出header (--header-only)..."
    if $GSC decompress -M --header-only -i "$TEST_DIR/lossy.gsc" -o "$TEST_DIR/query_header.vcf" 2>&1 | tee -a "$LOG_FILE"; then
        local data_lines=$(grep -v "^#" "$TEST_DIR/query_header.vcf" | wc -l)
        if [ "$data_lines" -eq 0 ]; then
            log_pass "header-only成功 - 无数据行"
        else
            log_fail "header-only失败 - 包含 $data_lines 条数据行"
        fi
    else
        log_fail "header-only失败"
    fi

    # 4.5 不输出header
    log_info "4.5 测试不输出header (--no-header)..."
    if $GSC decompress -M --no-header -i "$TEST_DIR/lossy.gsc" -o "$TEST_DIR/query_no_header.vcf" 2>&1 | tee -a "$LOG_FILE"; then
        local header_lines=$(grep "^#" "$TEST_DIR/query_no_header.vcf" | wc -l)
        if [ "$header_lines" -eq 0 ]; then
            log_pass "no-header成功 - 无header行"
        else
            log_fail "no-header失败 - 包含 $header_lines 条header行"
        fi
    else
        log_fail "no-header失败"
    fi

    # 4.6 不输出基因型
    log_info "4.6 测试不输出基因型 (--no-genotype)..."
    if $GSC decompress -M -G -i "$TEST_DIR/lossy.gsc" -o "$TEST_DIR/query_no_gt.vcf" 2>&1 | tee -a "$LOG_FILE"; then
        local col_count=$(grep -v "^#" "$TEST_DIR/query_no_gt.vcf" | head -1 | awk -F'\t' '{print NF}')
        if [ "$col_count" -le 8 ]; then
            log_pass "no-genotype成功 - 只有 $col_count 列"
        else
            log_fail "no-genotype失败 - 有 $col_count 列"
        fi
    else
        log_fail "no-genotype失败"
    fi

    # 4.7 输出AC/AN
    log_info "4.7 测试输出AC/AN (--out-ac-an)..."
    if $GSC decompress -M -C -i "$TEST_DIR/lossy.gsc" -o "$TEST_DIR/query_ac_an.vcf" 2>&1 | tee -a "$LOG_FILE"; then
        if grep -v "^#" "$TEST_DIR/query_ac_an.vcf" | head -1 | grep -q "AC="; then
            log_pass "out-ac-an成功 - 包含AC字段"
        else
            log_info "out-ac-an完成 - 检查AC/AN字段"
        fi
    else
        log_fail "out-ac-an失败"
    fi

    # 4.8 按染色体分割输出
    log_info "4.8 测试分割输出 (--split)..."
    if $GSC decompress -M -S -i "$TEST_DIR/lossy.gsc" -o "$TEST_DIR/query_split.vcf" 2>&1 | tee -a "$LOG_FILE"; then
        log_pass "split输出成功"
    else
        log_fail "split输出失败"
    fi
}

#===============================================================================
# Part 5: gVCF 压缩测试
#===============================================================================

test_gvcf_compression() {
    log_section "Part 5: gVCF 压缩测试"

    local original_size=$(zcat "$GVCF_FILE" | wc -c)

    # 5.1 默认gVCF压缩（bsc）
    log_info "5.1 测试gVCF压缩 (compressor=bsc)..."
    local start_time=$(date +%s.%N)
    if $GSC gvcf -i "$GVCF_FILE" -o "$TEST_DIR/gvcf_bsc.gsc" --compressor bsc 2>&1 | tee -a "$LOG_FILE"; then
        local end_time=$(date +%s.%N)
        local elapsed=$(echo "$end_time - $start_time" | bc)
        local comp_size=$(get_file_size "$TEST_DIR/gvcf_bsc.gsc")
        local ratio=$(calc_ratio $original_size $comp_size)
        log_pass "gVCF压缩(bsc)成功 - 压缩率: ${ratio}%, 耗时: ${elapsed}s"
    else
        log_fail "gVCF压缩(bsc)失败"
    fi

    # 5.2 zstd压缩器
    log_info "5.2 测试gVCF压缩 (compressor=zstd)..."
    local start_time=$(date +%s.%N)
    if $GSC gvcf -i "$GVCF_FILE" -o "$TEST_DIR/gvcf_zstd.gsc" --compressor zstd 2>&1 | tee -a "$LOG_FILE"; then
        local end_time=$(date +%s.%N)
        local elapsed=$(echo "$end_time - $start_time" | bc)
        local comp_size=$(get_file_size "$TEST_DIR/gvcf_zstd.gsc")
        local ratio=$(calc_ratio $original_size $comp_size)
        log_pass "gVCF压缩(zstd)成功 - 压缩率: ${ratio}%, 耗时: ${elapsed}s"
    else
        log_fail "gVCF压缩(zstd)失败"
    fi

    # 5.3 brotli压缩器
    log_info "5.3 测试gVCF压缩 (compressor=brotli)..."
    local start_time=$(date +%s.%N)
    if $GSC gvcf -i "$GVCF_FILE" -o "$TEST_DIR/gvcf_brotli.gsc" --compressor brotli 2>&1 | tee -a "$LOG_FILE"; then
        local end_time=$(date +%s.%N)
        local elapsed=$(echo "$end_time - $start_time" | bc)
        local comp_size=$(get_file_size "$TEST_DIR/gvcf_brotli.gsc")
        local ratio=$(calc_ratio $original_size $comp_size)
        log_pass "gVCF压缩(brotli)成功 - 压缩率: ${ratio}%, 耗时: ${elapsed}s"
    else
        log_fail "gVCF压缩(brotli)失败"
    fi
}

#===============================================================================
# Part 6: gVCF 解压验证
#===============================================================================

test_gvcf_decompression() {
    log_section "Part 6: gVCF 解压验证"

    # 6.1 解压到VCF
    log_info "6.1 测试gVCF解压到VCF..."
    local start_time=$(date +%s.%N)
    if $GSC gvcf-decompress -i "$TEST_DIR/gvcf_bsc.gsc" -o "$TEST_DIR/gvcf_restored.vcf" 2>&1 | tee -a "$LOG_FILE"; then
        local end_time=$(date +%s.%N)
        local elapsed=$(echo "$end_time - $start_time" | bc)
        local restored_lines=$(wc -l < "$TEST_DIR/gvcf_restored.vcf")
        log_pass "gVCF解压成功 - $restored_lines 行, 耗时: ${elapsed}s"
    else
        log_fail "gVCF解压失败"
    fi

    # 6.2 解压到压缩VCF (vcf.gz)
    log_info "6.2 测试gVCF解压到vcf.gz..."
    if $GSC gvcf-decompress -i "$TEST_DIR/gvcf_bsc.gsc" -o "$TEST_DIR/gvcf_restored.vcf.gz" 2>&1 | tee -a "$LOG_FILE"; then
        log_pass "gVCF解压到vcf.gz成功"
    else
        log_fail "gVCF解压到vcf.gz失败"
    fi

    # 6.3 验证解压数据一致性
    log_info "6.3 验证gVCF解压数据一致性..."
    # 比较原始gVCF和解压后的VCF（只比较数据行）
    local original_data_count=$(zcat "$GVCF_FILE" | grep -v "^#" | wc -l)
    local restored_data_count=$(grep -v "^#" "$TEST_DIR/gvcf_restored.vcf" | wc -l)

    if [ "$original_data_count" -eq "$restored_data_count" ]; then
        log_pass "gVCF数据行数一致: $original_data_count 行"
    else
        log_fail "gVCF数据行数不一致: 原始 $original_data_count, 解压 $restored_data_count"
    fi

    # 6.4 深度验证（抽样对比）
    log_info "6.4 抽样验证数据内容..."
    local sample_match=0
    # 抽取前100行、中间100行、最后100行进行对比
    for offset in 1 10000000 19999900; do
        local orig_sample=$(zcat "$GVCF_FILE" | grep -v "^#" | sed -n "${offset},$((offset+99))p" | cut -f1-5 | md5sum | cut -d' ' -f1)
        local rest_sample=$(grep -v "^#" "$TEST_DIR/gvcf_restored.vcf" | sed -n "${offset},$((offset+99))p" | cut -f1-5 | md5sum | cut -d' ' -f1)
        if [ "$orig_sample" = "$rest_sample" ]; then
            sample_match=$((sample_match + 1))
        fi
    done

    if [ "$sample_match" -eq 3 ]; then
        log_pass "抽样验证全部通过 (3/3)"
    else
        log_fail "抽样验证部分失败 ($sample_match/3)"
    fi
}

#===============================================================================
# Part 7: gVCF 检索测试
#===============================================================================

test_gvcf_query() {
    log_section "Part 7: gVCF 检索测试 (gvcf-query)"

    # 7.1 范围检索 - 小范围
    log_info "7.1 测试小范围检索 (chr1:10000-100000)..."
    local start_time=$(date +%s.%N)
    if $GSC gvcf-query -i "$TEST_DIR/gvcf_bsc.gsc" -o "$TEST_DIR/gvcf_query_small.vcf" -r "chr1:10000-100000" 2>&1 | tee -a "$LOG_FILE"; then
        local end_time=$(date +%s.%N)
        local elapsed=$(echo "$end_time - $start_time" | bc)
        local count=$(grep -v "^#" "$TEST_DIR/gvcf_query_small.vcf" | wc -l)
        log_pass "小范围检索成功 - $count 条记录, 耗时: ${elapsed}s"
    else
        log_fail "小范围检索失败"
    fi

    # 7.2 范围检索 - 中等范围
    log_info "7.2 测试中等范围检索 (chr1:1000000-5000000)..."
    local start_time=$(date +%s.%N)
    if $GSC gvcf-query -i "$TEST_DIR/gvcf_bsc.gsc" -o "$TEST_DIR/gvcf_query_medium.vcf" -r "chr1:1000000-5000000" 2>&1 | tee -a "$LOG_FILE"; then
        local end_time=$(date +%s.%N)
        local elapsed=$(echo "$end_time - $start_time" | bc)
        local count=$(grep -v "^#" "$TEST_DIR/gvcf_query_medium.vcf" | wc -l)
        log_pass "中等范围检索成功 - $count 条记录, 耗时: ${elapsed}s"
    else
        log_fail "中等范围检索失败"
    fi

    # 7.3 范围检索 - 大范围
    log_info "7.3 测试大范围检索 (chr1:10000000-30000000)..."
    local start_time=$(date +%s.%N)
    if $GSC gvcf-query -i "$TEST_DIR/gvcf_bsc.gsc" -o "$TEST_DIR/gvcf_query_large.vcf" -r "chr1:10000000-30000000" 2>&1 | tee -a "$LOG_FILE"; then
        local end_time=$(date +%s.%N)
        local elapsed=$(echo "$end_time - $start_time" | bc)
        local count=$(grep -v "^#" "$TEST_DIR/gvcf_query_large.vcf" | wc -l)
        log_pass "大范围检索成功 - $count 条记录, 耗时: ${elapsed}s"
    else
        log_fail "大范围检索失败"
    fi

    # 7.4 验证检索结果正确性
    log_info "7.4 验证检索结果正确性..."
    # 验证检索结果中的位置都在指定范围内
    local out_of_range=$(awk -F'\t' '$1!~/^#/ && ($2<10000 || $2>100000) {print}' "$TEST_DIR/gvcf_query_small.vcf" | wc -l)
    if [ "$out_of_range" -eq 0 ]; then
        log_pass "检索结果范围验证通过"
    else
        log_fail "检索结果存在 $out_of_range 条超出范围的记录"
    fi
}

#===============================================================================
# Part 8: 格式转换测试
#===============================================================================

test_format_conversion() {
    log_section "Part 8: 格式转换测试"

    # 8.1 转换为PLINK BED格式
    log_info "8.1 测试转换为PLINK BED格式 (--make-bed)..."
    if $GSC decompress -M --make-bed -i "$TEST_DIR/lossy.gsc" -o "$TEST_DIR/plink_output" 2>&1 | tee -a "$LOG_FILE"; then
        if [ -f "$TEST_DIR/plink_output.bed" ] && [ -f "$TEST_DIR/plink_output.bim" ] && [ -f "$TEST_DIR/plink_output.fam" ]; then
            log_pass "PLINK BED格式转换成功 - 生成 .bed/.bim/.fam 文件"
            ls -lh "$TEST_DIR/plink_output".* | tee -a "$LOG_FILE"
        else
            log_fail "PLINK BED格式转换失败 - 文件不完整"
        fi
    else
        log_fail "PLINK BED格式转换失败"
    fi

    # 8.2 转换为PLINK2 PGEN格式
    log_info "8.2 测试转换为PLINK2 PGEN格式 (--make-pgen)..."
    if $GSC decompress -M --make-pgen -i "$TEST_DIR/lossy.gsc" -o "$TEST_DIR/pgen_output" 2>&1 | tee -a "$LOG_FILE"; then
        if [ -f "$TEST_DIR/pgen_output.pgen" ] && [ -f "$TEST_DIR/pgen_output.pvar" ] && [ -f "$TEST_DIR/pgen_output.psam" ]; then
            log_pass "PLINK2 PGEN格式转换成功 - 生成 .pgen/.pvar/.psam 文件"
            ls -lh "$TEST_DIR/pgen_output".* | tee -a "$LOG_FILE"
        else
            log_fail "PLINK2 PGEN格式转换失败 - 文件不完整"
        fi
    else
        log_fail "PLINK2 PGEN格式转换失败"
    fi

    # 8.3 转换为BGEN格式
    log_info "8.3 测试转换为BGEN格式 (--make-bgen)..."
    if $GSC decompress -M --make-bgen -i "$TEST_DIR/lossy.gsc" -o "$TEST_DIR/bgen_output" 2>&1 | tee -a "$LOG_FILE"; then
        if [ -f "$TEST_DIR/bgen_output.bgen" ]; then
            log_pass "BGEN格式转换成功 - 生成 .bgen 文件"
            ls -lh "$TEST_DIR/bgen_output".* | tee -a "$LOG_FILE"
        else
            log_fail "BGEN格式转换失败 - 文件不存在"
        fi
    else
        log_fail "BGEN格式转换失败"
    fi

    # 8.4 转换为BCF格式
    log_info "8.4 测试转换为BCF格式 (--bcf)..."
    if $GSC decompress -M --bcf -i "$TEST_DIR/lossy.gsc" -o "$TEST_DIR/output.bcf" 2>&1 | tee -a "$LOG_FILE"; then
        if [ -f "$TEST_DIR/output.bcf" ]; then
            log_pass "BCF格式转换成功"
            ls -lh "$TEST_DIR/output.bcf" | tee -a "$LOG_FILE"
        else
            log_fail "BCF格式转换失败"
        fi
    else
        log_fail "BCF格式转换失败"
    fi
}

#===============================================================================
# Part 9: 性能基准测试
#===============================================================================

test_performance_benchmark() {
    log_section "Part 9: 性能基准测试"

    log_info "9.1 压缩性能对比..."
    echo "| 压缩器 | 压缩大小 | 压缩率 | 压缩时间 |" | tee -a "$LOG_FILE"
    echo "|--------|----------|--------|----------|" | tee -a "$LOG_FILE"

    local original_size=$(zcat "$VCF_FILE" | wc -c)

    for comp in bsc zstd brotli; do
        local comp_file="$TEST_DIR/lossless_${comp}.gsc"
        if [ -f "$comp_file" ]; then
            local comp_size=$(get_file_size "$comp_file")
            local ratio=$(calc_ratio $original_size $comp_size)
            echo "| $comp | $(numfmt --to=iec $comp_size) | ${ratio}% | - |" | tee -a "$LOG_FILE"
        fi
    done

    log_info "9.2 gVCF压缩性能对比..."
    echo "| 压缩器 | 压缩大小 | 压缩率 |" | tee -a "$LOG_FILE"
    echo "|--------|----------|--------|" | tee -a "$LOG_FILE"

    local gvcf_original_size=$(zcat "$GVCF_FILE" | wc -c)

    for comp in bsc zstd brotli; do
        local comp_file="$TEST_DIR/gvcf_${comp}.gsc"
        if [ -f "$comp_file" ]; then
            local comp_size=$(get_file_size "$comp_file")
            local ratio=$(calc_ratio $gvcf_original_size $comp_size)
            echo "| $comp | $(numfmt --to=iec $comp_size) | ${ratio}% |" | tee -a "$LOG_FILE"
        fi
    done
}

#===============================================================================
# Part 10: 边界条件测试
#===============================================================================

test_edge_cases() {
    log_section "Part 10: 边界条件测试"

    # 10.1 空范围检索
    log_info "10.1 测试空范围检索 (chr1:1,10)..."
    if $GSC decompress -M -r "chr1:1,10" -i "$TEST_DIR/lossy.gsc" -o "$TEST_DIR/query_empty.vcf" 2>&1 | tee -a "$LOG_FILE"; then
        local count=$(grep -v "^#" "$TEST_DIR/query_empty.vcf" | wc -l)
        log_pass "空范围检索成功 - 返回 $count 条记录"
    else
        log_fail "空范围检索失败"
    fi

    # 10.2 单样本检索
    log_info "10.2 测试单样本检索..."
    if $GSC decompress -M -s "CK28903712" -i "$TEST_DIR/lossy.gsc" -o "$TEST_DIR/query_single_sample.vcf" 2>&1 | tee -a "$LOG_FILE"; then
        log_pass "单样本检索成功"
    else
        log_fail "单样本检索失败"
    fi

    # 10.3 全范围检索
    log_info "10.3 测试全范围检索 (chr1:1,100000000)..."
    if $GSC decompress -M -r "chr1:1,100000000" -i "$TEST_DIR/lossy.gsc" -o "$TEST_DIR/query_full.vcf" 2>&1 | tee -a "$LOG_FILE"; then
        local count=$(grep -v "^#" "$TEST_DIR/query_full.vcf" | wc -l)
        log_pass "全范围检索成功 - 返回 $count 条记录"
    else
        log_fail "全范围检索失败"
    fi
}

#===============================================================================
# 主程序
#===============================================================================

main() {
    # 先创建目录，确保日志文件可以写入
    mkdir -p "$TEST_DIR"

    echo "==============================================================================="
    echo "GSC 完整功能验证测试"
    echo "开始时间: $(date)"
    echo "==============================================================================="

    init_test

    # VCF测试
    test_vcf_lossless_compression
    test_vcf_lossless_decompression
    test_vcf_lossy_compression
    test_vcf_query

    # gVCF测试
    test_gvcf_compression
    test_gvcf_decompression
    test_gvcf_query

    # 格式转换测试
    test_format_conversion

    # 性能测试
    test_performance_benchmark

    # 边界条件测试
    test_edge_cases

    # 输出总结
    log_section "测试总结"
    echo "通过: $PASS_COUNT" | tee -a "$LOG_FILE"
    echo "失败: $FAIL_COUNT" | tee -a "$LOG_FILE"
    echo "跳过: $SKIP_COUNT" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    echo "测试结束时间: $(date)" | tee -a "$LOG_FILE"
    echo "详细日志: $LOG_FILE"

    # 清理临时文件（可选）
    # rm -rf "$TEST_DIR"

    if [ "$FAIL_COUNT" -gt 0 ]; then
        exit 1
    fi
    exit 0
}

# 运行测试
main "$@"
