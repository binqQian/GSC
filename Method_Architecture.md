# GSC算法方法详解 - 面向代码理解

## 1. 核心设计思想

### 1.1 问题定义
**输入**：VCF/BCF文件（包含基因型数据和注释信息）  
**输出**：压缩的.gsc文件，支持快速随机访问  
**目标**：高压缩率 + 快速查询

### 1.2 设计策略
```
VCF文件 = 基因型数据 + 注释数据
           ↓              ↓
    层次化块压缩    字段分类压缩
           ↓              ↓
         BSC后端压缩
           ↓
       .gsc文件
```

**关键洞察**：
- 基因型数据经聚类后呈现**稀疏性**
- 相邻单倍体具有**相似性**（小汉明距离）
- 注释字段有**结构化特征**

---

## 2. 整体算法架构

```
┌─────────────────────────────────────────────────────────┐
│                    VCF/BCF 输入                         │
└──────────────────────┬──────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────┐
│                   预处理阶段                            │
│  ├─ 多等位基因拆分（保持无损）                         │
│  └─ 基因型二进制编码（2bit per genotype）             │
└──────────────────────┬──────────────────────────────────┘
                       ↓
        ┌──────────────┴──────────────┐
        ↓                              ↓
┌──────────────────┐          ┌──────────────────┐
│ 基因型数据压缩    │          │ 注释数据压缩      │
│                  │          │                  │
│ Step 1: 分块     │          │ 固定字段流:       │
│   └→ 块大小=s    │          │  CHROM, POS...    │
│                  │          │                  │
│ Step 2: 聚类     │          │ 可变字段:         │
│   └→ 汉明距离    │          │  INFO/FORMAT      │
│                  │          │  └→ DAG结构       │
│ Step 3: 稀疏化   │          │                  │
│   └→ XOR操作     │          │ 每个字段独立     │
│                  │          │ 8MB块压缩        │
│ Step 4: 编码     │          │                  │
│   └→ 位置索引    │          │                  │
│                  │          │                  │
│ Step 5: 合并     │          │                  │
│   └→ Chunk       │          │                  │
└────────┬─────────┘          └────────┬─────────┘
         └───────────┬─────────────────┘
                     ↓
         ┌───────────────────────┐
         │   BSC 后端压缩        │
         └───────────┬───────────┘
                     ↓
         ┌───────────────────────┐
         │   .gsc 压缩文件       │
         │   + B树索引           │
         └───────────────────────┘
```

---

## 3. 预处理阶段

### 3.1 多等位基因拆分
**输入**: 变异行 (POS, REF, ALT字段可能包含多个等位基因)
**输出**: 多个双等位基因变异行

**处理逻辑**:
1. 检测ALT字段中的逗号分隔符
2. 对每个备选等位基因创建独立变异:
   - 第一个: ALT = "<N>,第一个等位基因"
   - 后续: ALT = "<M>,对应等位基因"
3. 在REF前添加序号前缀保证同位置变异的顺序

**标记含义**:
- `<N>`: 标识第一个拆分的变异
- `<M>`: 标识后续拆分的变异
- 数字前缀: 保持原始顺序用于恢复

**代码实现要点**:
- 遍历每个变异行
- 解析ALT字段
- 生成新的变异行并保持元数据一致

### 2.2 基因型二进制编码

**编码映射表**:
```
VCF表示  →  2-bit编码  →  十进制
  0|0    →    00      →    0
  0|1    →    01      →    1  
  1|1    →    11      →    3
  ./.    →    10      →    2
```

**数据组织**:
- 每个变异 → 2个长度为h的比特向量
- 比特向量1: 存储第一个bit
- 比特向量2: 存储第二个bit

**存储结构**:
```
对于每个块 (s个变异, h个单倍体):
  bit_matrix[2s][h]
  - 行索引: 2s (每个变异2行)
  - 列索引: h (每个单倍体1列)
```

**代码实现要点**:
- 字节对齐: 使用uint8_t存储，每字节4个基因型
- 位操作: 使用位移和掩码提取/设置基因型

**关键点**：
- `<N>` 标记第一个拆分变异
- `<M>` 标记后续拆分变异  
- REF前的数字保证同位置变异的顺序

### 3.2 基因型二进制编码

**编码方案**：每个基因型 → 2个比特
```
基因型    →    2-bit编码
  0|0     →      00
  0|1     →      01
  1|1     →      11
  ./.     →      10
```

**数据结构**：
```
对于 h 个单倍体，s 个变异：
- 每个变异 → 2个长度为h的比特向量
- 每个块 → 2s个比特向量
- 比特矩阵大小 = 2s × h
```

---

## 4. 基因型数据压缩（核心方法）

### 4.1 Step 1: 数据分块

**块大小s的确定**：
```python
if h < 65536:
    s = h
else:
    s = 65536
```

**块的组成**：
```
Block = {
    bit_vectors: [2s][h],  # 2s个比特向量，每个长度h
    P: [h],                # 排列顺序数组
}
```

**分块策略**：
- 每个染色体独立分块
- 块之间可以并行处理
- 最后一个块可能不满

### 4.2 Step 2: 单倍体聚类

**目标**：将相似的列（单倍体）排列在一起

**算法**：基于汉明距离的最近邻聚类
```python
def haplotype_clustering(block):
    # 初始化
    P = []  # 排列顺序
    used = [False] * h
    
    # 选择起始列（如第一列）
    current = 0
    P.append(current)
    used[current] = True
    
    # 迭代选择最近邻
    for i in range(1, h):
        min_distance = infinity
        nearest = -1
        
        for j in range(h):
            if not used[j]:
                dist = hamming_distance(block[:, current], block[:, j])
                if dist < min_distance:
                    min_distance = dist
                    nearest = j
        
        P.append(nearest)
        used[nearest] = True
        current = nearest
    
    # 根据P重排列
    return reorder_columns(block, P), P
```

**示意图**：
```
原始块:                聚类后:
Col0 Col1 Col2        Col1 Col0 Col2
 0    1    0           1    0    0
 1    1    0    →      1    1    0
 0    0    1           0    0    1
 1    1    0           1    1    0

P = [1, 0, 2]  # 新顺序
```

### 4.3 Step 3: 稀疏化（XOR操作）

**核心思想**：利用相邻列的相似性，通过XOR减少1的数量

**判断条件**：
```python
# 总汉明距离
D = sum(hamming_distance(col[i], col[i+1]) for i in range(h-1))
D += hamming_weight(col[h-1])  # 最后一列的汉明权重

# 稀疏度判断
num_ones = count_ones(block)
if D < num_ones:
    apply_sparsification()
```

**稀疏化算法**：
```python
def sparsify_block(block):
    # 每8列一组（字节对齐）
    for group_start in range(0, h, 8):
        group_end = min(group_start + 8, h)
        
        for col in range(group_start + 1, group_end):
            prev_col = col - 1
            
            # 计算汉明距离和权重
            hamming_dist = hamming_distance(block[:, col], block[:, prev_col])
            hamming_weight = count_ones(block[:, col])
            
            # 如果XOR后更稀疏，则替换
            if hamming_dist < hamming_weight:
                block[:, col] = block[:, col] XOR block[:, prev_col]
```

**示意图**：
```
组内列:   Col0  Col1  Col2  Col3
原始:      1     1     0     1
           0     1     1     0
           1     0     1     1

操作:    保持  XOR   XOR   XOR
               ↓     ↓     ↓
结果:      1     0     1     0
           0     1     0     1
           1     1     1     0

减少了1的数量！
```

### 4.4 Step 4: 稀疏编码

**子步骤1：标记特殊向量**
```python
def mark_special_vectors(block):
    Vzero = []  # 全零向量标记
    Vcopy = []  # 重复向量标记
    Icopy = []  # 重复向量的原始位置
    
    for i in range(2s):
        if is_all_zero(block[i]):
            Vzero[i] = 1
        elif is_duplicate(block[i]):
            Vcopy[i] = 1
            Icopy.append(find_original_index(block[i]))
        else:
            Vzero[i] = 0
            Vcopy[i] = 0
    
    return Vzero, Vcopy, Icopy
```

**子步骤2：提取稀疏位置**
```python
def extract_sparse_positions(block, Vzero, Vcopy):
    Vindex = []
    
    for i in range(2s):
        if Vzero[i] == 1 or Vcopy[i] == 1:
            continue  # 跳过全零和重复向量
        
        # 记录该行中1的位置
        positions = [j for j in range(h) if block[i][j] == 1]
        Vindex.extend(positions)
        Vindex.append(0)  # 0作为行分隔符
    
    return Vindex
```

**示意图**：
```
稀疏矩阵:       位置索引:
Row0: 1 0 1 0   → [0, 2, 0,
Row1: 0 0 0 0   → (全零，跳过)
Row2: 0 1 1 0   → 1, 2, 0,
Row3: 1 0 1 0   → (与Row0重复，跳过)
                   Icopy记录: Row3→Row0]
```

**子步骤3：差分和变长编码**
```python
def encode_vindex(Vindex):
    # 差分编码
    delta_encoded = [Vindex[0]]
    for i in range(1, len(Vindex)):
        delta_encoded.append(Vindex[i] - Vindex[i-1])
    
    # 变长编码（Variable Byte Encoding）
    return variable_length_encode(delta_encoded)
```

### 4.5 Step 5: 排列顺序P的存储

**问题**：P数组可能占用大量空间

**解决方案：重排映射**

**方法**：
```python
def encode_permutation(P, POS, REF):
    # 根据P重排POS和REF
    reordered_POS = [POS[P[i]] for i in range(h)]
    reordered_REF = [REF[P[i]] for i in range(h)]
    
    # 只存储重排后的POS和REF
    # P可以通过对reordered_POS排序恢复
    return reordered_POS, reordered_REF

def decode_permutation(reordered_POS):
    # 恢复P
    indexed_POS = [(pos, i) for i, pos in enumerate(reordered_POS)]
    sorted_indexed = sorted(indexed_POS, key=lambda x: x[0])
    P = [i for _, i in sorted_indexed]
    return P
```

**示意图**：
```
原始:
P = [2, 0, 1]
POS = [100, 200, 300]
REF = [A, C, G]

重排后:
reordered_POS = [300, 100, 200]  # 按P[0]=2, P[1]=0, P[2]=1
reordered_REF = [G, A, C]

恢复P:
对reordered_POS排序 → [(100,1), (200,2), (300,0)]
提取索引 → P = [1, 2, 0] （实际上是逆排列）
```

**特殊情况**：
- 最后一个块（变异数 ≠ h）：直接存储P（变长编码）
- 完整块（变异数 = h）：使用重排映射

### 4.6 块合并为Chunk

**Chunk的定义**：
```python
chunk_size = 65536  # 约65k个变异
blocks_per_chunk = floor(chunk_size / s)

Chunk = {
    blocks: [blocks_per_chunk],
    compressed_data: bytes  # BSC压缩后的数据
}
```

**为什么需要Chunk？**
- 平衡压缩率和查询速度
- 查询时只需解压相关chunk
- 减少随机访问开销

---

## 5. 注释数据压缩

### 5.1 固定字段流

**字段列表**：CHROM, POS, ID, REF, ALT, QUAL, FILTER

**处理方式**：
```python
def compress_fixed_fields(vcf_file):
    for field in ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']:
        stream = extract_field(vcf_file, field)
        
        # 按基因型块大小分块
        blocks = split_into_blocks(stream, block_size_from_genotype)
        
        # 每块独立压缩
        compressed_blocks = [BSC.compress(block) for block in blocks]
        
        save_stream(field, compressed_blocks)
```

### 5.2 INFO和FORMAT字段的DAG处理

**问题**：子字段顺序在不同变异行中可能不同

**解决方案**：构建DAG记录顺序关系

**算法步骤**：

**Step 1: 构建keys表**
```python
def build_keys_table(vcf_header):
    keys_table = []
    for field in ['INFO', 'FORMAT']:
        for id, metadata in parse_header(field):
            key_id = get_key_id_from_htllib(id)
            keys_table.append({
                'Field': field,
                'ID': id,
                'key_id': key_id
            })
    return keys_table
```

**Step 2: 映射每行的子字段**
```python
def map_variant_to_keys(variant_line, keys_table):
    info_ids = extract_info_ids(variant_line)
    format_ids = extract_format_ids(variant_line)
    
    key_sequence = []
    for id in info_ids + format_ids:
        key_id = lookup_key_id(id, keys_table)
        key_sequence.append(key_id)
    
    return key_sequence
```

**Step 3: 构建DAG**
```python
class DAG:
    def __init__(self):
        self.graph = {}  # {node: [children]}
    
    def add_sequence(self, sequence):
        for i in range(len(sequence) - 1):
            current = sequence[i]
            next_node = sequence[i + 1]
            
            if current not in self.graph:
                self.graph[current] = set()
            self.graph[current].add(next_node)
    
    def topological_sort(self):
        # 拓扑排序恢复顺序
        visited = set()
        result = []
        
        def dfs(node):
            if node in visited:
                return
            visited.add(node)
            if node in self.graph:
                for child in self.graph[node]:
                    dfs(child)
            result.append(node)
        
        for node in self.graph:
            dfs(node)
        
        return reversed(result)
```

**示意图**：
```
变异行序列:
Row1: NS → DP → AF
Row2: NS → AF → DB
Row3: DP → AF

构建DAG:
    NS → DP → AF → DB
         ↓     ↑
         └─────┘

拓扑排序恢复keys表顺序
```

**Step 4: 字段值压缩**
```python
def compress_field_values(field_values):
    # 按数据类型组织
    streams = group_by_type(field_values)
    
    # 每个流分8MB块
    for stream in streams:
        blocks = split_into_blocks(stream, 8 * 1024 * 1024)
        compressed = [BSC.compress(block) for block in blocks]
        save_stream(compressed)
```

---

## 6. 解压与查询机制

### 6.1 数据结构：B树索引

**索引层次**：
```
B-Tree Index {
    Chromosome Level: {
        chr1 → Chunk[0..N1]
        chr2 → Chunk[0..N2]
        ...
    }
    
    Chunk Level: {
        Chunk_i → {
            start_position: int
            end_position: int
            blocks: Block[0..M]
            file_offset: long
        }
    }
    
    Block Level: {
        Block_j → {
            start_variant: int
            end_variant: int
            file_offset: long
        }
    }
}
```

### 6.2 变异查询（Variant Query）

**输入参数**：
- 染色体ID
- 位置范围 [start, end]
- 其他过滤条件（可选）

**查询流程**：
```python
def query_variants(chr_id, start_pos, end_pos):
    # Step 1: 定位相关chunk
    chunks = index.find_chunks(chr_id, start_pos, end_pos)
    
    # Step 2: 对每个chunk
    results = []
    for chunk in chunks:
        # 定位相关block
        blocks = chunk.find_blocks(start_pos, end_pos)
        
        # Step 3: 解压block
        for block in blocks:
            # 解压chunk（如果未缓存）
            decompressed_chunk = BSC.decompress(chunk.data)
            
            # 提取block数据
            block_data = extract_block(decompressed_chunk, block.offset)
            
            # Step 4: 解码block
            decoded_block = decode_block(block_data)
            # 包括：
            #   - 恢复P
            #   - 处理Vzero, Vcopy, Vindex
            #   - 反XOR
            #   - 反聚类
            
            # Step 5: 筛选变异
            variants = filter_variants(decoded_block, start_pos, end_pos)
            results.extend(variants)
    
    return results
```

**解码block的详细步骤**：
```python
def decode_block(block_data):
    # 1. 解析元数据
    Vzero = block_data.vzero
    Vcopy = block_data.vcopy
    Icopy = block_data.icopy
    Vindex_encoded = block_data.vindex
    reordered_POS = block_data.pos
    
    # 2. 恢复P
    P = decode_permutation(reordered_POS)
    
    # 3. 解码Vindex
    Vindex = variable_length_decode(Vindex_encoded)
    delta_decoded = inverse_delta_encoding(Vindex)
    
    # 4. 重建稀疏矩阵
    sparse_matrix = reconstruct_from_positions(delta_decoded)
    
    # 5. 恢复全零和重复向量
    full_matrix = []
    sparse_idx = 0
    for i in range(2s):
        if Vzero[i] == 1:
            full_matrix.append([0] * h)
        elif Vcopy[i] == 1:
            original_idx = Icopy[copy_counter]
            full_matrix.append(full_matrix[original_idx].copy())
            copy_counter += 1
        else:
            full_matrix.append(sparse_matrix[sparse_idx])
            sparse_idx += 1
    
    # 6. 反XOR（逆稀疏化）
    for group_start in range(0, h, 8):
        for col in range(group_start + 1, min(group_start + 8, h)):
            if was_xored(col):  # 需要记录哪些列被XOR了
                full_matrix[:, col] = full_matrix[:, col] XOR full_matrix[:, col-1]
    
    # 7. 反聚类（逆排列）
    inverse_P = compute_inverse_permutation(P)
    original_order = reorder_columns(full_matrix, inverse_P)
    
    return original_order
```

### 6.3 样本查询（Sample Query）

**输入参数**：
- 样本名称列表
- 位置范围（可选）

**查询流程**：
```python
def query_samples(sample_names, chr_id=None, start_pos=None, end_pos=None):
    # Step 1: 确定样本对应的单倍体列
    sample_haplotypes = []
    for name in sample_names:
        hap_indices = get_haplotype_indices(name)  # 通常每个样本2个单倍体
        sample_haplotypes.extend(hap_indices)
    
    # Step 2: 定位相关chunk和block（如果指定位置）
    if chr_id:
        chunks = index.find_chunks(chr_id, start_pos, end_pos)
    else:
        chunks = index.get_all_chunks()
    
    # Step 3: 部分解压
    results = []
    for chunk in chunks:
        decompressed = BSC.decompress(chunk.data)
        
        for block in chunk.blocks:
            block_data = extract_block(decompressed, block.offset)
            
            # 获取P
            P = decode_permutation(block_data.pos)
            
            # 确定需要解压的字节位置
            byte_positions = []
            for hap_idx in sample_haplotypes:
                # 找到该单倍体在排列后的位置
                permuted_idx = P.index(hap_idx)
                byte_pos = permuted_idx // 8
                byte_positions.append(byte_pos)
            
            # 只解压相关字节
            partial_data = extract_bytes(block_data, byte_positions)
            decoded = decode_partial_block(partial_data, sample_haplotypes, P)
            
            results.append(decoded)
    
    return results
```

**字节级查找表**：
```python
# 预计算查找表加速解码
LOOKUP_TABLE = {}
for byte_value in range(256):
    for bit_pos in range(8):
        genotype = extract_2bits(byte_value, bit_pos)
        LOOKUP_TABLE[(byte_value, bit_pos)] = genotype
```

### 6.4 转换为PLINK格式

**PLINK BED格式要求**：
- 每个基因型用2bit表示
- 0|0 → 纯合主要等位基因
- 1|1 → 纯合次要等位基因
- 0|1 或 1|0 → 杂合

**转换映射**：
```python
def convert_to_plink_bed(gsc_file):
    # VCF中: 0=主要, 1=次要, 2=其他
    # PLINK中: 0=次要, 1=主要
    
    mapping = {
        (0, 0): 0b00,  # 纯合主要 → 00
        (0, 1): 0b10,  # 杂合 → 10
        (1, 0): 0b10,  # 杂合 → 10
        (1, 1): 0b11,  # 纯合次要 → 11
        (2, 2): 0b00,  # 其他等位基因视为主要
        'missing': 0b01  # 缺失
    }
    
    for chunk in gsc_file.chunks:
        decoded = decode_genotype_chunk(chunk)
        
        # 处理多等位基因：2→0
        normalized = normalize_multiallelic(decoded)
        
        # 应用映射
        bed_data = apply_mapping(normalized, mapping)
        
        write_bed_chunk(bed_data)
```

**直接转换的优势**：
- 无需完全解压VCF
- 在解压genotype时直接应用映射
- 避免BCFtools预处理步骤

---

## 7. 关键数据结构设计

### 7.1 Block结构
```cpp
struct Block {
    // 原始数据
    BitVector bit_vectors[2*s][h];
    
    // 聚类信息
    int P[h];  // 或通过POS恢复
    
    // 稀疏编码数据
    BitVector Vzero[2*s];    // 全零标记
    BitVector Vcopy[2*s];    // 重复标记
    vector<int> Icopy;       // 重复向量原始位置
    vector<int> Vindex;      // 稀疏位置索引（编码后）
    
    // 元数据
    int num_variants;
    int start_position;
    int end_position;
};
```

### 7.2 Chunk结构
```cpp
struct Chunk {
    int chunk_id;
    vector<Block> blocks;
    
    // 压缩数据
    ByteBuffer compressed_data;
    size_t compressed_size;
    
    // 索引信息
    int start_position;
    int end_position;
    long file_offset;
    
    // 元数据
    string chromosome;
    int num_blocks;
};
```

### 7.3 索引结构
```cpp
struct Index {
    // 染色体级索引
    map<string, vector<Chunk*>> chromosome_chunks;
    
    // B树节点
    struct BTreeNode {
        int start_pos;
        int end_pos;
        union {
            vector<BTreeNode*> children;  // 内部节点
            Chunk* chunk;                  // 叶子节点
        };
    };
    
    BTreeNode* root;
    
    // 查询接口
    vector<Chunk*> find_chunks(string chr, int start, int end);
    vector<Block*> find_blocks(Chunk* chunk, int start, int end);
};
```

### 7.4 DAG结构（INFO/FORMAT）
```cpp
struct DAG {
    // 图表示
    map<int, set<int>> adjacency_list;  // key_id → {后继key_ids}
    
    // keys表
    struct KeyEntry {
        string field;      // "INFO" or "FORMAT"
        string id;         // 子字段ID
        int key_id;        // HTSlib的key_id
    };
    vector<KeyEntry> keys_table;
    
    // 操作
    void add_sequence(vector<int> key_sequence);
    vector<int> topological_sort();
    vector<KeyEntry> recover_keys_table();
};
```

---

## 8. 算法复杂度分析

### 8.1 时间复杂度

**压缩阶段**：
```
预处理: O(V × h)  // V=变异数, h=单倍体数
分块: O(V)
聚类: O(h² × s)   // 每个块
稀疏化: O(h × s)
编码: O(非零元素数)
BSC压缩: O(数据量 × log(数据量))

总计: O(V × h + 块数 × h² × s)
```

**查询阶段**：
```
定位chunk: O(log(chunk数))
解压chunk: O(chunk大小)
解码block: O(h × s)
筛选变异: O(结果数)

总计: O(log(chunk数) + h × s + 结果数)
```

### 8.2 空间复杂度

**内存使用**：
```
块缓冲: O(h × s)
索引结构: O(V / chunk_size)
解压缓冲: O(chunk大小)
查找表: O(1)  // 预计算固定大小

峰值: O(h × s + chunk大小)
```

---

## 9. 代码实现映射

### 9.1 核心模块结构
```
GSC/
├── src/
│   ├── preprocess/
│   │   ├── multiallelic.cpp      # 多等位基因拆分
│   │   └── binary_encoding.cpp   # 2-bit编码
│   │
│   ├── genotype/
│   │   ├── block.cpp              # 块管理
│   │   ├── clustering.cpp         # 单倍体聚类
│   │   ├── sparsification.cpp     # XOR稀疏化
│   │   ├── sparse_encoding.cpp    # 位置索引编码
│   │   └── permutation.cpp        # P的存储与恢复
│   │
│   ├── annotation/
│   │   ├── fixed_fields.cpp       # 固定字段流
│   │   └── dag.cpp                # INFO/FORMAT的DAG
│   │
│   ├── compression/
│   │   ├── bsc_wrapper.cpp        # BSC接口
│   │   └── chunk.cpp              # Chunk管理
│   │
│   ├── index/
│   │   └── btree.cpp              # B树索引
│   │
│   ├── query/
│   │   ├── variant_query.cpp      # 变异查询
│   │   ├── sample_query.cpp       # 样本查询
│   │   └── decoder.cpp            # 解码器
│   │
│   └── format/
│       └── plink_converter.cpp    # PLINK转换
│
└── include/
    └── gsc/
        ├── structures.h           # 数据结构定义
        └── algorithms.h           # 算法接口
```

### 9.2 关键函数映射

**压缩流程**：
```cpp
// 主压缩函数
void compress_vcf(const string& input_vcf, const string& output_gsc) {
    VCFReader reader(input_vcf);
    GSCWriter writer(output_gsc);
    
    // 1. 预处理
    PreprocessedData data = preprocess(reader);
    
    // 2. 压缩genotype
    for (auto& chromosome : data.chromosomes) {
        vector<Block> blocks = split_into_blocks(chromosome, s);
        
        for (auto& block : blocks) {
            // 聚类
            auto [clustered, P] = haplotype_clustering(block);
            
            // 稀疏化
            auto sparsified = sparsify_block(clustered);
            
            // 编码
            auto encoded = sparse_encode(sparsified, P);
            
            blocks_encoded.push_back(encoded);
        }
        
        // 合并chunk
        vector<Chunk> chunks = merge_into_chunks(blocks_encoded);
        
        // BSC压缩
        for (auto& chunk : chunks) {
            chunk.data = BSC::compress(chunk.raw_data);
            writer.write_chunk(chunk);
        }
    }
    
    // 3. 压缩annotation
    compress_annotations(data.annotations, writer);
    
    // 4. 写入索引
    writer.write_index(build_index(chunks));
}
```

**查询流程**：
```cpp
// 变异查询
vector<Variant> query_variants(const GSCFile& gsc, 
                               const string& chr,
                               int start, int end) {
    // 1. 查找相关chunk
    auto chunks = gsc.index.find_chunks(chr, start, end);
    
    vector<Variant> results;
    for (auto* chunk : chunks) {
        // 2. 解压chunk
        auto decompressed = BSC::decompress(chunk->data);
        
        // 3. 查找相关block
        auto blocks = chunk->find_blocks(start, end);
        
        for (auto* block : blocks) {
            // 4. 解码block
            auto decoded = decode_block(decompressed, block);
            
            // 5. 筛选变异
            auto variants = filter_variants(decoded, start, end);
            results.insert(results.end(), variants.begin(), variants.end());
        }
    }
    
    return results;
}

// 解码block
Block decode_block(const ByteBuffer& data, const BlockMeta* meta) {
    // 恢复P
    auto P = decode_permutation(meta->reordered_pos);
    
    // 解码Vindex
    auto Vindex = decode_vindex(data.vindex_encoded);
    
    // 重建稀疏矩阵
    auto sparse = reconstruct_sparse_matrix(Vindex, data.vzero, data.vcopy, data.icopy);
    
    // 反XOR
    auto unsparsified = reverse_xor(sparse);
    
    // 反聚类
    auto original = reverse_clustering(unsparsified, P);
    
    return original;
}
```

### 9.3 数据流转示意
```
压缩:
VCF → preprocess() → split_blocks() → haplotype_clustering() 
    → sparsify_block() → sparse_encode() → merge_chunks() 
    → BSC::compress() → write_gsc()

解压:
GSC → read_index() → find_chunks() → BSC::decompress() 
    → extract_blocks() → decode_vindex() → reconstruct_sparse()
    → reverse_xor() → reverse_clustering() → decode_permutation()
    → output_vcf()

查询:
Query → find_chunks() → find_blocks() → partial_decompress()
    → decode_relevant_data() → filter_results() → return
```

---

## 10. 实现要点

### 10.1 位操作优化
```cpp
// 快速汉明距离计算
inline int hamming_distance(uint64_t a, uint64_t b) {
    return __builtin_popcountll(a ^ b);
}

// 快速提取2-bit genotype
inline int extract_genotype(uint8_t byte, int pos) {
    return (byte >> (pos * 2)) & 0b11;
}
```

### 10.2 内存管理
```cpp
// 使用内存池避免频繁分配
class BlockPool {
    vector<Block*> free_blocks;
    
public:
    Block* acquire() {
        if (free_blocks.empty()) {
            return new Block();
        }
        auto block = free_blocks.back();
        free_blocks.pop_back();
        return block;
    }
    
    void release(Block* block) {
        block->clear();
        free_blocks.push_back(block);
    }
};
```

### 10.3 并行处理
```cpp
// 块级并行压缩
void parallel_compress_blocks(vector<Block>& blocks) {
    #pragma omp parallel for
    for (size_t i = 0; i < blocks.size(); i++) {
        blocks[i] = process_block(blocks[i]);
    }
}
```

### 10.4 流式IO
```cpp
// 避免一次性加载整个文件
class StreamingVCFReader {
    ifstream file;
    buffer_queue<Variant> buffer;
    
public:
    optional<Variant> next() {
        if (buffer.empty()) {
            fill_buffer();
        }
        return buffer.pop();
    }
};
```

---

## 11. 总结

GSC通过以下技术实现高效压缩和快速查询：

1. **层次化设计**：块 → Chunk → BSC，平衡压缩率和访问速度
2. **聚类利用相似性**：减少相邻列差异
3. **XOR稀疏化**：进一步增强稀疏性
4. **智能编码**：位置索引 + 差分 + 变长编码
5. **巧妙的P存储**：重排映射避免额外开销
6. **DAG记录结构**：灵活处理可变顺序
7. **B树索引**：支持快速定位
8. **部分解压**：样本查询只解压相关部分

这些方法的组合使GSC在保持高压缩率的同时支持高效的随机访问。