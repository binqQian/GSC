# LLM提示词：VCF.gz压缩加速优化任务

## 任务目标
优化VCF.gz文件的压缩和解压缩性能，参考`./ref_code/vcf_Parser.h`中的实现，确保无损压缩/解压，并通过性能测试验证优化效果。

## 具体要求

### 1. 项目理解阶段
- 仔细阅读现有项目代码结构，理解当前VCF处理流程
- 分析`./ref_code/vcf_Parser.h`的优化思路和关键技术点
- 识别当前实现的性能瓶颈

### 2. 优化实施
- 参考`vcf_Parser.h`的优化方法，改进压缩/解压算法
- 可能的优化方向：
  - 使用更高效的zlib参数配置
  - 实现多线程并行处理
  - 优化内存缓冲区管理
  - 改进I/O操作策略
- **关键约束**：必须保证无损压缩和解压，数据完整性优先

### 3. 测试验证
使用两个测试文件进行全面测试：
- **大文件测试**：`./toy/final_subset.vcf.gz`
- **小文件测试**：`./toy/toy.vcf.gz`

测试内容包括：
```bash
# 对每个测试文件执行：
1. 解压缩测试（记录时间）
2. 压缩测试（记录时间）
3. 数据完整性验证（MD5/SHA256校验）
4. 与原始实现对比性能提升百分比
```

### 4. 性能对比报告
生成测试报告，包含：
- 优化前后的处理时间对比表格
- 速度提升百分比
- 内存使用情况（如有变化）
- 压缩率对比（确保一致）

### 5. Git提交
完成优化和测试后：
```bash
git add .
git commit -m "feat: optimize vcf.gz compression/decompression performance

- Reference vcf_Parser.h implementation
- Improve compression speed by X%
- Improve decompression speed by Y%
- Verified lossless compression with test files
- Tested on final_subset.vcf.gz and toy.vcf.gz"
git push
```

## 验证清单
- [ ] 代码能正确处理vcf和vcf.gz格式
- [ ] 压缩后文件可以完整解压
- [ ] 解压后数据与原始数据完全一致（校验和验证）
- [ ] 大文件和小文件测试均通过
- [ ] 性能提升已量化并记录
- [ ] 代码已提交到git仓库

## 输出要求
请提供：
1. 优化代码的关键修改说明
2. 性能测试对比数据表格
3. Git提交确认信息