# GSC (Genotype Sparse Compression)
Genotype Sparse Compression (GSC) is an advanced tool for lossless compression of VCF files, designed to efficiently store and manage VCF files in a compressed format. It accepts VCF/BCF files as input and utilizes advanced compression techniques to significantly reduce storage requirements while ensuring fast query capabilities. In our study, we successfully compressed the VCF files from the 1000 Genomes Project (1000Gpip3), consisting of 2504 samples and 80 million variants, from an uncompressed VCF file of 803.70GB to approximately 1GB.

## Requirements
### GSC requires:

- **Compiler Compatibility**: GSC requires a modern C++14-ready compiler, such as:
  - g++ version 10.1.0 or higher

- **Build System**: CMake 3.16+ is required for building GSC.

- **Operating System**: GSC supports 64-bit operating systems, including:
  - Linux (Ubuntu 18.04+)

## Build
### Dockerfile
Dockerfile can be used to build a Docker image with all necessary dependencies and GSC compressor. The image is based on Ubuntu 18.04. To build a Docker image and run a Docker container, you need Docker Desktop (https://www.docker.com). Example commands (run it within a directory with Dockerfile):
```bash
# build
docker build -t gsc_project .

# run
docker run -it gsc_project
```
### Building the GSC command line tool

```bash
# Clone the repository
git clone https://github.com/luo-xiaolong/GSC.git
cd GSC

# Build using CMake
cd build
cmake ..
make -j4
```
To clean the GSC build use:
```bash
cd build && rm -rf *
```
## Usage
```bash
Usage: gsc [option] [arguments]
Available options:
        compress        - compress VCF/BCF file (multi-sample)
        decompress      - query and decompress to VCF/BCF file
        gvcf            - compress single-sample gVCF file (optimized)
        gvcf-decompress - decompress gVCF compressed file
```
### Compress (multi-sample VCF/BCF)
```bash
Usage of gsc compress:

    gsc compress [options] [--in [in_file]] [--out [out_file]]

Where:
    [options]              Optional flags and parameters for compression.
    -i,  --in [in_file]    Specify the input file (default: VCF or VCF.GZ). If omitted, input is taken from standard input (stdin).
    -o,  --out [out_file]  Specify the output file. If omitted, output is sent to standard output (stdout).
    --compressor [name]    Select compressor: bsc (default), zstd, brotli.

Options:
    -M,  --mode_lossly     Choose lossy compression mode (lossless by default).
    -b,  --bcf             Input is a BCF file (default: VCF or VCF.GZ).
    -p,  --ploidy [X]      Set ploidy of samples in input VCF to [X] (default: 2).
        --max-block-rows [X]  Max variants per GT block (default: 10000).
        --max-block-cols [X]  Max haplotypes (samples * ploidy) per GT column block (default: 10000).
    -t,  --threads [X]     Set number of threads to [X] (default: 1).
    -d,  --depth [X]       Set maximum replication depth to [X] (default: 100, 0 means no matches).
    -m,  --merge [X]       Specify files to merge, separated by commas (e.g., -m chr1.vcf,chr2.vcf), or '@' followed by a file containing a list of VCF files (e.g., -m @file_with_IDs.txt). By default, all VCF files are compressed.
```
### Decompress / Query
```bash
Usage of gsc decompress and query:

    gsc decompress [options] --in [in_file] --out [out_file]

Where:
    [options]              Optional flags and parameters for compression.
    -i,  --in [in_file]    Specify the input file. If omitted, input is taken from standard input (stdin).
    -o,  --out [out_file]  Specify the output file (default: VCF). If omitted, output is sent to standard output (stdout).
    --compressor [name]    Select compressor: bsc (default), zstd, brotli.

Options:
    General Options:
        -M,  --mode_lossly    Choose lossy compression mode (default: lossless).
        -b,  --bcf            Output a BCF file (default: VCF).

    Filter options (applicable in lossy compression mode only):
        -r,  --range [X]      Specify range in format [start],[end] (e.g., -r 4999756,4999852).
        -s,  --samples [X]    Samples separated by comms (e.g., -s HG03861,NA18639) OR '@' sign followed by the name of a file with sample name(s) separated by whitespaces (for exaple: -s @file_with_IDs.txt). By default all samples/individuals are decompressed.
        --header-only         Output only the header of the VCF/BCF.
        --no-header           Output without the VCF/BCF header (only genotypes).
        -G,  --no-genotype    Don't output sample genotypes (only #CHROM, POS, ID, REF, ALT, QUAL, FILTER, and INFO columns).
        -C,  --out-ac-an      Write AC/AN to the INFO field.
        -S,  --split          Split output into multiple files (one per chromosome).
        -I, [ID=^]            Include only sites with specified ID (e.g., -I "ID=rs6040355").
        --minAC [X]           Include only sites with AC <= X.
        --maxAC [X]           Include only sites with AC >= X.
        --minAF [X]           Include only sites with AF >= X (X: 0 to 1).
        --maxAF [X]           Include only sites with AF <= X (X: 0 to 1).
        --min-qual [X]        Include only sites with QUAL >= X.
        --max-qual [X]        Include only sites with QUAL <= X.
```

### gVCF Compression (single-sample optimized)
For single-sample gVCF files (e.g., from GATK HaplotypeCaller), use the optimized gVCF compression mode which achieves ~3.7% compression ratio.

```bash
Usage of gsc gvcf:

    gsc gvcf -i [in_file] -o [out_file] [--compressor name]

Where:
    -i,  --in [in_file]    Input gVCF file (required, supports .gvcf, .gvcf.gz, .vcf, .vcf.gz)
    -o,  --out [out_file]  Output compressed file (required, .gsc extension)
    --compressor [name]    Select compressor: bsc (default, best ratio), zstd (fast), brotli

Example:
    ./gsc gvcf -i sample.gvcf.gz -o sample.gsc --compressor bsc
```

### gVCF Decompression
```bash
Usage of gsc gvcf-decompress:

    gsc gvcf-decompress -i [in_file] -o [out_file]

Where:
    -i,  --in [in_file]    Input compressed file (.gsc)
    -o,  --out [out_file]  Output file (.vcf, .vcf.gz, or .bcf - auto-detected by extension)

Example:
    ./gsc gvcf-decompress -i sample.gsc -o restored.vcf.gz
```
## Example
There is an example VCF/VCF.gz/BCF file, `toy.vcf`/`toy.vcf.gz`/`toy.bcf`, in the toy folder, which can be used to test GSC
### Compress

#### Lossless compression:
The input file format is VCF. You can compress a VCF file in lossless mode using one of the following methods:
1. **Explicit input and output file parameters**:
   
   Use the `--in` option to specify the input VCF file and the `--out` option for the output compressed file.
   ```bash
   ./gsc compress --in toy/toy.vcf --out toy/toy_lossless.gsc
   ```
2. **Input file parameter and output redirection**:
   
   Use the `--out` option for the output compressed file and redirect the input VCF file into the command.
   ```bash
   ./gsc compress --out toy/toy_lossless.gsc < toy/toy.vcf
   ```
3. **Output file redirection and input file parameter**:
   
   Specify the input VCF file with the `--in` option and redirect the output to create the compressed file.
   ```bash
   ./gsc compress --in toy/toy.vcf > toy/toy_lossless.gsc
   ```
4. **Input and output redirection**:
   
   Use shell redirection for both input and output. This method does not use the `--in` and `--out` options.
   ```bash
   ./gsc compress < toy/toy.vcf > toy/toy_lossless.gsc
   ```
This will create a file:
* `toy_lossless.gsc` - The compressed archive of the entire VCF file.

#### Lossy compression:

The input file format is VCF. The commands are similar to those used for lossless compression, with the addition of the `-M` parameter to enable lossy compression.

   For example, to compress a VCF file in lossy mode:

   ```bash
   ./gsc compress -M --in toy/toy.vcf --out toy/toy_lossy.gsc
   ```
   or 
  
   Using redirection:
   ```bash
   ./gsc compress -M --out toy/toy_lossy.gsc < toy/toy.vcf
   ``` 
   This will create a file:
   * `toy_lossy.gsc` - The compressed archive of the entire VCF file is implemented with lossy compression. It only retains the 'GT' subfield within the INFO and FORMAT fields, and excludes all other subfields..
    
### Decompress   (The commands are similar to those used for compression)
#### Lossless decompression:

To decompress the compressed toy_lossless.gsc into a VCF file named toy_lossless.vcf:
```bash
./gsc decompress --in toy/toy_lossless.gsc --out toy/toy_lossless.vcf
```
#### Lossy decompression:

To decompress the compressed toy_lossy.gsc into a VCF file named toy_lossy.vcf:
```bash
./gsc decompress -M --in toy/toy_lossy.gsc --out toy/toy_lossy.vcf
```
### Query
#### Variant-based query
Retrieve entries for chromosome 20 with POS ranging from 1 to 1,000,000, and output to the toy/query_toy_r_20_3_1000000.vcf file.

```bash
./gsc decompress -M --range 20:1,1000000 --in toy/toy_lossy.gsc --out toy/query_toy_20_3_1000000.vcf
```
Retrieve entries for chromosome 20 with POS ranging from 1 to 1,000,000, and output to the terminal interface.
```bash
./gsc decompress -M --range 20:1,1000000 --in toy/toy_lossy.gsc
```
#### Sample-based query
Retrieve genotype columns for samples named NA00001 and NA00002, and output to the toy/query_toy_s_NA00001_NA00002.vcf file.
```bash
./gsc decompress -M --samples NA00001,NA00002 --in toy/toy_lossy.gsc --out toy/query_toy_s_NA00001_NA00002.vcf
```
or

The names NA00001 and NA00002 are stored in the toy/samples_name_file.
```bash
./gsc decompress -M --samples @toy/samples_name_file --in toy/toy_lossy.gsc --out toy/query_toy_s_NA00001_NA00002.vcf
```
#### Note
You can also perform mixed queries based on sample names and variants.

### gVCF Compression Example
For single-sample gVCF files (e.g., from GATK HaplotypeCaller):

#### Compress a gVCF file:
```bash
./gsc gvcf -i toy/HG002_first10k_fixed.gvcf -o toy/HG002.gsc --compressor bsc
```

#### Decompress to VCF:
```bash
./gsc gvcf-decompress -i toy/HG002.gsc -o toy/HG002_restored.vcf
```

#### Decompress to compressed VCF (gzip):
```bash
./gsc gvcf-decompress -i toy/HG002.gsc -o toy/HG002_restored.vcf.gz
```

## Compression Backend Comparison
GSC supports three compression backends:

| Backend | Compression Ratio | Speed | Recommended Use |
|---------|------------------|-------|-----------------|
| **bsc** | Best (~3.7%) | Medium | Default, best for archival |
| **zstd** | Good (~5.7%) | Fast | When speed matters |
| **brotli** | Good (~4.8%) | Slow | Not recommended |

## Citations
- **bio.tools ID**: `gsc_genotype_sparse_compression`
- **Research Resource Identifier (RRID)**: `SCR_025071`
- **Doi**:`https://doi.org/10.48546/WORKFLOWHUB.WORKFLOW.887.1`
