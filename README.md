# GSC (Genotype Sparse Compression)
Genotype Sparse Compression (GSC) is an advanced tool for compression of VCF files, designed to efficiently store and manage VCF files in a compressed format. It accepts VCF/BCF files as input and utilizes advanced compression techniques to significantly reduce storage requirements while ensuring fast query capabilities. In our study, we successfully compressed the VCF files from the 1000 Genomes Project (1000G phase 3), consisting of 2504 samples and 80 million variants, from an uncompressed VCF file of 803.70GB to approximately 1GB.

## Key Features
- Lossless compression of VCF/BCF into `.gsc` with fast decompression and range/sample queries.
- Pluggable compression backend for most streams: `bsc` (default), `zstd`, `brotli` via `--compressor`.
- Optional **adaptive FORMAT** stream (`adaptive_format_data`) for seen-in-practice FORMAT tags (currently most effective on AD/PL/DP-like data).
- Optional **integrity verification** using XXHash64 to ensure data consistency after compression/decompression.

## Requirements
GSC requires:

- **Compiler Compatibility**: GSC requires a modern C++14-ready compiler, such as:
  - g++ version 10.1.0 or higher

- **Build System**: `make` (Makefile build) or `cmake` >= 3.16 (CMake build).

- **Operating System**: GSC supports 64-bit operating systems, including:
  - Linux (Ubuntu 18.04)

- **Libraries**:
  - System: `zlib`, `libbz2`, `liblzma`
  - Logging/formatting: `spdlog`, `fmt`
  - Compression backend support: `brotli` (required to build; used when selecting `--compressor brotli`)

On Ubuntu you can usually install these with `apt` (package names vary by distro), or use `vcpkg`.

## Build
### Option A: CMake (recommended)
Clone:
```bash
git clone https://github.com/luo-xiaolong/GSC.git
cd GSC
```
If you use `vcpkg`, set `VCPKG_ROOT` so CMake can locate `spdlog`, `fmt`, and `brotli`:
```bash
export VCPKG_ROOT=/path/to/vcpkg
```
Build:
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/gsc
```

### Option B: Make
The Makefile can use `VCPKG_ROOT` to locate `spdlog`, `fmt`, and `brotli`:
```bash
export VCPKG_ROOT=/path/to/vcpkg/installed/x64-linux
make clean
make
./gsc
```

### Dockerfile
Dockerfile can be used to build a Docker image. The image is based on Ubuntu 18.04. To build a Docker image and run a Docker container, you need Docker Desktop (https://www.docker.com). Example commands (run it within a directory with Dockerfile):
```bash
#build
docker build -t gsc_project .

#run
docker run -it gsc_project
```

## Usage
Quick reference (run `./gsc compress` / `./gsc decompress` for the full list):

- `--compressor bsc|zstd|brotli`: select backend for the majority of streams (choice is stored in the `.gsc` file).
- Adaptive FORMAT (lossless mode only):
  - `compress --adaptive-format off|shadow|primary`
  - `compress --adaptive-format-compressor auto|follow|raw|bsc|zstd|brotli` (controls `adaptive_format_data` part compression; default `auto`)
  - `decompress --use-adaptive-format` (required when decompressing a `.gsc` created with `--adaptive-format primary`)

```bash
Usage: gsc [option] [arguments] 
Available options: 
        compress - compress VCF/BCF file
        decompress     - query and decompress to VCF/BCF file
```
- Compress the input VCF/BCF file
```bash
Usage of gsc compress:

        gsc compress [options] [--in [in_file]] [--out [out_file]]

Where:

        [options]              Optional flags and parameters for compression.
        -i,  --in [in_file]    Specify the input file (default: VCF or VCF.GZ). If omitted, input is taken from standard input (stdin).
        -o,  --out [out_file]  Specify the output file. If omitted, output is sent to standard output (stdout).
        --compressor [name]    Select compressor: bsc (default), zstd, brotli.
        --adaptive-format [off|shadow|primary]  Adaptive FORMAT mode (lossless only).
        --adaptive-format-compressor [name]  Compressor for `adaptive_format_data` parts: auto (default), follow, raw, bsc, zstd, brotli.

Options:

        -M,  --mode_lossly     Choose lossy compression mode (lossless by default).
        -b,  --bcf             Input is a BCF file (default: VCF or VCF.GZ).
        -p,  --ploidy [X]      Set ploidy of samples in input VCF to [X] (default: 2).
        --max-block-rows [X]   Max variants per GT block (default: 10000).
        --max-block-cols [X]   Max haplotypes (samples * ploidy) per GT column block (default: 10000).
        -t,  --threads [X]     Set number of threads to [X] (default: 5).
        -d,  --depth [X]       Set maximum replication depth to [X] (default: 100, 0 means no matches).
        -m,  --merge [X]       Specify files to merge, separated by commas (e.g., -m chr1.vcf,chr2.vcf), or '@' followed by a file containing a list of VCF files (e.g., -m @file_with_IDs.txt). By default, all VCF files are compressed.
        -I,  --integrity [mode]  Enable integrity hash: none (default), xxhash64/xxhash/fast.
```
- Decompress / Query
```bash
Usage of gsc decompress and query:

        gsc decompress [options] --in [in_file] --out [out_file]

Where:
        [options]              Optional flags and parameters for compression.
        -i,  --in [in_file]    Specify the input file . If omitted, input is taken from standard input (stdin).
        -o,  --out [out_file]  Specify the output file (default: VCF). If omitted, output is sent to standard output (stdout).
        --compressor [name]    Select compressor: bsc (default), zstd, brotli.
        --use-adaptive-format  Prefer `adaptive_format_data` stream for non-GT FORMAT reconstruction (recommended with VCF output).

Options:

    General Options:

        -M,  --mode_lossly      Choose lossy compression mode (default: lossless).
        -b,  --bcf              Output a BCF file (default: VCF).
        -t,  --threads [X]      Set number of threads to [X] (default: 5).
        -l,  --level [0-9]      When writing to `.vcf.gz`, set bgzip level (0 = uncompressed).
        --make-bed              Output PLINK BED/BIM/FAM instead of VCF/BCF.

    Integrity Verification:

        --verify                Force integrity verification if hash is present (default: enabled).
        --no-verify             Skip integrity verification even if hash is stored.

    Filter / Query options:

        -r,  --range [X]        Specify range in format `CHR[:START[,END]]` (e.g., `-r chr1:4999756,4999852` or `-r 20:1,1000000`).
        -s,  --samples [X]      Samples separated by commas (e.g., `-s HG03861,NA18639`) OR `@file` with one sample per whitespace. By default all samples/individuals are decompressed.
        --header-only           Output only the header of the VCF/BCF.
        --no-header             Output without the VCF/BCF header (only genotypes).
        -G,  --no-genotype      Don't output sample genotypes (only #CHROM, POS, ID, REF, ALT, QUAL, FILTER, and INFO columns).
        -C,  --out-ac-an        Write AC/AN to the INFO field.
        -S,  --split            Split output into multiple files (one per chromosome).
        -I, [ID=^]              Include only sites with specified ID (e.g., -I "ID=rs6040355").
        --minAC [X]             Include only sites with AC <= X.
        --maxAC [X]             Include only sites with AC >= X.
        --minAF [X]             Include only sites with AF >= X (X: 0 to 1).
        --maxAF [X]             Include only sites with AF <= X (X: 0 to 1).
        --min-qual [X]          Include only sites with QUAL >= X.
        --max-qual [X]          Include only sites with QUAL <= X.
```

## Adaptive FORMAT Compression (Experimental)

GSC has an optional **adaptive FORMAT** stream (`adaptive_format_data`) for lossless compression of non-GT FORMAT subfields. This is designed to:
- provide field-aware encoding for common tags (currently most effective on AD/PL/DP-like patterns),
- allow correct sample-subset output for non-GT FORMAT when using `-s/--samples` (legacy non-GT FORMAT cannot be safely subsetted),
- remain backward-compatible (old `.gsc` without the adaptive stream still decompress as before).

### Modes
Enable on compression with:
```bash
./gsc compress --adaptive-format off|shadow|primary ...
```
- `off` (default): legacy behavior; no adaptive stream is written.
- `shadow`: writes both legacy non-GT FORMAT (other_fields) and adaptive stream (larger file; useful for regression).
- `primary`: writes adaptive stream and **suppresses legacy only for tags actually emitted to adaptive**; tags not covered by adaptive remain in legacy streams.

Important: If you compress with `--adaptive-format primary`, you should decompress with `--use-adaptive-format`, otherwise suppressed tags will not be reconstructed.

### Part compression (adaptive stream only)
Adaptive stream parts use a small framing header (`AFD1`) and can be independently compressed:
```bash
./gsc compress --adaptive-format primary \
  --adaptive-format-compressor auto|follow|raw|bsc|zstd|brotli \
  --in input.vcf --out output.gsc
```
- `auto` (default): tries `zstd/bsc/brotli` and picks the smallest.
- `follow`: uses the same backend as `--compressor` (legacy streams).
- `raw`: stores uncompressed payloads (debug only; can be much larger).

Note: this option affects **only** `adaptive_format_data`; other streams still follow `--compressor`.

### Current limitations
- Adaptive decoding is currently intended for **VCF output** (`--out *.vcf` / `*.vcf.gz`). For BCF output, legacy typed FORMAT is preferred.
- PL “AB pattern” is only enabled for biallelic diploid PL vectors (len==3) for correctness; multi-allelic reconstruction needs further work.

### Practical guidance / FAQ

#### “Why is `--adaptive-format primary` sometimes larger than `off`?”
This can happen when legacy other_fields streams compress extremely well (chunk-level compression) while adaptive data is stored row-wise with additional metadata.

Recommendations:
- Start with `--adaptive-format-compressor auto` (default) or `--adaptive-format-compressor follow` (match `--compressor`), then compare sizes:
  ```bash
  ./gsc compress --adaptive-format off     -i input.vcf -o tmp/off.gsc
  ./gsc compress --adaptive-format primary -i input.vcf -o tmp/primary.gsc
  ls -lh tmp/off.gsc tmp/primary.gsc
  ```
- Use `--adaptive-format shadow` to validate correctness while keeping the legacy stream intact:
  ```bash
  ./gsc compress --adaptive-format shadow -i input.vcf -o tmp/shadow.gsc
  ./gsc decompress --in tmp/shadow.gsc --out tmp/shadow.vcf
  ```
- If you see many non-GT FORMAT tags that are effectively “free text” (codec falls back to `RawString`), adaptive will intentionally omit them in `primary` to avoid size regressions; in that case, overall gains may be limited to AD/PL/DP-like tags.

#### “Do I need `--use-adaptive-format` on decompression?”
- If the `.gsc` was created with `--adaptive-format off` or `shadow`: decompression without `--use-adaptive-format` is fine (legacy path includes all tags).
- If the `.gsc` was created with `--adaptive-format primary`: you should use `--use-adaptive-format`, otherwise suppressed tags cannot be reconstructed.

#### “Why are non-GT FORMAT fields missing when using `-s/--samples`?”
Legacy non-GT FORMAT streams cannot be safely subsetted without re-encoding, so in sample-subset mode GSC only outputs non-GT FORMAT fields that can be reconstructed safely. Enabling adaptive FORMAT (`--adaptive-format primary` on compression + `--use-adaptive-format` on decompression) is the intended solution.

## Integrity Verification (Hash)

GSC supports optional **integrity verification** using XXHash64 to ensure data consistency after compression/decompression.

### How It Works
- **Compression**: When enabled with `--integrity xxhash64`, GSC computes a hash of the VCF header and sample names, then writes this hash to a footer at the end of the `.gsc` file.
- **Decompression**: GSC automatically detects the hash footer and verifies the decompressed data matches the original. If verification fails, a warning is logged.

### Usage
```bash
# Compress with integrity hash
./gsc compress --integrity xxhash64 -i input.vcf -o output.gsc

# Decompress with automatic verification (default behavior)
./gsc decompress -i output.gsc -o restored.vcf

# Skip verification if needed
./gsc decompress --no-verify -i output.gsc -o restored.vcf
```

### Default Behavior
- **Compression**: Integrity checking is **disabled by default**. Use `--integrity xxhash64` to enable.
- **Decompression**: If a hash is present in the file, verification is **automatic**. Use `--no-verify` to skip.

### Technical Details
- Uses **XXHash64** algorithm (speed: ~13.8 GB/s, much faster than MD5/SHA)
- Hash is stored as a 14-byte footer at the end of the `.gsc` file
- Backward compatible: old GSC versions ignore the footer; new versions can read old files without hash

## Logging
GSC logs to stderr and to a rotating file `gsc.log` (up to ~10MB, keep 3 rotated files). Use `GSC_LOG_LEVEL` to control verbosity:
```bash
export GSC_LOG_LEVEL=debug   # trace|debug|info|warn|error|critical|off
```

## Validation / Smoke Tests
The `toy/` directory contains small fixtures for manual verification.

Lossless round-trip (legacy):
```bash
./gsc compress --in toy/toy.vcf --out tmp/toy_off.gsc
./gsc decompress --in tmp/toy_off.gsc --out tmp/toy_off.vcf
```

Lossless round-trip (adaptive primary, recommended flags):
```bash
./gsc compress --adaptive-format primary --in toy/toy.vcf --out tmp/toy_primary.gsc
./gsc decompress --use-adaptive-format --in tmp/toy_primary.gsc --out tmp/toy_primary.vcf
diff -q tmp/toy_off.vcf tmp/toy_primary.vcf
```

Compare file sizes:
```bash
ls -lh tmp/toy_off.gsc tmp/toy_primary.gsc
```

## Developer Notes
- Adaptive FORMAT implementation/plan:
  - `ADAPTIVE_FORMAT_HANDOFF.md`
  - `ADAPTIVE_FORMAT_IMPLEMENTATION_NOTES.md`
  - `ADAPTIVE_FORMAT_COMPRESSION_PLAN.md`

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
   Select a different compressor backend (the choice is stored in the `.gsc` file and auto-detected during decompression):
   ```bash
   ./gsc compress --compressor zstd --in toy/toy.vcf --out toy/toy_lossless_zstd.gsc
   ```
   Adaptive FORMAT (shadow for regression; primary for production):
   ```bash
   ./gsc compress --adaptive-format shadow  --in toy/toy.vcf --out tmp/toy_shadow.gsc
   ./gsc compress --adaptive-format primary --in toy/toy.vcf --out tmp/toy_primary.gsc
   ./gsc decompress --use-adaptive-format --in tmp/toy_primary.gsc --out tmp/toy_primary.vcf
   ```
   Control adaptive stream part compressor independently:
   ```bash
   ./gsc compress --adaptive-format primary --adaptive-format-compressor follow --compressor brotli \
     --in toy/toy.vcf --out tmp/toy_primary_follow_brotli.gsc
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
   * `toy_lossy.gsc` - The compressed archive of the entire VCF file is implemented with lossy compression. It keeps only the genotype (`GT`) and drops other INFO/FORMAT subfields.
    
### Decompress   (The commands are similar to those used for compression)
#### Lossless decompression:

To decompress the compressed toy_lossless.gsc into a VCF file named toy_lossless.vcf:
```bash
./gsc decompress --in toy/toy_lossless.gsc --out toy/toy_lossless.vcf
```
Write a gzipped VCF (bgzip) with a chosen compression level:
```bash
./gsc decompress --in toy/toy_lossless.gsc --out toy/toy_lossless.vcf.gz --level 6
```
Decompress a `.gsc` produced with `--adaptive-format primary`:
```bash
./gsc decompress --use-adaptive-format --in tmp/toy_primary.gsc --out tmp/toy_primary.vcf
```

Note: `decompress --compressor ...` is generally not needed. The backend used to compress the `.gsc` is stored in the file and auto-detected during decompression.
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

## Citations
- **bio.tools ID**: `gsc_genotype_sparse_compression`
- **Research Resource Identifier (RRID)**: `SCR_025071`
- **Doi**:`https://doi.org/10.48546/WORKFLOWHUB.WORKFLOW.887.1`
