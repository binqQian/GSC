# Agent Notes (GSC)
#
# Scope: this file applies to the whole repo.
#
# Build
- Configure: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`
- Build: `cmake --build build -j 8`
#
# Useful runtime knobs
- Enable field stats table: `GSC_LOG_LEVEL=debug`
- PL codec selector: `GSC_PL_CODEC_ID` (`4` default, `5` experimental, `6` recommended)
- AD codec selector: `GSC_AD_CODEC_ID` (`4` default, `5` recommended)
- PL block size: `GSC_PL_BLOCK_SIZE` (default `16384`)
- AD baseline block size: `GSC_AD_BASELINE_BLOCK` (default `16384`)
#
# Quick lossless check
- Compress: `GSC_PL_CODEC_ID=6 GSC_AD_CODEC_ID=5 ./build/gsc compress --in toy/final_subset.vcf.gz --out tmp/test.gsc --compressor bsc`
- Decompress: `./build/gsc decompress --in tmp/test.gsc --out tmp/test.vcf`
- Verify: `zcat toy/final_subset.vcf.gz | md5sum && md5sum tmp/test.vcf`
#
# Repo quirk
- `.gitignore` ignores `*.md` by default; `AGENTS.md` and `docs/optimization/AD_PL_Optimization.md` are explicitly un-ignored.
