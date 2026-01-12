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
- AD/PL codecs: built-in (PL=`6`, AD=`5`)
#
# Quick lossless check
- Compress: `./build/gsc compress --in toy/final_subset.vcf.gz --out tmp/test.gsc --compressor bsc`
- Decompress: `./build/gsc decompress --in tmp/test.gsc --out tmp/test.vcf`
- Verify: `zcat toy/final_subset.vcf.gz | md5sum && md5sum tmp/test.vcf`
#
# Repo quirk
- `.gitignore` ignores `*.md` by default; `AGENTS.md` and `docs/optimization/AD_PL_Optimization.md` are explicitly un-ignored.
