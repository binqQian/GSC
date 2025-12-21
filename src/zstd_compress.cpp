#include "zstd_compress.h"

namespace zstd {


 bool zstd_compress(const std::vector<uint8_t>& srcContent, std::vector<uint8_t>& cBuff) {
    size_t cSizeActual = 0;

    size_t cSize = ZSTD_compressBound(srcContent.size());

 
    cBuff.resize(cSize);


    cSizeActual = ZSTD_compress(cBuff.data(), cSize, srcContent.data(), srcContent.size(), 10);
    cBuff.resize(cSizeActual);

    return !ZSTD_isError(cSizeActual);
  }

  bool zstd_decompress(const std::vector<uint8_t>& cBuff, std::vector<uint8_t>& dBuff) {
    return zstd_decompress_ptr(cBuff.data(), cBuff.size(), dBuff);
  }

  bool zstd_decompress_ptr(const uint8_t* src, size_t src_size, std::vector<uint8_t>& dBuff) {
    const size_t hinted_size = ZSTD_getFrameContentSize(src, src_size);
    if (hinted_size != ZSTD_CONTENTSIZE_ERROR && hinted_size != ZSTD_CONTENTSIZE_UNKNOWN) {
      dBuff.resize(hinted_size);
      const size_t actual = ZSTD_decompress(dBuff.data(), hinted_size, src, src_size);
      if (ZSTD_isError(actual)) {
        return false;
      }
      dBuff.resize(actual);
      return true;
    }

    ZSTD_DCtx* dctx = ZSTD_createDCtx();
    if (!dctx) return false;
    size_t offset = 0;
    const size_t chunk = 1 << 16;
    dBuff.clear();
    dBuff.resize(chunk);
    ZSTD_inBuffer input = {src, src_size, 0};
    bool ok = true;
    while (input.pos < input.size) {
      ZSTD_outBuffer output = {dBuff.data() + offset, dBuff.size() - offset, 0};
      const size_t ret = ZSTD_decompressStream(dctx, &output, &input);
      if (ZSTD_isError(ret)) {
        ok = false;
        break;
      }
      offset += output.pos;
      if (input.pos < input.size) {
        dBuff.resize(dBuff.size() + chunk);
      }
    }
    dBuff.resize(offset);
    ZSTD_freeDCtx(dctx);
    return ok;
  }
}
