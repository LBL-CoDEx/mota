#ifndef _442187d7_0979_4aad_85d7_194ab1174c2e
#define _442187d7_0979_4aad_85d7_194ab1174c2e

#include <cstdint>

namespace mota {

uint32_t EncodeMorton2(uint32_t x, uint32_t y);
uint32_t EncodeMorton3(uint32_t x, uint32_t y, uint32_t z);
uint32_t EncodeMorton4(uint32_t x, uint32_t y, uint32_t z, uint32_t w);
uint32_t EncodeMorton5(uint32_t x, uint32_t y, uint32_t z, uint32_t w, uint32_t v);

}

#endif
