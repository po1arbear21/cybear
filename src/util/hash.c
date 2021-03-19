#include <inttypes.h>
#include <stddef.h>

// 32-bit integer hash function (nullprogram.com/blog/2018/07/31/)
uint32_t c_hash_int32(uint32_t i) {
  uint32_t h;

  h = i + 0x7f4a7c15; // avoid i == 0 -> h == 0
  h = h ^ (h >> 17);
  h = h * 0xed5ad4bb;
  h = h ^ (h >> 11);
  h = h * 0xac4c1b51;
  h = h ^ (h >> 15);
  h = h * 0x31848bab;
  h = h ^ (h >> 14);

  return h;
}

// hash 32-bit integer array
uint32_t c_hash_int32_array(const uint32_t *a, size_t n) {
  uint32_t h;
  size_t i;

  h = 0;
  for (i = 0; i < n; ++i) h = c_hash_int32(h ^ a[i]);

  return h;
}

// 64-bit integer hash function (splitmix64)
uint64_t c_hash_int64(uint64_t i) {
  uint64_t h;

  h = i + 0x9e3779b97f4a7c15; // avoid i == 0 -> h == 0
  h = h ^ (h >> 30);
  h = h * 0xbf58476d1ce4e5b9;
  h = h ^ (h >> 27);
  h = h * 0x94d049bb133111eb;
  h = h ^ (h >> 31);

  return h;
}

// hash 64-bit integer array
uint64_t c_hash_int64_array(const uint64_t *a, size_t n) {
  uint64_t h;
  size_t i;

  h = 0;
  for (i = 0; i < n; ++i) h = c_hash_int64(h ^ a[i]);

  return h;
}
