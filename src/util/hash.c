#include <inttypes.h>
#include <stddef.h>

union uint32 {
  int32_t i;
  uint32_t u;
};
union uint64 {
  int64_t i;
  uint64_t u;
};

// 32-bit integer hash function (nullprogram.com/blog/2018/07/31/)
uint32_t c_hash_uint32(uint32_t i) {
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
int32_t c_hash_int32(int32_t i) {
  union uint32 ui, uh;

  ui.i = i;
  uh.u = c_hash_uint32(ui.u);

  return uh.i;
}

// hash 32-bit integer array
uint32_t c_hash_uint32_array(const uint32_t *a, size_t n) {
  uint32_t h;
  size_t i;

  h = 0;
  for (i = 0; i < n; ++i) h = c_hash_int32(h ^ a[i]);

  return h;
}
int32_t c_hash_int32_array(const int32_t *a, size_t n) {
  union uint32 *ua, uh;

  ua   = (union uint32 *)a;
  uh.u = c_hash_uint32_array((uint32_t *)ua, n);

  return uh.i;
}

// 64-bit integer hash function (splitmix64)
uint64_t c_hash_uint64(uint64_t i) {
  uint64_t h;

  h = i + 0x9e3779b97f4a7c15; // avoid i == 0 -> h == 0
  h = h ^ (h >> 30);
  h = h * 0xbf58476d1ce4e5b9;
  h = h ^ (h >> 27);
  h = h * 0x94d049bb133111eb;
  h = h ^ (h >> 31);

  return h;
}
int64_t c_hash_int64(int64_t i) {
  union uint64 ui, uh;

  ui.i = i;
  uh.u = c_hash_uint64(ui.u);

  return uh.i;
}

// hash 64-bit integer array
uint64_t c_hash_uint64_array(const uint64_t *a, size_t n) {
  uint64_t h;
  size_t i;

  h = 0;
  for (i = 0; i < n; ++i) h = c_hash_int64(h ^ a[i]);

  return h;
}
int64_t c_hash_int64_array(const int64_t *a, size_t n) {
  union uint64 *ua, uh;

  ua = (union uint64 *)a;
  uh.u = c_hash_uint64_array((uint64_t *)ua, n);

  return uh.i;
}

// string hash function (FNV-1a for 32 and 64 bits)
int32_t c_hash_string32(const unsigned char *str) {
  union uint32 hash;
  uint32_t c;

  hash.u = 2166136261u;
  while ((c = *str++)) {
    hash.u ^= c;
    hash.u *= 16777619u;
  }

  return hash.i;
}
int64_t c_hash_string64(const unsigned char *str) {
  union uint64 hash;
  uint64_t c;

  hash.u = 14695981039346656037ull;
  while ((c = *str++)) {
    hash.u ^= c;
    hash.u *= 1099511628211ull;
  }

  return hash.i;
}
