// build_k21_prefix13_suffix8_sharded_vectors.cpp
//
// Reads UNSORTED `jellyfish dump -c` from stdin for k=21.
// Builds:
//   prefix_starts.bin : uint32_t starts[2^26+1]
//   suffixes.bin      : uint16_t suffix[N]   (sorted within prefix bucket)
//   counts.bin        : uint16_t counts[N]
//
// Approach: "a bunch of STL vectors" via sharding.
// - Shard by top SHARD_BITS of 26-bit prefix into 2^SHARD_BITS vectors.
// - One pass: push (prefix,suffix,count) into shard vector.
// - Then sort each shard by (prefix,suffix).
// - Emit shards in increasing shard order => globally sorted by prefix then suffix.
// - Build starts[] as we emit.

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

static inline bool is_space(unsigned char c) {
  return c==' ' || c=='\t' || c=='\n' || c=='\r';
}

static inline uint8_t base2bits(unsigned char c) {
  switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default:            return 255;
  }
}

static inline uint16_t clamp_u16(uint64_t x) {
  return (x > std::numeric_limits<uint16_t>::max())
           ? std::numeric_limits<uint16_t>::max()
           : static_cast<uint16_t>(x);
}

static void write_all(int fd, const void* buf, size_t nbytes) {
  const uint8_t* p = static_cast<const uint8_t*>(buf);
  while (nbytes) {
    ssize_t w = ::write(fd, p, nbytes);
    if (w < 0) throw std::runtime_error("write() failed");
    p += static_cast<size_t>(w);
    nbytes -= static_cast<size_t>(w);
  }
}

#pragma pack(push, 1)
struct Rec {
  uint32_t prefix; // lower 26 bits used
  uint16_t suffix; // last 8 bases
  uint16_t count;  // capped at 64k
};
#pragma pack(pop)
static_assert(sizeof(Rec) == 8, "Rec should be 8 bytes");

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "Usage: jellyfish dump -c <db.jf> | " << argv[0]
              << " <prefix_starts.bin> <suffixes.bin> <counts.bin> "
              << "[--shard_bits=12] [--reserve=N] [--max_n=N]\n";
    return 2;
  }

  const std::string out_prefix = argv[1];
  const std::string out_suf    = argv[2];
  const std::string out_cnt    = argv[3];

  int shard_bits = 12;     // 4k shards by default
  uint64_t reserve_n = 0;  // total reserve across all shards (rough)
  uint64_t max_n = 0;

  for (int i = 4; i < argc; ++i) {
    std::string a = argv[i];
    if (a.rfind("--shard_bits=", 0) == 0) shard_bits = std::atoi(a.c_str() + 12);
    else if (a.rfind("--reserve=", 0) == 0) reserve_n = std::strtoull(a.c_str() + 10, nullptr, 10);
    else if (a.rfind("--max_n=", 0) == 0) max_n = std::strtoull(a.c_str() + 8, nullptr, 10);
    else { std::cerr << "ERROR: unknown option: " << a << "\n"; return 2; }
  }
  if (shard_bits < 1 || shard_bits > 20) {
    std::cerr << "ERROR: shard_bits should be in [1,20]\n";
    return 2;
  }

  constexpr uint32_t P = (1u << 26);
  const uint32_t NSHARDS = (1u << shard_bits);
  const uint32_t SHIFT = 26 - shard_bits;

  std::vector<std::vector<Rec>> shards(NSHARDS);

  if (reserve_n) {
    // rough even reserve: reserve_n / NSHARDS per shard
    size_t per = static_cast<size_t>(reserve_n / NSHARDS);
    for (auto& v : shards) v.reserve(per);
  }

  // Fast stdin buffer
  static constexpr size_t INBUF_SZ = 32 * 1024 * 1024;
  std::vector<uint8_t> buf(INBUF_SZ);
  size_t fill = 0, pos = 0;

  auto refill = [&]() -> bool {
    if (pos < fill) {
      size_t rem = fill - pos;
      std::memmove(buf.data(), buf.data() + pos, rem);
      fill = rem;
      pos = 0;
    } else {
      fill = 0;
      pos = 0;
    }
    ssize_t r = ::read(0, buf.data() + fill, buf.size() - fill);
    if (r < 0) throw std::runtime_error("read() failed");
    if (r == 0) return false;
    fill += static_cast<size_t>(r);
    return true;
  };

  refill();

  uint64_t lines = 0, kept = 0, skipped = 0, clamped = 0;
  const uint64_t report_every = 10'000'000ULL;

  while (true) {
    if (max_n && kept >= max_n) break;

    if (pos >= fill) {
      if (!refill()) break;
      if (pos >= fill) break;
    }

    while (pos < fill && is_space(buf[pos])) ++pos;
    if (pos >= fill) continue;

    // Parse exactly k=21: prefix13 + suffix8
    uint32_t prefix = 0;
    uint16_t suffix = 0;
    bool bad = false;

    for (int i = 0; i < 21; ++i) {
      while (pos >= fill) {
        if (!refill()) { bad = true; break; }
      }
      if (bad) break;

      uint8_t c = buf[pos++];
      if (is_space(c)) { bad = true; break; }
      uint8_t b = base2bits(c);
      if (b == 255) { bad = true; break; }

      if (i < 13) prefix = (prefix << 2) | static_cast<uint32_t>(b);
      else        suffix = static_cast<uint16_t>((suffix << 2) | b);
    }

    // If token longer than 21, consume to whitespace and mark bad
    if (!bad) {
      while (pos >= fill) { if (!refill()) break; }
      if (pos < fill && !is_space(buf[pos])) {
        bad = true;
        while (pos < fill && !is_space(buf[pos])) ++pos;
      }
    }

    // Skip whitespace before count
    while (!bad) {
      while (pos >= fill) {
        if (!refill()) { bad = true; break; }
      }
      if (bad) break;
      if (!is_space(buf[pos])) break;
      ++pos;
    }

    // Parse count
    uint64_t cnt64 = 0;
    bool have_digit = false;
    while (!bad) {
      while (pos >= fill) { if (!refill()) break; }
      if (pos >= fill) break;
      uint8_t c = buf[pos];
      if (c < '0' || c > '9') break;
      have_digit = true;
      cnt64 = cnt64 * 10 + static_cast<uint64_t>(c - '0');
      ++pos;
    }

    // Consume to EOL
    while (pos < fill && buf[pos] != '\n') ++pos;
    if (pos < fill && buf[pos] == '\n') ++pos;

    ++lines;

    if (bad || !have_digit || prefix >= P) { ++skipped; continue; }
    uint16_t cnt16 = clamp_u16(cnt64);
    if (cnt64 > std::numeric_limits<uint16_t>::max()) ++clamped;

    uint32_t shard = prefix >> SHIFT; // top shard_bits of prefix
    shards[shard].push_back(Rec{prefix, suffix, cnt16});
    ++kept;

    if ((kept % report_every) == 0) {
      std::cerr << "kept=" << kept << " lines=" << lines
                << " skipped=" << skipped << " clamped=" << clamped << "\n";
    }
  }

  // Sort shards by (prefix, suffix)
  std::cerr << "Sorting " << NSHARDS << " shards...\n";
  for (uint32_t s = 0; s < NSHARDS; ++s) {
    auto& v = shards[s];
    std::sort(v.begin(), v.end(), [](const Rec& a, const Rec& b) {
      if (a.prefix != b.prefix) return a.prefix < b.prefix;
      return a.suffix < b.suffix;
    });
    if ((s % 256) == 0) std::cerr << "  sorted shard " << s << " / " << NSHARDS << "\n";
  }

  // Build starts while streaming out suffix/count
  constexpr uint32_t UNSET = 0xFFFFFFFFu;
  std::vector<uint32_t> starts((size_t)P + 1, UNSET);
  starts[0] = 0;

  int fd_suf = ::open(out_suf.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0644);
  if (fd_suf < 0) throw std::runtime_error("Cannot open " + out_suf);
  int fd_cnt = ::open(out_cnt.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0644);
  if (fd_cnt < 0) throw std::runtime_error("Cannot open " + out_cnt);

  static constexpr size_t OUT_CHUNK = 8'000'000;
  std::vector<uint16_t> out_suf_buf;
  std::vector<uint16_t> out_cnt_buf;
  out_suf_buf.reserve(OUT_CHUNK);
  out_cnt_buf.reserve(OUT_CHUNK);

  auto flush_out = [&]() {
    if (out_suf_buf.empty()) return;
    write_all(fd_suf, out_suf_buf.data(), out_suf_buf.size() * sizeof(uint16_t));
    write_all(fd_cnt, out_cnt_buf.data(), out_cnt_buf.size() * sizeof(uint16_t));
    out_suf_buf.clear();
    out_cnt_buf.clear();
  };

  uint32_t last_prefix = 0;
  bool have_last = false;
  uint64_t out_pos64 = 0;

  std::cerr << "Emitting outputs and building starts...\n";
  for (uint32_t s = 0; s < NSHARDS; ++s) {
    const auto& v = shards[s];
    for (const auto& r : v) {
      if (!have_last || r.prefix != last_prefix) {
        if (starts[r.prefix] == UNSET) starts[r.prefix] = static_cast<uint32_t>(out_pos64);
        last_prefix = r.prefix;
        have_last = true;
      }
      out_suf_buf.push_back(r.suffix);
      out_cnt_buf.push_back(r.count);
      ++out_pos64;
      if (out_suf_buf.size() == OUT_CHUNK) flush_out();
    }
  }
  flush_out();
  ::close(fd_suf);
  ::close(fd_cnt);

  if (out_pos64 > std::numeric_limits<uint32_t>::max()) {
    throw std::runtime_error("Too many records for uint32 offsets (>2^32-1)");
  }
  const uint32_t N = static_cast<uint32_t>(out_pos64);

  // Finalize starts: fill UNSET by carrying forward last value; starts[P]=N
  uint32_t cur = 0;
  for (uint32_t p = 0; p < P; ++p) {
    if (starts[p] == UNSET) starts[p] = cur;
    else cur = starts[p];
  }
  starts[P] = N;

  // Write prefix_starts
  int fd_pref = ::open(out_prefix.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0644);
  if (fd_pref < 0) throw std::runtime_error("Cannot open " + out_prefix);
  write_all(fd_pref, starts.data(), ((size_t)P + 1) * sizeof(uint32_t));
  ::close(fd_pref);

  std::cerr << "DONE\n"
            << "  N=" << N << "\n"
            << "  hit-rate should now be sane with binary search.\n";
  return 0;
}
