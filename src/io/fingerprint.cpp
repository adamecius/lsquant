#include "io/fingerprint.hpp"

#include <fstream>
#include <sstream>
#include <cstdint>
#include <cstdio>
#include <cstring>

// Self-contained SHA-256 (FIPS 180-4) + system fingerprint. No external deps.

namespace lsquant { namespace io {

namespace {

inline uint32_t rotr(uint32_t x, uint32_t n) { return (x >> n) | (x << (32 - n)); }

const uint32_t K[64] = {
	0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
	0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
	0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
	0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
	0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
	0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
	0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
	0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

struct Sha256 {
	uint32_t h[8];
	uint64_t len;
	uint8_t  buf[64];
	size_t   buflen;

	Sha256() : len(0), buflen(0) {
		h[0]=0x6a09e667; h[1]=0xbb67ae85; h[2]=0x3c6ef372; h[3]=0xa54ff53a;
		h[4]=0x510e527f; h[5]=0x9b05688c; h[6]=0x1f83d9ab; h[7]=0x5be0cd19;
	}

	void block(const uint8_t* p) {
		uint32_t w[64];
		for (int i = 0; i < 16; ++i)
			w[i] = (uint32_t(p[i*4])<<24)|(uint32_t(p[i*4+1])<<16)|(uint32_t(p[i*4+2])<<8)|uint32_t(p[i*4+3]);
		for (int i = 16; i < 64; ++i) {
			uint32_t s0 = rotr(w[i-15],7) ^ rotr(w[i-15],18) ^ (w[i-15]>>3);
			uint32_t s1 = rotr(w[i-2],17) ^ rotr(w[i-2],19) ^ (w[i-2]>>10);
			w[i] = w[i-16] + s0 + w[i-7] + s1;
		}
		uint32_t a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],hh=h[7];
		for (int i = 0; i < 64; ++i) {
			uint32_t S1 = rotr(e,6) ^ rotr(e,11) ^ rotr(e,25);
			uint32_t ch = (e & f) ^ ((~e) & g);
			uint32_t t1 = hh + S1 + ch + K[i] + w[i];
			uint32_t S0 = rotr(a,2) ^ rotr(a,13) ^ rotr(a,22);
			uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
			uint32_t t2 = S0 + maj;
			hh=g; g=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
		}
		h[0]+=a; h[1]+=b; h[2]+=c; h[3]+=d; h[4]+=e; h[5]+=f; h[6]+=g; h[7]+=hh;
	}

	void update(const uint8_t* data, size_t n) {
		len += n;
		while (n) {
			size_t take = 64 - buflen; if (take > n) take = n;
			std::memcpy(buf + buflen, data, take);
			buflen += take; data += take; n -= take;
			if (buflen == 64) { block(buf); buflen = 0; }
		}
	}

	std::string hex() {
		uint64_t bits = len * 8;
		uint8_t pad = 0x80; update(&pad, 1);
		uint8_t zero = 0x00;
		while (buflen != 56) update(&zero, 1);
		uint8_t lenbuf[8];
		for (int i = 0; i < 8; ++i) lenbuf[i] = uint8_t(bits >> (56 - i*8));
		update(lenbuf, 8);
		char out[65];
		for (int i = 0; i < 8; ++i) std::snprintf(out + i*8, 9, "%08x", h[i]);
		return std::string(out, 64);
	}
};

bool file_exists(const std::string& p) { std::ifstream f(p.c_str()); return f.good(); }

} // namespace

std::string sha256_hex(const std::string& data) {
	Sha256 s;
	s.update(reinterpret_cast<const uint8_t*>(data.data()), data.size());
	return s.hex();
}

std::string fingerprint_files(const std::vector<std::string>& paths) {
	Sha256 s;
	for (size_t k = 0; k < paths.size(); ++k) {
		const std::string& p = paths[k];
		// name + ':' + length + ':' frames each file so concatenation is unambiguous.
		std::ifstream f(p.c_str(), std::ios::binary);
		std::ostringstream header;
		if (f.is_open()) {
			std::stringstream ss; ss << f.rdbuf();
			const std::string content = ss.str();
			header << p << ":" << content.size() << ":";
			const std::string hs = header.str();
			s.update(reinterpret_cast<const uint8_t*>(hs.data()), hs.size());
			s.update(reinterpret_cast<const uint8_t*>(content.data()), content.size());
		} else {
			header << p << ":ABSENT:";
			const std::string hs = header.str();
			s.update(reinterpret_cast<const uint8_t*>(hs.data()), hs.size());
		}
	}
	return s.hex();
}

std::string fingerprint_system(const std::string& label, const std::string& operators_dir) {
	const std::string base = operators_dir + "/" + label + ".";
	// Hamiltonian first, then a fixed canonical order of common operator IDs, then
	// their descriptors. Unknown IDs are not enumerated (CSR set is open-ended), but
	// the present ones below pin the physics: H, velocities, spin, current.
	const char* ids[] = { "HAM", "VX", "VY", "VZ", "SX", "SY", "SZ", "JX", "JY", "JZ" };
	std::vector<std::string> paths;
	for (size_t k = 0; k < sizeof(ids)/sizeof(ids[0]); ++k) {
		const std::string csr = base + ids[k] + ".CSR";
		if (file_exists(csr)) paths.push_back(csr);
	}
	for (size_t k = 0; k < sizeof(ids)/sizeof(ids[0]); ++k) {
		const std::string desc = base + ids[k] + ".desc";
		if (file_exists(desc)) paths.push_back(desc);
	}
	return fingerprint_files(paths);
}

}} // namespace lsquant::io
