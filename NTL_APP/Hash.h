#pragma once

#include "Polynomial.h"
#include "picosha2.h"

namespace Hash {

	std::vector<unsigned char> messageToHash(std::string path) {
		std::ifstream f(path, std::ios::binary);
		std::vector<unsigned char> s(picosha2::k_digest_size);
		picosha2::hash256(f, s.begin(), s.end());
		return s;
	}

	Vec<GF2> fileToGF2(std::string path, long a, long modulus_deg) {
		long block_size = modulus_deg - a;
		auto hash = Hash::messageToHash(path);
		Vec<GF2> message_values;
		for (auto a : hash) {
			std::bitset<8> bit_vec(a);
			for (int i = 0; i < 8; i++) {
				if (message_values.length() == block_size) {
					return message_values;
				}
				message_values.append(GF2(bit_vec[i]));
			}
		}
		while (message_values.length() < block_size) {
			message_values.append(GF2(0));
		}
		return message_values;
	}


	
}
