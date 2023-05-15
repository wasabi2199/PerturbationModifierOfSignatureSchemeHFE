#define NOMINMAX
#include<NTL/GF2X.h>
#include<NTL/GF2E.h>
#include<NTL/GF2XFactoring.h>
#include<NTL/vec_GF2E.h>
#include<NTL/mat_GF2E.h>
#include<NTL/GF2.h>
#include<NTL/vec_GF2.h>
#include<NTL/mat_GF2.h>
#include<NTL/GF2EX.h>
#include"Polynomial.h"
#include"AffineTransformation.h"
#include"HFE.h"
#include "NTL/GF2EXFactoring.h"
#include "Signature.h"
#include "PublicKey.h"
#include <bitset>
#include <chrono>
#include <ctime>
#include <fstream>
#include <stdio.h>
#include<Windows.h>
#include "Hash.h"
#include "popl.h"
#include "KeyGeneration.h"
#include "SignatureGeneration.h"
#include "SignatureVerification.h"

NTL_CLIENT

int main(int argc, char* argv[])
{
	popl::OptionParser parser("options");
	auto help_option = parser.add<popl::Switch>("h", "help", "produce help message");
	auto keygen_option = parser.add<popl::Switch>("k", "key_generation", "generates keys");
	auto sign_option = parser.add<popl::Switch>("s", "sign", "generates digital signature");
	auto verify_option = parser.add<popl::Switch>("v", "verify", "verifies digital signature");
	auto modulus_deg_option = parser.add<popl::Value<long>>("n", "modulus_degree", "sets parameter n");
	auto hfe_deg_option = parser.add<popl::Value<long>>("d", "hfe_degree", "sets parameter d");
	auto t_option = parser.add<popl::Value<long>>("t", "perturbation", "sets parameter t");
	auto a_option = parser.add<popl::Value<long>>("a", "modifier", "sets parameter a");
	auto key_directory_option = parser.add<popl::Value<std::string>>("i", "key_directory", "path to key directory");
	auto file_option = parser.add<popl::Value<std::string>>("f", "file", "path to file that should be signed");
	auto signature_file_option = parser.add<popl::Value<std::string>>("g", "signature", "path to file with digital signature");
	auto perturbation_option = parser.add<popl::Switch>("p", "perturbation_mode", "ads perturbation to the scheme");

	parser.parse(argc, argv);

	if (keygen_option->is_set()) {
		KeyGeneration::generate(modulus_deg_option, t_option, hfe_deg_option, a_option, key_directory_option, perturbation_option);
	}
	else if (sign_option->is_set()) {
		SignatureGeneration::generate(modulus_deg_option, t_option, hfe_deg_option, a_option, key_directory_option, perturbation_option, file_option, signature_file_option);
	}

	else if (verify_option->is_set()) {
		bool verified = SignatureVerification::verify(modulus_deg_option, t_option, hfe_deg_option, a_option, key_directory_option, perturbation_option, file_option, signature_file_option);
		if (verified) {
			cout << "valid signature";
		}
		else {
			cout << "invalid signature";
		}
	}
	else if (help_option->is_set()) {
		cout << parser << endl;
	}
	
	return 0;
}


