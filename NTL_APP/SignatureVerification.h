#pragma once

#include "Polynomial.h"
#include "picosha2.h"
#include "Polynomial.h"
#include "HFE.h"
#include "PublicKey.h"
#include <fstream>
#include "popl.h"
#include "Check.h"

namespace SignatureVerification {


	bool verify(std::shared_ptr<popl::Value<long>> modulus_deg_option, std::shared_ptr<popl::Value<long>> t_option,
		std::shared_ptr<popl::Value<long>> hfe_deg_option, std::shared_ptr<popl::Value<long>> a_option,
		std::shared_ptr<popl::Value<string>> key_directory_option, std::shared_ptr<popl::Switch> perturbation_option,
		std::shared_ptr<popl::Value<string>> file_option, std::shared_ptr<popl::Value<string>> signature_file_option) {

		long modulus_deg = 263; 
		//long hfe_deg = 65;
		//long t = 6;
		long a = 7;
		std::string path = ".\\";
		std::string file_path = ".\\";
		std::string signature_path = ".\\";

		if (modulus_deg_option->is_set()) {
			modulus_deg = modulus_deg_option->value();
		}
		/*if (hfe_deg_option->is_set()) {
			hfe_deg = hfe_deg_option->value();
		}
		if (t_option->is_set()) {
			t = t_option->value();
		}*/
		if (a_option->is_set()) {
			a = a_option->value();
		}
		if (key_directory_option->is_set()) {
			path = key_directory_option->value() + "\\";
			Check::isDirectory(path);
		}
		else {
			cout << "parameter key directory not set";
			exit(1);
		}
		if (signature_file_option->is_set()) {
			signature_path = signature_file_option->value() + "\\";;
		}
		if (file_option->is_set()) {
			file_path = file_option->value();
		}
		else {
			cout << "parameter signed file not set";
			exit(1);
		}

		if (perturbation_option->is_set()) {
			path += "perturbed_";
			signature_path += "perturbed_";
		}


		path += "public_key.key";
		signature_path += "signature.sign";

		Check::isPathToFile(path);
		Check::isPathToFile(file_path);
		Check::isPathToFile(signature_path);
		
		std::ifstream file(path, std::ios::binary);
		Vec<Polynomial<GF2>> public_key;
		file >> public_key;
		file.close();
		
		std::ifstream file1(signature_path, std::ios::binary);
		Vec<GF2> signature;
		bool verifiedSignaturePerturbed;
		file1 >> signature;
		file1.close();

		Vec<GF2> message = Hash::fileToGF2(file_path, a, modulus_deg);
		try {
			verifiedSignaturePerturbed = Signature::verifySignature(a, signature, message, public_key, modulus_deg);
		}
		catch (std::exception &e) {
			cout << " chyba pri overeni >> " << e.what();
			exit(1);
		}
	}
}
