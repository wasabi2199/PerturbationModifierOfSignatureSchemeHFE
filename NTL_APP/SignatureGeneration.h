#pragma once

#include "Polynomial.h"
#include "HFE.h"
#include "PublicKey.h"
#include <fstream>
#include "popl.h"
#include "Hash.h"
#include "Signature.h"

namespace SignatureGeneration {
	void generate(std::shared_ptr<popl::Value<long>> modulus_deg_option, std::shared_ptr<popl::Value<long>> t_option,
		std::shared_ptr<popl::Value<long>> hfe_deg_option, std::shared_ptr<popl::Value<long>> a_option,
		std::shared_ptr<popl::Value<string>> key_directory_option, std::shared_ptr<popl::Switch> perturbation_option,
		std::shared_ptr<popl::Value<string>> file_option, std::shared_ptr<popl::Value<string>> signature_file_option) {

		long modulus_deg = 263; 
		//long hfe_deg = 65;
		long t = 6;
		long a = 7;
		std::string path = ".\\";
		std::string file_path = ".\\";
		std::string signature_path = ".\\";

		if (modulus_deg_option->is_set()) {
			modulus_deg = modulus_deg_option->value();
		}
		/*if (hfe_deg_option->is_set()) {
			hfe_deg = hfe_deg_option->value();
		}*/
		if (t_option->is_set()) {
			t = t_option->value();
		}
		if (a_option->is_set()) {
			a = a_option->value();
		}
		if (key_directory_option->is_set()) {
			path = key_directory_option->value() + "\\";
		}
		if (file_option->is_set()) {
			file_path = file_option->value();
			Check::isPathToFile(file_path);
		}
		else {
			cout << "parameter file not set";
			exit(1);
		}
		if (signature_file_option->is_set()) {
			signature_path = signature_file_option->value() + "\\";
			Check::isDirectory(signature_path);
		}

		GF2X modulus;
		GF2X modulus2;
		Mat<GF2> matrix_T;
		Vec<GF2> vector_T;
		Mat<GF2> matrix_S;
		Vec<GF2> vector_S;

		if (perturbation_option->is_set()) { 
			path += "perturbed_"; 
			signature_path += "perturbed_";
		}

		Check::isPathToFile(path + "private_key.key");
		std::ifstream file(path + "private_key.key", std::ios::binary);
		file >> matrix_T >> vector_T >> matrix_S >> vector_S >> modulus;
		file.close();

		GF2E::init(modulus);
		GF2EX hfe;

		Vec<GF2> message = Hash::fileToGF2(file_path, a, modulus_deg);
		Vec<GF2> signature;

		bool valid = false;
		int counter = 0;

		file.open(path + "private_key.key", std::ios::binary);
		file >> matrix_T >> vector_T >> matrix_S >> vector_S >> modulus2;
		if (perturbation_option->is_set()) {
			Vec<GF2E> betas;
			Vec<Mat<GF2>> perturbation_polynomials;
			file >> hfe >> betas >> perturbation_polynomials;
			while (valid == false) {
				valid = Signature::generateSignaturePerturbed(signature, a, matrix_T, vector_T, matrix_S, vector_S, hfe, perturbation_polynomials, betas, message, modulus_deg, t);
				counter++;
			}
		}
		else {
			file >> hfe;

			while (valid == false) {
				valid = Signature::generateSignature(signature, a, hfe, matrix_T, vector_T, matrix_S, vector_S, message, modulus_deg);
				counter++;
			}
		}
		file.close();
		
		std::ofstream file1(signature_path + "signature.sign", std::ios::binary);
		file1 << signature;
		file1.close();
	}
}