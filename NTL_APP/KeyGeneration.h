#pragma once


#include "Polynomial.h"
#include "HFE.h"
#include "PublicKey.h"
#include <fstream>
#include "popl.h"
#include "Check.h"


namespace KeyGeneration{

	void generate(std::shared_ptr<popl::Value<long>> modulus_deg_option, std::shared_ptr<popl::Value<long>> t_option,
		std::shared_ptr<popl::Value<long>> hfe_deg_option, std::shared_ptr<popl::Value<long>> a_option,
		std::shared_ptr<popl::Value<string>> key_directory_option, std::shared_ptr<popl::Switch> perturbation_option) {
		
		long modulus_deg = 263; 
		long hfe_deg = 65;
		long t = 6;
		//long a = 7;
		std::string path = ".\\";

		if (modulus_deg_option->is_set()) {
			modulus_deg = modulus_deg_option->value();
		}
		if (hfe_deg_option->is_set()) {
			hfe_deg = hfe_deg_option->value();
		}
		if (t_option->is_set()) {
			t = t_option->value();
		}
		//if (a_option->is_set()) {
			//a = a_option->value();
		//}
		if (key_directory_option->is_set()) {
			path = key_directory_option->value() + "\\";
			Check::isDirectory(path);
		}

		GF2EX hfe;
		GF2X modulus;
		Vec<Polynomial<GF2>> system_of_polynomials;
		Mat<GF2> matrix_T;
		Vec<GF2> vector_T;
		Mat<GF2> matrix_S;
		Vec<GF2> vector_S;
		Vec<Polynomial<GF2>> public_key;

		hfe = HFE::generateHFEPolynomial(modulus_deg, hfe_deg, modulus);
		system_of_polynomials = HFE::hfeToSystemOfPolynomials(modulus_deg, hfe);

		if (perturbation_option->is_set()) {
			path += "perturbed_";

			Vec<GF2E> betas;
			Vec<Mat<GF2>> perturbation_polynomials;
			Vec<Polynomial<GF2>> pert_sys_of_polynomials;

			pert_sys_of_polynomials = HFE::perturbation(t, betas, perturbation_polynomials, modulus_deg, system_of_polynomials);
			public_key = PublicKey::getPublicKey(matrix_T, vector_T, matrix_S, vector_S, pert_sys_of_polynomials);

			std::ofstream file1(path + "private_key.key", std::ios::binary);
			file1 << matrix_T << vector_T << matrix_S << vector_S << modulus << hfe << betas << perturbation_polynomials;
			file1.close();
		}
		else {
			public_key = PublicKey::getPublicKey(matrix_T, vector_T, matrix_S, vector_S, system_of_polynomials);

			std::ofstream file1(path + "private_key.key", std::ios::binary);
			file1 << matrix_T << vector_T << matrix_S << vector_S << modulus << hfe;
			file1.close();
		}

		std::ofstream file(path + "public_key.key", std::ios::binary);
		file << public_key;
		file.close();
	}
}