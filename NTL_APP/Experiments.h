#pragma once

#include "Polynomial.h"
#include "HFE.h"
#include "PublicKey.h"
#include <fstream>
#include "Signature.h"
#include"AffineTransformation.h"
#include <bitset>
#include <chrono>
#include <ctime>
#include <fstream>
#include <stdio.h>
#include<Windows.h>

namespace Experiments {

#define MeasuringStart auto start = get_cpu_time();
#define MeasuringEnd auto end = get_cpu_time(); auto CPU_time = end - start;

	double get_cpu_time() {
		FILETIME a, b, c, d;
		if (GetProcessTimes(GetCurrentProcess(), &a, &b, &c, &d) != 0) {
			return
				(double)(d.dwLowDateTime |
					((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
		}
		else {
			return 0;
		}
	}

	void experiment(){
		long modulus_deg = 263; //= deg(modulus);
		long hfe_deg = 130;
		long t = 5;
		long a = 7;

		GF2EX hfe;
		Vec<Polynomial<GF2>> sustava_polynomov;
		Mat<GF2> matrix_T;
		Vec<GF2> vector_T;
		Mat<GF2> matrix_S;
		Vec<GF2> vector_S;
		Vec<Polynomial<GF2>> public_key;
		GF2X modulus;

		hfe = HFE::generateHFEPolynomial(modulus_deg, hfe_deg, modulus);
		//cout << "hfe:" << hfe << endl;
		sustava_polynomov = HFE::hfeToSystemOfPolynomials(modulus_deg, hfe);
		cout << "sustava:" << sustava_polynomov << endl;

		public_key = PublicKey::getPublicKey(matrix_T, vector_T, matrix_S, vector_S, sustava_polynomov);
		cout << public_key;


		const std::string file_name = "wwwvystup_t_" + std::to_string(t) + "hfe" + std::to_string(hfe_deg) + "_3.csv";
		std::ofstream file(file_name);
		file << "hfe, gen public key, hfe+key gen, gen signature cpu, attempts, \
				verify signature, gen pub key pert,  hfe+pert key gen, \
				gen signature pert +, attempts +, verify signature pert +,  \
				gen signature pert proj, attempts proj, verify signature pert proj\n";
		file.flush();

		for (int i = 0; i < 1; i++) {
			std::string line = "";
			GF2EX hfe;
			Vec<Polynomial<GF2>> sustava_polynomov;
			Mat<GF2> matrix_T;
			Vec<GF2> vector_T;
			Mat<GF2> matrix_S;
			Vec<GF2> vector_S;
			Vec<Polynomial<GF2>> public_key;
			cout << "i = " << i << endl;
			{
				MeasuringStart;
				hfe = HFE::generateHFEPolynomial(modulus_deg, hfe_deg);
				//cout << "hfe deg: " << deg(hfe) << endl;
				sustava_polynomov = HFE::hfeToSystemOfPolynomials(modulus_deg, hfe);
				MeasuringEnd;
				file << CPU_time << ", ";
				file.flush();
			}

			{
				MeasuringStart;
				public_key = PublicKey::getPublicKey(matrix_T, vector_T, matrix_S, vector_S, sustava_polynomov);
				MeasuringEnd;
				file << CPU_time << ", , ";
				file.flush();
			}

			Vec<GF2> message = random_vec_GF2(modulus_deg - a);
			Vec<GF2> signature;
			bool verifiedSignature;

			{
				MeasuringStart;
				bool valid = false;
				int counter = 0;
				while (valid == false) {
					valid = Signature::generateSignature(signature, a, hfe, matrix_T, vector_T, matrix_S, vector_S, message, modulus_deg);
					counter++;
				}
				MeasuringEnd;
				file << CPU_time << ", " << counter << ", ";
				file.flush();
			}

			{
				MeasuringStart;
				verifiedSignature = Signature::verifySignature(a, signature, message, public_key, modulus_deg);
				MeasuringEnd;
				file << CPU_time << ", ";
				file.flush();
			}

			std::cout << endl << "verifiedSignature = " << verifiedSignature << endl;

			Vec<GF2E> betas;
			Vec<Mat<GF2>> perturbation_polynomials;
			Vec<Polynomial<GF2>> pert_sys_of_polynomials;
			Vec<Polynomial<GF2>> public_key_perturbed;

			{
				MeasuringStart;
				pert_sys_of_polynomials = HFE::perturbation(t, betas, perturbation_polynomials, modulus_deg, sustava_polynomov);
				public_key_perturbed = PublicKey::getPublicKey(matrix_T, vector_T, matrix_S, vector_S, pert_sys_of_polynomials);
				MeasuringEnd;
				file << CPU_time << ", , ";
				file.flush();
			}

			Vec<GF2> signature_perturbed;
			bool verifiedSignaturePerturbed;

			{
				MeasuringStart;
				bool valid = false;
				int counter = 0;
				while (valid == false) {
					valid = Signature::generateSignaturePerturbed(signature_perturbed, a, matrix_T, vector_T, matrix_S, vector_S, hfe, perturbation_polynomials, betas, message, modulus_deg, t);
					counter++;
				}
				MeasuringEnd;
				file << CPU_time << ", " << counter << ", ";
				file.flush();
			}

			{
				MeasuringStart;
				verifiedSignaturePerturbed = Signature::verifySignature(a, signature_perturbed, message, public_key_perturbed, modulus_deg);
				MeasuringEnd;
				//file << CPU_time << ", ";
				file << CPU_time << "\n";
				file.flush();
			}

			/*
			std::cout << "verifiedSignaturePerturbed = " << verifiedSignaturePerturbed << endl;


			Vec<GF2> signature_perturbed_projection;
			bool verifiedSignaturePerturbedProjection;

			{
				MeasuringStart;
				bool valid = false;
				int counter = 0;
				while (valid == false) {
					valid = Signature::generateSignaturePerturbedProjection(signature_perturbed_projection, a, pert_sys_of_polynomials, matrix_T, vector_T, matrix_S, vector_S, hfe, perturbation_polynomials, betas, message, modulus_deg, t);
					counter++;
				}
				MeasuringEnd;
				file << CPU_time << ", " << counter << ", ";
				file.flush();
			}

			{
				MeasuringStart;
				verifiedSignaturePerturbedProjection = Signature::verifySignature(a, signature_perturbed_projection, message, public_key_perturbed, modulus_deg);
				MeasuringEnd;
				file << CPU_time << "\n";
				file.flush();
			}
			std::cout << "verifiedSignaturePerturbedProjection = " << verifiedSignaturePerturbedProjection << endl;
			*/
		}
		file.close();
	}
}
