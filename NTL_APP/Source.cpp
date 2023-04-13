#include<NTL/ZZ_p.h>
#include<NTL/vec_ZZ_p.h>
#include<NTL/mat_ZZ_p.h>
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

#define MeasuringStart std::clock_t c_start = std::clock(); auto t_start = std::chrono::high_resolution_clock::now();
#define MeasuringEnd std::clock_t c_end = std::clock(); auto t_end = std::chrono::high_resolution_clock::now(); \
					auto CPU_time = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC; \
					auto WallClock_time = std::chrono::duration <double, std::milli> (t_end - t_start).count();

NTL_CLIENT

// generovanie verejneho kluca
// generovanie podpisu - inverzna afinna tranformacia ku T, riesenie sustavy P', inverzna afinna tranformacia ku S
// overenie podpisu - dosadenie hodnoty hodnot podpisu do verejneho kluca
// porovnanie metod invertovania

int main()
{
	long modulus_deg = 263; //= deg(modulus);
	//long hfe_deg = 7;
	long hfe_deg = 66;
	const long t = 3;
	long a = 7;
	/*
	long modulus_deg = 263; //= deg(modulus);
	long hfe_deg = 65;
	const long t = 6;
	long a = 7;
	*/

	/*
	GF2X modulus; BuildIrred(modulus, modulus_deg); GF2E::init(modulus);
	GF2EX hfe;
	Vec<Pair<GF2EX, long>> korene;
	clock_t clock1, clock2;

	for (long stupen = 64; stupen < 1000000; stupen *= 2)
	{
		cerr << stupen << " ";
		random(hfe, stupen); MakeMonic(hfe);
		clock1 = clock();
		berlekamp(korene, hfe);
		clock2 = clock();
		cerr << " Berlekamp done in " << (double)(((double)clock2 - clock1) / (double)CLOCKS_PER_SEC) << " s " ;
		clock1 = clock();
		CanZass(korene, hfe);
		clock2 = clock();
		cerr << " CanZass done in " << (double)(((double)clock2 - clock1) / (double)CLOCKS_PER_SEC) << " s" << endl;

	}
	return 0;
	*/

	const std::string file_name = "vystup_" + std::to_string(modulus_deg) + "_"
		+ std::to_string(hfe_deg) + "_"
		+ std::to_string(t) + "_"
		+ std::to_string(a) + ".csv";
	std::ofstream file(file_name);
	file << "hfe system cpu, hfe system wall, gen public key cpu, gen public key wall, gen signature cpu, gen signature wall, \
			verify signature cpu, verify signature wall, gen pub key pert cpu, gen pub key pert wall, \
			gen signature pert cpu, gen signature pert wall, verify signature pert cpu, verify signature pert wall, \
			gen signature pert proj cpu, gen signature pert proj wall, verify signature pert proj cpu, verify signature pert proj wall\n";
	file.flush();

	for (int i = 0; i < 10; i++) {
		std::string line = "";
		GF2EX hfe;
		Vec<Polynomial<GF2>> sustava_polynomov;
		Mat<GF2> matrix_T;
		Vec<GF2> vector_T;
		Mat<GF2> matrix_S;
		Vec<GF2> vector_S;
		Vec<Polynomial<GF2>> public_key;

		{
			MeasuringStart;
			hfe = HFE::generateHFEPolynomial(modulus_deg, hfe_deg);
			cout << "hfe deg: " << deg(hfe) << endl;
			sustava_polynomov = HFE::hfeToSystemOfPolynomials(modulus_deg, hfe);
			MeasuringEnd;
			file << CPU_time << ", " <<  WallClock_time <<", ";
			file.flush();
		}

		{
			MeasuringStart;
			public_key = PublicKey::getPublicKey(matrix_T, vector_T, matrix_S, vector_S, sustava_polynomov);
			MeasuringEnd;
			file << CPU_time << ", " << WallClock_time << ", ";
			file.flush();
		}

		Vec<GF2> message = random_vec_GF2(modulus_deg - a);
		Vec<GF2> signature;
		bool verifiedSignature;

		{
			MeasuringStart;
			bool valid = false;
			while (valid == false) {
				valid = Signature::generateSignature(signature, a, hfe, matrix_T, vector_T, matrix_S, vector_S, message, modulus_deg);
			}
			MeasuringEnd;
			file << CPU_time << ", " << WallClock_time << ", ";
			file.flush();
		}

		{
			MeasuringStart;
			verifiedSignature = Signature::verifySignature(a, signature, message, public_key, modulus_deg);
			MeasuringEnd;
			file << CPU_time << ", " << WallClock_time << ", ";
			file.flush();
		}

		std::cout << "i = " << i << endl << "verifiedSignature = " << verifiedSignature << endl;

		Vec<GF2E> betas;
		Vec<Mat<GF2>> perturbation_polynomials;
		Vec<Polynomial<GF2>> pert_sys_of_polynomials;
		Vec<Polynomial<GF2>> public_key_perturbed;

		{
			MeasuringStart;
			pert_sys_of_polynomials = HFE::perturbation(t, betas, perturbation_polynomials, modulus_deg, sustava_polynomov);
			public_key_perturbed = PublicKey::getPublicKey(matrix_T, vector_T, matrix_S, vector_S, pert_sys_of_polynomials);
			MeasuringEnd;
			file << CPU_time << ", " << WallClock_time << ", ";
			file.flush();
		}

		Vec<GF2> signature_perturbed;
		bool verifiedSignaturePerturbed;

		{
			MeasuringStart;
			bool valid = false;
			while (valid == false) {
				valid = Signature::generateSignaturePerturbed(signature_perturbed, a, matrix_T, vector_T, matrix_S, vector_S, hfe, perturbation_polynomials, betas, message, modulus_deg, t);
			}
			MeasuringEnd;
			file << CPU_time << ", " << WallClock_time << ", ";
			file.flush();
		}

		{
			MeasuringStart;
			verifiedSignaturePerturbed = Signature::verifySignature(a, signature_perturbed, message, public_key_perturbed, modulus_deg);
			MeasuringEnd;
			file << CPU_time << ", " << WallClock_time << ", ";
			file.flush();
		}

		std::cout << "verifiedSignaturePerturbed = " << verifiedSignaturePerturbed << endl;

		Vec<GF2> signature_perturbed_projection;
		bool verifiedSignaturePerturbedProjection;

		{
			MeasuringStart;
			bool valid = false;
			while (valid == false) {
				valid = Signature::generateSignaturePerturbedProjection(signature_perturbed_projection, a, pert_sys_of_polynomials, matrix_T, vector_T, matrix_S, vector_S, hfe, perturbation_polynomials, betas, message, modulus_deg, t);
			}
			MeasuringEnd;
			file << CPU_time << ", " << WallClock_time << ", ";
			file.flush();
		}

		{
			MeasuringStart;
			verifiedSignaturePerturbedProjection = Signature::verifySignature(a, signature_perturbed_projection, message, public_key_perturbed, modulus_deg);
			MeasuringEnd;
			file << CPU_time << ", " << WallClock_time << "\n";
			file.flush();
		}
		std::cout << "verifiedSignaturePerturbedProjection = " << verifiedSignaturePerturbedProjection << endl;
	}
	file.close();
	

	return 0;
}


