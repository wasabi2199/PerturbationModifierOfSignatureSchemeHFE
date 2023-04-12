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

NTL_CLIENT

// generovanie sukromneho kluca
// generovanie verejneho kluca
// generovanie podpisu - inverzna afinna tranformacia ku T, riesenie sustavy P', inverzna afinna tranformacia ku S
// overenie podpisu - dosadenie hodnoty hodnot podpisu do verejneho kluca

int main()
{
	long modulus_deg = 263; //= deg(modulus);
	long hfe_deg = 65;
	const long t = 6;
	long a = 7;
	/*long modulus_deg = 263; //= deg(modulus);
	long hfe_deg = 65;
	const long t = 6;
	long a = 7;*/
	GF2EX hfe = HFE::generateHFEPolynomial(modulus_deg, hfe_deg);
	cout << "hfe" << endl;

	Vec<Polynomial<GF2>> sustava_polynomov = HFE::hfeToSystemOfPolynomials(modulus_deg, hfe);
	cout << "sustava_polynomov" << endl;	

	Mat<GF2> matrix_T;
	Vec<GF2> vector_T;
	Mat<GF2> matrix_S;
	Vec<GF2> vector_S;

	Vec<Polynomial<GF2>> public_key = PublicKey::getPublicKey(matrix_T, vector_T, matrix_S, vector_S, sustava_polynomov);
	cout << "public_key" << endl;

	Vec<GF2> message = random_vec_GF2(modulus_deg - a);
	cout << "message " << message << endl;
	Vec<GF2> signature;
	bool valid = false;
	int counter = 0;
	while (valid == false) {
		valid = Signature::generateSignature(signature, a, hfe, matrix_T, vector_T, matrix_S, vector_S, message, modulus_deg);
		cout << "counter = " << counter << endl;
		counter += 1;
		if (counter >= 500) {
			break;
		}
	}
	if (counter < 500) {
		//cout << "message " << message << endl;
		cout << "signature " << signature << endl;
		bool verifiedSignature = Signature::verifySignature(a, signature, message, public_key, modulus_deg);
		if (verifiedSignature == true) {
			cout << "podpis platny" << verifiedSignature << endl;
		}
		else {
			cout << "podpis neplatny = " << verifiedSignature << endl;
		}
		
	}


	Vec<GF2E> betas;
	Vec<Mat<GF2>> perturbation_polynomials;
	Vec<Polynomial<GF2>> pert_sys_of_polynomials = HFE::perturbation(t, betas, perturbation_polynomials, modulus_deg, sustava_polynomov);
	//cout << "betas " << betas << endl;
	//cout << "perturbation_polynomials " << endl << perturbation_polynomials << endl;
	//cout << "pert_sys_of_polynomials " << pert_sys_of_polynomials << endl;

	Vec<Polynomial<GF2>> public_key_perturbed = PublicKey::getPublicKey(matrix_T, vector_T, matrix_S, vector_S, pert_sys_of_polynomials);
	Vec<GF2> signature_perturbed;
	valid = false;
	counter = 0;
	while (valid == false) {
		valid = Signature::generateSignaturePerturbed(signature_perturbed, a, matrix_T, vector_T, matrix_S, vector_S, hfe, perturbation_polynomials, betas, message, modulus_deg, t);
		counter += 1;
		if (counter >= 20) {
			break;
		}
	}
	if (counter < 20) {
		cout << "signature_perturbed " << endl;
		cout << "podpis platny = " << Signature::verifySignature(a, signature_perturbed, message, public_key_perturbed, modulus_deg) << endl;
	}
	else {
		cout << "//////////////////////////////////////////" << endl;
	}
	

	Vec<GF2> signature_perturbed_projection;
	valid = false;
	counter = 0;
	while (valid == false) {
		valid = Signature::generateSignaturePerturbedProjection(signature_perturbed_projection, a, pert_sys_of_polynomials, matrix_T, vector_T, matrix_S, vector_S, hfe, perturbation_polynomials, betas, message, modulus_deg, t);
		counter += 1;
		if (counter >= 20) {
			break;
		}
	}
	if (counter < 20) {
		cout << "signature_perturbed_projection " << signature_perturbed_projection << endl;
		cout << "podpis platny = " << Signature::verifySignature(a, signature_perturbed_projection, message, public_key_perturbed, modulus_deg) << endl;
	}
	else {
		cout << "//////////////////////////////////////////" << endl;
	}

	return 0;
}


