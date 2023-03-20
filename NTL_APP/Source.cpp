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
	//GF2X modulus;
	//cout << "zadaj modulus polynom> ";
	//cin >> modulus;
	//GF2E::init(modulus);
	//SetSeed(conv<ZZ>(1));

	long modulus_deg = 3; //= deg(modulus);
	long hfe_deg = 7;
	long const t = 2;
	long a = 1;
	GF2EX hfe = HFE::generateHFEPolynomial(modulus_deg, hfe_deg);

	Vec<Polynomial<GF2>> sustava_polynomov = HFE :: hfeToSystemOfPolynomials(modulus_deg, hfe);
	cout << sustava_polynomov<<endl;

	Mat<GF2> matrix_T;
	Vec<GF2> vector_T;
	Mat<GF2> matrix_S;
	Vec<GF2> vector_S;
	
	Vec<Polynomial<GF2>> public_key = PublicKey::getPublicKey(matrix_T, vector_T, matrix_S, vector_S, sustava_polynomov);

	Vec<GF2> message = random_vec_GF2(modulus_deg - a);
	Vec<GF2> signature = Signature::generateSignature(a, hfe, matrix_T, vector_T, matrix_S, vector_S, message, modulus_deg);

	cout << "message " << message << endl;
	cout << "signature " << signature << endl;

	cout << Signature::verifySignature(a, signature, message, public_key, modulus_deg) << endl;

	Vec<GF2E> betas;
	Vec<Mat<GF2>> perturbation_polynomials;
	Vec<Polynomial<GF2>> pert_sys_of_polynomials = HFE::perturbation(t, betas, perturbation_polynomials, modulus_deg, sustava_polynomov);
	cout << "//////" << endl << "betas " << betas << endl;
	cout << "perturbation_polynomials " << endl << perturbation_polynomials << endl;
	cout << "pert_sys_of_polynomials " << pert_sys_of_polynomials << endl;

	Vec<Polynomial<GF2>> public_key_perturbed = PublicKey::getPublicKey(matrix_T, vector_T, matrix_S, vector_S, pert_sys_of_polynomials);
	Vec<GF2> signature_perturbed = Signature :: generateSignaturePerturbed(a, matrix_T, vector_T, matrix_S, vector_S, hfe, perturbation_polynomials,betas, message, modulus_deg, t);
	cout << signature_perturbed << endl;
	cout << Signature::verifySignature(a, signature_perturbed, message, public_key_perturbed, modulus_deg) << endl;
	
	Vec<GF2> signature_perturbedP = Signature::generateSignaturePerturbedP(a, pert_sys_of_polynomials, matrix_T, vector_T, matrix_S, vector_S, hfe, perturbation_polynomials, betas, message, modulus_deg, t);
	cout << signature_perturbedP << endl;
	cout << Signature::verifySignature(a, signature_perturbedP, message, public_key_perturbed, modulus_deg) << endl;

	return 0;
}


