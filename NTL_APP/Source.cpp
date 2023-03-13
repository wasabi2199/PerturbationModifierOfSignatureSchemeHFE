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
	GF2EX hfe = HFE::generateHFEPolynomial(modulus_deg, hfe_deg);

	Vec<Polynomial<GF2>> sustava_polynomov = HFE :: hfeToSystemOfPolynomials(modulus_deg, hfe);

	Mat<GF2> matrix_T;
	Vec<GF2> vector_T;
	Mat<GF2> matrix_S;
	Vec<GF2> vector_S;
	
	Vec<Polynomial<GF2>> public_key = PublicKey::getPublicKey(matrix_T, vector_T, matrix_S, vector_S, sustava_polynomov);

	Vec<GF2> message = random_vec_GF2(modulus_deg);
	Vec<GF2> signature = Signature::generateSignature(hfe, matrix_T, vector_T, matrix_S, vector_S, message, modulus_deg);

	cout << "message " << message << endl;
	cout << "signature " << signature << endl;

	cout << Signature::verifySignature(signature, message, public_key, modulus_deg);

	return 0;
}


