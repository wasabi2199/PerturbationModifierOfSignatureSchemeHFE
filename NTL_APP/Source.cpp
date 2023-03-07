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

NTL_CLIENT

// generovanie sukromneho kluca
// generovanie verejneho kluca
// generovanie podpisu - inverzna afinna tranformacia ku T, riesenie sustavy P', inverzna afinna tranformacia ku S
// overenie podpisu - dosadenie hodnoty hodnot podpisu do verejneho kluca

int main()
{
	GF2X modulus;
	cout << "zadaj modulus polynom> ";
	cin >> modulus;
	GF2E::init(modulus);
	//GF2EX hfe;
	//cout << "zadaj hfe polynom> ";
	//cin >> hfe;
	//SetSeed(conv<ZZ>(1));
	long modulus_deg = deg(modulus);
	long hfe_deg = 7;
	GF2EX hfe = HFE::generateHFEPolynomial(modulus_deg, hfe_deg);
	cout << "hfe " << hfe << endl;
	
	cout << "modulus_deg je :" << modulus_deg << "  a hfe_deg je: " << hfe_deg << endl;
	
	//long modulus_deg = 3;
	//long hfe_deg = 6;
	
	Vec<GF2> abs_coefficients = HFE::getAbsoluteCoefficients(modulus_deg, hfe);
	Vec<Vec<GF2>> linear_coefficients = HFE::getLinearCoefficients(modulus_deg, hfe);
	Vec<Mat<GF2>> quadratic_coefficients = HFE::getQuadraticCoefficients(modulus_deg, hfe);

	Vec<Polynomial<GF2>> sustava_polynomov;
	for (int i = 0; i < modulus_deg; i++) {
		Polynomial<GF2> temp(modulus_deg, quadratic_coefficients[i], linear_coefficients[i], abs_coefficients[i]);
		temp.fixMatrix();
		sustava_polynomov.append(temp);
	}

	cout << "abs_coefficients: " << abs_coefficients << endl;
	cout << "linear_coefficients: " << linear_coefficients << endl;
	cout << "quadratic_coefficients: " << endl;
	for (int i = 0; i < quadratic_coefficients.length();i++) {
		cout<< quadratic_coefficients[i]<<endl;
	}
	
	cout << endl << "pred transformaciou: " << endl;
	for (int i = 0; i < sustava_polynomov.length(); i++) {
		cout << endl << "quadr coeff: " << endl << sustava_polynomov[i].getQuadraticCoefficient() << endl;
		cout << "lin coeff: " << sustava_polynomov[i].getLinearCoefficient() << endl;
		cout << "abs coeff: " << sustava_polynomov[i].getConstant() << endl;
	}

	Mat<GF2> matrix_T;
	Vec<GF2> vector_T;
	Vec<Polynomial<GF2>> sustava2 = AffineTransformation::affineTransformation(sustava_polynomov, matrix_T, vector_T);
	cout << endl<<"prva transformacia: " << endl;
	for (int i = 0; i < sustava2.length(); i++) {
		cout << endl<<"quadr coeff: "<<endl << sustava2[i].getQuadraticCoefficient() << endl;
		cout << "lin coeff: "<<sustava2[i].getLinearCoefficient() << endl;
		cout << "abs coeff: " << sustava2[i].getConstant() << endl;
	}

	Mat<GF2> matrix_S;
	Vec<GF2> vector_S;
	Vec<Polynomial<GF2>> sustava3 = AffineTransformation::affineTransformationS(sustava2, matrix_S, vector_S);
	cout <<endl<< "druha transformacia: " << endl;
	for (int i = 0; i < sustava3.length(); i++) {
		cout << endl<<"quadr coeff: " << endl << sustava3[i].getQuadraticCoefficient() << endl;
		cout << "lin coeff: " << sustava3[i].getLinearCoefficient() << endl;
		cout << "abs coeff: " <<  sustava3[i].getConstant() << endl;
	}

	
	for (long i = 0; i < 5; i++) {
		Vec<GF2> hodnoty = random_vec_GF2(modulus_deg);
		//cout << "hodnoty "<<hodnoty << endl;
		GF2E hodnoty_alfa = conv<GF2E>(conv<GF2X> (hodnoty));
		GF2E hfe_eval = eval(hfe, hodnoty_alfa);
		Vec<GF2> hodnoty_PK;
		for (int i = 0; i < modulus_deg; i++) {
			hodnoty_PK.append(sustava3[i].getConstant());
			hodnoty_PK[i] += sustava3[i].getLinearCoefficient() * hodnoty + hodnoty * sustava3[i].getQuadraticCoefficient() * hodnoty;

		}
		
		cout << "hodnoty pk "<<hodnoty_PK << endl;

		Vec<GF2> hodnoty_ = hodnoty * matrix_T + vector_T;
		Vec<GF2> hodnoty_T;
		for (int i = 0; i < modulus_deg; i++) {
			hodnoty_T.append(sustava_polynomov[i].getConstant());
			hodnoty_T[i] += sustava_polynomov[i].getLinearCoefficient() * hodnoty_ + hodnoty_ * sustava_polynomov[i].getQuadraticCoefficient() * hodnoty_;

		}
		//cout <<"hodnoty T "<< hodnoty_ << endl;

		//cout <<"hodnoty TF "<< hodnoty_T << endl;
		hodnoty_ = hodnoty_T * matrix_S + vector_S;
		cout << "hodnoty TFS "<<hodnoty_ << endl;
		cout << "/////////////"<<endl;
	}

	
	Vec<GF2> hodnoty_T;
	Vec<GF2> hodnoty_TP;
	Vec<GF2> hodnoty_TPS;
	
	//cout <<"hodnoty T "<< hodnoty_ << endl;

	//cout <<"hodnoty TF "<< hodnoty_T << endl;
	hodnoty_TPS = random_vec_GF2(modulus_deg);
	hodnoty_TP = AffineTransformation::applyInverseAffineTransformation(hodnoty_TPS, matrix_S, vector_S);
	//cout << x<<endl<<hodnoty_T<<endl<<hodnoty_TP<<endl<<hodnoty_TPS;
	GF2E Y = conv<GF2E>(conv<GF2X>(hodnoty_TP));
	GF2EX hfe_y = hfe - Y;
	MakeMonic(hfe_y);
	Vec<Pair<GF2EX, long>> korene;
	korene = berlekamp(hfe_y);
	GF2E X;
	bool existuje_koren = false;
	for (auto c : korene) {
		if (deg(c.a) == 1) {
			X = ConstTerm(c.a);
			existuje_koren = true;
			break;
		}
	}
	if (existuje_koren == false) {
		cout << "nema riesenie"<<endl;
	}
	Vec<GF2> x = conv<Vec<GF2>>(conv<GF2X>(X));
	x.SetLength(modulus_deg);
	x = AffineTransformation::applyInverseAffineTransformation(x, matrix_T, vector_T);

	Vec<GF2> hodnoty_PK;
	for (int i = 0; i < modulus_deg; i++) {
		hodnoty_PK.append(sustava3[i].getConstant());
		hodnoty_PK[i] += sustava3[i].getLinearCoefficient() * x + x * sustava3[i].getQuadraticCoefficient() * x;

	}
	cout << hodnoty_PK << " " << hodnoty_TPS;
	return 0;
}


