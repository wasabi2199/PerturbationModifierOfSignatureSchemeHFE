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

NTL_CLIENT

// generovanie sukromneho kluca
// generovanie verejneho kluca
// generovanie podpisu - inverzna afinna tranformacia ku T, riesenie sustavy P', inverzna afinna tranformacia ku S
// overenie podpisu - dosadenie hodnoty hodnot podpisu do verejneho kluca

int main()
{
	//Vec<Polynomial<GF2>> sustava;
	//Polynomial<GF2> f1;
	//f1.fixMatrix();
	//cout << endl << f1.getQuadraticCoefficient();

	//Polynomial<GF2> f2;
	//sustava.append(f1);
	//sustava.append(f2);
	//Vec<Polynomial<GF2>> sustava2 = affineTransformation<GF2>(sustava);
	/*cout << endl << "Q' 1: " << endl << sustava2[0].getQuadraticCoefficient() << endl << "L' 1: " << endl << sustava2[0].getLinearCoefficient();
	cout <<endl<< "Q' 2: " << endl << sustava2[1].getQuadraticCoefficient() << endl << "L' 2: " << endl << sustava2[1].getLinearCoefficient();*/
	//Vec<Polynomial<GF2>> sustava3 = affineTransformation<GF2>(sustava);
	//cout << endl << "Q'' 1: " << endl << sustava3[0].getQuadraticCoefficient() << endl << "L'' 1: " << endl << sustava3[0].getLinearCoefficient() << endl << "A'' 1: " << endl << sustava3[0].getConstant();
	//cout << endl << "Q'' 2: " << endl << sustava3[1].getQuadraticCoefficient() << endl << "L'' 2: " << endl << sustava3[1].getLinearCoefficient() << endl << "A'' 2: " << endl << sustava3[1].getConstant();
	
	
	//long modulus_deg = 3;
	//long hfe_deg = 7;

	//GF2EX hfe = HFE::generateHFEPolynomial(modulus_deg, hfe_deg);

	GF2X modulus;
	cout << "zadaj modulus polynom> ";
	cin >> modulus;
	GF2E::init(modulus);
	GF2EX hfe;
	cout << "zadaj hfe polynom> ";
	cin >> hfe;
	long modulus_deg = deg(modulus);
	long hfe_deg = deg(hfe);
	cout << "modulus_deg je :" << modulus_deg << "  a hfe_deg je: " << hfe_deg << endl;
	

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
		//cout << sustava_polynomov[0].getQuadraticCoefficient() << endl << sustava_polynomov[1].getQuadraticCoefficient()<<endl<<sustava_polynomov[2].getQuadraticCoefficient();
	
	/*
	for (long i = 0; i < 5; i++) {
		Vec<GF2> hodnoty = random_vec_GF2(modulus_deg);
		GF2E hodnoty_alfa = conv<GF2E>(conv<GF2X> (hodnoty));
		GF2E hfe_eval = eval(hfe, hodnoty_alfa);
		Vec<GF2> hodnoty_eval = abs_coefficients;
		hodnoty_eval[0] += linear_coefficients[0] * hodnoty + hodnoty * quadratic_coefficients[0] * hodnoty;
		hodnoty_eval[1] += linear_coefficients[1] * hodnoty + hodnoty * quadratic_coefficients[1] * hodnoty;
		hodnoty_eval[2] += linear_coefficients[2] * hodnoty + hodnoty * quadratic_coefficients[2] * hodnoty;
		cout << hfe_eval << " " << hodnoty_eval<<endl;
	}
	*/

	return 0;
}


