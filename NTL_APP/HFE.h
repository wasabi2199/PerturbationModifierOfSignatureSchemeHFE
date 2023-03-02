#pragma once

#include "Polynomial.h"
#include<NTL/GF2.h>
#include<NTL/vec_GF2.h>
#include<NTL/mat_GF2.h>
#include<NTL/GF2X.h>
#include<NTL/GF2XFactoring.h>
#include<NTL/GF2E.h>
#include<NTL/vec_GF2E.h>
#include<NTL/mat_GF2E.h>
#include<NTL/GF2EX.h>

namespace HFE {

	GF2EX generateHFEPolynomial(long modulus_deg, long hfe_deg) 
	{
		GF2X modulus;
		BuildIrred(modulus, modulus_deg); //nahodny ireduc polynom modulus_deg stupna
		GF2E::init(modulus);
		GF2EX hfe;
		SetCoeff(hfe, 0, random_GF2E());
		for (long i = 0; (1 << i) < hfe_deg; i++) {
			SetCoeff(hfe, (1 << i), random_GF2E());
		}
		for (long i = 0; (i << 1) < hfe_deg;i++) {
			for (long j = i + 1; ((i << 1) + (j << 1)) < hfe_deg; j++) {
				SetCoeff(hfe, ((i << 1) + (j << 1)), random_GF2E());
			}
		}
		return hfe;
	}

	GF2E getAlpha() 
	{
		GF2X alfa_temp;
		SetCoeff(alfa_temp, 1);
		return conv<GF2E>(alfa_temp);
	}

	Vec<GF2> getAbsoluteCoefficients(long modulus_deg, GF2EX hfe)
	{
		auto d = deg(hfe[0]._GF2E__rep);
		Vec<GF2> abs_coefficients;
		for (int i = 0; i <= d; i++)
		{
			abs_coefficients.append(coeff(hfe[0]._GF2E__rep, i));
		}

		while (abs_coefficients.length() < modulus_deg) {
			abs_coefficients.append(GF2(0));
		}

		return abs_coefficients;
	}

	Vec<Vec<GF2>> getLinearCoefficients(long modulus_deg, GF2EX hfe)
	{
		Vec<Vec<GF2>> linear_coefficients;
		for (int i = 0; i < modulus_deg; i++) {
			Vec<GF2> vector;
			vector.SetLength(modulus_deg);
			linear_coefficients.append(vector);
		}

		long hfe_deg = deg(hfe);
		GF2E B;  //current hfe (2^i) coefficient 
		GF2E alpha = getAlpha();
		cout << "alfa je: "<< alpha << endl;

		for (int i = 0; i < hfe_deg; i++)
		{
			if (pow(2, i) > deg(hfe)) {
				break;
			}
			B = hfe[pow(2, i)];
			//cout << "B je: " << B << endl;
			for (int j = 0; j < 3; j++) {
				//cout << "alfa^" << (j * pow(2, i)) << " je: " << power(alpha, (j * pow(2, i))) << endl;
				//cout << "B * alfa je: " << B * power(alpha, (j * pow(2, i))) << endl;
				GF2E aaa = B * power(alpha, (j * pow(2, i)));
				for (int k = 0; k < modulus_deg; k++) {
					GF2 pozicia = coeff(aaa._GF2E__rep, k);
					if (pozicia == 1) {
						linear_coefficients[k][j] += 1;
					}
				}
				//cout << "lin coeffs: " << endl << linear_coefficients << endl << endl;
			}
		}

		return linear_coefficients;
	}

	Vec<Mat<GF2>> getQuadraticCoefficients(long modulus_deg, GF2EX hfe) 
	{
		Vec<Mat<GF2>> quadratic_coefficients;
		for (int i = 0; i < modulus_deg; i++) 
		{
			Mat<GF2> matrix;
			matrix.SetDims(modulus_deg, modulus_deg);
			quadratic_coefficients.append(matrix);
		}

		long hfe_deg = deg(hfe);
		GF2E alpha = getAlpha();
		GF2E A;
		GF2E temp;
		long index;
		for (int i = 0; i < hfe_deg; i++) {
			for (int j = i; j < hfe_deg; j++) {
				if (i != j) {
					index = (pow(2, i) + pow(2, j));
					if (index <= deg(hfe)) {
						A = hfe[index];
						//cout << "A" << i << "," << j << " je: "<<A<<endl;
						for (int r = 1; r <= modulus_deg; r++) {
							for (int s = 1; s <= modulus_deg; s++) {
								//cout << "r=" << r << ",s=" << s << "   alfa^"<< (((r - 1) * pow(2, i)) + ((s - 1) * pow(2, j))) <<endl;
								temp = A * power(alpha, (((r - 1) * pow(2, i)) + ((s - 1) * pow(2, j))));
								//cout << "temp = " << temp << endl<<endl;
								for (int k = 0; k < modulus_deg; k++) {
									GF2 pozicia = coeff(temp._GF2E__rep, k);
									if (pozicia == 1) {
										//cout << "coeff je 1 pre i=" << i << ", j=" << j << ", r=" << r << ", s=" << s << endl;
										quadratic_coefficients[k][r - 1][s - 1] += 1;
									}
								}
							}
						}
					}
				}
			}
		}

		/*for (int i = 0; i < quadratic_coefficients.length(); i++) {
			cout << quadratic_coefficients[i] << endl;
		}*/

		return quadratic_coefficients;
	}



}