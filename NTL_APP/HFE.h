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

		for (long i = 0; (1 << i) < hfe_deg;i++) {
			for (long j = i + 1; ((1 << i) + (1 << j)) < hfe_deg; j++) {
				SetCoeff(hfe, ((1 << i) + (1 << j)), random_GF2E());
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

		for (int i = 0; i < hfe_deg; i++)
		{
			if (pow(2, i) > deg(hfe)) {
				break;
			}
			B = hfe[pow(2, i)];
			for (int j = 0; j < modulus_deg; j++) {
				GF2E aaa = B * power(alpha, (j * pow(2, i)));
				for (int k = 0; k < modulus_deg; k++) {
					GF2 pozicia = coeff(aaa._GF2E__rep, k);
					if (pozicia == 1) {
						linear_coefficients[k][j] += 1;
					}
				}
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
		long long index;
		for (int i = 0; (1<<i) < hfe_deg; i++) {
			for (int j = i; (1<<j) < hfe_deg; j++) {
				if (i != j) {
					//index = (pow(2, i) + pow(2, j));
					index = (1 << i) + (1 << j);
					//todo deg+1 alebo deg?
					if (index > (deg(hfe))) {
						//cout << "break " << "i = " << i << " j = " << j << endl;
						break;
					}
					if (index <= deg(hfe)) {
						A = hfe[index];
						for (int r = 1; r <= modulus_deg; r++) {
							for (int s = 1; s <= modulus_deg; s++) {
								temp = A * power(alpha, (((r - 1) * pow(2, i)) + ((s - 1) * pow(2, j))));
								for (int k = 0; k < modulus_deg; k++) {
									GF2 position = coeff(temp._GF2E__rep, k);
									if (position == 1) {
										quadratic_coefficients[k][r - 1][s - 1] += 1;
									}
								}
							}
						}
					}
				}
			}
		}

		return quadratic_coefficients;
	}

	Vec<Polynomial<GF2>> hfeToSystemOfPolynomials(long modulus_deg, GF2EX hfe) 
	{
		Vec<GF2> abs_coefficients = getAbsoluteCoefficients(modulus_deg, hfe);
		Vec<Vec<GF2>> linear_coefficients = getLinearCoefficients(modulus_deg, hfe);
		Vec<Mat<GF2>> quadratic_coefficients = getQuadraticCoefficients(modulus_deg, hfe);
		Vec<Polynomial<GF2>> system_of_polynomials;

		for (int i = 0; i < modulus_deg; i++) {
			Polynomial<GF2> temp(modulus_deg, quadratic_coefficients[i], linear_coefficients[i], abs_coefficients[i]);
			temp.fixMatrix();
			system_of_polynomials.append(temp);
		}

		return system_of_polynomials;
	}

	Vec<Polynomial<GF2>> perturbation(long t, Vec<GF2E>& betas, Vec<Mat<GF2>>& perturbation_polynomials, long modulus_deg, const Vec<Polynomial<GF2>>& system_of_polynomials) {
		GF2E alpha = getAlpha();
		Vec<Polynomial<GF2>> pert_sys_of_poly = system_of_polynomials;
		
		for (int i = 0; i < t; i++) {
			perturbation_polynomials.append(random_mat_GF2(modulus_deg, modulus_deg));
			betas.append(random_GF2E());
		}

		for (int i = 0; i < t; i++) {
			GF2E beta_alpha = betas[i];
			Vec<GF2> beta_alpha_vec = conv<Vec<GF2>>(conv<GF2X>(beta_alpha));
			for (int k = 0; k < beta_alpha_vec.length(); k++) {
				if (IsOne(beta_alpha_vec[k])) {
					pert_sys_of_poly[k].setQuadraticCoefficient(pert_sys_of_poly[k].getQuadraticCoefficient() + perturbation_polynomials[i]);
				}
			}
		}

		return pert_sys_of_poly;
	}



}