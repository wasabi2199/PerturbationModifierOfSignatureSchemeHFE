#pragma once

#include "Polynomial.h"
#include "AffineTransformation.h"
#include <bitset>

namespace Signature {

	bool generateSignature(Vec<GF2>& signature, long a,  const GF2EX& hfe, const Mat<GF2>& matrix_T, const Vec<GF2>& vector_T, const Mat<GF2>& matrix_S, const Vec<GF2>& vector_S, const Vec<GF2>& y_without_a, long modulus_deg) {
		Vec<GF2> hodnoty_TP;
		Vec<GF2> y = y_without_a;
		for (int i = 0; i < a; i++) {
			y.append(random_GF2());
		}

		//inverzna transformacia ku T na vektor pre spravu (vektor y)
		hodnoty_TP = AffineTransformation::applyInverseAffineTransformation(y, matrix_S, vector_S);

		//invertovanie podla P', vektor hodnoty_TP zobrazime na prvok pola GF(2^n) a faktorizaciou zistime korene 
		GF2E Y = conv<GF2E>(conv<GF2X>(hodnoty_TP));
		GF2EX hfe_y = hfe - Y;
		MakeMonic(hfe_y);
		Vec<Pair<GF2EX, long>> roots;
		roots = CanZass(hfe_y); 
		GF2E X;
		bool root_exists = false;
		for (auto c : roots) {
			if (deg(c.a) == 1) {
				X = ConstTerm(c.a); //X je koren
				root_exists = true;
				break;
			}
		}
		if (root_exists == false) {
			return false;
		}

		//invertovanie vektora x podla transformacie T (x = platny podpis)
		Vec<GF2> x = conv<Vec<GF2>>(conv<GF2X>(X));
		x.SetLength(modulus_deg);
		signature = AffineTransformation::applyInverseAffineTransformation(x, matrix_T, vector_T);

		return true;
	}

	bool verifySignature(long a, const Vec<GF2>& signature, const Vec<GF2>& message, const Vec<Polynomial<GF2>>& public_key, long modulus_deg) {
		Vec<GF2> PK_values;
		for (int i = 0; i < modulus_deg - a; i++) {
			PK_values.append(public_key[i].getConstant());
			PK_values[i] += public_key[i].getLinearCoefficient() * signature + signature * public_key[i].getQuadraticCoefficient() * signature;
		}
		if (message == PK_values) {
			return true;
		}
		return false;
	}

	bool generateSignaturePerturbed(Vec<GF2>& signature, long a, const Mat<GF2>& matrix_T, const Vec<GF2>& vector_T, const Mat<GF2>& matrix_S, const Vec<GF2>& vector_S, const GF2EX& hfe, const Vec<Mat<GF2>>& perturbation_polynomials, const Vec<GF2E>& betas, const Vec<GF2>& y_without_a, long modulus_deg, long const t) {
		Vec<GF2> hodnoty_TP;
		Vec<GF2> y = y_without_a;
		for (int i = 0; i < a; i++) {
			y.append(random_GF2());
		}

		//inverzna transformacia ku T na vektor pre spravu (vektor y)
		hodnoty_TP = AffineTransformation::applyInverseAffineTransformation(y, matrix_S, vector_S);
		bool root_exists = false;
		GF2E X;
		Vec<GF2> pert_root;
		Vec<GF2> X_vec;
		Vec<Vec<GF2>> vektory2;
		for (int i = 0; i < pow(2, t); i++) {
			Vec<GF2> bit_vec_GF2;
			std::bitset<10> bit_vec(i);
			for (int j = 0; j < t; j++) {
				bit_vec_GF2.append(GF2(bit_vec[j]));
			}
			vektory2.append(bit_vec_GF2);
		}
		for (Vec<GF2> t_values : vektory2) {
			GF2EX hfe_betas = hfe;
			for (int i = 0; i < t; i++) {
				hfe_betas += betas[i] * t_values[i];
			}
			GF2E Y = conv<GF2E>(conv<GF2X>(hodnoty_TP));
			GF2EX hfe_y = hfe_betas - Y;
			MakeMonic(hfe_y);
			Vec<Pair<GF2EX, long>> roots;
			roots = CanZass(hfe_y);
			pert_root.SetLength(t);
			for (auto c : roots) {
				if (deg(c.a) == 1) {
					X = ConstTerm(c.a); //X je koren
					//overenie ci je koren platne riesenie perturbacnych polynomov
					X_vec = conv<Vec<GF2>>(conv<GF2X>(X));
					X_vec.SetLength(modulus_deg);
					for (int i = 0; i < t; i++) {
						pert_root[i] = X_vec * perturbation_polynomials[i] * X_vec;
					}
					if (pert_root == t_values) {
						root_exists = true;
						break;
					}
				}
			}
			if (root_exists == true) {
				break;
			}
		}
		if (root_exists == false) {
			return false;
		}
		//invertovanie vektora x podla transformacie T (x = platny podpis)
		signature = AffineTransformation::applyInverseAffineTransformation(X_vec, matrix_T, vector_T);

		return true;
	}

	bool generateSignaturePerturbedProjection(Vec<GF2>& signature, long a, Vec<Polynomial<GF2>>pert_sys_of_polynomials, const Mat<GF2>& matrix_T, const Vec<GF2>& vector_T, const Mat<GF2>& matrix_S, const Vec<GF2>& vector_S, const GF2EX& hfe, const Vec<Mat<GF2>>& perturbation_polynomials, const Vec<GF2E>& betas, const Vec<GF2>& y_without_a, long modulus_deg, long const t) {
		Vec<GF2> hodnoty_TP;
		Vec<GF2> y = y_without_a;
		for (int i = 0; i < a; i++) {
			y.append(random_GF2());
		}
		//inverzna transformacia ku T na vektor pre spravu (vektor y)
		hodnoty_TP = AffineTransformation::applyInverseAffineTransformation(y, matrix_S, vector_S);
		bool root_exists = false;
		GF2E X;
		Vec<GF2> pert_root;
		Vec<GF2> X_vec;
		//ak by t bolo const
		Vec<Vec<GF2>> vektory2;
		Vec<GF2E> betas_linear;
		for (int i = 0; i < pow(2, t); i++) {
			Vec<GF2> bit_vec_GF2;
			std::bitset<10> bit_vec(i);
			GF2E temp = GF2E();
			for (int j = 0; j < t; j++) {
				bit_vec_GF2.append(GF2(bit_vec[j]));
				temp += bit_vec_GF2[j] * betas[j];
			}
			vektory2.append(bit_vec_GF2);
			betas_linear.append(temp);
		}
		GF2EX projection_poly;
		BuildFromRoots(projection_poly, betas_linear);
		GF2E hodnoty_TP_projection = eval(projection_poly, conv<GF2E>(conv<GF2X>(hodnoty_TP)));
		GF2EX hfe_projection = GF2EX();
		for (long i = 0; i <= deg(projection_poly); i++) {
			if (!IsZero(coeff(projection_poly, i))) {
				hfe_projection += coeff(projection_poly, i) * power(hfe, i);
			}
		}
		hfe_projection -= hodnoty_TP_projection;
		MakeMonic(hfe_projection);
		Vec<Pair<GF2EX, long>> roots;
		roots = CanZass(hfe_projection);
		pert_root.SetLength(modulus_deg);
		for (auto c : roots) {
			if (deg(c.a) == 1) {
				X = ConstTerm(c.a); //X je koren
				//overenie ci je koren platne riesenie perturbacnych polynomov
				X_vec = conv<Vec<GF2>>(conv<GF2X>(X));
				X_vec.SetLength(modulus_deg);
				for (int i = 0; i < modulus_deg; i++) {
					pert_root[i] = X_vec * pert_sys_of_polynomials[i].getQuadraticCoefficient() * X_vec;
					pert_root[i] += X_vec * pert_sys_of_polynomials[i].getLinearCoefficient();
					pert_root[i] += pert_sys_of_polynomials[i].getConstant();
				}
				if (pert_root == hodnoty_TP) {
					root_exists = true;
					break;
				}
			}
		}
		if (root_exists == false) {
			return false;
		}

		//invertovanie vektora x podla transformacie T (x = platny podpis)
		signature = AffineTransformation::applyInverseAffineTransformation(X_vec, matrix_T, vector_T);
		return true;
	}

}