#pragma once

#include "Polynomial.h"
#include "AffineTransformation.h"

namespace Signature {
	template <class T>
	Vec<T> generateSignature(const GF2EX& hfe,const Mat<T>& matrix_T, const Vec<T>& vector_T, const Mat<T>& matrix_S, const Vec<T>& vector_S, const Vec<GF2>& y, long modulus_deg) {
		Vec<GF2> hodnoty_TP;

		//inverzna transformacia ku T na vektor pre spravu (vektor y)
		hodnoty_TP = AffineTransformation::applyInverseAffineTransformation(y, matrix_S, vector_S);

		//invertovanie podla P', vektor hodnoty_TP zobrazime na prvok pola GF(2^n) a faktorizaciou zistime korene 
		GF2E Y = conv<GF2E>(conv<GF2X>(hodnoty_TP));
		GF2EX hfe_y = hfe - Y;
		MakeMonic(hfe_y);
		Vec<Pair<GF2EX, long>> roots;
		roots = berlekamp(hfe_y); 
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
			cout << "nema riesenie" << endl;
		}

		//invertovanie vektora x podla transformacie T (x = platny podpis)
		Vec<GF2> x = conv<Vec<GF2>>(conv<GF2X>(X));
		x.SetLength(modulus_deg);
		x = AffineTransformation::applyInverseAffineTransformation(x, matrix_T, vector_T);

		return x;
	}

	bool verifySignature(Vec<GF2> signature, Vec<GF2> message, Vec<Polynomial<GF2>> public_key, long modulus_deg) {
		Vec<GF2> PK_values;
		for (int i = 0; i < modulus_deg; i++) {
			PK_values.append(public_key[i].getConstant());
			PK_values[i] += public_key[i].getLinearCoefficient() * signature + signature * public_key[i].getQuadraticCoefficient() * signature;
		}

		cout << "PK_values " << PK_values << endl;

		if (message == PK_values) {
			return true;
		}
		return false;
	}

	
}