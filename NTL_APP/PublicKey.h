#pragma once

#include "Polynomial.h"
#include "AffineTransformation.h"

namespace PublicKey {

	template <class T>
	Vec<Polynomial<T>> getPublicKey(Mat<T>& matrix_T, Vec<T>& vector_T, Mat<T>& matrix_S, Vec<T>& vector_S,const Vec<Polynomial<GF2>>& P) {
		
		Vec<Polynomial<T>> TP = AffineTransformation::affineTransformation(P, matrix_T, vector_T);
		Vec<Polynomial<T>> TPS = AffineTransformation::affineTransformationS(TP, matrix_S, vector_S);

		return TPS;
	}

}
