#pragma once

#include "Polynomial.h"

namespace AffineTransformation {
	template <class T>
	Polynomial<T> applyAffineTransformation(const Polynomial<T>& polynomial, const Mat<T>& affineTransformationMatrix, const Vec<T>& affineTransformationVector) {
		int n = polynomial.getQuadraticCoefficient().NumRows();
		Polynomial<T> newPolynomial(n);

		newPolynomial.setQuadraticCoefficient(affineTransformationMatrix * polynomial.getQuadraticCoefficient() * transpose(affineTransformationMatrix));

		newPolynomial.setLinearCoefficient(affineTransformationVector * polynomial.getQuadraticCoefficient() * transpose(affineTransformationMatrix) +
			affineTransformationVector * transpose(polynomial.getQuadraticCoefficient()) * transpose(affineTransformationMatrix) +
			polynomial.getLinearCoefficient() * transpose(affineTransformationMatrix));

		newPolynomial.setConstant(affineTransformationVector * polynomial.getQuadraticCoefficient() * affineTransformationVector +
			polynomial.getLinearCoefficient() * affineTransformationVector + polynomial.getConstant());

		return newPolynomial;
	}

	template <class T>
	Vec<T> applyInverseAffineTransformation(const Vec<T>& y, const Mat<T>& affineTransformationMatrix, const Vec<T>& affineTransformationVector) {
		Vec<T> x;
		Mat<T> inverseMatrix;
		long m = affineTransformationMatrix.NumCols();
		inverseMatrix.SetDims(m, m);
		inv(inverseMatrix, affineTransformationMatrix);
		x = (y - affineTransformationVector) * inverseMatrix ;
		return x;
	}

	template <class T>
	Mat<T> generateAffineTransformationMatrix(long n) {
		Mat<T> transformationMatrix;
		Mat<T> temp_matrix;

		transformationMatrix.SetDims(n, n);
		temp_matrix.SetDims(n, n);

		int c = 0;
		while (1)
		{
			random(transformationMatrix, n, n);
			temp_matrix = transformationMatrix;
			if (gauss(temp_matrix) == n) {
				break;
			}
		}

		return transformationMatrix;
	}

	template <class T>
	Vec<T> generateAffineTransformationVector(long n) {
		Vec<T> transformationVector;
		transformationVector.SetLength(n);
		random(transformationVector, n);
		return transformationVector;
	}

	template <class T>
	Vec<Polynomial<T>> affineTransformation(Vec<Polynomial<T>> polynomials, Mat<T>& matrix_T, Vec<T>& vector_T) {
		long n = polynomials[0].getQuadraticCoefficient().NumRows();
		matrix_T = generateAffineTransformationMatrix<T>(n); //musi byt invertovatelne nad GF2
		vector_T = generateAffineTransformationVector<T>(n); //nahodny vektor hodnot GF2
		//cout << endl << "matrix_T: " << endl << matrix_T<<endl;
		//cout << "vector_T: " << endl << vector_T<<endl;
		Vec<Polynomial<T>> newPolynomials;

		for (int i = 0; i < polynomials.length(); i++) {
			newPolynomials.append(applyAffineTransformation(polynomials[i], matrix_T, vector_T));
		}

		return newPolynomials;
	}

	template <class T>
	Vec<Polynomial<T>> affineTransformationS(Vec<Polynomial<T>> polynomials, Mat<T>& matrix_T, Vec<T>& vector_T) {
		long m = polynomials.length();
		long n = polynomials[0].getQuadraticCoefficient().NumRows();
		matrix_T = generateAffineTransformationMatrix<T>(m);
		vector_T = generateAffineTransformationVector<T>(m);
		//cout << endl << "matrix_S: " << endl << matrix_T << endl;
		//cout << "vector_S: " << endl << vector_T << endl;
		Vec<Polynomial<T>> newPolynomials;
		Mat<T> tempCuadraticCoefficient;
		Vec<T> tempLinearCoefficient;
		T tempAbsolutCoefficient;

		for (int i = 0; i < m; i++) {
			tempCuadraticCoefficient.kill();
			tempLinearCoefficient.kill();
			tempAbsolutCoefficient = 0;
			tempCuadraticCoefficient.SetDims(n, n);
			tempLinearCoefficient.SetLength(n);

			for (int j = 0; j < m; j++) {
				tempCuadraticCoefficient += matrix_T[j][i] * polynomials[j].getQuadraticCoefficient();
				tempLinearCoefficient += matrix_T[j][i] * polynomials[j].getLinearCoefficient();
				tempAbsolutCoefficient += matrix_T[j][i] * polynomials[j].getConstant();
			}

			tempAbsolutCoefficient += vector_T[i];
			Polynomial<T> tempPolynomial(n, tempCuadraticCoefficient, tempLinearCoefficient, tempAbsolutCoefficient);
			newPolynomials.append(tempPolynomial);
		}

		return newPolynomials;
	}
}
