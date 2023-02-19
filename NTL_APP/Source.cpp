#include<NTL/ZZ_p.h>
#include<NTL/vec_ZZ_p.h>
#include<NTL/mat_ZZ_p.h>
#include<NTL/GF2X.h>
#include<NTL/GF2E.h>
#include"Polynomial.h"
#include<NTL/GF2XFactoring.h>
#include<NTL/GF2E.h>
#include<NTL/vec_GF2E.h>
#include<NTL/mat_GF2E.h>
#include<NTL/GF2.h>
#include<NTL/vec_GF2.h>
#include<NTL/mat_GF2.h>

NTL_CLIENT

// generovanie sukromneho kluca
// generovanie verejneho kluca
// generovanie podpisu - inverzna afinna tranformacia ku T, riesenie sustavy P', inverzna afinna tranformacia ku S
// overenie podpisu - dosadenie hodnoty hodnot podpisu do verejneho kluca

template <class T>
Polynomial<T> applyAffineTransformation(const Polynomial<T> &polynomial, const Mat<T> &affineTransformationMatrix, const Vec<T> &affineTransformationVector) {
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
Vec<Polynomial<T>> affineTransformation(Vec<Polynomial<T>> polynomials) {
	long n = polynomials[0].getQuadraticCoefficient().NumRows();
	/*Mat<T> matrix_T = generateAffineTransformationMatrix<T>(n); //musi byt invertovatelne nad GF2
	Vec<T> vector_T = generateAffineTransformationVector<T>(n); //nahodny vektor hodnot GF2*/
	Mat<T> matrix_T;
	Vec<T> vector_T;
	cout << endl << "transformacna matica: ";
	cin >> matrix_T;
	cout << endl << "transformacny vektor: ";
	cin >> vector_T;
	Vec<Polynomial<T>> newPolynomials;

	for (int i = 0; i < polynomials.length(); i++) {
		newPolynomials.append(applyAffineTransformation(polynomials[i], matrix_T, vector_T));
	}

	return newPolynomials;
}

template <class T>
Vec<Polynomial<T>> affineTransformationS(Vec<Polynomial<T>> polynomials) {
	long m = polynomials.length();
	long n = polynomials[0].getQuadraticCoefficient().NumRows();
	//Mat<T> matrix_T = generateAffineTransformationMatrix<T>(m); 
	//Vec<T> vector_T = generateAffineTransformationVector<T>(m);
	Mat<T> matrix_T;
	Vec<T> vector_T;
	cout << endl << "transformacna matica: ";
	cin >> matrix_T;
	cout << endl << "transformacny vektor: ";
	cin >> vector_T;
	Vec<Polynomial<T>> newPolynomials;
	Mat<T> tempCuadraticCoefficient;
	Vec<T> tempLinearCoefficient;
	T tempAbsolutCoefficient;
	
	for (int i = 0; i < m; i++) {
		tempCuadraticCoefficient.kill();
		tempLinearCoefficient.kill();
		tempAbsolutCoefficient =0 ;
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


int main()
{
	Vec<Polynomial<GF2>> sustava;
	Polynomial<GF2> f1;
	Polynomial<GF2> f2;
	sustava.append(f1);
	sustava.append(f2);
	//Vec<Polynomial<GF2>> sustava2 = affineTransformation<GF2>(sustava);
	/*cout << endl << "Q' 1: " << endl << sustava2[0].getQuadraticCoefficient() << endl << "L' 1: " << endl << sustava2[0].getLinearCoefficient();
	cout <<endl<< "Q' 2: " << endl << sustava2[1].getQuadraticCoefficient() << endl << "L' 2: " << endl << sustava2[1].getLinearCoefficient();*/
	Vec<Polynomial<GF2>> sustava3 = affineTransformation<GF2>(sustava);
	cout << endl << "Q'' 1: " << endl << sustava3[0].getQuadraticCoefficient() << endl << "L'' 1: " << endl << sustava3[0].getLinearCoefficient() << endl << "A'' 1: " << endl << sustava3[0].getConstant();
	cout << endl << "Q'' 2: " << endl << sustava3[1].getQuadraticCoefficient() << endl << "L'' 2: " << endl << sustava3[1].getLinearCoefficient() << endl << "A'' 2: " << endl << sustava3[1].getConstant();
	return 0;
}

