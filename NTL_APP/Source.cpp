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
#include<NTL/GF2EX.h>

NTL_CLIENT

// generovanie sukromneho kluca
// generovanie verejneho kluca
// generovanie podpisu - inverzna afinna tranformacia ku T, riesenie sustavy P', inverzna afinna tranformacia ku S
// overenie podpisu - dosadenie hodnoty hodnot podpisu do verejneho kluca
//aaa

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

//template <class T>
GF2EX generateHFEPolynomial(GF2X modulus, long n) {
	//BuildIrred(modulus, 3); //nahodny ireduc polynom 3. stupna
	//GF2E::init(modulus);
	cout << "zadaj modulus polynom> ";
	cin >> modulus;
	long deg_n = deg(modulus);
	GF2E::init(modulus);
	GF2EX hfe;
	//todo zadat stupen poly 7 pre 0-7?
	random(hfe, deg_n);
	cout << "zadaj hfe polynom> ";
	cin >> hfe;
	return hfe;
}

Vec<GF2> getAbsoluteCoefficients(GF2EX hfe) {
	GF2E coefficient = hfe[0];
	auto dd = deg(coefficient._GF2E__rep);
	Vec<GF2> abs_coefficients;
	for (int i = 0; i <= dd; i++)
	{
		abs_coefficients.append(coeff(coefficient._GF2E__rep, i));
	}

	//todo preco 3? ako nastavit?
	while (abs_coefficients.length() < 3) {
		abs_coefficients.append(GF2(0));
	}
	//cout << abs_coeffs;
	return abs_coefficients;
}

Vec<Vec<GF2>> getLinearCoefficients(long n, GF2X modulus, GF2EX hfe) {
	long deg_n = deg(modulus);
	Vec<Vec<GF2>> linear_coefficients;
	for (int i = 0; i < deg_n; i++) {
		Vec<GF2> vector;
		vector.SetLength(deg_n);
		linear_coefficients.append(vector);
	}

	GF2E B;
	//todo ako nastavit alfa inak ako zo vstupu [0 1]?
	GF2E alfa;
	GF2X alfa_temp;
	SetCoeff(alfa_temp, 1);
	alfa = conv<GF2E>(alfa_temp);

	cout << alfa << endl;
	//cout << "zadaj alfu ";
	//cin >> alfa;
	for (int i = 0; i < n; i++) {
		cout << endl <<"i je: " << i << endl;
		if (pow(2, i) > deg(hfe)) {
			break;
		}
		B = hfe[pow(2, i)];
		cout << "B je: " << B<<endl;
		for (int j = 0; j < 3; j++) {
			//cout << "alfa^" << (j * pow(2, i)) << " je: " << power(alfa, (j * pow(2, i))) << endl;
			cout << "B * alfa je: "<<B * power(alfa, (j * pow(2, i)))<<endl;
			GF2E aaa = B * power(alfa, (j * pow(2, i)));
			for (int k = 0; k <= n; k++) {
				GF2 pozicia = coeff(aaa._GF2E__rep, k);
				if (pozicia == 1) {
					linear_coefficients[k][j] += 1;
				}
			}
			cout << "lin coeffs: " << endl << linear_coefficients<<endl<<endl;
		}
	}
	//cout << deg(hfe);
	return linear_coefficients;
}

Vec<Mat<GF2>> getQuadraticCofficient(long deg_n, long n, GF2EX hfe) {
	Vec<Mat<GF2>> quadratic_coefficients;
	for (int i = 0; i < deg_n; i++) {
		Mat<GF2> matrix;
		matrix.SetDims(deg_n, deg_n);
		quadratic_coefficients.append(matrix);
	}
	
	GF2E alfa = hfe[0];
	GF2E A;
	GF2E temp;
	long index;
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			if (i != j) {
				index = (pow(2, i) + pow(2, j));
				if (index <= deg(hfe)) {
					A = hfe[index];
					//cout << "A" << i << "," << j << " je: "<<A<<endl;
					for (int r = 1; r <= deg_n; r++) {
						for (int s = 1; s <= deg_n; s++) {
							//cout << "r=" << r << ",s=" << s << "   alfa^"<< (((r - 1) * pow(2, i)) + ((s - 1) * pow(2, j))) <<endl;
							temp = A * power(alfa, (((r - 1) * pow(2, i)) + ((s - 1) * pow(2, j))));
							//cout << "temp = " << temp << endl<<endl;
							for (int k = 0; k <= n; k++) {
								GF2 pozicia = coeff(temp._GF2E__rep, k);
								if (pozicia == 1) {
									cout << "coeff je 1 pre i=" << i << ", j=" << j << ", r=" << r << ", s=" << s << endl;
									quadratic_coefficients[k][r - 1][s - 1] += 1;
								}
							}
						}
					}
				}
			}
			//cout << endl;
		}
	}
	return quadratic_coefficients;
}


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
	
	//stupen poly
	
	GF2X modulus;
	cout << "zadaj modulus polynom> ";
	cin >> modulus;
	long deg_n = deg(modulus);
	GF2E::init(modulus);
	GF2EX hfe;
	long n;

	GF2E alfa;
	GF2X alfa_temp;
	SetCoeff(alfa_temp, 1);
	alfa = conv<GF2E>(alfa_temp);
	//random(hfe, n);
	//for (int i = 0; i < n; i++) {
	//	cout << "koeficient X^"<< i <<": " << endl;
	//	cin >> hfe[i];
	//}
	//long deg_d = 6;
	
	//random(hfe, deg_d); 
	cout << "zadaj hfe polynom> ";
	cin >> hfe;
	n = deg(hfe);
	GF2E coefficient = hfe[0];

	auto dd = deg(coefficient._GF2E__rep);
	Vec<GF2> abs_coefficients;
	for (int i = 0; i <= dd; i++)
	{
		abs_coefficients.append(coeff(coefficient._GF2E__rep, i));
	}

	//todo preco 3? ako nastavit?
	while (abs_coefficients.length() < deg(modulus)) {
		abs_coefficients.append(GF2(0));
	}
	for (int i = 0; i < abs_coefficients.length(); i++) {
		cout<<abs_coefficients[i]<<endl;
	}
	//cout << abs_coefficients;

	Vec<Vec<GF2>> linear_coefficients;
	for (int i = 0; i < deg_n; i++) {
		Vec<GF2> vector;
		vector.SetLength(deg_n);
		linear_coefficients.append(vector);
	}
	
	GF2E B;
	//todo ako nastavit alfa inak ako zo vstupu [0 1]?
	//GF2E alfa = hfe[1];
	cout << alfa << endl;
	//cout << "zadaj alfu ";
	//cin >> alfa;
	for (int i = 0; i < n; i++) {
		cout << endl <<"i je: " << i << endl;
		if (pow(2, i) > deg(hfe)) {
			break;
		}
		B = hfe[pow(2, i)];
		cout << "B je: " << B<<endl;
		for (int j = 0; j < 3; j++) {
			//cout << "alfa^" << (j * pow(2, i)) << " je: " << power(alfa, (j * pow(2, i))) << endl;
			cout << "B * alfa je: "<<B * power(alfa, (j * pow(2, i)))<<endl;
			GF2E aaa = B * power(alfa, (j * pow(2, i)));
			for (int k = 0; k < deg_n; k++) {
				GF2 pozicia = coeff(aaa._GF2E__rep, k);
				if (pozicia == 1) {
					linear_coefficients[k][j] += 1;
				}
			}
			cout << "lin coeffs: " << endl << linear_coefficients<<endl<<endl;
		}
	}

	//cout << deg(hfe);

	Vec<Mat<GF2>> quadratic_coefficients;
	for (int i = 0; i < deg_n; i++) {
		Mat<GF2> matrix;
		matrix.SetDims(deg_n, deg_n);
		quadratic_coefficients.append(matrix);
	}

	GF2E A;
	GF2E temp;
	long index;
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			if (i != j) {
				index = (pow(2, i) + pow(2, j));
				if (index <= deg(hfe)) {
					A = hfe[index];
					//cout << "A" << i << "," << j << " je: "<<A<<endl;
					for (int r = 1; r <= deg_n; r++) {
						for (int s = 1; s <= deg_n; s++) {
							//cout << "r=" << r << ",s=" << s << "   alfa^"<< (((r - 1) * pow(2, i)) + ((s - 1) * pow(2, j))) <<endl;
							temp = A * power(alfa, (((r - 1) * pow(2, i)) + ((s - 1) * pow(2, j))));
							//cout << "temp = " << temp << endl<<endl;
							for (int k = 0; k < deg_n; k++) {
								GF2 pozicia = coeff(temp._GF2E__rep, k);
								if (pozicia == 1) {
									cout << "coeff je 1 pre i=" << i << ", j=" << j << ", r=" << r << ", s=" << s<<endl;
									quadratic_coefficients[k][r-1][s-1] += 1;
								}
							}
						}
					}
				}
			}
			//cout << endl;
		}
	}
	for (int i = 0; i < quadratic_coefficients.length(); i++) {
		cout << quadratic_coefficients[i]<<endl;
	}
	/*/Vec<Polynomial<GF2>> sustava_polynomov;
	for (int i = 0; i < deg_n; i++) {
		Polynomial<GF2> temp(deg_n,quadratic_coefficients[i], linear_coefficients[i], abs_coefficients[i]);
		temp.fixMatrix();
		sustava_polynomov.append(temp);
	}
	cout << quadratic_coefficients[0] << endl << quadratic_coefficients[1] << endl << quadratic_coefficients[2];
	//cout << sustava_polynomov[0].getQuadraticCoefficient() << endl << sustava_polynomov[1].getQuadraticCoefficient()<<endl<<sustava_polynomov[2].getQuadraticCoefficient();
	*/
	for (long i = 0; i < 5; i++) {
		Vec<GF2> hodnoty = random_vec_GF2(deg_n);
		GF2E hodnoty_alfa = conv<GF2E>(conv<GF2X> (hodnoty));
		GF2E hfe_eval = eval(hfe, hodnoty_alfa);
		Vec<GF2> hodnoty_eval = abs_coefficients;
		hodnoty_eval[0] += linear_coefficients[0] * hodnoty + hodnoty * quadratic_coefficients[0] * hodnoty;
		hodnoty_eval[1] += linear_coefficients[1] * hodnoty + hodnoty * quadratic_coefficients[1] * hodnoty;
		hodnoty_eval[2] += linear_coefficients[2] * hodnoty + hodnoty * quadratic_coefficients[2] * hodnoty;
		cout << hfe_eval << " " << hodnoty_eval<<endl;
	}
	
	return 0;
}


