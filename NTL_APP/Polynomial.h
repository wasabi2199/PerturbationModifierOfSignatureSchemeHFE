#pragma once

#include<NTL/ZZ_p.h>
#include<NTL/vec_ZZ_p.h>
#include<NTL/mat_ZZ_p.h>
#include<type_traits>
#include<NTL/GF2.h>
#include<NTL/vec_GF2.h>
#include<NTL/mat_GF2.h>
#include <iostream>

NTL_CLIENT

template <class T> 
class Polynomial {
public:
	/*
	Polynomial() {
		cout << "matica kvadratickych koeficientov: ";
		cin >> m_quadratic_coefficient;
		cout << "vektor linearnych koeficientov: ";
		cin >> m_linear_coefficient;
		cout << "konstanta: ";
		cin >> m_constant;
	}
	*/
	Polynomial() {
	}
	Polynomial(int n) {
		m_quadratic_coefficient.SetDims(n, n);
		m_linear_coefficient.SetLength(n);
		generateRandom();
	}
	Polynomial(int n, const Mat<T> &quadraticCoefficient, const Vec<T> &linearCoefficientl, const T &constant) {
		m_quadratic_coefficient.SetDims(n, n);
		m_linear_coefficient.SetLength(n);
		this->m_constant = constant;
		this->m_linear_coefficient = linearCoefficientl;
		this->m_quadratic_coefficient = quadraticCoefficient;
	}
	Mat<T> getQuadraticCoefficient() const {
		return m_quadratic_coefficient;
	}
	Vec<T> getLinearCoefficient() const {
		return m_linear_coefficient;
	}
	T getConstant() const {
		return m_constant;
	}
	void setQuadraticCoefficient(const Mat<T>& input) {
		m_quadratic_coefficient = input;
	}
	void setLinearCoefficient(const Vec<T>& input) {
		m_linear_coefficient = input;
	}
	void setConstant(const T& input) {
		m_constant = input;
	}
	void fixMatrix() {
		for (int i = 0; i < m_quadratic_coefficient.NumCols(); i++) {
			for (int j = i; j < m_quadratic_coefficient.NumCols(); j++) {
				if (i < j) {
					m_quadratic_coefficient[i][j] += m_quadratic_coefficient[j][i];
					m_quadratic_coefficient[j][i] = 0;
				}
			}
		}
	}
	
	friend istream& operator>>(istream& os, Polynomial<T>& poly)
	{
		//os << endl << "quadr coeff: " << endl << poly.getQuadraticCoefficient() << endl;
		//os << "lin coeff: " << poly.getLinearCoefficient() << endl;
		//os << "abs coeff: " << poly.getConstant() << endl;
		os >> poly.m_quadratic_coefficient >> poly.m_linear_coefficient >> poly.m_constant;
		return os;
	}

	friend istream& operator>>(istream& os, Vec<Polynomial<T>>& polynomials)
	{
		int counter = 256;
		while(counter>0) {
			Polynomial poly;
			os >> poly;
			polynomials.append(poly);
			//poly.print();
			counter--;
		}
		return os;
	}
	void print() {
		cout << "lin coeff: " << this.getLinearCoefficient() << endl;
	}
private:
	void generateRandom(){
		for (int i = 0; i < m_quadratic_coefficient.NumRows(); i++) {
			for (int j = i; j < m_quadratic_coefficient.NumCols(); j++) {
				/*if (std::is_same<T, ZZ_p>::value) {
					m_quadratic_coefficient[i][j] = random_ZZ_p();
				}*/
				if (std::is_same<T, GF2>::value) {
					m_quadratic_coefficient[i][j] = random_GF2();
				}
			}
		}
		for (int i = 0; i < m_linear_coefficient.length(); i++) {
			/*if (std::is_same<T, ZZ_p>::value) {
				m_linear_coefficient[i] = random_ZZ_p();
			}*/
			if (std::is_same<T, GF2>::value) {
				m_linear_coefficient[i] = random_GF2();
			}
		}
		/*if (std::is_same<T, ZZ_p>::value) {
			m_constant = random_ZZ_p();
		}*/
		if (std::is_same<T, GF2>::value) {
			m_constant = random_GF2();
		}
	}
	 
public:
	Mat<T> m_quadratic_coefficient;
	Vec<T> m_linear_coefficient;
	T m_constant;

};

template <class T>
ostream& operator<<(ostream& os, const Polynomial<T>& poly)
{
	//os << endl << "quadr coeff: " << endl << poly.getQuadraticCoefficient() << endl;
	//os << "lin coeff: " << poly.getLinearCoefficient() << endl;
	//os << "abs coeff: " << poly.getConstant() << endl;
	os << poly.getQuadraticCoefficient() << poly.getLinearCoefficient() << poly.getConstant();
	return os;
}

template <class T>
ostream& operator<<(ostream& os, const Vec<Polynomial<T>>& polynomials)
{
	for (const auto& poly : polynomials) {
		os << poly;
	}
	return os;
}

