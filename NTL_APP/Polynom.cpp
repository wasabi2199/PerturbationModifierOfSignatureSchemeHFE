#include "polynom.h"

Polynom::Polynom() {
	cin >> m_quadratic_coeficient;
	cin >> m_linear_coeficient;
	cin >> m_constant;
}

mat_ZZ_p Polynom::getQuadraticCoeficient() const {
	return m_quadratic_coeficient;
}

vec_ZZ_p Polynom::getLinearCoeficient() const {
	return m_linear_coeficient;
}

ZZ_p Polynom::getConstant() const {
	return m_constant;
}

void Polynom::setQuadraticCoeficient(const mat_ZZ_p &input) {
	m_quadratic_coeficient = input;
}

void Polynom::setLinearCoeficient(const vec_ZZ_p &input) {
	 m_linear_coeficient = input;
}

void Polynom::setConstant(const ZZ_p &input) {
	 m_constant = input;
}