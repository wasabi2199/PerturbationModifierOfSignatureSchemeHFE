#pragma once

#include<NTL/ZZ_p.h>
#include<NTL/vec_ZZ_p.h>
#include<NTL/mat_ZZ_p.h>

NTL_CLIENT

class Polynom {
public:
	Polynom();
	mat_ZZ_p getQuadraticCoeficient() const;
	vec_ZZ_p getLinearCoeficient() const;
	ZZ_p getConstant() const;
	void setQuadraticCoeficient(const mat_ZZ_p &input);
	void setLinearCoeficient(const vec_ZZ_p &input);
	void setConstant(const ZZ_p &input);

private:
	mat_ZZ_p m_quadratic_coeficient;
	vec_ZZ_p m_linear_coeficient;
	ZZ_p m_constant;
};