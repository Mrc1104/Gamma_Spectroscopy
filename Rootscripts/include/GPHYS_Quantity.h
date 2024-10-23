#ifndef GPHYS_QUANTITY_H
#define GPHYS_QUANTITY_H
/***
 * File: GPHYS_Quantity.h
 * Author: Ryan Conaway
 * Email: mc321015@ohio.edu
 * Modified: 09/18/2024 
 *
 * Discription:
 *	This is a container class that is used to hold 
 *	physical quantities that have an associated errors.
 *	Operators are overloaded to automatically propogate
 *	the errors.
 *
 *	If there are covariances between the two quantities,
 *	the individual mult(...), div(...), add(...), and sub(...)
 *  functions can be called to propogate the covariances
 * 
 * Operations implemented (w/ error prop.) : operator overloaded ?
 * Addition/ Subtraction 	 : Y
 * Multiplication / Division : Y
 * Logorithm                 : N
 **/
#include <iostream>
#include <cmath>

class GPHYS_Quantity
{
public:
	// Constructor
	GPHYS_Quantity(); 
	GPHYS_Quantity(const GPHYS_Quantity&);
	GPHYS_Quantity(double q, double e = 0);
	// Copy Constructor
	GPHYS_Quantity(GPHYS_Quantity&);

public:
	// Operator overloading
	GPHYS_Quantity operator*(const GPHYS_Quantity&) const;
	GPHYS_Quantity operator/(const GPHYS_Quantity&) const;
	GPHYS_Quantity operator+(const GPHYS_Quantity&) const;
	GPHYS_Quantity operator-(const GPHYS_Quantity&) const;
	GPHYS_Quantity operator=(const GPHYS_Quantity&);
	bool operator==(const GPHYS_Quantity&) const;
	bool operator!=(const GPHYS_Quantity&) const;

public:
	// Error Propogration with covariance
	void mult(const GPHYS_Quantity&, const double cov);
	void div(const GPHYS_Quantity&, const double cov);
	void add(const GPHYS_Quantity&, const double cov);
	void sub(const GPHYS_Quantity&, const double cov);
public:
	// Error Propogation w/o covariance
	void log(const GPHYS_Quantity&);

private:
	// Operator overloading (<< and >>)
	friend std::ostream& operator<<(std::ostream&, const GPHYS_Quantity&);
	friend std::istream& operator<<(std::istream&, GPHYS_Quantity&);

public:
	void set(double q, double e);
	double getQuantity() const;
	double getError() const;
private:
	// Functions used for error propogation
	GPHYS_Quantity _propogate_multiplication_error(const double A, const double SA, const double B, const double SB, const double covAB, const int sign) const;
	GPHYS_Quantity _propogate_addition_error(const double A, const double SA, const double B, const double SB, const double covAB, const int sign) const;

private:
	double _quantity;
	double _error;

};

///////////////////////////////////////////////////////////////////////////////////
////////////////////// IMPLEMENTATION /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

GPHYS_Quantity::GPHYS_Quantity() : _quantity(0), _error(0) { }
GPHYS_Quantity::GPHYS_Quantity(double q, double e)
{
	_quantity = q;
	_error    = e;
}
GPHYS_Quantity::GPHYS_Quantity(const GPHYS_Quantity& other)
{
	*this = other;
}
GPHYS_Quantity::GPHYS_Quantity(GPHYS_Quantity& other)
{
	_quantity = other._quantity;
	_error    = other._error;
}

////////////////////////////////////////////////
/* Operator Overloading */
GPHYS_Quantity GPHYS_Quantity::operator*(const GPHYS_Quantity& other) const
{
	return 	_propogate_multiplication_error(_quantity, _error, other._quantity, other._error, 0, 1);
}

GPHYS_Quantity GPHYS_Quantity::operator/(const GPHYS_Quantity& other) const
{
	return 	_propogate_multiplication_error(_quantity, _error, other._quantity, other._error, 0, -1);
}

GPHYS_Quantity GPHYS_Quantity::operator+(const GPHYS_Quantity& other) const
{
	return 	_propogate_addition_error(_quantity, _error, other._quantity, other._error, 0, 1);
}

GPHYS_Quantity GPHYS_Quantity::operator-(const GPHYS_Quantity& other) const
{
	return 	_propogate_addition_error(_quantity, _error, other._quantity, other._error, 0, -1);
}
GPHYS_Quantity GPHYS_Quantity::operator=(const GPHYS_Quantity& other)
{
	// Check for self assignment
	if (this != &other) {
		_quantity = other._quantity;
		_error    = other._error;
	}
	return *this;
}

bool GPHYS_Quantity::operator==(const GPHYS_Quantity& other) const
{
	if( abs(_quantity - other._quantity) < 1e-10 ) return true;
	return false;
}

bool GPHYS_Quantity::operator!=(const GPHYS_Quantity& other) const
{
	return !(*this==other);
}
std::ostream& operator<<(std::ostream &os, const GPHYS_Quantity& tmp)
{
	return os << tmp._quantity << " +- " << tmp._error;
}

std::istream& operator<<(std::istream& is, GPHYS_Quantity& tmp)
{
	is >> tmp._quantity >> tmp._error;
	return is;
}
/* Operator Overloading */
////////////////////////////////////////////////




////////////////////////////////////////////////
/* Error Propogation Functions */

/*
 * Func: _propogate_multiplication_error
 * Brief:
 *	Calculates the error resulting from multiplying or dividing two quantities 
 *  A +- sa and B +- sb
 * Formula:
 *	F = A*B OR F = A/B
 *  sF = F^2 * { ( sA / A )^2 + ( sb / B )^2 +- 2 * cov(A,B) / (A*B) }
 * Comment:
 *	cov(A,B) is 0 for two independent variables.
 *	Multiplying has [ + cov(A,B)/ (A*B) ]
 *	Dividing    has [ - cov(A,B)/ (A*B) ]
 */
GPHYS_Quantity GPHYS_Quantity::_propogate_multiplication_error(const double A, const double SA, const double B, const double SB, const double covAB, const int sign) const
{
	GPHYS_Quantity tmp;
	if( sign == 1) tmp._quantity = A * B;
	else           tmp._quantity = A / B;
	double term1 = pow( (SA / A) , 2);	
	double term2 = pow( (SB / B) , 2);
	double term3 = sign*2 * covAB / (A*B) ;
	double e_sq  = pow(tmp._quantity , 2) * (term1 + term2 + term3);
	tmp._error   = sqrt(   e_sq     );
	return tmp;
}

/*
 * Func: _propogate_addition_error
 * Brief:
 *	Calculates the error resulting from adding or subtracting two quantities 
 *  A +- sa and B +- sb
 * Formula:
 *	F = A+-B
 *  sF =  ( sA )^2 + ( sb )^2 +- 2 * cov(A,B) 
 * Comment:
 *	cov(A,B) is 0 for two independent variables.
 *	Adding      has [ + cov(A,B)/ (A*B) ]
 *	Subtracting has [ - cov(A,B)/ (A*B) ]
 */
GPHYS_Quantity GPHYS_Quantity::_propogate_addition_error(const double A, const double SA, const double B, const double SB, const double covAB, const int sign) const
{
	GPHYS_Quantity tmp;
	if (sign == 1) tmp._quantity = A + B;
	else 		   tmp._quantity = A - B;
	double term1 = pow( SA , 2);
	double term2 = pow( SB,  2);
	double term3 = sign*2*covAB;
	double e_sq  = term1 + term2 + term3;
	tmp._error   = sqrt( e_sq );
	return tmp;
}
/* Error Propogation Functions */
////////////////////////////////////////////////


////////////////////////////////////////////////
/* Error Propogation w/ Covariance*/

void GPHYS_Quantity::mult(const GPHYS_Quantity& other, const double cov)
{
	GPHYS_Quantity tmp = _propogate_multiplication_error(_quantity, _error, other._quantity, other._error, cov, 1);
	_quantity = tmp._quantity;
	_error    = tmp._error;
}
void GPHYS_Quantity::div(const GPHYS_Quantity& other, const double cov)
{
	GPHYS_Quantity tmp = _propogate_multiplication_error(_quantity, _error, other._quantity, other._error, cov, -1);
	_quantity = tmp._quantity;
	_error    = tmp._error;
}
void GPHYS_Quantity::add(const GPHYS_Quantity& other, const double cov)
{
	GPHYS_Quantity tmp = _propogate_addition_error(_quantity, _error, other._quantity, other._error, cov, 1);
	_quantity = tmp._quantity;
	_error    = tmp._error;
}

void GPHYS_Quantity::sub(const GPHYS_Quantity& other, const double cov)
{
	GPHYS_Quantity tmp = _propogate_addition_error(_quantity, _error, other._quantity, other._error, cov, -1);
	_quantity = tmp._quantity;
	_error    = tmp._error;
}
/* Error Propogation w/ Covariance*/
////////////////////////////////////////////////

////////////////////////////////////////////////
/* Error Propogation w/o Covariance*/
/*
 * Func: log()
 * Brief:
 *	Calculates the error taking the log of a quantity w/ error 
 * Formula:
 *	F = Ln(A)
 *  sF = sA / A
 */
void GPHYS_Quantity::log(const GPHYS_Quantity& other)
{
	GPHYS_Quantity tmp;
	tmp = other;
	_quantity = std::log(tmp._quantity);
	_error    = tmp._error / _quantity;
}
/* Error Propogation w/o Covariance*/
////////////////////////////////////////////////


void GPHYS_Quantity::set(double q, double e)
{
	_quantity = q;
	_error    = e;
}

double GPHYS_Quantity::getQuantity() const
{
	return _quantity;
}

double GPHYS_Quantity::getError() const
{
	return _error;
}

#endif
