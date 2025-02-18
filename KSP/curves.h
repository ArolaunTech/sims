#include <vector>
#include <cmath>
#include <iostream>

#include "random.h"
#include "consts.h"

struct TaylorCurve {
	std::vector<double> coeffs;

	double evaluate(double x) {
		int numCoeffs = (int)coeffs.size();

		double y = 0;
		for (int i = numCoeffs - 1; i >= 0; i--) {
			y *= x;
			y += coeffs[i];
		}
		return y;
	}

	double derivative(double x) {
		std::size_t numCoeffs = coeffs.size();

		double y = 0;
		for (std::size_t i = numCoeffs - 1; i > 0; i--) {
			y *= x;
			y += coeffs[i]*i;
		}
		return y;
	}

	TaylorCurve derivativeFunction() {
		TaylorCurve out;
		out.coeffs.clear();

		std::size_t numCoeffs = coeffs.size();
		for (std::size_t i = 1; i < numCoeffs; i++) {
			out.coeffs.push_back(coeffs[i] * i);
		}
		return out;
	}

	void mutate(double bias, double deviation) {
		std::size_t numCoeffs = coeffs.size();

		for (std::size_t i = 0; i < numCoeffs; i++) {
			coeffs[i] *= bias + deviation * uniform(generator) * normal(generator);
		}
	}

	void multiply(double a) {
		std::size_t numCoeffs = coeffs.size();

		for (std::size_t i = 0; i < numCoeffs; i++) {
			coeffs[i] *= a;
		}
	}

	void initRandom(int numCoeffs, double length, double bound) {
		coeffs.clear();

		for (int i = 0; i < numCoeffs; i++) {
			coeffs.push_back(
				bound * (normal(generator)) / 
				std::pow(length, 1 * (double)i)
			);
		}
	}

	TaylorCurve operator+(const TaylorCurve& other) {
		std::size_t resultSize = coeffs.size();
		std::size_t otherSize = other.coeffs.size();
		if (otherSize < resultSize) resultSize = otherSize;

		TaylorCurve out;
		out.coeffs.clear();

		for (std::size_t i = 0; i < resultSize; i++) {
			out.coeffs.push_back(coeffs[i] + other.coeffs[i]);
		}
		return out;
	}

	TaylorCurve operator-(const TaylorCurve& other) {
		std::size_t resultSize = coeffs.size();
		std::size_t otherSize = other.coeffs.size();
		if (otherSize < resultSize) resultSize = otherSize;

		TaylorCurve out;
		out.coeffs.clear();

		for (std::size_t i = 0; i < resultSize; i++) {
			out.coeffs.push_back(coeffs[i] - other.coeffs[i]);
		}
		return out;
	}

	TaylorCurve operator*(const double other) {
		std::size_t numCoeffs = coeffs.size();

		TaylorCurve out;
		out.coeffs.clear();

		for (std::size_t i = 0; i < numCoeffs; i++) {
			out.coeffs.push_back(coeffs[i] * other);
		}
		return out;
	}

	void print() {
		std::size_t numCoeffs = coeffs.size();
		for (std::size_t i = 0; i < numCoeffs; i++) {
			std::cout << coeffs[i] << "x^" << i << "+\n";
		}
	}
};

TaylorCurve operator*(const double lhs, const TaylorCurve& rhs) {
	std::size_t numCoeffs = rhs.coeffs.size();

	TaylorCurve out;
	out.coeffs.clear();

	for (std::size_t i = 0; i < numCoeffs; i++) {
		out.coeffs.push_back(rhs.coeffs[i] * lhs);
	}
	return out;
}

struct Complex {
	double a;
	double b;

	Complex& operator+=(const Complex& rhs);
	Complex& operator*=(const Complex& rhs);
};

Complex operator+(const Complex& lhs, const Complex& rhs) {
	Complex out;
	out.a = lhs.a + rhs.a;
	out.b = lhs.b + rhs.b;
	return out;
}

Complex& Complex::operator+=(const Complex& rhs) {
	*this = *this + rhs;
	return *this;
}

Complex operator*(const Complex& lhs, const Complex& rhs) {
	Complex out;
	out.a = lhs.a * rhs.a - lhs.b * rhs.b;
	out.b = lhs.a * rhs.b + lhs.b * rhs.a;
	return out;
}

Complex operator*(const Complex& lhs, const double rhs) {
	Complex out;
	out.a = lhs.a * rhs;
	out.b = lhs.b * rhs;
	return out;
}

Complex& Complex::operator*=(const Complex& rhs) {
	*this = *this * rhs;
	return *this;
}

Complex exp(Complex x) {
	Complex out;
	out.a = std::exp(x.a) * std::cos(x.b);
	out.b = std::exp(x.a) * std::sin(x.b);
	return out;
}

Complex imaginaryNumber = {0, 1};

struct FourierCurve {
	double period;
	std::vector<Complex> coeffs;

	double evaluate(double x) {
		Complex y;
		y.a = 0;
		y.b = 0;

		std::size_t numCoeffs = coeffs.size();

		for (std::size_t i = 0; i < numCoeffs; i++) {
			Complex exponent;
			exponent.a = 0;
			exponent.b = 2 * pi * i * x/period;
			y += coeffs[i] * exp(exponent);
		}

		return y.a;
	}

	FourierCurve derivativeFunction() {
		FourierCurve out;
		out.period = period;
		out.coeffs.clear();
		std::size_t numCoeffs = coeffs.size();

		for (std::size_t i = 0; i < numCoeffs; i++) {
			out.coeffs.push_back(coeffs[i] * imaginaryNumber * (2 * pi * i/period));
		}
		return out;
	}

	double derivative(double x) {
		return derivativeFunction().evaluate(x);
	}

	void initRandom(int numCoeffs, double scale, double shrink, double newPeriod) {
		period = newPeriod;
		coeffs.clear();

		for (int i = 0; i < numCoeffs; i++) {
			double waveScale = scale * std::pow(shrink, (double)i);
			Complex coeff;
			coeff.a = waveScale * normal(generator);
			coeff.b = waveScale * normal(generator);
			coeffs.push_back(coeff);
		}
	}

	void mutate(double bias, double deviation) {
		std::size_t numCoeffs = coeffs.size();

		for (std::size_t i = 0; i < numCoeffs; i++) {
			Complex angleShift;
			angleShift.a = 0;
			angleShift.b = deviation * normal(generator);
			coeffs[i] *= exp(angleShift) * (bias + deviation * normal(generator));
		}
	}

	void setInitialPositionDerivative(double newPosition, double newDerivative) {
		double currDerivative = derivative(0);
		coeffs[1].b -= (newDerivative - currDerivative) * period/(2 * pi);

		double currPosition = evaluate(0);
		coeffs[0].a += newPosition - currPosition;
	}

	void print() {
		std::cout << 
			"Period: " << period << "s\n" <<
			"Coeffs:\n";

		std::size_t numCoeffs = coeffs.size();
		for (std::size_t i = 0; i < numCoeffs; i++) {
			std::cout << " - Coeff " << i << ": " << coeffs[i].a << " + " << coeffs[i].b << "i\n";
		}
		std::cout << "\n";
	}
};

FourierCurve operator+(const FourierCurve& lhs, const FourierCurve& rhs) {
	std::size_t resultSize = lhs.coeffs.size();
	std::size_t rhsSize = rhs.coeffs.size();
	if (rhsSize < resultSize) resultSize = rhsSize;

	FourierCurve out;
	out.period = lhs.period;
	out.coeffs.clear();

	for (std::size_t i = 0; i < resultSize; i++) {
		out.coeffs.push_back(lhs.coeffs[i] + rhs.coeffs[i]);
	}
	return out;
}

FourierCurve operator*(const FourierCurve& lhs, const double rhs) {
	std::size_t numCoeffs = lhs.coeffs.size();

	FourierCurve out;
	out.period = lhs.period;
	out.coeffs.clear();

	for (std::size_t i = 0; i < numCoeffs; i++) {
		out.coeffs.push_back(lhs.coeffs[i] * rhs);
	}
	return out;
}

FourierCurve operator*(const double lhs, const FourierCurve& rhs) {
	return rhs * lhs;
}

FourierCurve operator-(const FourierCurve& lhs, const FourierCurve& rhs) {
	return lhs + rhs * -1.0;
}