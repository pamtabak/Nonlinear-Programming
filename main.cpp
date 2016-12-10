#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <strings.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

double function1 (vector<double> params)
{
	double x1 = params[0];
	double x2 = params[1];
	double y  = sqrt(pow(x1,2) + pow((exp(x1) - x2), 2));
	return y;
}

vector<double> derivedFunction1 (vector<double> params)
{
	double x1 = params[0];
	double x2 = params[1];

	double derivedX1 = (exp(2*x1) - x2*exp(x1) + x1) / sqrt(pow(x1, 2) + pow((exp(x1) - x2), 2));
	double derivedX2 = (x2 - exp(x1)) / sqrt(pow(x1, 2) + pow((exp(x1) - x2), 2));

	vector<double> result;
	result.push_back(derivedX1);
	result.push_back(derivedX2);

	return result;
}

vector<double> derivedFunction2 (vector<double> params)
{
	double x1 = params[0];
	double x2 = params[1];
	double x3 = params[2];

	double derived = 1 / (2 * sqrt(x1 + x2 + x3));
	
	vector<double> result;
	result.push_back(derived);
	result.push_back(derived);
	result.push_back(derived);

	return result;
}

double function2 (vector<double> params)
{
	double x1 = params[0];
	double x2 = params[1];
	double x3 = params[2];
	double y  = sqrt(x1 + x2 + x3);
	return y;
}

// phi(t) = f(x + t*d)
double phi (double t, vector<double> x, vector<double> d, double (*function)(vector<double>))
{
	// O phi(t) será sua f(x + t*d)
	// y = x + t*d
	vector<double> y;
	for (int i = 0; i < x.size(); i++)
	{
		y.push_back(x[i] + t*d[i]);
	}

	return function(y);
}

bool vectorIsZero (vector<double> x, double eps)
{
	for (int i = 0; i < x.size(); i++)
	{
		if (fabs(x[i]) > eps)
		{
			return false;
		}
	}
	return true;
}

bool vectorHasChanged (vector<double> x0, vector<double> x1, double eps)
{
	double result = 0.0;
	for (int i = 0; i < x0.size(); i++)
	{
		result += pow (x1[i] - x0[i], 2);
	}
	result = sqrt(result);

	if (result < eps)
	{
		return false;
	}
	return true;
}

// CHECAR SE A FUNCAO E UNIMODAL
double goldenSectionSearch (double eps, double ro, vector<double> x0, vector<double> d, double (*function)(vector<double>))
{
	double teta1 = (3.0 - sqrt(5))/2;
	double teta2 = 1.0  - teta1;

	// Initializing a, b and s
	double a = 0;
	double s = ro;
	double b = 2*ro;

	while (phi(b, x0, d, function) < phi(s, x0, d, function))
	{
		a = s;
		s = b;
		b = 2*b;
	}

	double u = a + (teta1 * (b - a));
	double v = a + (teta2 * (b - a));

	while ((b - a) >  eps)
	{
		if (phi(u, x0, d, function) < phi(v, x0, d, function))
		{
			b = v;
			v = u;
			u = a + (teta1 * (b - a));
		}
		else
		{
			a = u;
			u = v;
			v = a + (teta2 * (b - a));
		}
	}

	return (u+v)/2;
}

vector<double> gradientMethod (vector<double> x0,int iterationLimit, double (*function)(vector<double>), vector<double> (*derivedFunction)(vector<double>))
{
	cout << "started gradient method" << endl;
	int k = 0;
	vector<double> lastXk = x0;
	vector<double> xk = x0;

	while (!vectorIsZero(derivedFunction(xk), 0.001) && k <= iterationLimit)
	{
		cout << "current interaction: ";
		cout << k << endl;

		vector<double> dk = derivedFunction(xk);
		for (int i = 0; i < dk.size(); i++)
		{
			dk[i] = -1.0*dk[i];
		}
		
		double tk = goldenSectionSearch(0.00001, 5, xk, dk, function);

		for (int i = 0; i < xk.size(); i++)
		{
			lastXk[i] = xk[i];
			xk[i]     = xk[i] + tk*dk[i];
		}
		k = k + 1;

		if (!vectorHasChanged(lastXk, xk, 0.00001))
		{
			break;
		}
	}

	return xk;
}

int main(int argc, char const *argv[])
{
	cout << "Hello World" << endl;
	vector<double> x0;
	x0.push_back(1.0);
	x0.push_back(1.0);
	vector<double> min = gradientMethod(x0, 100, function1, derivedFunction1);
	for (int i = 0; i < min.size(); i++)
	{
		cout << min[i] << endl;
	}
	return 0;
}