#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include "functions.cpp"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// g++ main.cpp -o main.out -I/usr/local/Cellar/eigen/3.2.4/include/eigen3/

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

double goldenSectionSearch (double eps, double ro, vector<double> x0, vector<double> d, double (*function)(vector<double>))
{
	double teta1 = (3.0 - sqrt(5))/2;
	double teta2 = 1.0  - teta1;
	double a     = 0;
	double s     = ro;
	double b     = 2*ro;

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

vector<double> gradientMethod (vector<double> x0, int iterationLimit, double ro, double (*function)(vector<double>), vector<double> (*derivedFunction)(vector<double>))
{
	cout << "started gradient method" << endl;
	int k                 = 0;
	vector<double> lastXk = x0;
	vector<double> xk     = x0;

	while (!vectorIsZero(derivedFunction(xk), 0.001) && k <= iterationLimit)
	{
		cout << "current interaction: ";
		cout << k << endl;

		vector<double> dk = derivedFunction(xk);
		for (int i = 0; i < dk.size(); i++)
		{
			dk[i] = -1.0*dk[i];
		}
		
		double tk = goldenSectionSearch(0.00001, ro, xk, dk, function);

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

vector<double> newtonMethod (vector<double> x0,int iterationLimit, double ro, double (*function)(vector<double>), vector<double> (*derivedFunction)(vector<double>), MatrixXd (*secondDerivateFunction)(vector<double>))
{
	cout << "started newton method" << endl;
	int k                 = 0;
	vector<double> lastXk = x0;
	vector<double> xk     = x0;

	while (!vectorIsZero(derivedFunction(xk), 0.001) && k <= iterationLimit)
	{
		cout << "current interaction: ";
		cout << k << endl;

		vector<double> firstDerivate         = derivedFunction(xk);
		MatrixXd secondDerivateMatrix        = secondDerivateFunction(xk);	
		MatrixXd secondDerivateMatrixInverse = secondDerivateMatrix.inverse();

		vector<double> dk;

		for (int i = 0; i < secondDerivateMatrixInverse.rows(); i++)
		{
			double value = 0.0;
			for (int j = 0; j < firstDerivate.size(); j++)
			{
				value += secondDerivateMatrixInverse(i,j)*firstDerivate[j];
			}
			dk.push_back(-1.0*value);
		}
		
		double tk = goldenSectionSearch(0.00001, ro, xk, dk, function);

		for (int i = 0; i < xk.size(); i++)
		{
			lastXk[i] = xk[i];
			xk[i]     = xk[i] + tk*dk[i];
		}
		k = k + 1;

		if (!vectorHasChanged(lastXk, xk, 0.0001))
		{
			break;
		}
	}

	return xk;
}

vector<double> quasiNewtonMethod (vector<double> x0, int iterationLimit, double ro, double (*function)(vector<double>), vector<double> (*derivedFunction)(vector<double>), MatrixXd (*secondDerivateFunction)(vector<double>))
{
	cout << "started quasi-newton method" << endl;
	int k                 = 0;
	vector<double> lastXk = x0;
	vector<double> xk     = x0;

	while (!vectorIsZero(derivedFunction(xk), 0.001) && k <= iterationLimit)
	{
		cout << "current interaction: ";
		cout << k << endl;

		// TO DO

		k = k + 1;
		if (!vectorHasChanged(lastXk, xk, 0.0001))
		{
			break;
		}
	}

	return xk;
}

int main(int argc, char const *argv[])
{

	Functions functions;
	double (*function)(vector<double>);
	vector<double> (*functionFirstDerivative)(vector<double>);
	MatrixXd (*functionSecondDerivative)(vector<double>);

	vector<double> x0;

	if (argc < 2)
	{
		cout << "wrong number of parameters" << endl;
		exit(0);
	}

	if (argc - 2 == 2)
	{
		// Function 1
		function                  = functions.function1;
		functionFirstDerivative   = functions.function1FirstDerivative;
		functionSecondDerivative  = functions.function1SecondDerivative;
		x0.push_back(atof(argv[2]));
		x0.push_back(atof(argv[3]));
	}
	else if (argc - 2 == 3)
	{
		// Function 2
		function                  = functions.function2;
		functionFirstDerivative   = functions.function2FirstDerivative;
		functionSecondDerivative  = functions.function2SecondDerivative;
		x0.push_back(atof(argv[2]));
		x0.push_back(atof(argv[3]));
		x0.push_back(atof(argv[4]));
	}

	vector<double> min;

	std::string method(argv[1]);
	
	if (method == "Gradient")
	{
		min = gradientMethod(x0, 100, 5, function, functionFirstDerivative);
	}
	else if (method == "Newton")
	{
		min = newtonMethod(x0, 100, 1.0, function, functionFirstDerivative, functionSecondDerivative);
	}
	else if (method == "Quasi-Newtown")
	{
		
	}

	for (int i = 0; i < min.size(); i++)
	{
		cout << min[i] << endl;
	}	
	
	return 0;
}