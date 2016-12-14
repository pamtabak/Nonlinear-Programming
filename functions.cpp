#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class Functions
{
public:
	Functions()
	{

	}

	~Functions()
	{

	}

	static double function1 (vector<double> params)
	{	
		double x1 = params[0];
		double x2 = params[1];
		double y  = sqrt(pow(x1,2) + pow((exp(x1) - x2), 2));
		return y;
	}

	static vector<double> function1FirstDerivative (vector<double> params)
	{
		double x1 = params[0];
		double x2 = params[1];
	
		double derivedX1 = (exp(2*x1) - x2*exp(x1) + x1) / (sqrt(pow(x1, 2) - 2*exp(x1)*x2 + exp(2*x1) + pow(x2,2)));
		double derivedX2 = (x2 - exp(x1)) / (sqrt(pow(x1, 2) + pow((exp(x1) - x2), 2)));
	
		vector<double> result;
		result.push_back(derivedX1);
		result.push_back(derivedX2);
	
		return result;
	}

	static MatrixXd function1SecondDerivative (vector<double> params)
	{
		double x1 = params[0];
		double x2 = params[1];

		double derivedX1X1 = (2*exp(x1)*(exp(x1) - x2) + 2*exp(2*x1) + 2)/(2*sqrt(pow(x1,2) + pow((exp(x1) - x2),2)));
		derivedX1X1        = derivedX1X1 - (pow((2*exp(x1)*(exp(x1) - x2) + 2*x1), 2)/(4*pow(pow(x1,2) + pow(exp(x1) - x2,2), 3.0/2)));

		double derivedX1X2 = ((2*exp(x1)*(exp(x1) - x2)+2*x1) * (exp(x1)  - x2))/(2*pow(pow(x1, 2) + pow(exp(x1) - x2, 2), 3.0/2));
		derivedX1X2        = derivedX1X2 - exp(x1)/(sqrt(pow(x1,2) + pow(exp(x1)- x2, 2)));

		double derivedX2X1 = ((2*exp(x1)*(exp(x1) - x2) + 2*x1) * (exp(x1) -x2))/(2*pow(pow(x1,2) + pow(exp(x1) -x2, 2), 3.0/2));
		derivedX2X1        = derivedX2X1 - exp(x1)/ sqrt(pow(x1,2) + pow(exp(x1) - x2, 2));
		
		
		double derivedX2X2 = 1/(sqrt(pow(x1,2) + pow(exp(x1) - x2,2)));
		derivedX2X2        = derivedX2X2 - (pow(exp(x1) - x2, 2)/(pow((pow(x1,2) + pow(exp(x1) - x2, 2)), 3.0/2)));

		MatrixXd result(2,2);

		result(0,0) = derivedX1X1;
		result(0,1) = derivedX1X2;
		result(1,0) = derivedX2X1;
		result(1,1) = derivedX2X2;

	
		return result;
	}

	static double function2 (vector<double> params)
	{
		double x1 = params[0];
		double x2 = params[1];
		double x3 = params[2];
		double y  = sqrt(x1 + x2 + x3);
		return y;
	}

	static vector<double> function2FirstDerivative (vector<double> params)
	{
		double x1 = params[0];
		double x2 = params[1];
		double x3 = params[2];

		double value = 1.0/(2.0*sqrt(x1 + x2 + x3));
		
		vector<double> result;
		result.push_back(value);
		result.push_back(value);
		result.push_back(value);

		return result;
	}

	static MatrixXd function2SecondDerivative (vector<double> params)
	{
		double x1 = params[0];
		double x2 = params[1];
		double x3 = params[2];

		double derivate = -1.0 * (1 / (4 * pow (x1 + x2 + x3, 3.0/2)));

		MatrixXd result(3,3);
		
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				result(i,j) = derivate;
			}
		}
	
		return result;
	}
};