#include <iostream>

using namespace std;

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
	
		double derivedX1 = (exp(2*x1) - x2*exp(x1) + x1) / sqrt(pow(x1, 2) + pow((exp(x1) - x2), 2));
		double derivedX2 = (x2 - exp(x1)) / sqrt(pow(x1, 2) + pow((exp(x1) - x2), 2));
	
		vector<double> result;
		result.push_back(derivedX1);
		result.push_back(derivedX2);
	
		return result;
	}


	// TO DO!!!!!!!
	static vector<double> function1SecondDerivate (vector<double> params)
	{
		double x1 = params[0];
		double x2 = params[1];

		// TO DO!!!!!!!!
		double derivedX1 = 0.0;
		double derivedX2 = 0.0;

		vector<double> result;
		result.push_back(derivedX1);
		result.push_back(derivedX2);
	
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

		double derived = 1 / (2 * sqrt(x1 + x2 + x3));
		
		vector<double> result;
		result.push_back(derived);
		result.push_back(derived);
		result.push_back(derived);

		return result;
	}

	static vector<double> function2SecondDerivate (vector<double> params)
	{
		double x1 = params[0];
		double x2 = params[1];
		double x3 = params[2];

		double derived = -1.0 * (1 / (4 * pow (x1 + x2 + x3, 3.0/2)));

		vector<double> result;
		result.push_back(derived);
		result.push_back(derived);
		result.push_back(derived);
	
		return result;
	}
};