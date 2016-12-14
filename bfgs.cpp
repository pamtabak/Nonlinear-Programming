#include <iostream>

using namespace std;

class BFGS
{
public:
	BFGS()
	{

	}

	~BFGS()
	{

	}

	static vector<vector<double> > method(vector<vector<double> > hk, vector<double> pk, vector<double> qk)
	{
		vector<vector<double> > newHk;

		double denominator = vectorTransposeXvector(pk, qk);

		double firstComponent = vectorTransposeXvector(vectorTranposeXmatrix(qk, hk), qk)/denominator;
		
		vector<vector<double> > secondComponent = vectorXvectorTranspose(pk, pk);
		for (int i = 0; i < secondComponent.size(); i++)
		{
			for (int j = 0; j < secondComponent[0].size(); j++)
			{
				secondComponent[i][j] = (1 + firstComponent)*(secondComponent[i][j]/denominator);
			}
		}

		vector<vector<double> > thirdComponent1 = matrixXmatrix(vectorXvectorTranspose(pk, qk), hk);
		vector<vector<double> > thirdComponent2 = vectorXvectorTranspose(matrixXvector(hk, qk), pk);
		vector<vector<double> > thirdComponent;
		for (int i = 0; i < thirdComponent1.size(); i++)
		{
			vector<double> row;
			for (int j = 0; j < thirdComponent1[0].size(); j++)
			{
				row.push_back((thirdComponent1[i][j] + thirdComponent2[i][j])/denominator);
			}
			thirdComponent.push_back(row);
		}

		for (int i = 0; i < hk.size(); i++)
		{
			vector<double> row;
			for (int j = 0; j < hk[0].size(); j++)
			{
				row.push_back(hk[i][j] + secondComponent[i][j] - thirdComponent[i][j]);
			}
			newHk.push_back(row);
		}

		return newHk;
	}

	static double vectorTransposeXvector(vector<double> vector1, vector<double> vector2)
	{
		double result = 0.0;
		for (int i = 0; i < vector1.size(); i++)
		{
			result += vector1[i]*vector2[i];
		}
		return result;
	}

	static vector<vector<double> > vectorXvectorTranspose(vector<double> vector1, vector<double> vector2)
	{
		vector<vector<double> > result;
		for (int i = 0; i < vector1.size(); i++)
		{
			vector<double> row;
			for (int j = 0; j < vector2.size(); j++)
			{
				row.push_back(vector1[i]*vector2[j]);
			}
			result.push_back(row);
		}
		return result;
	}

	static vector<double> vectorTranposeXmatrix(vector<double> vector1, vector<vector<double> > matrix)
	{
		vector<double> result;
		for (int j = 0; j < vector1.size(); j++)
		{
			double value = 0.0;
			for (int i = 0; i < matrix[0].size(); i++)
			{
				value += matrix[i][j]*vector1[i];
			}
			result.push_back(value);
		}
		return result;
	}

	static vector<vector<double> > matrixXmatrix (vector<vector<double> > matrix1, vector<vector<double> > matrix2)
	{
		vector<vector<double> > result;
		for (int i = 0; i < matrix1.size(); i++)
		{
			vector<double> row;
			for (int j = 0; j < matrix2[0].size(); j++)
			{
				double value = 0.0;
				for (int k = 0; k < matrix1[0].size(); k++)
				{
					value += matrix1[i][k] * matrix2[k][j];
				}
				row.push_back(value);
			}
			result.push_back(row);
		}
		return result;
	}

	static vector<double> matrixXvector (vector<vector<double> > matrix, vector<double> vector1)
	{
		vector<double> result;
		for (int i = 0; i < matrix.size(); i++)
		{
			double value = 0.0;
			for (int j = 0; j < matrix[0].size(); j++)
			{
				value += matrix[i][j]*vector1[j];
			}
			result.push_back(value);
		}
		return result;
	}
};