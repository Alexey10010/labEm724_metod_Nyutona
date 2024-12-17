#define _USE_MATH_DEFINES
#include <iostream>
#include <clocale>
#include <cmath>

using namespace std;

const double eps = 0.0001;

double f2pr(double x)
{
	double y = (pow(M_E, 2 * x) - 1) / pow(M_E, x);
	return y;
}

double f1pr(double x)
{
	double y = (pow(M_E, 2 * x) + 1) / pow(M_E, x);
	return y;
}

double f(double x)
{
	double y = pow(M_E, x) - pow(M_E, -x) - 2;
	return y;
}

double solve_equation(double(*f)(double), double x0)
{
	return f(x0);
}

int main()
{
	setlocale(LC_ALL, "Rus");
	float x0, x1, x, delta;
	int a = 0, b = 1, n = 0; // n - количество итераций
	bool iskX = 0; //solve_equation(f, a или b) * solve_equation(f2pr, a или b) == 0 -- iskomiy X 
	//найден
	if (solve_equation(f, a) * solve_equation(f2pr, a) > 0)
		x0 = a;
	if (solve_equation(f, b) * solve_equation(f2pr, b) > 0)
		x0 = b;
	if (solve_equation(f2pr, a) != 0)
	{
		if (solve_equation(f, a) * solve_equation(f2pr, a) == 0)
		{
			x = a;
			iskX = 1;
		}
	}
	if (solve_equation(f2pr, b) != 0)
	{
		if (solve_equation(f, b) * solve_equation(f2pr, b) == 0)
		{
			x = b;
			iskX = 1;
		}
	}
	if (!iskX)
	{
		do
		{
			x1 = x0 - solve_equation(f, x0) / solve_equation(f1pr, x0);
			//Вывод промежуточных результатов
			cout << "x" << n << " = " << x0 << ", f(x" << n << ") = " << solve_equation(f, x0)
				<< endl;
			delta = abs(x0 - x1);
			//cout << "Delta = " << delta << endl;
			x0 = x1;
			n = n + 1;
		} while (delta > eps);
		cout << "x" << n << " = " << x1 << ", f(x" << n << ") = " << solve_equation(f, x1)
			<< endl;
		x = x1;
		cout << "Ответ: x" << n << " = " << x << ", f(x" << n << ") = " << solve_equation(f, x)
			<< ", eps = " << eps;
	}
	else
		cout << "Ответ: x" << n << " = " << x << ", f(x" << n << ") = " << solve_equation(f, x)
		<< ", eps = " << eps;
	

	return 0;
}









