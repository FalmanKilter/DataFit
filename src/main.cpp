/*	the program takes in two vectors of data - x and y - that make up points on a 2d plane,
	and returns coefficients for the chosen fit law - linear or non-linear (quadratic, cubic, etc).
	
	the coefficients are determined using gradient descent method to minimize mse of the predicted output.

	the user is prompted to input:
	x and y - input and output data vectors;
	
	fit_type - linear or non-linear fit law;
	alpha - learning rate for the gradient descent;
	tol - absolute tolerance for the solution;

	max_iter - maximum number of iterations.
*/

#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <iomanip>
#include <math.h>

using namespace std;

// inserts a blank line into the console
void insert_line()
{
	cout << endl;
}

// returns a square of the input
double square(double a)
{
	return a * a;
}

// reads a vector from the console
vector<double> read_array()
{
	vector<double> a;
	double el;

	while (cin >> el)
		a.push_back(el);

	cin.clear();
	cin.ignore(numeric_limits<streamsize>::max(), '\n');


	return a;
}

// writes string s into the console
void display_message(const string& s)
{
	cout << s << endl;
}

// displays results of theta calculation
void display_theta(const vector<double>& theta, const int &fit)
{
	// leave one blank line before output
	insert_line();

	switch (fit)
	{
	case 1:
		display_message("Linear fit coefficients for THETA_0 * X + THETA_1 are:");
		cout << "THETA_O = " << theta[0] << endl;
		cout << "THETA_1 = " << theta[1] << endl;
		break;
	case 2:
		display_message("Quadratic fit coefficients for THETA_0 * X^2 + THETA_1 * X + THETA_2 are:");
		cout << "THETA_O = " << theta[0] << endl;
		cout << "THETA_1 = " << theta[1] << endl;
		cout << "THETA_2 = " << theta[2] << endl;
		break;
	}

}

// displays goodness of fit in terms of RMSE
void display_fit_results(vector<double>& X, vector<double>& y, size_t &l,  const vector<double>& theta, int &fit)
{

	// using fancy iterators
	vector<double>::iterator x_iter = X.begin();
	vector<double>::iterator y_iter = y.begin();

	double pred = 0;

	double sum = 0;			// sum of squared errors
	double RMSE = 0;		// RMSE

	while (x_iter!=X.end())
	{
		// calculate current prediction 
		switch (fit)
		{
		case 1:
			pred = theta[0] * (*x_iter) + theta[1];
			break;
		case 2:
			pred = theta[0] * square(*x_iter) + theta[1]*(*x_iter) + theta[2];
			break;
		// ...

		default:
			pred = 0;
			break;
		}
		
		sum += square((*y_iter) - pred);

		++x_iter;
		++y_iter;
	}

	RMSE = sqrt(sum / l);

	insert_line();

	display_message("Goodness of fit (RMSE):");
	cout << RMSE << endl;

}

// makes sure vectors a and b meet certain requirements
void validate_input(vector<double>&a, vector<double> &b)
{
	vector<double>::size_type la = a.size();
	vector<double>::size_type lb = b.size();

	assert(la != 0 && lb != 0);
	assert(la == lb);
}

// makes sure the chosen fit type exists
void validate_fit_type(int& a)
{
	assert(a == 1 || a == 2);	// or some other fit types
}

// displays whether maximum iteration or tolerance level were reached
void assess_fit(int& it, double& max_it)
{
	if (it == max_it)
		display_message("Maximum iteration reached!");
	else 
		display_message("The best fit is found within the given tolerance");
}


int main()
{
	//================================================ READ FIT PARAMETERS ================================================//

	// read vector X from console
	display_message("Enter vector X: ");
	vector<double> X = read_array();

	// read vector y from console
	display_message("Enter vector y: ");
	vector<double> y = read_array();
	
	// make sure the input vectors are valid: they are not empty and have equal length
	validate_input(X, y);

	vector<double>::size_type l_data = X.size();

	// read fit type
	display_message("Choose fit type: \n 1) linear; \n 2) quadratic; \n ...");
	int fit_type;
	cin >> fit_type;

	// make sure the selected type of fit exists
	validate_fit_type(fit_type);

	// read learning rate
	display_message("Enter learning rate: ");
	double alpha;
	cin >> alpha;

	assert(alpha > 0);		// alpha can only be a positive floating point

	// read absolute tolerance
	display_message("Enter absolute tolerance: ");
	double tol;
	cin >> tol;

	assert(tol > 0);		// tol can only be a positive floating point

	// read maximum number of iterations
	display_message("Enter maximum number of iterations: ");
	double max_iter;		// type double to allow *number*e*power of 10*
	cin >> max_iter;

	assert(max_iter > 0);	// number of iterations must be a positive integer

	insert_line();
	//================================================ CALCULATE FIT ================================================//

	int it = 0;				// iteration counter	

	// generate initial guess vector based on the fit type
	vector<double> theta; 

	switch (fit_type)
	{
	case 1: 
		theta = { 1, 1 };		// linear fit
		break;
	case 2: 
		theta = { 1, 1, 1 };	// quadratic fit
		break;

	// ...
	
	}

	double h = 0;						// predicted output

	vector<double> J;				// cost function values
	double J_current;

	double temp;					// auxillary temporary variable for theta updates

	while (it != max_iter)
	{
		// reinitilize sum variables at each iteration

		double sum = 0;					// sum of squared errors

		double sum_der_lin_1 = 0;		// sum of errors used in cost function derivative 
		double sum_der_lin_2 = 0;		// for linear fit

		double sum_der_quad_1 = 0;		// sum of errors used in cost function derivative
		double sum_der_quad_2 = 0;		// for quadratic fit
		double sum_der_quad_3 = 0;		// 


		// ...							// sum of errors for other types of fit
		// calculate predicted output
		for (vector<double>::size_type i = 0; i != l_data; ++i)
		{
			switch (fit_type)
			{
			case 1:
				h = theta[0] * X[i] + theta[1];	
				sum_der_lin_1 += (h - y[i]) * X[i];
				sum_der_lin_2 += (h - y[i]);
				break;
			case 2:
				h = theta[0] * square(X[i]) + theta[1] * X[i] + theta[2];
				sum_der_quad_1 += (h - y[i]) * square(X[i]);
				sum_der_quad_2 += (h - y[i]) * X[i];
				sum_der_quad_3 += (h - y[i]);
				break;

			// ...
			default:
				h = 0;
				break;
			}

			sum += square(h - y[i]);
		}


		// append the current cost function value to the vector J
		J_current = (1/ double(l_data)) * sum;
		J.push_back(J_current);

		// check exit condition
		if (it >= 2)
		{
			if (abs((J_current - J[it - 1]) / J[it - 1]) <= tol)
				break;
		}

		// update theta
		switch (fit_type)
		{
		case 1:
			temp = theta[0] - alpha * (2 / double(l_data)) * sum_der_lin_1;
			theta[0] = temp;
			temp = theta[1] - alpha * (2 / double(l_data)) * sum_der_lin_2;
			theta[1] = temp;
			break;

		case 2:
			temp = theta[0] - alpha * (2 / double(l_data)) * sum_der_quad_1;
			theta[0] = temp;
			temp = theta[1] - alpha * (2 / double(l_data)) * sum_der_quad_2;
			theta[1] = temp;
			temp = theta[2] - alpha * (2 / double(l_data)) * sum_der_quad_3;
			theta[2] = temp;
			break;

		}

		it++;
	}

	assess_fit(it, max_iter);

	//================================================ GENERATE OUTPUT ================================================//

	display_theta(theta, fit_type);

	display_fit_results(X, y, l_data, theta, fit_type);

	return 0;
}