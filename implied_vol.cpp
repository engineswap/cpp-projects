#include <iostream>
#include <cmath>
#include <vector>

/*
	This file calculates implied volatility of an option using black scholes + Newton-Raphson method
	Doesnt account for dividends
	Run command: g++ -std=c++11 implied_vol.cpp -o iv && ./iv
*/

const double e = 2.7182818284;
const double tolerance = 0.001;

// Option specifics
double S = 369.67;
double K = 290; 
double T = 1.01917808219;
double r = 0.0482;
double price = 110;

double normal_cdf(double value){
	return 0.5 * erfc(-value * M_SQRT1_2);
}

// Price of a call and put option using the black scholes model
std::vector<double> price_option(double S, double K, double T, double r, double vol){
	const double d1 = (log(S/K) + (r + pow(vol, 2)/2) * T) / (vol * pow(T, 0.5));
	const double d2 = d1 - vol*(pow(T,0.5));

	double call_price = (S * normal_cdf(d1)) - (K * pow(e, (-1)*r*T) * normal_cdf(d2));
	double put_price = K * pow(e, (-1)*r*T) * normal_cdf((-1)*d2) - (S * normal_cdf((-1)*d1));

	std::vector<double> res = {call_price,put_price};
	return res;
}

// Estimated price - actual price
double f(double vol_){
	std::vector<double> option_prices_predicted = price_option(S,K,T,r,vol_);
	return option_prices_predicted[0] - price;
}

// Finite difference derivative
double df(double x){
	double dx = 0.0000001;
	return (f(x+dx) - f(x)) / dx;
}

// x0 is initial guess
double newton(double x0){
	for(int i=0; i<100; i++){
		// std::cout << "Guess " << i+1 << ": " << x0 << "\n";
		// std::cout << x0 << "\n" << f(x0) << "\n" << df(x0) << "\n\n";
		double tmp_x0 = x0 - f(x0) / df(x0);
		if (std::abs(tmp_x0-x0) <= tolerance){
			break;
		}
		x0 = tmp_x0;
	}
	return x0;
}

int main(){
	std::cout << "Call price: $" << price << "\n";


	std::cout << "Estimated implied vol: " << newton(0.47) << "\n";
	return 0;
}