#include <iostream>
#include <cmath>
#include <vector>

/*
	Creates a vol surface using black scholes model
	Run command: g++ -std=c++11 price_option.cpp -o price_option && ./price_option
*/

const double e  = 2.7182818284;

double normal_cdf(double value){
	return 0.5 * erfc(-value * M_SQRT1_2);
}

std::vector<double> price_option(double S, double K, double T, double r, double vol){
	const double d1 = (log(S/K) + (r + pow(vol, 2)/2) * T) / (vol * pow(T, 0.5));
	//std::cout << "d1: " << d1 << "\n";
	const double d2 = d1 - vol*(pow(T,0.5));
	//std::cout << "d2: " << d2 << "\n";

	double call_price = (S * normal_cdf(d1)) - (K * pow(e, (-1)*r*T) * normal_cdf(d2));
	double put_price = K * pow(e, (-1)*r*T) * normal_cdf((-1)*d2) - (S * normal_cdf((-1)*d1));

	std::vector<double> res = {call_price,put_price};
	return res;
}

int main(){
	std::cout << "Please enter the stock price (S): $";
	double S;
	std::cin >> S;
	
	std::cout << "Please enter the strike price (K): $";
	double K;
	std::cin >> K;

	std::cout << "Please enter the time until expiration time in years (T): ";
	double T;
	std::cin >> T;

	std::cout << "Please enter the annualized risk-free rate as a decimal (r): ";
	double r;
	std::cin >> r;

	std::cout << "Please enter the volatility of the stock (vol): ";
	double vol;
	std::cin >> vol;
	
	std::vector<double> prices = price_option(S,K,T,r,vol);
	double call_price = prices[0];
	double put_price = prices[1];
	std::cout << "Call price: $" << call_price << "\n";
	std::cout << "Put price: $" << put_price << "\n";


	return 0;
}