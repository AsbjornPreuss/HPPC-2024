#include <iostream>
using namespace std;

double ask_for_value(string text){
    cout << text;
    double value;
    cin >> value;
    return value;
}

int main() {
    double beta_factor = ask_for_value("Please enter your beta factor\n");
    double gamma_factor = ask_for_value("Please enter your gamma factor\n");
    double population = ask_for_value("Please enter your population factor\n");
    
    cout << "You have chosen a beta factor of " << beta_factor <<
    ", and a gamma factor of "<< gamma_factor << " and a population of " << population << "\n\n";



    return 0;
}