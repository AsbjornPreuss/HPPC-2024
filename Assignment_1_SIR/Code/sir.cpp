/* 
* THIS FILE WAS WRITTEN BY:
* ASBJOERN BONEFELD PREUSS,
* DANIEL LOMHOLT CHRISTENSEN
* ELIE CUETO
*
* DURING THEIR COURSE IN HIGH PERFORMANCE PARALLEL COMPUTING,
* AT THE UNIVERSITY OF COPENHAGEN, NIELS BOHR INSTITUTE.
* SPRING 2024
* */

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double ask_for_value(string text){
    // This function takes a text, and prompts the user to give a value. A float must be given
    cout << text;
    double value;
    cin >> value;
    return value;
}

double change_in_susceptiple_people(double& beta_factor, double& infected_people, double& susceptible_people, double& population){
    return 0 - beta_factor*infected_people*susceptible_people/population;
}

double change_in_infected_people(double& beta_factor, double& infected_people, double& susceptible_people, double& population, double& gamma_factor){
    return beta_factor*infected_people*susceptible_people/population - gamma_factor*infected_people;
}

double change_in_recovered_people(double& gamma_factor, double& infected_people){
    return gamma_factor*infected_people;
}
int main() {
    // Declare variables used in SIR model. S is susceptible_people, I is infected people, R is recovered people
    double susceptible_people = 0;
    double infected_people = 0;
    double recovered_people = 0;
    // Declare variables that will be the change in S, I and R. Work as temp variables during the simulation
    double susceptible_people_differential = 0;
    double infected_people_differential = 0;
    double recovered_people_differential = 0;

    // Declare and retrieve the beta gamma and population factors from the users,
    // as well as how long the simulation will run, and the timesteps at which it is evaluated.
    double beta_factor = ask_for_value("Please enter your beta factor in days^-1\n");
    double gamma_factor = ask_for_value("Please enter your gamma factor in days^-1\n");
    double population = ask_for_value("Please enter your population in an integer number of people\n");
    double modelled_time= ask_for_value("Please enter the amount of days you would like to simulate the infection\n");
    double dt = ask_for_value("Please enter the size of the timesteps\n");
    // Tell the user what they have chosen
    cout << "You have chosen a beta factor of " << beta_factor <<
    ", and a gamma factor of "<< gamma_factor << " and a population of " << population << "\n\n";

    // Open output file
    string filename = "sir_out.txt";
    ofstream fout(filename);
    fout << "File generated from sir.cpp. Contains: I, R, S, t\n";

    // Now the infection starts:
    double time = 0; // The model starts at day zero.
    infected_people = 1; // The initial number of infected people are 1. TODO: Allow this to be user determined
    susceptible_people = population - infected_people;


    // Run the simulation for the required time
    for (time; time<=modelled_time; time += dt){
        // Print the values
        fout << floor(infected_people) << " "
                << floor(recovered_people)<< " "
                << floor(susceptible_people) << " "
                << time << "\n";

        // Calculate the differential values
        susceptible_people_differential = change_in_susceptiple_people(beta_factor, infected_people, 
                                                    susceptible_people, population);
        infected_people_differential = change_in_infected_people(beta_factor, infected_people,
                                                    susceptible_people, population, gamma_factor);
        recovered_people_differential = change_in_recovered_people(gamma_factor, infected_people);
        // Add the differential value to SIR values
        susceptible_people += susceptible_people_differential*dt;
        infected_people += infected_people_differential*dt;
        recovered_people += recovered_people_differential*dt;

        cout << "The amount of infected people are " << floor(infected_people) << "\n";
        cout << "The amount of recovered people are " << floor(recovered_people) << "\n";
        cout << "The amount of susceptible people are " << floor(susceptible_people) << "\n";
        cout << "We have a total population of " << infected_people + recovered_people + susceptible_people << "\n";
        cout << "The time is " << time << "\n\n";   
    }
    // Print the values one last time
    fout << floor(infected_people) << " "
        << floor(recovered_people)<< " "
        << floor(susceptible_people) << " "
        << time << "\n";

    fout.close();
    return 0;
}
