/* 
* THIS FILE WAS WRITTEN BY:
* ASBJOERN BONEFELD PREUSS,
* DANIEL LOMHOLT CHRISTENSEN
* ELIE CUETO
*
* DURING THEIR COURSE IN HIGH PERFORMANCE PARALLAL COMPUTING,
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
    double susceptible_people;
    double infected_people;
    double recovered_people;
    // Declare variables that will be the change in S, I and R. Work as temp variables during the simulation
    double susceptible_people_differential;
    double infected_people_differential;
    double recovered_people_differential;

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

    // Declare a Class that will contain S, I and R
    
    double SIR_output[int(modelled_time/dt)][3];


    // Now the infection starts:
    double time = 0; // The model starts at day zero.
    infected_people = 1; // The initial number of infected people are 1. TODO: Allow this to be user determined
    susceptible_people = population - infected_people;

    // Run the simulation for the required time
    for (time; time<=modelled_time; time += dt){
        // Calculate the differential values
        susceptible_people_differential = change_in_susceptiple_people(beta_factor, infected_people, 
                                                    susceptible_people, population);
        infected_people_differential = change_in_infected_people(beta_factor, infected_people,
                                                    susceptible_people, population, gamma_factor);
        recovered_people_differential = change_in_recovered_people(gamma_factor, infected_people);
        
        // Add the differential value to SIR values
        susceptible_people += susceptible_people_differential;
        infected_people += infected_people_differential;
        recovered_people += recovered_people_differential;

        // Print the values
        cout << "The amount of infected people are " << floor(infected_people) << "\n";
        cout << "The amount of recovered people are " << floor(recovered_people) << "\n";
        cout << "The amount of susceptible people are " << floor(susceptible_people) << "\n";
        cout << "We have a total population of " << infected_people + recovered_people + susceptible_people << "\n";
        cout << "The time is " << time << "\n\n";   

        // Add the values to the SIR_output object
        SIR_output[int(dt*time)][0] = susceptible_people;
        SIR_output[int(dt*time)][1] = infected_people;
        SIR_output[int(dt*time)][2] = recovered_people;
    }
    // Now we write the result to file
    string filename = "sir_out.txt";
    ofstream fout(filename);
    fout << "File generated from sir.cpp. Contains: S, I, R\n";
    for (int i = 0; i <int(modelled_time/dt); i++){
        fout << SIR_output[i][0]<< " " <<
          SIR_output[i][1]<< "  "<<
          SIR_output[i][2] << "\n";
    }
 
    return 0;
}
