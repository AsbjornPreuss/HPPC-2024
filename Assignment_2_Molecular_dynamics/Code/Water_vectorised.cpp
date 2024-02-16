#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cassert>
#include <math.h>
#include <chrono>

const double deg2rad = acos(-1)/180.0; // pi/180 for changing degs to radians
double accumulated_forces_bond  = 0.;     // Checksum: accumulated size of forces
double accumulated_forces_angle = 0.;     // Checksum: accumulated size of forces
double accumulated_forces_non_bond = 0.;  // Checksum: accumulated size of forces

class Vec3 {
public:
    double x, y, z;
    // initialization of vector
    Vec3(double x, double y, double z): x(x), y(y), z(z) {}
    // size of vector
    double mag() const{
        return sqrt(x*x+y*y+z*z);
    }
    Vec3 operator-(const Vec3& other) const{
        return {x - other.x, y - other.y, z - other.z};
    }
    Vec3 operator+(const Vec3& other) const{
        return {x + other.x, y + other.y, z + other.z};
    }
    Vec3 operator*(double scalar) const{
        return {scalar*x, scalar*y, scalar*z};
    }
    Vec3 operator/(double scalar) const{
        return {x/scalar, y/scalar, z/scalar};
    }
    Vec3& operator+=(const Vec3& other){
        x += other.x; y += other.y; z += other.z;
        return *this;
    }
    Vec3& operator-=(const Vec3& other){
        x -= other.x; y -= other.y; z -= other.z;
        return *this;
    }
    Vec3& operator*=(double scalar){
        x *= scalar; y *= scalar; z *= scalar;
        return *this;
    }
    Vec3& operator/=(double scalar){
        x /= scalar; y /= scalar; z /= scalar;
        return *this;
    }
};
Vec3 operator*(double scalar, const Vec3& y){
    return y*scalar;
}
Vec3 cross(const Vec3& a, const Vec3& b){
    return { a.y*b.z-a.z*b.y,
             a.z*b.x-a.x*b.z,
             a.x*b.y-a.y*b.x };
}
double dot(const Vec3& a, const Vec3& b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/* a class for the bond between two atoms U = 0.5k(r12-L0)^2 */
class Bond {
public:
    double K;    // force constant
    double L0;   // relaxed length
    int a1, a2;  // the indexes of the atoms at either end
};

/* a class for the angle between three atoms  U=0.5K(phi123-phi0)^2 */
class Angle {
public:
    double K;
    double Phi0;
    int a1, a2, a3; // the indexes of the three atoms, with a2 being the centre atom
};


// ===============================================================================
// Two new classes arranging Atoms in a Structure-of-Array data structure
// ===============================================================================

/* atom class, represent N instances of identical atoms */
class Atoms {
public:
    // The mass of the atom in (U)
    double mass;
    double ep;            // epsilon for LJ potential
    double sigma;         // Sigma, somehow the size of the atom
    double charge;        // charge of the atom (partial charge)
    std::string name;     // Name of the atom
    // the position in (nm), velocity (nm/ps) and forces (k_BT/nm) of the atom
    std::vector<Vec3> v,f;
    std::vector<double> px,py,pz;
    // constructor, takes parameters and allocates p, v and f properly to have N_identical elements
    Atoms(double mass, double ep, double sigma, double charge, std::string name, size_t N_identical) 
    : mass{mass}, ep{ep}, sigma{sigma}, charge{charge}, name{name}, 
      px(N_identical, 0), py(N_identical, 0), pz(N_identical, 0), v{N_identical, {0,0,0}}, f{N_identical, {0,0,0}}
    {}
};

/* molecule class for no_mol identical molecules */
class Molecules {
public:
    std::vector<Atoms> atoms;         // list of atoms in the N identical molecule
    std::vector<Bond> bonds;          // the bond potentials, eg for water the left and right bonds
    std::vector<Angle> angles;        // the angle potentials, for water just the single one, but keep it a list for generality
    int no_mol;
};

// ===============================================================================


/* system class */
class System {
public:
    Molecules molecules; // The molecules in the system. Changed from seq.
    double time = 0;     // current simulation time
};

class Sim_Configuration {
public:
    int steps = 10000;     // number of steps
    int no_mol = 4;        // number of molecules
    double dt = 0.0005;    // integrator time step
    int data_period = 100; // how often to save coordinate to trajectory
    std::string filename = "trajectory.txt";   // name of the output file with trajectory
    // system box size. for this code these values are only used for vmd, but in general md codes, period boundary conditions exist

    // simulation configurations: number of step, number of the molecules in the system, 
    // IO frequency, time step and file name
    Sim_Configuration(std::vector <std::string> argument){
        for (long unsigned int i = 1; i<argument.size() ; i += 2){
            std::string arg = argument.at(i);
            if(arg=="-h"){ // Write help
                std::cout << "MD -steps <number of steps> -no_mol <number of molecules>"
                          << " -fwrite <io frequency> -dt <size of timestep> -ofile <filename> \n";
                exit(0);
                break;
            } else if(arg=="-steps"){
                steps = std::stoi(argument[i+1]);
            } else if(arg=="-no_mol"){
                no_mol = std::stoi(argument[i+1]);
            } else if(arg=="-fwrite"){
                data_period = std::stoi(argument[i+1]);
            } else if(arg=="-dt"){
                dt = std::stof(argument[i+1]);
            } else if(arg=="-ofile"){
                filename = argument[i+1];
            } else{
                std::cout << "---> error: the argument type is not recognized \n";
            }
        }

        dt /= 1.57350; /// convert to ps based on having energy in k_BT, and length in nm
    }
};

// Given a bond, updates the force on all atoms correspondingly
void UpdateBondForces(System& sys){
    Molecules& molecules = sys.molecules;
    // Loops over the (2 for water) bond constraints
    for (Bond& bond : molecules.bonds){
        auto& atoms1=molecules.atoms[bond.a1];
        auto& atoms2=molecules.atoms[bond.a2];
        #pragma omp simd reduction(+:accumulated_forces_bond)
        for (int i=0; i<molecules.no_mol; i++){
	    Vec3 p1 = {atoms1.px[i],atoms1.py[i],atoms1.pz[i]};
	    Vec3 p2 = {atoms2.px[i],atoms2.py[i],atoms2.pz[i]};
            Vec3 dp  = p1-p2;
            Vec3 f   = -bond.K*(1-bond.L0/dp.mag())*dp;
            atoms1.f[i] += f;
            atoms2.f[i] -= f; 
            accumulated_forces_bond += f.mag();
        }
    }
}

// Iterates over all bonds in molecules (for water only 2: the left and right)
// And updates forces on atoms correpondingly
void UpdateAngleForces(System& sys){
    Molecules& molecules = sys.molecules;
    for (Angle& angles : molecules.angles){
        //====  angle forces  (H--O---H bonds) U_angle = 0.5*k_a(phi-phi_0)^2
        //f_H1 =  K(phi-ph0)/|H1O|*Ta
        // f_H2 =  K(phi-ph0)/|H2O|*Tc
        // f_O = -f1 - f2
        // Ta = norm(H1O x (H1O x H2O))
        // Tc = norm(H2O x (H2O x H1O))
        //=============================================================
        auto& atoms1=molecules.atoms[angles.a1];
        auto& atoms2=molecules.atoms[angles.a2];
        auto& atoms3=molecules.atoms[angles.a3];
        for (int i=0; i<molecules.no_mol; i++){
	    Vec3 p1 = {atoms1.px[i],atoms1.py[i],atoms1.pz[i]};
	    Vec3 p2 = {atoms2.px[i],atoms2.py[i],atoms2.pz[i]};
	    Vec3 p3 = {atoms3.px[i],atoms3.py[i],atoms3.pz[i]};
            Vec3 d21 = p2-p1;   // Oxygen is index 0, but atom 2 in molecule
            Vec3 d23 = p2-p3;    

            // phi = d21 dot d23 / |d21| |d23|
            double norm_d21 = d21.mag();
            double norm_d23 = d23.mag();
            double phi = acos(dot(d21, d23) / (norm_d21*norm_d23));

            // d21 cross (d21 cross d23)
            Vec3 c21_23 = cross(d21, d23);
            Vec3 Ta = cross(d21, c21_23);
            Ta /= Ta.mag();

            // d23 cross (d23 cross d21) = - d23 cross (d21 cross d23) = c21_23 cross d23
            Vec3 Tc = cross(c21_23, d23);
            Tc /= Tc.mag();

            Vec3 f1 = Ta*(angles.K*(phi-angles.Phi0)/norm_d21);
            Vec3 f3 = Tc*(angles.K*(phi-angles.Phi0)/norm_d23);

            atoms1.f[i] += f1;
            atoms2.f[i] -= f1+f3;
            atoms3.f[i] += f3;

            accumulated_forces_angle += f1.mag() + f3.mag();
        }
    }
}

// Iterates over all atoms in both molecules
// And updates forces on atoms correpondingly
void UpdateNonBondedForces(System& sys){
    /* nonbonded forces: only a force between atoms in different molecules
       The total non-bonded forces come from Lennard Jones (LJ) and coulomb interactions
       U = ep[(sigma/r)^12-(sigma/r)^6] + C*q1*q2/r */
    for (auto& atoms1 : sys.molecules.atoms)
    for (auto& atoms2 : sys.molecules.atoms){
            double ep = sqrt(atoms1.ep*atoms2.ep); // ep = sqrt(ep1*ep2)
            double sigma = 0.5*(atoms1.sigma+atoms2.sigma);  // sigma = (sigma1+sigma2)/2
	    double sigma2 = sigma*sigma;
            double q1q2 = atoms1.charge*atoms2.charge;

            double KC = 80*0.7;          // Coulomb prefactor
		for (int i = 0; i< sys.molecules.no_mol; i++){ // Note, this way of counting, while better for cache
							  // use might lead to floating point errors
		Vec3 af = {0,0,0};
		#pragma omp declare reduction(+:Vec3:omp_out += omp_in) initializer(omp_priv = omp_orig)
		#pragma omp simd reduction(+:accumulated_forces_non_bond,af)
		for (int j = i+1; j<sys.molecules.no_mol; j++){
			Vec3 p1 = {atoms1.px[i],atoms1.py[i],atoms1.pz[i]};
			Vec3 p2 = {atoms2.px[j],atoms2.py[j],atoms2.pz[j]};
            		Vec3 dp = p1-p2;
            		double r  = dp.mag();                   
            		double r2 = r*r;

            		double sir = sigma2/r2; // crossection**2 times inverse squared distance
			double sir3 = sir*sir*sir;
            		Vec3 f = ep*(12*sir3*sir3-6*sir3)*sir*dp + KC*q1q2/(r*r2)*dp; // LJ + Coulomb forces
            		af += f;
            		atoms2.f[j] -= f;

            		accumulated_forces_non_bond += f.mag();
		}
		atoms1.f[i] += af;
	}
	}
}

// integrating the system for one time step using Leapfrog symplectic integration
void Evolve(System &sys, Sim_Configuration &sc){

    // Kick velocities and zero forces for next update
    // Drift positions: Loop over molecules and atoms inside the molecules
    Molecules& molecules = sys.molecules;
    for (auto& atoms : molecules.atoms){ 		// Loop over all the different types of atoms O, H1, H2
        for (int i=0; i<molecules.no_mol; i++){
	    Vec3 p = {atoms.px[i],atoms.py[i],atoms.pz[i]}; // Load position into a vector
            atoms.v[i] += sc.dt/atoms.mass*atoms.f[i];    // Update the velocities
            atoms.f[i]  = {0,0,0};                      // set the forces zero to prepare for next potential calculation
            p += sc.dt* atoms.v[i];             // update position
	    atoms.px[i] = p.x; atoms.py[i] = p.y; atoms.pz[i] = p.z;
        }
    }

    // Update the forces on each particle based on the particles positions
    // Calculate the intermolecular forces in all molecules
    UpdateBondForces(sys);
    UpdateAngleForces(sys);
    // Calculate the intramolecular LJ and Coulomb potential forces between all molecules
    UpdateNonBondedForces(sys);

    sys.time += sc.dt; // update time
}

// Setup one water molecule
System MakeWater(int N_molecules){
    //===========================================================
    // creating water molecules at position X0,Y0,Z0. 3 atoms
    //                        H---O---H
    // The angle is 104.45 degrees and bond length is 0.09584 nm
    //===========================================================
    // mass units of dalton
    // initial velocity and force is set to zero for all the atoms by the constructor
    const double L0 = 0.09584;
    const double angle = 104.45*deg2rad;    

    //         mass    ep    sigma charge name
    Atoms Oatoms(16, 0.65,    0.31, -0.82, "O", N_molecules);  // Oxygen atom
    Atoms Hatoms1( 1, 0.18828, 0.238, 0.41, "H", N_molecules); // Hydrogen atom
    Atoms Hatoms2( 1, 0.18828, 0.238, 0.41, "H", N_molecules); // Hydrogen atom

    // bonds beetween first H-O and second H-O respectively
    std::vector<Bond> waterbonds = {
        { .K = 20000, .L0 = L0, .a1 = 0, .a2 = 1},
        { .K = 20000, .L0 = L0, .a1 = 0, .a2 = 2}
    };

    // angle between H-O-H
    std::vector<Angle> waterangles = {
        { .K = 1000, .Phi0 = angle, .a1 = 1, .a2 = 0, .a3 = 2 }
    };   

    System sys;
    for (int i = 0; i < N_molecules; i++){
        Vec3 P0{sin(2*M_PI/N_molecules*i), cos(2*M_PI/N_molecules*i), 2.0*i/N_molecules}; //Spiral configuration
        Oatoms.px[i]  = P0.x; Oatoms.py[i] = P0.y; Oatoms.pz[i] = P0.z;
	Hatoms1.px[i] = P0.x+L0*sin(angle/2); Hatoms1.py[i] = P0.y+L0*cos(angle/2); Hatoms1.pz[i] = P0.z;
	Hatoms2.px[i] = P0.x-L0*sin(angle/2); Hatoms2.py[i] = P0.y+L0*cos(angle/2); Hatoms2.pz[i] = P0.z;
    }
        std::vector<Atoms> atoms {Oatoms, Hatoms1, Hatoms2};
        sys.molecules ={atoms, waterbonds, waterangles, N_molecules};
        // Above we are passing an extra argument to sys.molecules, compared to seq.
    
    
    // Store atoms, bonds and angles in Water class and return
    return sys;
}

// Write the system configurations in the trajectory file.
void WriteOutput(System& sys, std::ofstream& file){  
    // Loop over all atoms in model one molecule at a time and write out position

    Molecules& molecules = sys.molecules; // First make an alias for the molecules in the system
    for (auto& atoms : molecules.atoms){ //Loop over all the different types of atoms
        for (int i = 0; i<molecules.no_mol; i++){ //Write the position of each different atom, that is the same type.
            file << sys.time << " " << atoms.name << " " 
                << atoms.px[i] << " " 
                << atoms.py[i] << " " 
                << atoms.pz[i] << '\n';
        }
    }
}

//======================================================================================================
//======================== Main function ===============================================================
//======================================================================================================
int main(int argc, char* argv[]){
    Sim_Configuration sc({argv, argv+argc}); // Load the system configuration from command line data
    
    System sys  = MakeWater(sc.no_mol);   // this will create a system containing sc.no_mol water molecules
    std::ofstream file(sc.filename); // open file

    WriteOutput(sys, file);    // writing the initial configuration in the trajectory file
    auto tstart = std::chrono::high_resolution_clock::now(); // start time (nano-seconds)
    
    // All forces in the simulation are commented out, to check that makewater is vectorised.
    // Molecular dynamics simulation
    for (int step = 0;step<sc.steps ; step++){

        Evolve(sys, sc); // evolving the system by one step
        if (step % sc.data_period == 0){
            //writing the configuration in the trajectory file
            WriteOutput(sys, file);
        }
    }

    auto tend = std::chrono::high_resolution_clock::now(); // end time (nano-seconds)

    std::cout <<  "Elapsed time:" << std::setw(9) << std::setprecision(4)
              << (tend - tstart).count()*1e-9 << "\n";
    std::cout <<  "Accumulated forces Bonds   : "  << std::setw(9) << std::setprecision(5) 
              << accumulated_forces_bond << "\n";
    std::cout <<  "Accumulated forces Angles  : "  << std::setw(9) << std::setprecision(5)
              << accumulated_forces_angle << "\n";
    std::cout <<  "Accumulated forces Non-bond: "  << std::setw(9) << std::setprecision(5)
              << accumulated_forces_non_bond << "\n";
}
