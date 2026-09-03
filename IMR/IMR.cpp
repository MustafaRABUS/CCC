// Ions in Matter – Bethe CSDA Range (IMR)
//
// Author:"Mustafa RABUS";
//
// This program requires the user to manually enter certain input parameters through the terminal.
// IMR is a nuclear-interaction simulation code that calculates the energy
// deposited by charged particles as they traverse matter and estimates their
// range in the target material following energetic charged-particle bombardment.
//
// The input data are read from a CSV file, and the simulation results,
// including the calculated particle ranges, are written to CSV files in the "RangeOut" folder.
//
// Videos demonstrating the simulation of three nuclear reactions are
// available on the following channel:https://www.youtube.com/@MstfRbs/videos


#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
using namespace std;

//CONSTANTS:

double c=299792458;//(m/s)light speed.
double me=0.51099895;//(MeV)electron mass energy equivalent.
double Na=6.022140*pow(10,26);//(nucleon/kg)avagadro number.
double re=2.817940*pow(10,-15);//(meter)classical electron radius.
double pi=3.1416;
double Range=0;

//PARAMETERS:
double dE=0.01;//Step Energy (MeV)
double E;
int zi;
int Ai;
double mi;//(MeV) projectile particle mass energy value.(from nist_raw_data.txt)
double A;//(amu) Target Material Relative Atomic Mass. (from PubChemElements.csv)
double d;//(kg/m3) Target material Density. (from PubChemElements.csv)
int Z;//Target material atomic number.(from PubChemElements.csv)
string Target;
string IonName;
string InterAction;
string Outputs="RangeOut/";
string Elements="PubChemElements.csv";
string Ions="nist_atomic_data.txt";
bool tryCalculate(const string& z_str, const string& a_str, const string& mass_str, int zi, int Ai, double& mi) {
    if (z_str.empty() || a_str.empty() || mass_str.empty()) return false;
    if (stoi(z_str) == zi && stoi(a_str) == Ai) {
        string pure_mass = mass_str.substr(0, mass_str.find('('));
        mi = 931.49410242 * stod(pure_mass)-(zi*me)+(14.4381 * pow(zi, 2.39) + 1.55468 * pow(10, -6) * pow(zi, 5.35)) / 1000000.0;
        return true;
    }
    return false;
}

//FUNCTIONS:
int main (){
    cout<<"Input Ion Charge/Atomic Number = ";
    cin >> zi ;
    cout<<"Input Ion Mass Number = ";
    cin >> Ai;
    cout<<"Input Ion Initial Energy (MeV) = ";
    cin >> E ;
    cout<<"Input Target Periodic Table Symbol = ";
    cin >> Target ;

    bool found = false;
    string lines1, key, value, cur_zi, cur_Ai, cur_mass, cur_symbol;
    stringstream ss;
    ifstream file1(Ions);
    if (!file1.is_open()) {cerr << "Can't open file\n";return 1;}
    while (getline(file1 >> ws, lines1)) {
        ss.clear();
        ss.str(lines1);

        if (getline(ss, key, '=') && getline(ss >> ws, value)) {
            value = value.substr(0, value.find_last_not_of(" \t\r\n") + 1);
            key = key.substr(0, key.find_last_not_of(" \t\r\n") + 1);
            if (key == "Atomic Number") cur_zi = value;
            else if (key == "Atomic Symbol") cur_symbol = value;
            else if (key == "Mass Number") cur_Ai = value;
            else if (key == "Relative Atomic Mass") {
                cur_mass = value;
                if (tryCalculate(cur_zi, cur_Ai, cur_mass, zi, Ai, mi)) {IonName = cur_symbol;found = true;break;}
            }
        }
    }
    file1.close();

    int kamma=0;
    int search=0;
    string lines2;
    ifstream file2(Elements);
    if (!file2.is_open()) {cerr << "Can't open file\n"; return 1;}
    while (getline(file2, lines2)) {
        stringstream ss(lines2);
        string item;
        while (getline(ss,item,',')) {
            if (item==Target){search=1;}
            if (search==1){kamma=kamma+1;}
            if (search==1 and kamma==3){Z=stoi(item);}
            if (search==1 and kamma==4){A=stod(item);}
            if (search==1 and kamma==14){d=1000*stod(item);search=0;}
            }
        }
    file2.close();

    //CALCULATE PENETRATION:

    double Ei=E;
    double Nv=d*Na/A;//(atom/m3) the number of atoms per unit volume of material penetrated by the particle.
    double Ih=(9.76+58.8*pow(Z,(-1.19)))*Z*pow(10,(-6));//(MeV) Z>=13 Target material ionization potential.
    double Il=(12*Z+7)*pow(10,(-6));//(MeV) Z<13 Target material ionization potential.
    double bethe=0;
    double gamma1=0;
    double gamma2=0;
    double beta1=0;
    double beta2=0;
    double dEdx=0;
    double dx=0; //Step (m).
    double R=0; //Range (m).

    vector<double>depth;
    vector<double>energy;
    InterAction=IonName+" in "+Target;
    ofstream file3(Outputs+InterAction+".csv");
    if (!file3.is_open()) {cerr << "Can't open file\n";return 1;}
    file3 << InterAction <<endl<<"ProjectileMass(MeV):,"<<mi<<endl<<"TargetAtomicNumber:,"<< Z <<endl<< "TargetAtomicMass(amu):,"<< A <<endl<< "TargetDensity(Kg/m3):," << d <<endl <<"\n";
    file3 << "ENERGY(MeV),RANGE(m)\n";
    for (int i=1;Ei>=0;i++) {
        gamma1=(Ei+ mi)/mi;
        gamma2=pow(((Ei+ mi)/mi),2);
        beta2=1-(1/gamma2);
        beta1=sqrt(beta2);
        bethe=(4*pi*pow(re,2)*pow(zi,2)*me*Nv*Z/beta2);
        if (Z>=13){
            if (mi > 0.511){dEdx=bethe*(log(2*me*beta2*gamma2/Ih)-beta2);}
            else if (zi<0) {dEdx=bethe*(log(beta1*gamma1*sqrt(gamma1-1)*me/Ih)+((pow(gamma1-1,2)/8)+1-log(2)*(gamma2+2*gamma1-1))/(2*gamma2));}
            else {dEdx=bethe*(log(beta1*gamma1*sqrt(gamma1-1)*me/Ih)-(23+14/(gamma1+1)+10/pow(gamma1+1,2)+4/pow(gamma1+1,3)*(beta2/24)+log(2)/2));}
            }
        else if (Z<13){
            if (mi > 0.511){dEdx=bethe*(log(2*me*beta2*gamma2/Il)-beta2);}
            else if (zi<0) {dEdx=bethe*(log(beta1*gamma1*sqrt(gamma1-1)*me/Il)+((pow(gamma1-1,2)/8)+1-log(2)*(gamma2+2*gamma1-1))/(2*gamma2));}
            else {dEdx=bethe*(log(beta1*gamma1*sqrt(gamma1-1)*me/Il)-(23+14/(gamma1+1)+10/pow(gamma1+1,2)+4/pow(gamma1+1,3)*(beta2/24)+log(2)/2));}
            }
        dx=(1/dEdx)*dE;
        depth.push_back(R);
        if (dx>=0) R=R+dx;
        energy.push_back(Ei);
        Ei=Ei-dE;
        }
    Range=R;
    for (int i=0;i<energy.size();i++) file3 << energy[i] << "," << (Range-depth[i]) << "\n";
    file3.close();
    depth.clear();
    energy.clear();
    cout<<InterAction<<endl<<mi<<endl<<Z<<endl<<A<<endl<<d<<endl<<R<<endl;
    return 0;
}
