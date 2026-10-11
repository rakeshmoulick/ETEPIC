/*
 spic_v2: particle in cell for space plasma 
 1D-1V PIC Plasma Code with PBC:
 This is the second version of the previous spic.   
 This is a Normalized Code with two electron species
 The normalization scheme has been changed in this code. 
 The spatial lengths are normalized by LD=sqrt(eps0*Tec/(n0*e^2))
 Tec = cold electron temperature and n0 = Total plasma density.
 All densities has been normalized by the total electron plasma density
 and frequencies by wp.
 ****************************************************
 Developer:
 DR. RAKESH MOULICK, CPP-IPR, ASSAM, INDIA
 ****************************************************
 */

# include <iostream>
# include <cmath>
# include <cstdlib>
# include <vector>
# include <list>
# include <ctime>
# include <random>
# include <cstring>
# include <fstream>
# include <sstream>
# include <chrono>
extern "C"
{
	# include "iniparser.h"

}
using namespace std;

/***************** INIPARSER ***************/
int  parse_ini_file(char * ini_name);

/* Random Number Generator */
std::mt19937 mt_gen(0);
std::uniform_real_distribution<double> rnd_dist(0,1.0);
double rnd()
{
    return rnd_dist(mt_gen);
}

/* Define universal constants */
double EPS;
double K;
double massE;
double chargeE;
double eV;
double e;
double AMU;
double massI;
double EV_TO_K;
double pi;

/* Define Simulation Parameters*/
double density;       // Plasma Density
double stepSize;      // Cell Spacing
double DT;            // Time steps
double DT_coeff;

double tempEcold;     // Temperature of the cold electrons in eV units
double tempEhot;      // Temperature of the hot electrons in eV units
double tempEbeam;     // Temperature of the Beam electrons in eV units
double tempI;  		  // Temperature of the ion species in eV


/*Simulation Normalizing Parameters*/
double n0, nec0, neh0, neb0, ni0;
double LD, LDH, LDC;
double wp, wpec, wpeh;
double CS;
int vd; // Multiple of electron thermal velocity (for beam drift)

/* CHANGED TYPE FROM CONST TO VAR FOR INPUT DATA CONTROL  */
int nParticlesE;     // Number of simulation electrons
int nParticlesI;     // Number of simulation ions
double alp, beta;			 // The fraction of cold electrons
int NC;              // Total number of cells
int NUM_TS;          // Total Time steps (default)
int write_interval;

string output; // Open up the output file

/* Class Domain: Hold the domain paramassEters*/
class Domain
{
public:
    int ni;      // Number of nodes
    double x0;   // initial position
    double dx;   // cell spacing
    double xl;   // domain length
    double xmax; // domain maximum position

    /* Field Data structures */
    double *phi; // Electric Potential
    double *ef;  // Electric field
    double *rho; // Charge Density
};

/* Class Particle: Hold particle position, velocity and particle identity*/
class Particle
{
public:
    double pos;  // particle position
    double vel; // particle velocity
    int id;  // hold particle identity

    // Add a constructor
    Particle(double x, double v):pos(x), vel(v){};
};

/* Class Species: Hold species data*/
class Species
{
public:
	int nc = 0;
    // Use linked list for the particles
    list<Particle> part_list;
    double mass;
    double charge;
    double spwt;
    string name;

    int NUM;
    double Temp;
    double *den;
    double *vel;

    void add(Particle part)
    {
        part.id=part_id++;
        part_list.push_back(part);
    }

    // Add a constructor
    Species(string name, double mass, double charge, double spwt, int NUM, double Temp)
    {
        setName(name);
        setMass(mass);
        setCharge(charge);
        setSpwt(spwt);
        setNum(NUM);
        setTemp(Temp);
    }

    // Define the constructor functions
    void setName(string name){this->name = name;}
    void setMass(double mass){this->mass = mass;}
    void setCharge(double charge){this->charge = charge;}
    void setSpwt(double spwt){this->spwt = spwt;}
    void setNum (int NUM){this->NUM = NUM;}
    void setTemp(double Temp){this->Temp = Temp;}

private:
    int part_id = 0;
};

// Define Domain and File as the global variable
Domain domain;

// Open files to hold numerical data
FILE *file_res; // File to hold the data of overall data
FILE *file_ke; 	// File to hold the Kinetic energy data
FILE *file_pe; // File to hold the Potential energy data
FILE *f1; 	// File to hold the ion particle data
FILE *f2; 	// File to hold the cold electron particle data
FILE *f3; 	// File to hold the hot electron particle data
FILE *f4;	// File to hold the beam electron particle data
FILE *file_loc; // File to hold data at a specific location of simulation domain

FILE *file_denflucI; 	// File to hold the density fluctuations for ions
FILE *file_denflucEC;	// File to hold the density fluctuations for cold electrons
FILE *file_denflucEH;	// File to hold the density fluctuations for hot electrons

// Define Helper functions
void Init(Species *species, string flag);
void ScatterSpecies(Species *species);
void ScatterSpeciesVel(Species *species);
void ComputeRho(Species *ions, Species *electrons_cold, Species *electrons_hot, Species *electrons_beam);
void ComputeEF(double *phi, double *ef);
void PushSpecies(Species *species, double *ef);
void RewindSpecies(Species *species, double *ef);
void SanityCheck(double new_pos, double old_pos, double cell_len);

// [Write Outputs]
void Write_ts(int ts, Species *ions,Species *electrons_cold, Species *electrons_hot, Species *electrons_beam);
void Write_Particle(FILE *file, int ts, Species *species);
void WriteKE(double Time, Species *ions, Species *electrons_cold, Species *electrons_hot, Species *electrons_beam);
void WritePE(double Time, Species *electrons_cold);
void WriteLocation(double Time, double pos);
void WriteDenFluc(FILE *file, double Time, double pos, Species *species);

double ComputeKE(Species *species, Species *electrons_cold);
double ComputePE(Species *electrons_cold);
double XtoL(double pos);
double gather(double lc, double *field);

// [Sample Velocities]
double SampleVel(double T, double mass);
double SampleVelCold(double T, double mass);
double SampleVelHot(double T, double mass);
double SampleVelBeam(double T, double mass);

//[ Potential Solvers]
bool SolvePotential(double *phi, double *rho);
bool SolvePotentialDirect(double *phi, double *rho);

// [Ini Parser File]
int parse_ini_file(char * ini_name)
{
    dictionary  *   ini;

    ini = iniparser_load(ini_name);
    if (ini==NULL) {
        fprintf(stderr, "cannot parse file: %s\n", ini_name);
        return -1 ;
    }
    //iniparser_dump(ini, stderr); // Comment out to fix issues with iniparser

	/*Universal Constants*/
	EPS = iniparser_getdouble(ini,"constants:EPS",-1.0);
	K   = iniparser_getdouble(ini,"constants:K",-1.0);
	eV  = iniparser_getdouble(ini,"constants:eV",-1.0);
	e   = iniparser_getdouble(ini,"constants:e",-1.0);
	AMU = iniparser_getdouble(ini,"constants:AMU",-1.0);
	EV_TO_K  = iniparser_getdouble(ini,"constants:EV_TO_K",-1.0);
	pi  = iniparser_getdouble(ini,"constants:pi",-1.0);

    /* SPECIES INFO */
    nParticlesI    = iniparser_getint(ini,"population:nParticlesI",-1);
    nParticlesE    = iniparser_getint(ini,"population:nParticlesE",-1);
    massI          = iniparser_getdouble(ini,"population:massI",-1.0);
    massI		   = massI*AMU;
    massE          = iniparser_getdouble(ini,"population:massE",-1.0);

    chargeE        = iniparser_getdouble(ini,"population:chargeE",-1.0);
    density        = iniparser_getdouble(ini,"population:density",-1.0);

    tempEcold	  = iniparser_getdouble(ini,"population:tempEcold",-1.0);
    tempEhot 	  = iniparser_getdouble(ini,"population:tempEhot",-1.0);
    tempEbeam     = iniparser_getdouble(ini,"population:tempEbeam",-1.0);
    tempI		  = iniparser_getdouble(ini,"population:tempI",-1.0);
    alp 		  = iniparser_getdouble(ini,"population:alp",-1.0);
    beta		  = iniparser_getdouble(ini,"population:beta",-1.0);
    vd 			  = iniparser_getdouble(ini,"population:vd",-1);
    output		  = iniparser_getstring(ini,"file:output",NULL);

    /* Normalizing Parameters */
    n0 = density;
    nec0 = n0/(1+alp+beta);
    neh0 = alp*nec0;
    neb0 = beta*nec0;
    ni0 = n0; 

    LDC = sqrt((EPS*tempEcold*eV)/(nec0*chargeE*chargeE)); //cold electron Debye length
    LDH = sqrt((EPS*tempEhot*eV)/(neh0*chargeE*chargeE)); // Hot electron Debye length
    LD  = sqrt((EPS*tempEcold*eV)/(n0*chargeE*chargeE)); // Characteristic Debye Length (provides smallest spatial resolution)
    
    wp = sqrt((n0*chargeE*chargeE)/(massE*EPS)); // Total Electron Plasma Frequency
    wpec = sqrt((nec0*chargeE*chargeE)/(massE*EPS)); // cold electron plasma frquency
    wpeh = sqrt((neh0*chargeE*chargeE)/(massE*EPS)); // hot electron plasma frquency
    CS = sqrt(tempEcold*eV/massI); // Ion acoustic speed

    /*Get Simulation Parameters */
    DT_coeff = iniparser_getdouble(ini,"diagnostics:DT_coeff",-1.0);
    NUM_TS     = iniparser_getint(ini,"time:NUM_TS",-1);
    DT		   = DT_coeff*(1.0/wp); // This is the unnormalized time interval
    stepSize   = LD;
    NC         = iniparser_getint(ini,"grid:NC",-1);

    /* DIAGNOSTICS */
    write_interval = iniparser_getint(ini,"diagnostics:write_interval",-1);
   	
    //cout << massI << '\t' << massE << endl;
    cout << "*************** Input Sanity Check ***************" << '\n';
    bool SFLAG = true;
    if (stepSize >= 10) {
      cout<<"ERROR: stepSize is bigger than Debye length."<<endl;
      SFLAG = false;
    }
    if (DT > 0.01) {
      cout<<"ERROR: timeStep is too big. The recommended value: <"<<(0.01/wp)<<" s"<<endl;
      SFLAG = false;
    }

    if (SFLAG==true) {
      cout<<"STATUS: Input parameters are compatible."<<endl;
    }
    else {
      cout<<"ERROR: Input parameters are incompatible."<<endl;
      exit (EXIT_FAILURE);
    }


    iniparser_freedict(ini);
    return 0;
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/********************* MAIN FUNCTION ***************************/
int main(int argc, char *argv[])
{
    /************* INIPARSER ********************/
    if(argc<2) {
      cout<<"ERROR, at least one argument expected (the input file)."<<endl;
      exit (EXIT_FAILURE);
    }
    parse_ini_file(argv[1]);
    /********************************************/
    cout << "Debye Length: " << LD << '\t' << "Plasma Frequency: " << wp << '\t' << "Ion Acoustic Speed: " << CS << endl;
    cout << "DX: " << stepSize << '\t' << "DT: " << DT << endl;
    double Time = 0;
    /*Construct the domain paramassEters*/
    domain.ni = NC+1;
    DT = wp*DT; // This is the normalized time interval
    domain.dx = stepSize/LD;
    domain.x0 = 0;
    domain.xl = (domain.ni-1)*domain.dx;
    domain.xmax = domain.x0 + domain.xl;

    /*Allocate massEmory to the domain data structures (Field variables)*/
    domain.phi = new double[domain.ni];
    domain.ef = new double[domain.ni];
    domain.rho = new double[domain.ni];

    /*Redifine the field variables */
    double *phi = domain.phi;
    double *ef = domain.ef;
    double *rho = domain.rho;

    /* Clear the domain fields*/
    memset(phi,0,sizeof(double)*domain.ni);
    memset(ef, 0,sizeof(double)*domain.ni);
    memset(rho,0,sizeof(double)*domain.ni);

    /**************************************************/

    /*Species Info: Create vector to hold the data*/
    vector <Species> species_list;

    /*Calculate the specific weights of the ions and electrons*/
    double ion_spwt = (ni0*domain.xl*LD)/(nParticlesI);    // Normalized with LD
    double electron_cold_spwt = (nec0*domain.xl*LD)/(nParticlesE); // Normalized with LD
    double electron_hot_spwt = (neh0*domain.xl*LD)/(nParticlesE); // Normalized with LD
    double electron_beam_spwt = (neb0*domain.xl*LD)/(nParticlesE); // Normalized with LD

    /* Add singly charged Positive ions, electrons and Background Neutrals */
    /************************************************************************/
    /* Create the species lists*/
    species_list.emplace_back("Argon Ion",massI,chargeE,ion_spwt, nParticlesI, tempI);
    species_list.emplace_back("Cold Electrons",massE,-chargeE,electron_cold_spwt, nParticlesE, tempEcold);
    species_list.emplace_back("Hot Electrons",massE,-chargeE,electron_hot_spwt, nParticlesE, tempEhot);
    species_list.emplace_back("Beam Electrons",massE,-chargeE,electron_beam_spwt, nParticlesE, tempEbeam);

    /*Assign the species list as ions and electrons*/
    Species &ions = species_list[0];
    Species &electrons_cold = species_list[1];
    Species &electrons_hot = species_list[2];
    Species &electrons_beam = species_list[3];

    /*Initiate the species density and velocity fields*/
    ions.den = new double[domain.ni];
    electrons_cold.den = new double[domain.ni];
    electrons_hot.den = new double[domain.ni];
    electrons_beam.den = new double[domain.ni];

    ions.vel = new double[domain.ni];
    electrons_cold.vel = new double[domain.ni];
    electrons_hot.vel = new double[domain.ni];
    electrons_beam.vel = new double[domain.ni];

    /*Initialize electrons and ions */
    Init(&ions,"ion");
    Init(&electrons_cold,"cold");
    Init(&electrons_hot,"hot");
    Init(&electrons_beam,"beam");

    for(auto &p:species_list)
        cout<< p.name << '\n' << p.mass<< '\n' << p.charge << '\n' << p.spwt << '\n' << p.NUM << endl;
    /***************************************************************************/

    /*Compute Number Density*/
    ScatterSpecies(&ions);
    ScatterSpecies(&electrons_cold);
    ScatterSpecies(&electrons_hot);
    ScatterSpecies(&electrons_beam);

    /*Compute charge density, solve for potential
     and compute the electric field*/

    ComputeRho(&ions, &electrons_cold, &electrons_hot, &electrons_beam);
    SolvePotential(phi, rho);
    ComputeEF(phi,ef);

    RewindSpecies(&ions,ef);
    RewindSpecies(&electrons_cold,ef);
    RewindSpecies(&electrons_hot,ef);
    RewindSpecies(&electrons_beam,ef);

    /*------------- Print Output ---------------*/

    /*create a folder named output and
	delete the previous output folder: print statement is just to show*/
    printf("Deleting the output folder ... \n");
    system(("rm -rf "+ output).c_str());

    /*create an output folder*/
    //system("mkdir output");
    system(("mkdir "+ output).c_str());
    system(("mkdir "+ output + "/files").c_str());

    char NAME[50];

    sprintf(NAME,"%s/files/results_%d.txt",output.c_str(),vd);
    file_res = fopen(NAME,"w");

    sprintf(NAME,"%s/files/ke_%d.txt",output.c_str(),vd);
    file_ke = fopen(NAME,"w");
    
    sprintf(NAME,"%s/files/pe_%d.txt",output.c_str(),vd);
    file_pe = fopen(NAME,"w");

    sprintf(NAME,"%s/files/potloc_%d.txt",output.c_str(),vd);
    //file_loc = fopen(NAME,"w");

    sprintf(NAME,"%s/files/denlocI_%d.txt",output.c_str(),vd);
    //file_denflucI = fopen(NAME,"w");
	
    sprintf(NAME,"%s/files/denlocEC_%d.txt",output.c_str(),vd);
    //file_denflucEC = fopen(NAME,"w");

    sprintf(NAME,"%s/files/denlocEH_%d.txt",output.c_str(),vd);
    //file_denflucEH = fopen(NAME,"w");

    /*MAIN LOOP*/
    clock_t start = clock();
    for (int ts=0; ts<NUM_TS+1; ts++)
    {
        //Compute number density
        ScatterSpecies(&ions);
        ScatterSpecies(&electrons_cold);
 	    ScatterSpecies(&electrons_hot);
	    ScatterSpecies(&electrons_beam);

        //Compute velocities
        ScatterSpeciesVel(&ions);
        ScatterSpeciesVel(&electrons_cold);
	    ScatterSpeciesVel(&electrons_hot);
	    ScatterSpeciesVel(&electrons_beam);

        //Compute charge density
        ComputeRho(&ions, &electrons_cold, &electrons_hot, &electrons_beam);

        //SolvePotential(phi, rho);
        SolvePotentialDirect(phi, rho);
        ComputeEF(phi, ef);

	    //move particles
        //PushSpecies(&ions, ef);
        PushSpecies(&electrons_cold, ef);
	    PushSpecies(&electrons_hot, ef);
	    PushSpecies(&electrons_beam, ef);

        //Write diagnostics
        if(ts%write_interval == 0)
        {
            sprintf(NAME,"%s/i%d.txt",output.c_str(),ts);
            //f1 = fopen(NAME,"w");

            sprintf(NAME,"%s/ec%d.txt",output.c_str(),ts);
            //f2 = fopen(NAME,"w");

            sprintf(NAME,"%s/eh%d.txt",output.c_str(),ts);
            //f3 = fopen(NAME,"w");

	        sprintf(NAME,"%s/eb%d.txt",output.c_str(),ts);
	        //f4 = fopen(NAME,"w");
            //===========================================================================================
            double max_phi = phi[0];
            for(int i=0; i<domain.ni; i++)
                if (phi[i]>max_phi) max_phi=phi[i];

			double max_vel_ion = ions.vel[0];
			for(int i=0; i<domain.ni; i++)
				if(ions.vel[i]>max_vel_ion) max_vel_ion = ions.vel[i];

            /*print diagnostics to the screen*/
	        printf("TS: %i \t delta_phi: %.3g \t max_vel_ion:%.3g \t nI:%ld \t nEC:%ld \t nEH:%ld \t nEB:%ld\n",
				   ts, max_phi-phi[0], max_vel_ion,ions.part_list.size(),electrons_cold.part_list.size(),electrons_hot.part_list.size(),electrons_beam.part_list.size());

	        /*Write time evolution of plasma profiles and Kinetic energy*/
    	    Write_ts(ts, &ions, &electrons_cold, &electrons_hot,&electrons_beam);
	        WriteKE(Time, &ions, &electrons_cold, &electrons_hot, &electrons_beam);
            WritePE(Time, &electrons_cold);
            
	        /*Write Electric Field Data at the Mid Plane*/
	        //WriteLocation(Time,domain.xl/2);

	        /*Write the Density Fluctuation at the Mid Plane*/
	        //WriteDenFluc(file_denflucI, Time, domain.xl/2, &ions);
	        //WriteDenFluc(file_denflucEC, Time, domain.xl/2, &electrons_cold);
	        //WriteDenFluc(file_denflucEH, Time, domain.xl/2, &electrons_hot);

	        /*Write individual particle data to  the file*/
    	    //Write_Particle(f1,ts, &ions);
    	    //Write_Particle(f2,ts, &electrons_cold);
	        //Write_Particle(f3,ts, &electrons_hot);
	        //Write_Particle(f4,ts, &electrons_beam);
    }

        Time += DT;
    }
    clock_t end = clock();

    /*free up memory*/
    delete phi;
    delete rho;
    delete ef;
    cout << "Time = " << ((end-start)/(double)CLOCKS_PER_SEC)/60 << " minutes" << endl;
    return 0;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/********************* HELPER FUNCTIONS ***************************/

/*Initialize the particle data : initial positions and velocities of each particle of each species*/
void Init(Species *species, string flag)
{
    // sample particle positions and velocities
    for(int p=0; p<species->NUM; p++)
    {
	double v;
        double x = domain.x0 + rnd()*(domain.ni-1)*domain.dx;
	if	(flag == "cold")
	{
		v = SampleVelCold(species->Temp*EV_TO_K, species->mass);
	}
        else if (flag == "hot")
	{
		v = SampleVelHot(species->Temp*EV_TO_K, species->mass);
	}
	else if (flag == "beam")
	{
		v = SampleVelBeam(species->Temp*EV_TO_K, species->mass);
	}
	else
	{
		v = SampleVel(species->Temp*EV_TO_K, species->mass);
	}

	// Add to the particle list
        species->add(Particle(x,v));
    }
}

/*Sample Velocity (According to Birdsall)*/
double SampleVel(double T, double mass)
{
    double v_th = sqrt(2*K*T/mass);
	double vt = v_th*sqrt(2)*(rnd()+rnd()+rnd()-1.5);
    // Normalize particle velocity by the thermal velocity of cold electrons & return
	return vt/(wpec*LDC); 
}

double SampleVelCold(double T, double mass)
{
	double v_th = sqrt(2*K*T/mass);
    double vt = v_th*sqrt(2)*(rnd()+rnd()+rnd()-1.5);
    // Normalize particle velocity by the thermal velocity of cold electrons & return
    return vt/(wpec*LDC);
}

double SampleVelHot(double T, double mass)
{
    double v_th = sqrt(2*K*T/mass);
    double vt = v_th*sqrt(2)*(rnd()+rnd()+rnd()-1.5);
    // Normalize particle velocity by the thermal velocity of cold electrons & return
    return vt/(wpec*LDC);
}

double SampleVelBeam(double T, double mass)
{
	double v_th = sqrt(2*K*T/mass);
	double vt = v_th*sqrt(2)*(rnd()+rnd()+rnd()-1.5) + vd*wpec*LDC;
    // Normalize particle velocity by the thermal velocity of cold electrons & return
	return vt/(wpec*LDC); 
}

/*Covert the physical coordinate to the logical coordinate*/
double XtoL(double pos)
{
    double li = (pos-domain.x0)/domain.dx;
    return li;
}

/*scatter the particle data to the massEsh and collect the densities at the massEsh */
void scatter(double lc, double value, double *field)
{
    int i = (int)lc;
    double di = lc-i;
    field[i] += value*(1-di);
    field[i+1] += value*(di);
}

/* Gather field values at logical coordinates*/
double gather(double lc, double *field)
{
    int i=(int)lc;
    double di = lc-i;
    double val = field[i]*(1-di) + field[i+1]*(di);
    return val;
}

/*Scatter the particles to the massEsh for evaluating densities*/
void ScatterSpecies(Species *species)
{
    /*grab a pointer to the species density data and change
     the density field using the pointer*/

    double *field = species->den;

    /*clear the field*/
    memset(field,0,sizeof(double)*domain.ni);

    /*scatter particles to the massEsh*/
    for(auto &p:species->part_list)
    {
        double lc = XtoL(p.pos);
        scatter(lc,species->spwt,field);
    }

    /*divide by cell volume*/
	/*Hint: we divide by the cell volume because we have scattered
	the spwt to the grid, which is just a number. Again number per
	unit volume is density. Hence we further divide the field by the cell volume*/
    for(int i=0; i<domain.ni; i++){
    	field[i] /=(domain.dx*LD);}

	/*Normalize the field value*/
	for(int i=0; i<domain.ni; i++){
		field[i] /=density;}

    field[0] *=2.0;
    field[domain.ni-1] *= 2.0;
}

/*Scatter the particles to the massEsh for evaluating velocities*/
void ScatterSpeciesVel(Species *species)
{
    /*grab a pointer to the species velocity field and change
     the velocity field using the pointer*/
    double *field = species->vel;

    /*clear the field*/
    memset(field,0,sizeof(double)*domain.ni);

    /*scatter particles to the mesh*/
    for(auto &p:species->part_list)
    {
        double lc = XtoL(p.pos);
        scatter(lc,species->spwt*p.vel,field);
		//scatter(lc,p.vel,field);
    }

    /*divide by cell volume*/
    for(int i=0; i<domain.ni; i++){
        field[i] /=(species->den[i]*density*domain.dx*LD);}


    field[0] *=2.0;
    field[domain.ni-1] *= 2.0;
}

//*******************************************************
void PushSpecies(Species *species, double *ef)
{
    // compute charge to mass ratio
    double qm = species->charge/species->mass;
    list<Particle>::iterator it = species->part_list.begin();

    // loop over particles
    while (it!=species->part_list.end())
    {
        // grab a reference to the pointer
        Particle &part = *it;

        // compute particle node position
        double lc = XtoL(part.pos);

        // gather electric field onto particle position
        double part_ef = gather(lc,ef);

        // advance velocity
	    //double wl = LD*LD*wp*wpec;
	    
        // Grab the old positions of the particles
        double old_pos = part.pos;
        //printf("old_pos: %f\n",old_pos);
        
        //part.vel += (1/wl)*chargeE*(qm*tempEcold/chargeE)*part_ef*DT;
	    part.vel += qm*(massE/e)*(wp/wpec)*(LD/LDC)*part_ef*DT;
        
        // Advance particle position
        part.pos += (wpec/wp)*(LDC/LD)*part.vel*DT;
	    
        // Grab the new positions of the particles
        double new_pos = part.pos;
        //printf("new_pos: %f\n",new_pos);
        
        // Check whether a particle is crossing one full cell length
        SanityCheck(new_pos, old_pos, domain.dx);
        
        // Take care of the particle that left the Domain (PBC)
		if (part.pos < domain.x0)
		{
			part.pos = part.pos + domain.xl;
		}
		else if(part.pos>=(domain.x0+domain.xl))
		{
			part.pos = part.pos - domain.xl;
		}
			it++;
    }
}

void SanityCheck(double new_pos, double old_pos, double cell_len)
{
    if (abs(new_pos-old_pos) > cell_len)
    {
        printf("Alert! Particles are crossing one full cell!\n");
        exit(-1);
    }
        
}

//*********************************************************
/*Rewind particle velocities by -0.5*DT */
void RewindSpecies(Species *species, double *ef)
{
    // compute charge to mass ratio
    double qm = species->charge/species->mass;
    for(auto &p:species->part_list)
    {
        // compute particle node position
        double lc = XtoL(p.pos);
        // gather electric field onto the particle position
        double part_ef = gather(lc,ef);
        //advance velocity
		//double wl = LD*LD*wp*wpec;
        //p.vel -= 0.5*chargeE*(1/wl)*(qm*tempEcold/chargeE)*part_ef*DT;
        p.vel -= 0.5*qm*(massE/e)*(wp/wpec)*(LD/LDC)*part_ef*DT;
    }
}

/* Compute the charge densities */
void ComputeRho(Species *ions, Species *electrons_cold, Species *electrons_hot, Species *electrons_beam)
{
    double *rho = domain.rho;
    memset(rho,0,sizeof(double)*domain.ni);

    for(int i=0; i<domain.ni; i++)
        rho[i] = (ions->den[i] - electrons_cold->den[i] - electrons_hot->den[i] - electrons_beam->den[i]);
}

/* Potential Solver: 1. Gauss-Seidel 2. Direct-Solver*/
bool SolvePotential(double *phi, double *rho)
{
    double L2;
    double dx2 = domain.dx*domain.dx;

    // Initialize boundaries
    phi[0]=phi[domain.ni-1]=0;

    // Main Solver
    for(int it=0; it<2000000; it++)
    {
        for(int i=1; i<domain.ni-1; i++)
        {
            double g = 0.5*(phi[i-1] + phi[i+1] + dx2*rho[i]);
            phi[i]=phi[i] + 1.4*(g-phi[i]);
        }
        // Check for convergence
        if(it%25==0)
        {
            double sum = 0;
            for(int i=1; i<domain.ni-1; i++)
            {
                double R = - rho[i] - (phi[i-1]-2*phi[i]+phi[i+1])/dx2;
                sum += R*R;
            }
            L2 = sqrt(sum)/domain.ni;
            if(L2<1e-4){return true;}

        }
        //printf("GS-Converged! L2=%g\n",L2);
    }
    printf("Gauss-Siedel solver failed to converge, L2=%g\n",L2);
    return false;
}

/* Potential Direct Solver */

bool SolvePotentialDirect(double *x, double *rho)
{
    /* Set coefficients, precompute them*/
    int ni = domain.ni;
    double dx2 = domain.dx*domain.dx;
    double *a = new double[ni];
    double *b = new double[ni];
    double *c = new double[ni];

    /*Centtral difference on internal nodes*/
    for(int i=1; i<ni-1; i++)
    {
        a[i] = 1; b[i] = -2; c[i] = 1;
    }

    /*Apply dirichlet boundary conditions on boundaries*/
    a[0]=0; b[0]=1; c[0]=0;
    a[ni-1]=0; b[ni-1]=1; c[ni-1]=0;

    /*multiply R.H.S.*/
    for (int i=1; i<ni-1; i++)
        x[i]=-rho[i]*dx2;

    x[0] = 0;
    x[ni-1] = 0;

    /*Modify the coefficients*/
    c[0] /=b[0];
    x[0] /=b[0];

    for(int i=1; i<ni; i++)
    {
        double id = (b[i]-c[i-1]*a[i]);
        c[i] /= id;
        x[i] = (x[i]-x[i-1]*a[i])/id;
    }

    /* Now back substitute */
    for(int i=ni-2; i>=0; i--)
        x[i] = x[i] - c[i]*x[i+1];

    return true;
}

/*Compute electric field (differentiating potential)*/
void ComputeEF(double *phi, double *ef)
{
    /*Apply central difference to the inner nodes*/
    for(int i=1; i<domain.ni-1; i++)
        ef[i] = -(phi[i+1]-phi[i-1])/(2*domain.dx);

    /*Apply one sided difference at the boundary nodes*/
    ef[0] = -(phi[1]-phi[0])/domain.dx;
    ef[domain.ni-1] = -(phi[domain.ni-1]-phi[domain.ni-2])/domain.dx;
}


/*Write the output with Time*/
void Write_ts(int ts, Species *ions,Species *electrons_cold, Species *electrons_hot,Species *electrons_beam)
{
    for(int i=0; i<domain.ni; i++)
    {
        fprintf(file_res,"%g \t %g \t %g \t %g \t %g \t %g\n", i*domain.dx, electrons_cold->den[i],
		electrons_hot->den[i], electrons_beam->den[i], domain.phi[i], domain.ef[i]);

    }
    fflush(file_res);
}

/* Write the data of a particular location */
void WriteLocation(double Time, double pos)
{
	double lc = XtoL(pos);
	int i = (int) lc;
	fprintf(file_loc,"%g \t %g\n",Time,domain.ef[i]);
}

void WriteDenFluc(FILE *file, double Time, double pos, Species *species)
{
	double lc = XtoL(pos);
	int i = (int) lc;
	fprintf(file,"%g \t %g\n",Time,species->den[i] - density);
	fflush(file);
}

/* Write the Output results*/
void Write_Particle(FILE *file, int ts, Species *species)
{
    for(auto& p: species->part_list)
    {
        fprintf(file,"%g \t %g\n",p.pos, p.vel);
    }
    fflush(file_res);
}


void WriteKE(double Time, Species *ions, Species *electrons_cold, Species *electrons_hot, Species *electrons_beam)
{
    double ke_ions = ComputeKE(ions, electrons_cold);
    double ke_electrons_cold = ComputeKE(electrons_cold, electrons_cold);
    double ke_electrons_hot = ComputeKE(electrons_hot, electrons_cold);
    double ke_electrons_beam = ComputeKE(electrons_beam, electrons_cold);

    fprintf(file_ke,"%g \t %g \t %g \t %g \t %g\n",Time, ke_ions, ke_electrons_cold, ke_electrons_hot, ke_electrons_beam);

    fflush(file_ke);
}
void WritePE(double Time, Species *electrons_cold)
{
    double pe = ComputePE(electrons_cold);
    fprintf(file_pe, "%g \t %g\n", Time, pe);
}
double ComputeKE(Species *species, Species *electrons_cold)
{
    double ke = 0;
    for (auto &p:species->part_list)
    {
        // un-normalize the velocity by multiplying with the cold thermal velocity
        ke += (p.vel*p.vel)*(wpec*LDC)*(wpec*LDC);
    }
    /*Multiply 0.5*mass for all particles*/
    ke *= 0.5*(species->spwt*species->mass);
    
    // Calculate the total thermal energy of all the cold electrons
    double Th = (electrons_cold->Temp*eV)*(electrons_cold->spwt)*nParticlesE;

    // Normalize the kinetic energy by the total thermal energy of cold electrons    
    ke = ke/Th;
    return ke;
}

double ComputePE(Species *electrons_cold)
{
    double pe;
    // obtain the sum of electric field over 
    // the whole domain at a definite time step
    double ef = 0;
    for (int i=0; i<domain.ni; i++)
    {
        ef += domain.ef[i];
    }
    // un-normalize this electric field
    ef = ef*(massE*wp*wp*LD/e);
    // calculate the un-normalized electric potential energy
    pe = 0.5*EPS*(ef*ef);
    
    // Calculate the total thermal energy of all the cold electrons
    double Th = (electrons_cold->Temp*eV)*(electrons_cold->spwt)*nParticlesE;
    // normalize the potential energy by the total thermal 
    // energy of the cold electrons
    pe = pe/Th;
    return pe;
}
