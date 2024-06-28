#ifndef AHDC_HITPROCESS_H
#define AHDC_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"
#include <string>

class ahdcConstants
{
public:
	
	// Database parameters
	int    runNo;
	string date;
	string connection;
	char   database[80];
	
	// translation table
	TranslationTable TT;
};


// Class definition
/// \class ahdc_HitProcess
/// <b> Alert Drift Chamber Hit Process Routine</b>\n\n

class ahdc_HitProcess : public HitProcess
{
public:
	
	~ahdc_HitProcess(){;}
	
	// constants initialized with initWithRunNumber
	static ahdcConstants atc;
	
	void initWithRunNumber(int runno);
	
	// - integrateDgt: returns digitized information integrated over the hit
	map<string, double> integrateDgt(MHit*, int);
	
	// - multiDgt: returns multiple digitized information / hit
	map< string, vector <int> > multiDgt(MHit*, int);
	
	// - charge: returns charge/time digitized information / step
	virtual map< int, vector <double> > chargeTime(MHit*, int);
	
	// - voltage: returns a voltage value for a given time. The input are charge value, time
	virtual double voltage(double, double, double);
	
	// The pure virtual method processID returns a (new) identifier
	// containing hit sharing information
	vector<identifier> processID(vector<identifier>, G4Step*, detector);
	
	// creates the HitProcess
	static HitProcess *createHitClass() {return new ahdc_HitProcess;}
	
	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();
	
public:
	// AHDC geometry parameters
	float PAD_W, PAD_L, PAD_S, RTPC_L;
	float phi_per_pad;
	
	// parameters for drift and diffustion equations for drift time, 
	// drift angle, and diffusion in z
	float a_t, b_t, c_t, d_t;
	float a_phi, b_phi, c_phi, d_phi;
	float a_z, b_z;
	
	// variables for storing drift times and diffusion in time
	float t_2GEM2, t_2GEM3, t_2PAD, t_2END;
	float sigma_t_2GEM2, sigma_t_2GEM3, sigma_t_2PAD, sigma_t_gap;
	
	// variables for storing drift angle and diffusion in phi
	float phi_2GEM2, phi_2GEM3, phi_2PAD, phi_2END;
	float sigma_phi_2GEM2, sigma_phi_2GEM3, sigma_phi_2PAD, sigma_phi_gap;
	
	float z_cm;
	float TPC_TZERO;
	
	map<int, double> timeShift_map;
	double shift_t;

public:
	// added by Felix
	void 	ShowMeHitContent	(MHit* aHit, int hitn); // to be modified (I want to center the wire at the origin of the plan, get insired by H_abh or doca calculation
	
};

namespace futils {
	bool cart2polar3D(double x, double y, double z, double & rho, double & theta, double & phi);
}

#include "Math/PdfFuncMathCore.h"

class ahdcSignal {
	// MHit identifiers
	private : 
		int hitn;
		int sector;
		int layer;
		int component;
		int nsteps;
	// vectors
	private :
		std::vector<double> Edep; // keV
		std::vector<double> G4Time; // ns
		std::vector<double> Doca; //mm
		std::vector<double> DriftTime; // ns
		void ComputeDocaAndTime(MHit * aHit);
		std::vector<double> Dgtz; 
		std::vector<double> Noise;
	// setting parameters for digitization
	private : 
		double tmin = 0; 
		double tmax = 6000; 
		double delay = 1000;
		double samplingTime = 44;      // ns
		double electronYield = 9500;   // ADC_gain
		int    adc_max = 50000;
	// public methods
	public :
		ahdcSignal() = default;
		ahdcSignal(MHit * aHit, int hitn_){
			// read identifiers
			hitn = hitn_;
			vector<identifier> identity = aHit->GetId();
			sector = 0;
			layer = 10 * identity[0].id + identity[1].id ; // 10*superlayer + layer
			component = identity[2].id;
			// fill vectors
			Edep = aHit->GetEdep();
			nsteps = Edep.size();
			for (int s=0;s<nsteps;s++){Edep.at(s) = Edep.at(s)*1000;} // convert MeV to keV
			G4Time = aHit->GetTime();
			this->ComputeDocaAndTime(aHit); // fills Doca and DriftTime
		}

		~ahdcSignal(){;}
		double 					GetSamplingTime()	 {return samplingTime;}
		double                                  GetElectronYield()	 {return electronYield;}
		int	                                GetAdcMax()		 {return adc_max;}
		double 					GetTmin() 		{return tmin;}
		double                                  GetTmax()               {return tmax;}
		std::vector<double> 			GetNoise() 		{return Noise;}
		double 					GetDelay()		{return delay;}
		void SetSamplingTime(double samplingTime_)			{samplingTime = samplingTime_;}
		void SetElectronYield(double electronYield_)			{electronYield = electronYield_;}
		void SetAdcMax(int adc_)				{adc_max = adc_;}
		void SetTmin(double tmin_)				{tmin = tmin_;}
		void SetTmax(double tmax_)                              {tmax = tmax_;}
		void SetNoise(std::vector<double> Noise_)		{Noise = Noise_;}
		void SetDelay(double delay_) 				{delay = delay_;}

		double operator()(double t){
			double res = 0;
			for (int s=0; s<nsteps; s++){
				res += Edep.at(s)*ROOT::Math::landau_pdf(t-delay,600/2.5,DriftTime.at(s)); // landau width parameter = 600/2.5; think to create a new para.
			}
			return res;
		}

		void Digitize();
		void GenerateNoise(double mean, double stdev);

		void PrintBeforeProcessing();
		void PrintAllShapes();
		void PrintAfterProcessing();
		void PrintSignal(); // to be modified // change hist by graph ?
		void PrintNoise(); 
};




#endif




















