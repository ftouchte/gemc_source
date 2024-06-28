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
	void 	ComputeDocaAndTime	(MHit* aHit, std::vector<double> & Height, std::vector<double> & Time); // could it be better to add this in ahdcSignal ? Same for ShowMeHitContent ?
	
};

namespace futils {
	bool cart2polar3D(double x, double y, double z, double & rho, double & theta, double & phi);
}

#include "Math/PdfFuncMathCore.h"

class ahdcSignal {
	private :
		MHit * aHit;
		int hitn;
		std::vector<double> Location;  // ns   // Geant4 time associated to each step
		std::vector<double> Amplitude; // keV  // Edep associated to each step
		void ComputeDocaAndTime();
		std::vector<double> Height;    // mm   // doca of each step
		std::vector<double> Time;      // ns   // drift time computed with the doca
		
		std::vector<int> Shape;        // allow to have a specific shape associated to each step (ex: gaussians with different std_dev )
		std::vector<double> Width;     // ns   // width parameter of the Landau distribution
	private : // for digitization
		double samplingTime = 44;      // [ns]
		double electronYield = 9500;   // ADC_gain
		int adc_max = 50000;           // saturation for digitization
	private : // for plot
		double tmin; 
		double tmax; 
		double delay = 1000;
		std::vector<double> Dgtz;
		std::vector<double> Noise;
	
	public :
		ahdcSignal(MHit * aHit_, int hitn_){
			aHit = aHit_;
			hitn = hitn;
			// to be complete ...
			std::cout << "     >> Inside ahdcSignal : before this->ComputeDocaAndTime();" << std::endl;
			this->ComputeDocaAndTime(); // That fills the vectors Height
			vector<G4double> Edep = aHit->GetEdep(); // à améliorer
			int nsteps = Edep.size();
			for (int s=0;s<nsteps;s++){
				this->Add(Time.at(s),Edep.at(s)*1000,600/2.5,1); // Edep converted in keV
			}
		}
		ahdcSignal(std::vector<double> Location_, std::vector<double> Amplitude_, std::vector<double> Width_, std::vector<int> Shape_){
			Location = Location_;
			Amplitude = Amplitude_;
			Width = Width_;
			Shape = Shape_;
		}
		ahdcSignal() = default;
		ahdcSignal(const ahdcSignal & obj){
			Location = obj.Location;
			Amplitude = obj.Amplitude;
			Width = obj.Width;
			Shape = obj.Shape;
			samplingTime = obj.samplingTime;
			electronYield = obj.electronYield;
			adc_max = obj.electronYield;
		}
		~ahdcSignal(){;}
		void Add(double time, double amplitude, double width, int shape){
			Location.push_back(time);
			Amplitude.push_back(amplitude);
			Width.push_back(width);
			Shape.push_back(shape);
		}
		std::vector<double>                     GetAmplitude()		{return Amplitude;}
		std::vector<double>                     GetLocation() 		{return Location;}
		std::vector<double>                     GetWidth()		{return Width;}
		std::vector<int>                        GetShape()		{return Shape;}
		double 					GetSamplingTime()	 {return samplingTime;}
		double                                  GetElectronYield()	 {return electronYield;}
		int	                                GetAdcMax()		 {return adc_max;}
		double 					GetTmin() 		{return tmin;}
		double                                  GetTmax()               {return tmax;}
		std::vector<double> 			GetNoise() 		{return Noise;}
		double 					GetDelay()		{return delay;}
		void SetAmplitude(std::vector<double> Amplitude_) 	{Amplitude = Amplitude_;}
		void SetLocation(std::vector<double> Location_)		{Location = Location_;}
		void SetWidth(std::vector<double> Width_)		{Width = Width_;}
		void SetShape(std::vector<int> Shape_)			{Shape = Shape_;}
		void SetSamplingTime(double samplingTime_)			{samplingTime = samplingTime_;}
		void SetElectronYield(double electronYield_)			{electronYield = electronYield_;}
		void SetAdcMax(int adc_)				{adc_max = adc_;}
		void SetTmin(double tmin_)				{tmin = tmin_;}
		void SetTmax(double tmax_)                              {tmax = tmax_;}
		void SetNoise(std::vector<double> Noise_)		{Noise = Noise_;}
		void SetDelay(double delay_) 				{delay = delay_;}
		bool is_safe(){
			int n1 = Location.size();
			int n2 = Amplitude.size();
			int n3 = Width.size();
			int n4 = Shape.size();
			return (n1 == n2) && (n2 == n3) && (n3 == n4);
		}
		// The reason why this class has been created
		double operator()(double x){
			double res = 0;
			int nLoc = Location.size();
			for (int l=0;l<nLoc;l++){
				if (Shape.at(l) == 0) {
					res += Amplitude.at(l)*ROOT::Math::gaussian_pdf(x-delay,Width.at(l),Location.at(l));
				}
				else if (Shape.at(l) == 1){
					res += Amplitude.at(l)*ROOT::Math::landau_pdf(x-delay,Width.at(l),Location.at(l));
				}
			}
			return res; // in keV
		}
		void PrintBeforeProcessing(const char * filename);
		void PrintAllShapes(const char * filename);
		void PrintAfterProcessing(const char * filename);
		void PrintDgtz(const char * filename); // to be modified // change hist by graph ?
		void PrintNoise(const char * filename); 
		void PrintResult(const char * filename);
		
		void GenerateNoise(double mean, double stdev);
		void Digitize();
};




#endif




















