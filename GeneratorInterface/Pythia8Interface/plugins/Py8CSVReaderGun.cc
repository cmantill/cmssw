#include <fstream>
#include "GeneratorInterface/Core/interface/GeneratorFilter.h"
#include "GeneratorInterface/ExternalDecays/interface/ExternalDecayDriver.h"
#include "GeneratorInterface/Pythia8Interface/interface/Py8GunBase.h"
namespace gen {
class Py8CSVReaderGun : public Py8GunBase {
   
   public:
      
      Py8CSVReaderGun( edm::ParameterSet const& );
      ~Py8CSVReaderGun() override {}
      bool generatePartonsAndHadronize() override;
      const char* classname() const override;
	 
   private:
      
      // PtGun particle(s) characteristics
      double  fMinEta;
      double  fMaxEta;
      double  fMinPt ;
      double  fMaxPt ;
      double  fMinProdRadius;
      double  fMaxProdRadius;
      bool    fMakeDisplaced;
      int     fNumParticlesPerEvent;
      std::string fFilename;
      std::vector<float> all_ee, all_px, all_py, all_pz;
      std::vector<int> used_events;
};
// implementation 
//
Py8CSVReaderGun::Py8CSVReaderGun( edm::ParameterSet const& ps )
   : Py8GunBase(ps), used_events({}) {
   // ParameterSet defpset ;
   edm::ParameterSet pgun_params = 
      ps.getParameter<edm::ParameterSet>("PGunParameters"); // , defpset ) ;
   fFilename = pgun_params.getParameter<std::string>("Filename");
   fMinEta     = pgun_params.getParameter<double>("MinEta"); // ,-2.2);
   fMaxEta     = pgun_params.getParameter<double>("MaxEta"); // , 2.2);
   fMinPt      = pgun_params.getParameter<double>("MinPt"); // ,  0.);
   fMaxPt      = pgun_params.getParameter<double>("MaxPt"); // ,  0.);
   fMinProdRadius = pgun_params.getParameter<double>("MinProdRadius"); // , 0.);
   fMaxProdRadius = pgun_params.getParameter<double>("MaxProdRadius"); // , 0.);
   fMakeDisplaced = pgun_params.getParameter<bool>("MakeDisplaced"); //, true);
   fNumParticlesPerEvent = pgun_params.getParameter<int>("NumParticlesPerEvent"); // 5
   std::cout << "[Py8CSVReaderGun constructor] Begin reading CSV input file..." << std::endl;
   std::cout << "[Py8CSVReaderGun constructor] Filename: " << fFilename << std::endl;
   std::ifstream infile(fFilename);
   float ee, px, py, pz;
   while (infile >> ee >> px >> py >> pz) {
      all_ee.push_back(ee);
      all_px.push_back(px);
      all_py.push_back(py);
      all_pz.push_back(pz);
   }
   infile.close();
   std::cout << "[Py8CSVReaderGun constructor] Finished reading CSV input file!" << std::endl;
}
bool Py8CSVReaderGun::generatePartonsAndHadronize()
{
   fMasterGen->event.reset();
   // compute common radius (only used if fMakeDisplaced set to True in config)
   double radius = (fMaxProdRadius-fMinProdRadius) * randomEngine().flat() + fMinProdRadius;
   double phi_prod = (fMaxPhi-fMinPhi) * randomEngine().flat() + fMinPhi;
   double vx = radius * cos(phi_prod);
   double vy = radius * sin(phi_prod);
   double vz = (70 - (-70)) * randomEngine().flat() + (-70); // luminous region in Z: (-70, 70) mm
   // ensure event is unique within single node (accept repetition after 100 times though)
   // note: this of course does not apply for batch production, where probability of repetition exists
   // (this is minimized by randomly sampling pluto list of events -- birthday problem)
   int randomNumber, count = 0;
   do {
      randomNumber = (int)(30000 * randomEngine().flat()) * fNumParticlesPerEvent;
      count++;
   }
   while (std::find(used_events.begin(), used_events.end(), randomNumber) != used_events.end() && count < 100);
   used_events.push_back(randomNumber);
   std::cout << "Retrieving CSV random event number " << randomNumber/fNumParticlesPerEvent << "..." << std::endl;
   // Get the gamma (2nd in CSV) and electrons (4th and 5th in CSV) momenta
   float ee, px, py, pz, mass, pID; 
   // gamma
   ee = all_ee.at(randomNumber+1);
   px = all_px.at(randomNumber+1);
   py = all_py.at(randomNumber+1);
   pz = all_pz.at(randomNumber+1);
   pID = 22;
   mass = (fMasterGen->particleData).m0(pID);
   (fMasterGen->event).append(pID, 1, 0, 0, px, py, pz, ee, mass); 
   // ele plus
   ee = all_ee.at(randomNumber+3);
   px = all_px.at(randomNumber+3);
   py = all_py.at(randomNumber+3);
   pz = all_pz.at(randomNumber+3);
   pID = 11;
   mass = (fMasterGen->particleData).m0(pID);
   (fMasterGen->event).append(pID, 1, 0, 0, px, py, pz, ee, mass); 
   // ele minus
   ee = all_ee.at(randomNumber+4);
   px = all_px.at(randomNumber+4);
   py = all_py.at(randomNumber+4);
   pz = all_pz.at(randomNumber+4);
   pID = -11;
   mass = (fMasterGen->particleData).m0(pID);
   (fMasterGen->event).append(pID, 1, 0, 0, px, py, pz, ee, mass); 

   if ( !fMasterGen->next() ) return false;
   evtGenDecay();
   
   event().reset(new HepMC::GenEvent);
   return toHepMC.fill_next_event( fMasterGen->event, event().get() );
  
}
const char* Py8CSVReaderGun::classname() const
{
   return "Py8CSVReaderGun"; 
}
typedef edm::GeneratorFilter<gen::Py8CSVReaderGun, gen::ExternalDecayDriver> Pythia8CSVReaderGun;
} // end namespace
using gen::Pythia8CSVReaderGun;
DEFINE_FWK_MODULE(Pythia8CSVReaderGun);
