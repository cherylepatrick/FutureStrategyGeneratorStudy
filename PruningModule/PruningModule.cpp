#include "PruningModule.h"

using namespace std;


DPP_MODULE_REGISTRATION_IMPLEMENT(PruningModule,"PruningModule");
PruningModule::PruningModule() : dpp::base_module()
{
  filename_output_="PrunedSimulation.root";
}

PruningModule::~PruningModule() {
  if (is_initialized()) this->reset();
}

void PruningModule::initialize(const datatools::properties& myConfig,
                                   datatools::service_manager& flServices,
                                   dpp::module_handle_dict_type& /*moduleDict*/){

  // Look for services
  if (flServices.has("geometry")) {
    const geomtools::geometry_service& GS = flServices.get<geomtools::geometry_service> ("geometry");
    geometry_manager_ = &GS.get_geom_manager();
    DT_THROW_IF(!geometry_manager_,
                std::runtime_error,
                "Null pointer to geometry manager return by geometry_service");
  }
  // Extract the filename_out key from the supplied config, if
  // the key exists. datatools::properties throws an exception if
  // the key isn't in the config, so catch this if thrown and don't do
  // anything
  try {
    myConfig.fetch("filename_out",this->filename_output_);
  } catch (std::logic_error& e) {
  }

  // Use the method of PTD2ROOT to create a root file with just trueparticle.kinenergy

  hfile_ = new TFile(filename_output_.c_str(),"RECREATE","Output file of pruned simulation data");
  hfile_->cd();
  tree_ = new TTree("SimData","SimData");
  tree_->SetDirectory(hfile_);
 
  // Some basic counts
  tree_->Branch("trueparticle.kinenergy",&pruned_.energy_);
  this->_set_initialized(true);
}
//! [ValidationModule::Process]
dpp::base_module::process_status
PruningModule::process(datatools::things& workItem) {
  

  // We need to run this before we start populating vectors. Put all your vectors in this function to clear them
  ResetVars();
  
  // Grab calibrated data bank
  // Calibrated data will only be present in reconstructed files,
  // so wrap in a try block



  try {
    // Get (true) primary vertex position
      const mctools::simulated_data& simData = workItem.get<mctools::simulated_data>("SD");
      if (simData.has_data())
      {
        mctools::simulated_data::primary_event_type primaryEvent=simData.get_primary_event ();
        
        for (int i=0;i<primaryEvent.get_number_of_particles();i++)// should be 2 particles for 0nubb
        {
          genbb::primary_particle trueParticle= primaryEvent.get_particle(i);
          double energy=trueParticle.get_kinetic_energy();
          pruned_.energy_.push_back(energy);
        }
      }
  }// end try for SD bank
    catch (std::logic_error& e) {
      //std::cerr << "failed to grab SD bank : " << e.what() << std::endl;
      //return dpp::base_module::PROCESS_ERROR;
      // This is OK, if it's data there will be no SD bank
    }
  tree_->Fill();
  // MUST return a status, see ref dpp::processing_status_flags_type
  return dpp::base_module::PROCESS_OK;
}


void PruningModule::ResetVars()
{
  pruned_.energy_.clear();
  return;
}


//! [PruningModule::reset]
void PruningModule::reset() {
  hfile_->cd();
  tree_->Write();
  hfile_->Close(); //
  std::cout << "In reset: finished conversion, file closed " << std::endl;

  // clean up
  delete hfile_;
  filename_output_ = "PrunedSimulation.root";
  this->_set_initialized(false);

}

