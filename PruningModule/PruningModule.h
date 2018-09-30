//! \file    PruningModule.h
//! \brief   Prune a simulated data tree to just save particle energies
//! \details Output basic branches for the Validation parser
#ifndef PRUNEMODULE_HH
#define PRUNEMODULE_HH
// Standard Library
// Third Party
//#include <boost/foreach.hpp>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "TVector3.h"

// - Bayeux
#include "bayeux/dpp/base_module.h"
#include "bayeux/mctools/simulated_data.h"
#include "bayeux/genbb_help/primary_particle.h"
#include "bayeux/genbb_help/primary_event.h"
#include "bayeux/datatools/service_manager.h"
#include "bayeux/geomtools/manager.h"
#include "bayeux/geomtools/geometry_service.h"
#include "bayeux/geomtools/line_3d.h"
#include "bayeux/geomtools/helix_3d.h"
#include "bayeux/geomtools/geomtools.h"

// - Falaise
#include "falaise/snemo/datamodels/calibrated_data.h"
#include "falaise/snemo/datamodels/tracker_clustering_data.h"
#include "falaise/snemo/datamodels/tracker_clustering_solution.h"
#include "falaise/snemo/datamodels/particle_track_data.h"


typedef struct PrunedEventStorage{
  // Quantities to histogram (h_)
  std::vector<double> energy_;


}Prunedeventstorage;


// This Project
class PruningModule : public dpp::base_module {

 public:
  //! Construct module
  PruningModule();
  //! Destructor
  virtual ~PruningModule();
  //! Configure the module
  virtual void initialize(const datatools::properties& myConfig,
                          datatools::service_manager& flServices,
                          dpp::module_handle_dict_type& moduleDict);
  //! Process supplied data record
  virtual dpp::base_module::process_status process(datatools::things& workItem);
  //! Reset the module
  virtual void reset();
 private:
  TFile* hfile_;
  TTree* tree_;
  PrunedEventStorage pruned_;

  // configurable data member
  std::string filename_output_;

  // geometry service
  const geomtools::manager* geometry_manager_; //!< The geometry manager

  void ResetVars();

  // Macro which automatically creates the interface needed
  // to enable the module to be loaded at runtime
  DPP_MODULE_REGISTRATION_INTERFACE(PruningModule);
};

#endif // PRUNEMODULE_HH


