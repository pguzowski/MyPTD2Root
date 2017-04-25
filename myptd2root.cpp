// - Implementation of MyPTD2Root
// Example on access to data in 'brio' files from flsimulate:
// Uses the 'PTD' data bank
// 1) Access all available data in PTD data banks
// 2) Write data to a flat TTree ROOT file

// Ourselves
#include "myptd2root.h"

#include <mctools/simulated_data.h>

// Standard Library
// Third Party
// - A
// This Project
// Macro which automatically implements the interface needed
// to enable the module to be loaded at runtime
// The first argument is the typename
// The second is the string key used to access the module in pipeline
// scripts. This must be globally unique.

DPP_MODULE_REGISTRATION_IMPLEMENT(MyPTD2Root,"MyPTD2Root");

// Construct
MyPTD2Root::MyPTD2Root() : dpp::base_module()
{
  filename_output_="default.root";
}

// Destruct
MyPTD2Root::~MyPTD2Root() {
  // MUST reset module at destruction
  this->reset();
}

// Initialize
void MyPTD2Root::initialize(const datatools::properties& myConfig,
			  datatools::service_manager& flServices,
			  dpp::module_handle_dict_type& /*moduleDict*/) {

  // Throw logic exception if we've already initialized this instance
  DT_THROW_IF(this->is_initialized(),
	      std::logic_error,
	      "MyPTD2Root already initialized");
  // Extract the filename_out key from the supplied config, if
  // the key exists. datatools::properties throws an exception if
  // the key isn't in the config, so catch this if thrown and don't do
  // anything
  try {
    myConfig.fetch("filename_out",this->filename_output_);
  } catch (std::logic_error& e) {
  }

  // Look for services
  if (flServices.has("geometry"));
  {
    const geomtools::geometry_service& GS = flServices.get<geomtools::geometry_service> ("geometry");

    // initialize geometry manager
    //    std::cout << "Initialize geo manager " << std::endl;
    geometry_manager_ = &GS.get_geom_manager();
    DT_THROW_IF(!geometry_manager_,
                std::runtime_error,
                "Null pointer to geometry manager return by geometry_service");
  }

  std::cout << "In INIT: create TFile " << std::endl;
  // Next all root file output here

  hfile_ = new TFile(filename_output_.c_str(),"RECREATE","Output file of Simulation data");
  hfile_->cd();

  tree_ = new TTree("PTD","PTD");
  tree_->SetDirectory(hfile_);

  // header data
  tree_->Branch("header.runnumber",&header_.runnumber_);
  tree_->Branch("header.eventnumber",&header_.eventnumber_);
  tree_->Branch("header.date",&header_.date_);
  tree_->Branch("header.runtype",&header_.runtype_);
  tree_->Branch("header.simulated",&header_.simulated_);

  // generator data
  tree_->Branch("truth.vetex_x",&gen_.vertex_x_);
  tree_->Branch("truth.vetex_y",&gen_.vertex_y_);
  tree_->Branch("truth.vetex_z",&gen_.vertex_z_);

  // particle data
  tree_->Branch("particle.nofparticles",&particle_.nofparticles_);
  tree_->Branch("particle.nofgammaonly",&particle_.nofgammas_);
  tree_->Branch("particle.particleID",&particle_.particle_id_);
  tree_->Branch("particle.charge",&particle_.charge_);
  tree_->Branch("particle.vertex1_type",&particle_.vertex1_type_);
  tree_->Branch("particle.vertex1_x",&particle_.vertex1_x_);
  tree_->Branch("particle.vertex1_y",&particle_.vertex1_y_);
  tree_->Branch("particle.vertex1_z",&particle_.vertex1_z_);
  tree_->Branch("particle.foil1_dir_x",&particle_.foil1_dir_x_);
  tree_->Branch("particle.foil1_dir_y",&particle_.foil1_dir_y_);
  tree_->Branch("particle.foil1_dir_z",&particle_.foil1_dir_z_);
  tree_->Branch("particle.vertex2_type",&particle_.vertex2_type_);
  tree_->Branch("particle.vertex2_x",&particle_.vertex2_x_);
  tree_->Branch("particle.vertex2_y",&particle_.vertex2_y_);
  tree_->Branch("particle.vertex2_z",&particle_.vertex2_z_);
  tree_->Branch("particle.foil2_dir_x",&particle_.foil2_dir_x_);
  tree_->Branch("particle.foil2_dir_y",&particle_.foil2_dir_y_);
  tree_->Branch("particle.foil2_dir_z",&particle_.foil2_dir_z_);
  tree_->Branch("particle.traj_length",&particle_.traj_length_);
  tree_->Branch("particle.traj_cl_del",&particle_.traj_cluster_delayed_);
  tree_->Branch("particle.traj_cl_del_time",&particle_.traj_cluster_delayed_time_);
  tree_->Branch("particle.calo1_associated",&particle_.calo1_associated_);
  tree_->Branch("particle.calo1_type",&particle_.calo1_type_);
  tree_->Branch("particle.calo1_energy",&particle_.calo1_energy_);
  tree_->Branch("particle.calo1_sigma_energy",&particle_.calo1_sigma_energy_);
  tree_->Branch("particle.calo1_time",&particle_.calo1_time_);
  tree_->Branch("particle.calo1_sigma_time",&particle_.calo1_sigma_time_);
  tree_->Branch("particle.calo1_side",&particle_.calo1_side_);
  tree_->Branch("particle.calo1_column",&particle_.calo1_column_);
  tree_->Branch("particle.calo1_row",&particle_.calo1_row_);
  tree_->Branch("particle.calo1_wall",&particle_.calo1_wall_);
  tree_->Branch("particle.calo1_loc_x",&particle_.calo1_loc_x_);
  tree_->Branch("particle.calo1_loc_y",&particle_.calo1_loc_y_);
  tree_->Branch("particle.calo1_loc_z",&particle_.calo1_loc_z_);
  tree_->Branch("particle.calo2_associated",&particle_.calo2_associated_);
  tree_->Branch("particle.calo2_type",&particle_.calo2_type_);
  tree_->Branch("particle.calo2_energy",&particle_.calo2_energy_);
  tree_->Branch("particle.calo2_sigma_energy",&particle_.calo2_sigma_energy_);
  tree_->Branch("particle.calo2_time",&particle_.calo2_time_);
  tree_->Branch("particle.calo2_sigma_time",&particle_.calo2_sigma_time_);
  tree_->Branch("particle.calo2_side",&particle_.calo2_side_);
  tree_->Branch("particle.calo2_column",&particle_.calo2_column_);
  tree_->Branch("particle.calo2_row",&particle_.calo2_row_);
  tree_->Branch("particle.calo2_wall",&particle_.calo2_wall_);
  tree_->Branch("particle.calo2_loc_x",&particle_.calo2_loc_x_);
  tree_->Branch("particle.calo2_loc_y",&particle_.calo2_loc_y_);
  tree_->Branch("particle.calo2_loc_z",&particle_.calo2_loc_z_);

  this->_set_initialized(true);
}

// Process
dpp::base_module::process_status MyPTD2Root::process(datatools::things& workItem) {
  // Local variables

  // particle event data
  std::vector<int> particleid;
  std::vector<int> charge;
  std::vector<int> vertex1type;
  std::vector<double> vertex1_x;
  std::vector<double> vertex1_y;
  std::vector<double> vertex1_z;
  std::vector<double> foil1_dir_x;
  std::vector<double> foil1_dir_y;
  std::vector<double> foil1_dir_z;
  std::vector<int> vertex2type;
  std::vector<double> vertex2_x;
  std::vector<double> vertex2_y;
  std::vector<double> vertex2_z;
  std::vector<double> foil2_dir_x;
  std::vector<double> foil2_dir_y;
  std::vector<double> foil2_dir_z;
  std::vector<double> traj_length;
  std::vector<int> traj_cl_delayed;
  std::vector<double> traj_cl_delayed_time;
  std::vector<int> calo1associated;
  std::vector<int> calo1type;
  std::vector<double> calo1energy;
  std::vector<double> calo1sigmaenergy;
  std::vector<double> calo1time;
  std::vector<double> calo1sigmatime;
  std::vector<int> calo1side;
  std::vector<int> calo1column;
  std::vector<int> calo1row;
  std::vector<int> calo1wall;
  std::vector<double> calo1_loc_x;
  std::vector<double> calo1_loc_y;
  std::vector<double> calo1_loc_z;
  std::vector<int> calo2associated;
  std::vector<int> calo2type;
  std::vector<double> calo2energy;
  std::vector<double> calo2sigmaenergy;
  std::vector<double> calo2time;
  std::vector<double> calo2sigmatime;
  std::vector<int> calo2side;
  std::vector<int> calo2column;
  std::vector<int> calo2row;
  std::vector<int> calo2wall;
  std::vector<double> calo2_loc_x;
  std::vector<double> calo2_loc_y;
  std::vector<double> calo2_loc_z;

  // Access the workItem

  // look for reconstructed data
  if(workItem.has("PTD"))
    {
      const snemo::datamodel::particle_track_data & PTD = workItem.get<snemo::datamodel::particle_track_data>("PTD");

      // Extract particle track data
      if (PTD.has_particles()) {
	particle_.nofparticles_ = PTD.get_number_of_particles();

	for (size_t i=0; i<PTD.get_number_of_particles();++i) {
	  const snemo::datamodel::particle_track & the_particle = PTD.get_particle(i);

	  particleid.push_back(the_particle.get_track_id());
	  charge.push_back(the_particle.get_charge());

	  // first the vertices
	  if (the_particle.has_vertices()) {
            unsigned int i = 0;
	    if (i<the_particle.get_vertices().size()) {
	      const geomtools::blur_spot & vertex = the_particle.get_vertices().at(i).get();
	      const geomtools::vector_3d & translation  = vertex.get_placement().get_translation();
	      vertex1_x.push_back(translation.x());
	      vertex1_y.push_back(translation.y());
	      vertex1_z.push_back(translation.z());
	      if (vertex.has_geom_id()) { // vertex is on calorimeter (calo or xcalo or gveto)
                vertex1type.push_back(vertex.get_geom_id().get_type());
                foil1_dir_x.push_back(0.);
                foil1_dir_y.push_back(0.);
                foil1_dir_z.push_back(0.);
	      }
	      else if (translation.x() < 0.01 && translation.x() > -0.01) { // test for zero x foil coordinate
		vertex1type.push_back(0);
		// get the line direction at the foil vertex position
		if (the_particle.has_trajectory()) {
		  const snemo::datamodel::tracker_trajectory & the_trajectory = the_particle.get_trajectory();
		  const snemo::datamodel::base_trajectory_pattern & the_base_pattern = the_trajectory.get_pattern();
		  if (the_base_pattern.get_pattern_id()=="line") {
		    const geomtools::line_3d & the_shape = (const geomtools::line_3d&)the_base_pattern.get_shape();
		    geomtools::vector_3d direction = the_shape.get_direction_on_curve(the_shape.get_first());
		    foil1_dir_x.push_back(direction.x());
		    foil1_dir_y.push_back(direction.y());
		    foil1_dir_z.push_back(direction.z());
		  }
		  else {
		    const geomtools::helix_3d & the_shape = (const geomtools::helix_3d&)the_base_pattern.get_shape();
		    geomtools::vector_3d direction = the_shape.get_direction_on_curve(the_shape.get_first());
		    foil1_dir_x.push_back(direction.x());
		    foil1_dir_y.push_back(direction.y());
		    foil1_dir_z.push_back(direction.z());
		  }
		}
                else {
                  foil1_dir_x.push_back(0.);
                  foil1_dir_y.push_back(0.);
                  foil1_dir_z.push_back(0.);
                }
	      }
	      else {// vertex is on wire (i.e. neither calo, nor foil)
		vertex1type.push_back(1);
                foil1_dir_x.push_back(0.);
                foil1_dir_y.push_back(0.);
                foil1_dir_z.push_back(0.);
              }
	    }
            else {
              vertex1type.push_back(-1);
	      vertex1_x.push_back(0.);
	      vertex1_y.push_back(0.);
              vertex1_z.push_back(0.);
              foil1_dir_x.push_back(0.);
              foil1_dir_y.push_back(0.);
              foil1_dir_z.push_back(0.);
            }
            i = 1;
	    if (i<the_particle.get_vertices().size()) {
	      const geomtools::blur_spot & vertex = the_particle.get_vertices().at(i).get();
	      const geomtools::vector_3d & translation  = vertex.get_placement().get_translation();
	      vertex2_x.push_back(translation.x());
	      vertex2_y.push_back(translation.y());
	      vertex2_z.push_back(translation.z());
	      if (vertex.has_geom_id()) { // vertex is on calorimeter (calo or xcalo or gveto)
                vertex2type.push_back(vertex.get_geom_id().get_type());
                foil2_dir_x.push_back(0.);
                foil2_dir_y.push_back(0.);
                foil2_dir_z.push_back(0.);
	      }
	      else if (translation.x() < 0.01 && translation.x() > -0.01) { // test for zero x foil coordinate
		vertex2type.push_back(0);
		// get the line direction at the foil vertex position
		if (the_particle.has_trajectory()) {
		  const snemo::datamodel::tracker_trajectory & the_trajectory = the_particle.get_trajectory();
		  const snemo::datamodel::base_trajectory_pattern & the_base_pattern = the_trajectory.get_pattern();
		  if (the_base_pattern.get_pattern_id()=="line") {
		    const geomtools::line_3d & the_shape = (const geomtools::line_3d&)the_base_pattern.get_shape();
		    geomtools::vector_3d direction = the_shape.get_direction_on_curve(the_shape.get_first());
		    foil2_dir_x.push_back(direction.x());
		    foil2_dir_y.push_back(direction.y());
		    foil2_dir_z.push_back(direction.z());
		  }
		  else {
		    const geomtools::helix_3d & the_shape = (const geomtools::helix_3d&)the_base_pattern.get_shape();
		    geomtools::vector_3d direction = the_shape.get_direction_on_curve(the_shape.get_first());
		    foil2_dir_x.push_back(direction.x());
		    foil2_dir_y.push_back(direction.y());
		    foil2_dir_z.push_back(direction.z());
		  }
		}
                else {
                  foil2_dir_x.push_back(0.);
                  foil2_dir_y.push_back(0.);
                  foil2_dir_z.push_back(0.);
                }
	      }
	      else {// vertex is on wire (i.e. neither calo, nor foil)
		vertex2type.push_back(1);
                foil2_dir_x.push_back(0.);
                foil2_dir_y.push_back(0.);
                foil2_dir_z.push_back(0.);
              }
	    }
            else {
              vertex2type.push_back(-1);
	      vertex2_x.push_back(0.);
	      vertex2_y.push_back(0.);
              vertex2_z.push_back(0.);
              foil2_dir_x.push_back(0.);
              foil2_dir_y.push_back(0.);
              foil2_dir_z.push_back(0.);
            }
	  }
          else {
            vertex1type.push_back(-1);
            vertex1_x.push_back(0.);
            vertex1_y.push_back(0.);
            vertex1_z.push_back(0.);
            foil1_dir_x.push_back(0.);
            foil1_dir_y.push_back(0.);
            foil1_dir_z.push_back(0.);
            vertex2type.push_back(-1);
            vertex2_x.push_back(0.);
            vertex2_y.push_back(0.);
            vertex2_z.push_back(0.);
            foil2_dir_x.push_back(0.);
            foil2_dir_y.push_back(0.);
            foil2_dir_z.push_back(0.);
          }

          // then associated calorimeter hits
          if (the_particle.has_associated_calorimeter_hits()) {
            unsigned int i = 0;
            if (i<the_particle.get_associated_calorimeter_hits().size()) {
              const snemo::datamodel::calibrated_calorimeter_hit & calo_hit = the_particle.get_associated_calorimeter_hits().at(i).get();
	      const geomtools::mapping & the_mapping = geometry_manager_->get_mapping();
	      if (! the_mapping.validate_id(calo_hit.get_geom_id())) {
		std::vector<geomtools::geom_id> gids;
		the_mapping.compute_matching_geom_id(calo_hit.get_geom_id(), gids); // front calo block = last entry
		const geomtools::geom_info & info = the_mapping.get_geom_info(gids.back()); // in vector gids
		const geomtools::vector_3d & loc  = info.get_world_placement().get_translation();
		calo1_loc_x.push_back(loc.x());
		calo1_loc_y.push_back(loc.y());
		calo1_loc_z.push_back(loc.z());
	      }
	      else {
		const geomtools::geom_info & info = the_mapping.get_geom_info(calo_hit.get_geom_id());
		const geomtools::vector_3d & loc  = info.get_world_placement().get_translation();
		calo1_loc_x.push_back(loc.x());
		calo1_loc_y.push_back(loc.y());
		calo1_loc_z.push_back(loc.z());
	      }
	      calo1associated.push_back(1);
	      calo1energy.push_back(calo_hit.get_energy());
	      calo1sigmaenergy.push_back(calo_hit.get_sigma_energy());
	      calo1time.push_back(calo_hit.get_time());
	      calo1sigmatime.push_back(calo_hit.get_sigma_time());
	      calo1type.push_back(calo_hit.get_geom_id().get_type());

	      if (calo_hit.get_geom_id ().get_type () == 1302)
		{
		  // CALO
		  calo1side.push_back(calo_hit.get_geom_id().get(1));
		  calo1column.push_back(calo_hit.get_geom_id().get(2));
		  calo1row.push_back(calo_hit.get_geom_id().get(3));
		  calo1wall.push_back(-1);
		}
              else if (calo_hit.get_geom_id ().get_type () == 1232)
		{
		  // XCALO
		  calo1side.push_back(calo_hit.get_geom_id().get(1));
		  calo1column.push_back(calo_hit.get_geom_id().get(3));
		  calo1wall.push_back(calo_hit.get_geom_id().get(2));
		  calo1row.push_back(calo_hit.get_geom_id().get(0));
		}
              else if (calo_hit.get_geom_id ().get_type () == 1252)
		{
		  // GVETO
		  calo1side.push_back(calo_hit.get_geom_id().get(1));
		  calo1column.push_back(calo_hit.get_geom_id().get(3));
		  calo1wall.push_back(calo_hit.get_geom_id().get(2));
		  calo1row.push_back(-1);
		}
              else {
		  calo1side.push_back(-1);
		  calo1column.push_back(-1);
		  calo1wall.push_back(-1);
		  calo1row.push_back(-1);
              }
	    }
            else {
              calo1associated.push_back(-1); // non particle
              calo1energy.push_back(-1.0);
              calo1sigmaenergy.push_back(-1.0);
              calo1time.push_back(-1.0);
              calo1sigmatime.push_back(-1.0);
              calo1type.push_back(-1);
              calo1side.push_back(-1);
              calo1column.push_back(-1);
              calo1wall.push_back(-1);
              calo1row.push_back(-1);
              calo1_loc_x.push_back(1.0e3);
              calo1_loc_y.push_back(1.0e4);
              calo1_loc_z.push_back(1.0e4);
            }
            i = 1;
            if (i<the_particle.get_associated_calorimeter_hits().size()) {
              const snemo::datamodel::calibrated_calorimeter_hit & calo_hit = the_particle.get_associated_calorimeter_hits().at(i).get();
	      const geomtools::mapping & the_mapping = geometry_manager_->get_mapping();
	      if (! the_mapping.validate_id(calo_hit.get_geom_id())) {
		std::vector<geomtools::geom_id> gids;
		the_mapping.compute_matching_geom_id(calo_hit.get_geom_id(), gids); // front calo block = last entry
		const geomtools::geom_info & info = the_mapping.get_geom_info(gids.back()); // in vector gids
		const geomtools::vector_3d & loc  = info.get_world_placement().get_translation();
		calo2_loc_x.push_back(loc.x());
		calo2_loc_y.push_back(loc.y());
		calo2_loc_z.push_back(loc.z());
	      }
	      else {
		const geomtools::geom_info & info = the_mapping.get_geom_info(calo_hit.get_geom_id());
		const geomtools::vector_3d & loc  = info.get_world_placement().get_translation();
		calo2_loc_x.push_back(loc.x());
		calo2_loc_y.push_back(loc.y());
		calo2_loc_z.push_back(loc.z());
	      }
	      calo2associated.push_back(1);
	      calo2energy.push_back(calo_hit.get_energy());
	      calo2sigmaenergy.push_back(calo_hit.get_sigma_energy());
	      calo2time.push_back(calo_hit.get_time());
	      calo2sigmatime.push_back(calo_hit.get_sigma_time());
	      calo2type.push_back(calo_hit.get_geom_id().get_type());

	      if (calo_hit.get_geom_id ().get_type () == 1302)
		{
		  // CALO
		  calo2side.push_back(calo_hit.get_geom_id().get(1));
		  calo2column.push_back(calo_hit.get_geom_id().get(2));
		  calo2row.push_back(calo_hit.get_geom_id().get(3));
		  calo2wall.push_back(-1);
		}
              else if (calo_hit.get_geom_id ().get_type () == 1232)
		{
		  // XCALO
		  calo2side.push_back(calo_hit.get_geom_id().get(1));
		  calo2column.push_back(calo_hit.get_geom_id().get(3));
		  calo2wall.push_back(calo_hit.get_geom_id().get(2));
		  calo2row.push_back(calo_hit.get_geom_id().get(0));
		}
              else if (calo_hit.get_geom_id ().get_type () == 1252)
		{
		  // GVETO
		  calo2side.push_back(calo_hit.get_geom_id().get(1));
		  calo2column.push_back(calo_hit.get_geom_id().get(3));
		  calo2wall.push_back(calo_hit.get_geom_id().get(2));
		  calo2row.push_back(-1);
		}
              else {
		  calo2side.push_back(-1);
		  calo2column.push_back(-1);
		  calo2wall.push_back(-1);
		  calo2row.push_back(-1);
              }
	    }
            else {
              calo2associated.push_back(-1); // non particle
              calo2energy.push_back(-1.0);
              calo2sigmaenergy.push_back(-1.0);
              calo2time.push_back(-1.0);
              calo2sigmatime.push_back(-1.0);
              calo2type.push_back(-1);
              calo2side.push_back(-1);
              calo2column.push_back(-1);
              calo2wall.push_back(-1);
              calo2row.push_back(-1);
              calo2_loc_x.push_back(1.0e3);
              calo2_loc_y.push_back(1.0e4);
              calo2_loc_z.push_back(1.0e4);
            }
	  }
	  else { // not a calo vertex, fill blanks
	    calo1associated.push_back(-1); // non particle
	    calo1energy.push_back(-1.0);
	    calo1sigmaenergy.push_back(-1.0);
	    calo1time.push_back(-1.0);
	    calo1sigmatime.push_back(-1.0);
	    calo1type.push_back(-1);
	    calo1side.push_back(-1);
	    calo1column.push_back(-1);
	    calo1wall.push_back(-1);
	    calo1row.push_back(-1);
	    calo1_loc_x.push_back(1.0e3);
	    calo1_loc_y.push_back(1.0e4);
	    calo1_loc_z.push_back(1.0e4);
	    calo2associated.push_back(-1); // non particle
	    calo2energy.push_back(-1.0);
	    calo2sigmaenergy.push_back(-1.0);
	    calo2time.push_back(-1.0);
	    calo2sigmatime.push_back(-1.0);
	    calo2type.push_back(-1);
	    calo2side.push_back(-1);
	    calo2column.push_back(-1);
	    calo2wall.push_back(-1);
	    calo2row.push_back(-1);
	    calo2_loc_x.push_back(1.0e3);
	    calo2_loc_y.push_back(1.0e4);
	    calo2_loc_z.push_back(1.0e4);
	  }


	  // then the trajectory length, always
	  if (the_particle.has_trajectory()) {
	    const snemo::datamodel::tracker_trajectory & the_trajectory = the_particle.get_trajectory();
	    const snemo::datamodel::base_trajectory_pattern & the_base_pattern = the_trajectory.get_pattern();
	    const geomtools::i_shape_1d & the_shape = the_base_pattern.get_shape();
	    traj_length.push_back(the_shape.get_length());

	    const snemo::datamodel::tracker_cluster & the_cluster = the_trajectory.get_cluster();
	    traj_cl_delayed.push_back((int)the_cluster.is_delayed());
	    if (the_cluster.is_delayed()>0)
	      traj_cl_delayed_time.push_back(the_cluster.get_hit(0).get_delayed_time());
	    else
	      traj_cl_delayed_time.push_back(-1.0);
	  }
	  else {
	    traj_length.push_back(-1.0);
	    traj_cl_delayed.push_back(-1);
	    traj_cl_delayed_time.push_back(-1.0);
	  }
	}
      }
      else
	particle_.nofparticles_ = 0;

      // Extract un-associated calorimeter hits
      if (PTD.has_non_associated_calorimeters ()) {
	particle_.nofgammas_ = PTD.get_non_associated_calorimeters().size();

	for (unsigned int i=0; i<PTD.get_non_associated_calorimeters().size();++i) {
	  const snemo::datamodel::calibrated_calorimeter_hit & un_calo_hit = PTD.get_non_associated_calorimeters().at(i).get();

 	  const geomtools::mapping & the_mapping = geometry_manager_->get_mapping();
	  if (! the_mapping.validate_id(un_calo_hit.get_geom_id())) {
	    std::vector<geomtools::geom_id> gids;
	    the_mapping.compute_matching_geom_id(un_calo_hit.get_geom_id(), gids); // front calo block = last entry
	    const geomtools::geom_info & info = the_mapping.get_geom_info(gids.back()); // in vector gids
	    const geomtools::vector_3d & loc  = info.get_world_placement().get_translation();
	    calo1_loc_x.push_back(loc.x());
	    calo1_loc_y.push_back(loc.y());
	    calo1_loc_z.push_back(loc.z());
	  }
	  else {
	    const geomtools::geom_info & info = the_mapping.get_geom_info(un_calo_hit.get_geom_id());
	    const geomtools::vector_3d & loc  = info.get_world_placement().get_translation();
	    calo1_loc_x.push_back(loc.x());
	    calo1_loc_y.push_back(loc.y());
	    calo1_loc_z.push_back(loc.z());
	  }

	  calo1associated.push_back(0);
	  calo1energy.push_back(un_calo_hit.get_energy());
	  calo1sigmaenergy.push_back(un_calo_hit.get_sigma_energy());
	  calo1time.push_back(un_calo_hit.get_time());
	  calo1sigmatime.push_back(un_calo_hit.get_sigma_time());
	  calo1type.push_back(un_calo_hit.get_geom_id().get_type());

	  if (un_calo_hit.get_geom_id ().get_type () == 1302)
	    {
	      // CALO
	      calo1side.push_back(un_calo_hit.get_geom_id().get(1));
	      calo1column.push_back(un_calo_hit.get_geom_id().get(2));
	      calo1row.push_back(un_calo_hit.get_geom_id().get(3));
	      calo1wall.push_back(-1);
	    }
          else if (un_calo_hit.get_geom_id ().get_type () == 1232)
	    {
	      // XCALO
	      calo1side.push_back(un_calo_hit.get_geom_id().get(1));
	      calo1column.push_back(un_calo_hit.get_geom_id().get(3));
	      calo1wall.push_back(un_calo_hit.get_geom_id().get(2));
	      calo1row.push_back(un_calo_hit.get_geom_id().get(0));
	    }
          else if (un_calo_hit.get_geom_id ().get_type () == 1252)
	    {
	      // GVETO
	      calo1side.push_back(un_calo_hit.get_geom_id().get(1));
	      calo1column.push_back(un_calo_hit.get_geom_id().get(3));
	      calo1wall.push_back(un_calo_hit.get_geom_id().get(2));
	      calo1row.push_back(-1);
	    }
          else {
          }
	  // calo unassociated, append blanks
	  particleid.push_back(-1); // not a particle
	  charge.push_back(-1);
	  vertex1_x.push_back(0.); // calo only
	  vertex1_y.push_back(0.);
	  vertex1_z.push_back(0.);
	  foil1_dir_x.push_back(0.); // double fill for vertex
	  foil1_dir_y.push_back(0.);
	  foil1_dir_z.push_back(0.);
	  vertex1type.push_back(-1);
	  vertex2_x.push_back(0.); // calo only
	  vertex2_y.push_back(0.);
	  vertex2_z.push_back(0.);
	  foil2_dir_x.push_back(0.); // double fill for vertex
	  foil2_dir_y.push_back(0.);
	  foil2_dir_z.push_back(0.);
	  vertex2type.push_back(-1);
	  traj_length.push_back(-1.0);
	  traj_cl_delayed.push_back(-1);
          traj_cl_delayed_time.push_back(-1.0);
          calo2associated.push_back(-1); // non particle
          calo2energy.push_back(-1.0);
          calo2sigmaenergy.push_back(-1.0);
          calo2time.push_back(-1.0);
          calo2sigmatime.push_back(-1.0);
          calo2type.push_back(-1);
          calo2side.push_back(-1);
          calo2column.push_back(-1);
          calo2wall.push_back(-1);
          calo2row.push_back(-1);
          calo2_loc_x.push_back(1.0e3);
          calo2_loc_y.push_back(1.0e4);
          calo2_loc_z.push_back(1.0e4);
	}
      }
      else {
	particle_.nofgammas_ = 0;
      }

      particle_.particle_id_ = &particleid;
      particle_.charge_ = &charge;
      particle_.vertex1_type_ = &vertex1type;
      particle_.vertex1_x_ = &vertex1_x;
      particle_.vertex1_y_ = &vertex1_y;
      particle_.vertex1_z_ = &vertex1_z;
      particle_.foil1_dir_x_ = &foil1_dir_x;
      particle_.foil1_dir_y_ = &foil1_dir_y;
      particle_.foil1_dir_z_ = &foil1_dir_z;
      particle_.vertex2_type_ = &vertex2type;
      particle_.vertex2_x_ = &vertex2_x;
      particle_.vertex2_y_ = &vertex2_y;
      particle_.vertex2_z_ = &vertex2_z;
      particle_.foil2_dir_x_ = &foil2_dir_x;
      particle_.foil2_dir_y_ = &foil2_dir_y;
      particle_.foil2_dir_z_ = &foil2_dir_z;
      particle_.traj_length_ = &traj_length;
      particle_.traj_cluster_delayed_ = &traj_cl_delayed;
      particle_.traj_cluster_delayed_time_ = &traj_cl_delayed_time;
      particle_.calo1_associated_ = &calo1associated;
      particle_.calo1_type_ = &calo1type;
      particle_.calo1_energy_ = &calo1energy;
      particle_.calo1_sigma_energy_ = &calo1sigmaenergy;
      particle_.calo1_time_ = &calo1time;
      particle_.calo1_sigma_time_ = &calo1sigmatime;
      particle_.calo1_side_ = &calo1side;
      particle_.calo1_column_ = &calo1column;
      particle_.calo1_row_ = &calo1row;
      particle_.calo1_wall_ = &calo1wall;
      particle_.calo1_loc_x_ = &calo1_loc_x;
      particle_.calo1_loc_y_ = &calo1_loc_y;
      particle_.calo1_loc_z_ = &calo1_loc_z;
      particle_.calo2_associated_ = &calo2associated;
      particle_.calo2_type_ = &calo2type;
      particle_.calo2_energy_ = &calo2energy;
      particle_.calo2_sigma_energy_ = &calo2sigmaenergy;
      particle_.calo2_time_ = &calo2time;
      particle_.calo2_sigma_time_ = &calo2sigmatime;
      particle_.calo2_side_ = &calo2side;
      particle_.calo2_column_ = &calo2column;
      particle_.calo2_row_ = &calo2row;
      particle_.calo2_wall_ = &calo2wall;
      particle_.calo2_loc_x_ = &calo2_loc_x;
      particle_.calo2_loc_y_ = &calo2_loc_y;
      particle_.calo2_loc_z_ = &calo2_loc_z;

    }

  // look for event header
  if(workItem.has("EH"))
    {
      const snemo::datamodel::event_header & EH = workItem.get<snemo::datamodel::event_header>("EH");
      //      std::cout << "In process: found EH event header " << std::endl;
      header_.runnumber_ = EH.get_id ().get_run_number ();
      header_.eventnumber_ = EH.get_id ().get_event_number ();
      header_.date_ = 0;
      header_.runtype_ = 0;
      header_.simulated_ = (EH.is_simulated () ? true : false);
    }

  if(workItem.has("SD")) {
    const mctools::simulated_data& SD = workItem.get<mctools::simulated_data>("SD");
    gen_.vertex_x_ = SD.get_vertex().x();
    gen_.vertex_y_ = SD.get_vertex().y();
    gen_.vertex_z_ = SD.get_vertex().z();
  }

  tree_->Fill();

  // MUST return a status, see ref dpp::processing_status_flags_type
  return dpp::base_module::PROCESS_OK;
}

// Reset
void MyPTD2Root::reset() {
  // write the output, finished streaming
  hfile_->cd();
  tree_->Write();
  hfile_->Close(); //
  std::cout << "In reset: finished conversion, file closed " << std::endl;

  // clean up
  delete hfile_;
  filename_output_ = "default.root";
  this->_set_initialized(false);
}
