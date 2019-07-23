/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <limits>

#include "Tudat/JsonInterface/jsonInterface.h"
#include "Tudat/SimulationSetup/tudatEstimationHeader.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"

using namespace tudat;
using namespace tudat::json_interface;
using namespace tudat::ephemerides;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::propagators;
using namespace tudat::coordinate_conversions;
using namespace tudat::reference_frames;
using namespace tudat::observation_models;
using namespace tudat::estimatable_parameters;
using namespace tudat::orbit_determination;

//! Get path for output directory.
static inline std::string getOutputPath( const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) -
                                std::string( "dopplerDataSimulation.cpp" ).length( ) );
    std::string outputPath = reducedPath + "SimulationOutput/";
    if( extraDirectory != "" )
    {
        outputPath += extraDirectory;
    }

    if( outputPath.at( outputPath.size( ) - 1 ) != '/' )
    {
        outputPath += "/";
    }

    return outputPath;
}

//! This function reads the main.json file, and runs the Doppler data simulation based on teh settings provided there.
//! Simulated Doppler shifts are provided in an output file. Data is simulated for a single ground station tracking a single
//! spacecraft, at all epochs where the spacecraft is visible (defined by a user-provided elevation angle cutoff, and data cadence)
int main( )
{
    // Load settings file
    std::string inputFile = "/home/dominic/Software/numericalAstrodynamicsTudatBundle/tudatBundle/tudatApplications/Flaress/main.json";
    JsonSimulationManager< > jsonSimulation( inputFile );
    jsonSimulation.updateSettings( );

    // Retrieve and update physical environment.
    NamedBodyMap bodyMap = jsonSimulation.getBodyMap( );
    bodyMap[ "Spacecraft" ]->setEphemeris( std::make_shared< TabulatedCartesianEphemeris< > >(
                                               std::shared_ptr< interpolators::OneDimensionalInterpolator
                                               < double, Eigen::Vector6d > >( ), "Earth", "J2000" ) );
    {
        Eigen::Vector3d nominalPosition = getValue< Eigen::Vector3d >(
                    jsonSimulation.getOriginalJsonObject( ), "stationPosition" );

        createGroundStation( bodyMap.at( "Earth" ), "GroundStation", nominalPosition, cartesian_position );
    }

    // Retrieve propagation/integration settings and simulate dynamics
    std::shared_ptr< IntegratorSettings< > > integratorSettings = jsonSimulation.getIntegratorSettings( );
    std::shared_ptr< PropagatorSettings< > > propagatorSettings = jsonSimulation.getPropagatorSettings( );
    SingleArcDynamicsSimulator< > dynamicsSimulator =
            SingleArcDynamicsSimulator< >( bodyMap, integratorSettings, propagatorSettings, true, false, true );

    // Define one-way Doppler link
    LinkEnds oneWayLinkEnds;
    oneWayLinkEnds[ transmitter ] = std::make_pair( "Earth", "GroundStation" );
    oneWayLinkEnds[ receiver ] = std::make_pair( "Spacecraft", "" );
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsList;
    linkEndsList[ one_way_doppler ].push_back( oneWayLinkEnds );

    // Define nominal times to simulation observation
    typedef std::map< ObservableType, std::map< LinkEnds, std::shared_ptr< ObservationSettings > > > SortedObservationSettingsMap;
    SortedObservationSettingsMap observationSettingsList;
    std::map< ObservableType, std::map< LinkEnds, std::shared_ptr< ObservationSimulationTimeSettings< double > > > >
            observationTimeSettings;
    std::vector< double > unconstrainedObservationTimes;
    double initialTime = dynamicsSimulator.getEquationsOfMotionNumericalSolution( ).begin( )->first + 100.0,
            finalTime = dynamicsSimulator.getEquationsOfMotionNumericalSolution( ).rbegin( )->first - 100.0,
            timeStep = getValue< double >(
                jsonSimulation.getOriginalJsonObject( ), "dataCadence" );
    double currentTime = initialTime;
    while( currentTime <= finalTime )
    {
        unconstrainedObservationTimes.push_back( currentTime );
        currentTime += timeStep;
    }

    // Define settings for observable
    std::vector< std::string > perturbingBodies = { "Earth", "Sun" };
    observationSettingsList[ one_way_doppler ][ oneWayLinkEnds ] = std::make_shared<
            OneWayDopplerObservationSettings >
            (  std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ),
               std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
               std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ) );
    observationTimeSettings[ one_way_doppler ][ oneWayLinkEnds ] = std::make_shared<
            TabulatedObservationSimulationTimeSettings< double > >(
                transmitter, unconstrainedObservationTimes );

    // Define minimum elevation angles for Earth/Mars stations
    double earthMinimumElevationAngle = getValue< double >(
                jsonSimulation.getOriginalJsonObject( ), "elevationCutoff" ) * mathematical_constants::PI / 180.0;

    // Create observation viability calculators
    std::vector< std::shared_ptr< ObservationViabilitySettings > > observationViabilitySettings;
    observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                minimum_elevation_angle, std::make_pair( "Earth", "" ), "",
                                                earthMinimumElevationAngle ) );
    PerObservableObservationViabilityCalculatorList viabilityCalculators = createObservationViabilityCalculators(
                bodyMap, linkEndsList, observationViabilitySettings );

    // Create observation simulators
    std::map< ObservableType,  std::shared_ptr< ObservationSimulatorBase< double, double > > > observationSimulators =
            createObservationSimulators( observationSettingsList , bodyMap );

    // Simulate observations with viability constraints directly from simulateObservations function
    std::pair< Eigen::VectorXd, std::vector< double > >
            constrainedSimulatedObservables = removeLinkIdFromSimulatedObservations(
                simulateObservations( observationTimeSettings, observationSimulators, viabilityCalculators ) ).at(
                one_way_doppler ).at( oneWayLinkEnds );

    // Retrieve
    Eigen::MatrixXd simulationOutput = Eigen::MatrixXd( constrainedSimulatedObservables.first.rows( ), 2 );
    simulationOutput.block( 0, 1, constrainedSimulatedObservables.first.rows( ), 1 ) =
            constrainedSimulatedObservables.first;
    simulationOutput.block( 0, 0, constrainedSimulatedObservables.first.rows( ), 1 ) =
            utilities::convertStlVectorToEigenVector(
                constrainedSimulatedObservables.second );

    input_output::writeMatrixToFile( simulationOutput,
                                     getValue< std::string >( jsonSimulation.getOriginalJsonObject( ), "outputFile" ), 16,
                                     getOutputPath( ) );

}


