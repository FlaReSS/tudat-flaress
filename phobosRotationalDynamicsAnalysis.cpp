/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN


#include <limits>


#include "Tudat/Astrodynamics/Ephemerides/keplerEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/constantRotationalEphemeris.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/Astrodynamics/Propagators/getZeroProperModeRotationalInitialState.h"
#include "Tudat/SimulationSetup/tudatEstimationHeader.h"
#include "Tudat/InputOutput/readHistoryFromFile.h"

//Using declarations.
using namespace tudat;
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


NamedBodyMap getTestBodyMap( const double phobosSemiMajorAxis,
                             const bool useSymmetricEquator = 0,
                             double finalEphemerisTime = 0.0 )
{
    // Define time range of test.
    double initialEphemerisTime = 0.0;

    NamedBodyMap bodyMap = createBodies(
                getDefaultBodySettings( { "Sun", "Jupiter", "Earth" },
                                        initialEphemerisTime - 100.0 * 86400.0, finalEphemerisTime + 100.0 * 86400.0 ) );

    bodyMap[ "Mars" ] = std::make_shared< Body >( );
    bodyMap[ "Mars" ]->setEphemeris( createBodyEphemeris(
                                         getDefaultEphemerisSettings(
                                             "Mars", initialEphemerisTime - 86400.0, finalEphemerisTime + 86400.0 ), "Mars" ) );
    std::shared_ptr< SphericalHarmonicsGravityFieldSettings > marsGravityFieldSettings =
            std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >(
                getDefaultGravityFieldSettings( "Mars", TUDAT_NAN, TUDAT_NAN ) );
    marsGravityFieldSettings->resetGravitationalParameter( spice_interface::getBodyGravitationalParameter( "Mars" ) );
    bodyMap[ "Mars" ]->setGravityFieldModel(
                createGravityFieldModel( marsGravityFieldSettings, "Mars", bodyMap ) );
    bodyMap[ "Mars" ]->setShapeModel(
                createBodyShapeModel( getDefaultBodyShapeSettings( "Mars", TUDAT_NAN, TUDAT_NAN ), "Mars" ) );
    //                std::make_shared< gravitation::GravityFieldModel >(
    //                    spice_interface::getBodyGravitationalParameter( "Mars" ) ) );

    Eigen::Vector7d fixedMarsRotation = Eigen::Vector7d::Zero( );
    fixedMarsRotation( 0 ) = 1.0;
    bodyMap[ "Mars" ]->setRotationalEphemeris(
                //                std::make_shared< ConstantRotationalEphemeris >(
                //                    fixedMarsRotation, "ECLIPJ2000", "IAU_Mars " ) );
                createRotationModel( getDefaultRotationModelSettings(
                                         "Mars", TUDAT_NAN, TUDAT_NAN ), "Mars" ) );

    bodyMap[ "Phobos" ] = std::make_shared< Body >( );
    std::shared_ptr< EphemerisSettings > phobosEphemerisSettings =
            getDefaultEphemerisSettings( "Phobos" );
    phobosEphemerisSettings->resetFrameOrigin( "Mars" );
    bodyMap[ "Phobos" ]->setEphemeris(
                getTabulatedEphemeris( createBodyEphemeris( phobosEphemerisSettings, "Phobos" ),
                                       -3600.0, finalEphemerisTime + 3600.0, 120.0,
                                       std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ), "Mars" ) ) ;

    bodyMap[ "Phobos" ]->setRotationalEphemeris(
                //                std::make_shared< ConstantRotationalEphemeris >(
                //                    fixedMarsRotation, "ECLIPJ2000", "IAU_Mars " ) );
                getTabulatedRotationalEphemeris(
                    createRotationModel( getDefaultRotationModelSettings(
                                             "Phobos", TUDAT_NAN, TUDAT_NAN ), "Phobos" ),
                    -3600.0, finalEphemerisTime + 3600.0, 120.0,
                    std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) ) );
    Eigen::Matrix3d phobosInertiaTensor = Eigen::Matrix3d::Zero( );
    phobosInertiaTensor( 0, 0 ) = 0.3615;
    phobosInertiaTensor( 1, 1 ) = 0.4265;
    phobosInertiaTensor( 2, 2 ) = 0.5024;

    if( useSymmetricEquator )
    {
        phobosInertiaTensor( 0, 0 ) = phobosInertiaTensor( 1, 1 );
    }

    double phobosReferenceRadius = 11.27E3;
    double phobosMass = 1.0659E16;

    phobosInertiaTensor *= (phobosReferenceRadius * phobosReferenceRadius * phobosMass );
    bodyMap[ "Phobos" ]->setBodyInertiaTensor( phobosInertiaTensor );
    bodyMap[ "Phobos" ]->setShapeModel(
                std::make_shared< SphericalBodyShapeModel >( 15.0E3 ) );

    double phobosGravitationalParameter = phobosMass * physical_constants::GRAVITATIONAL_CONSTANT;

    Eigen::MatrixXd phobosCosineGravityFieldCoefficients = Eigen::Matrix4d::Zero( ),
            phobosSineGravityFieldCoefficients = Eigen::Matrix4d::Zero( );
    double phobosScaledMeanMomentOfInertia;
    gravitation::getDegreeTwoSphericalHarmonicCoefficients(
                phobosInertiaTensor, phobosGravitationalParameter, phobosReferenceRadius, true,
                phobosCosineGravityFieldCoefficients, phobosSineGravityFieldCoefficients, phobosScaledMeanMomentOfInertia );

    phobosCosineGravityFieldCoefficients( 3, 0 ) = 0.01480 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 3, 0 );
    phobosCosineGravityFieldCoefficients( 3, 1 ) = 0.00785 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 3, 1 );
    phobosCosineGravityFieldCoefficients( 3, 2 ) = -0.00272 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 3, 2 );
    phobosCosineGravityFieldCoefficients( 3, 3 ) = -0.00024 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 3, 3 );

    phobosSineGravityFieldCoefficients( 3, 1 ) = -0.00186 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 3, 1 );
    phobosSineGravityFieldCoefficients( 3, 2 ) = 0.00004 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 3, 2 );
    phobosSineGravityFieldCoefficients( 3, 3 ) = 0.0013 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 3, 3 );

    bodyMap[ "Phobos" ]->setGravityFieldModel(
                std::make_shared< gravitation::SphericalHarmonicsGravityField >(
                    phobosGravitationalParameter, phobosReferenceRadius, phobosCosineGravityFieldCoefficients,
                    phobosSineGravityFieldCoefficients, "IAU_Phobos",
                    std::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodyMap.at( "Phobos" ), true ) ) );

    std::cout<<"Cosine "<<std::endl<<phobosCosineGravityFieldCoefficients<<std::endl;
    std::cout<<"Sine "<<std::endl<<phobosSineGravityFieldCoefficients<<std::endl;
    std::cout<<"MOI "<<std::endl<<bodyMap.at( "Phobos" )->getBodyInertiaTensor( )<<std::endl;
    std::cout<<"MOI "<<std::endl<<bodyMap.at( "Phobos" )->getScaledMeanMomentOfInertia( )<<std::endl;

    bodyMap[ "Deimos" ] = std::make_shared< Body >( );
    std::shared_ptr< EphemerisSettings > deimosEphemerisSettings =
            getDefaultEphemerisSettings( "Deimos" );
    deimosEphemerisSettings->resetFrameOrigin( "Mars" );
    bodyMap[ "Deimos" ]->setEphemeris( createBodyEphemeris( deimosEphemerisSettings, "Deimos" ) );
    bodyMap[ "Deimos" ]->setGravityFieldModel(
                std::make_shared< gravitation::GravityFieldModel >( physical_constants::GRAVITATIONAL_CONSTANT * 1.4762E15 ) );

    return bodyMap;
}

int main( )
{
    bool propagateDynamicsOnly = false;
    bool propagateVariationalEquations = false;
    bool performPerturbationAnalysis = false;
    bool findDampedInitialState = false;
    bool performParameterPerturbationAnalysis = true;


    double printInterval = 365.25 / 4.0 * 86400.0;

    std::string dataFolderBase = "/home/dominic/Documents/Articles/PhobosCoupledDynamics/Data/BehaviourAnalysisParameters/";
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "meta_PhobosPCK.tm" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "NOE-4-2015-b.bsp" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de430.bsp" );

    // Define time range of test.
    int numberOfWeeks = 26;
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + numberOfWeeks * 7 * 86400.0;

    // Retrieve list of body objects.
    NamedBodyMap bodyMap = getTestBodyMap( 9376.0E3, 0, finalEphemerisTime );
    createGroundStation( bodyMap.at( "Phobos" ), "Lander", ( Eigen::Vector3d( ) << 0.1, 0.35, 0.0 ).finished( ), geodetic_position );
    setGlobalFrameBodyEphemerides( bodyMap, "Mars", "ECLIPJ2000" );

    for( unsigned int propagationType = 2; propagationType < 3; propagationType++ )
    {
        std::string propagationTypeFolder;

        if( propagationType == 0 )
        {
            propagationTypeFolder = "Translational/";
        }
        else if( propagationType == 1 )
        {
            propagationTypeFolder = "Rotational/";
        }
        else if( propagationType == 2 )
        {
            propagationTypeFolder = "Coupled/";
        }


        for( unsigned int dynamicsType = 2; dynamicsType < 3; dynamicsType++ )
        {
            std::string dynamicsTypeFolder;

            if( dynamicsType == 0 )
            {
                dynamicsTypeFolder = "Keplerian/";
            }
            else if( dynamicsType == 1 )
            {
                dynamicsTypeFolder = "Simplified/";
            }
            else if( dynamicsType == 2 )
            {
                dynamicsTypeFolder = "Full/";
            }

            std::string dataFolder = dataFolderBase + propagationTypeFolder + dynamicsTypeFolder;

            std::string propagationDataFilePrefixBase = "Phobos_" + std::to_string( numberOfWeeks ) + "_weeks";

            for( unsigned int integrationType = 0; integrationType < 1; integrationType++ )
            {
                std::cout<<"RUNNING ANALYSIS FOR ****************************************** "<<" "<<propagationType<<" "<<
                           dynamicsType<<" "<<integrationType<<std::endl;

                // Set torques between bodies that are to be taken into account.
                std::vector< std::string > bodiesToIntegrate;
                bodiesToIntegrate.push_back( "Phobos" );

                // Define mean motion (equal to rotation rate).
                Eigen::VectorXd initialRotationalState = //Eigen::VectorXd( 7 );
                        //initialRotationalState<<0.7591829940812994, 0.1566361286612512, 0.1622275306326094, 0.6105641104754064, -1.849456448817688e-09, 2.899615136761079e-08, 0.0002323023538292387;
                        propagators::getInitialRotationalStateOfBody(
                            "Phobos", "ECLIPJ2000", bodyMap, initialEphemerisTime );

                // Create torque models
                SelectedTorqueMap torqueMap;
                if( dynamicsType > 0 )
                {
                    torqueMap[ "Phobos" ][ "Mars" ].push_back( std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );
                }
                else if( dynamicsType > 1 )
                {
                    torqueMap[ "Phobos" ][ "Sun" ].push_back( std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );
                    torqueMap[ "Phobos" ][ "Jupiter" ].push_back( std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );
                    torqueMap[ "Phobos" ][ "Deimos" ].push_back( std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );
                }

                // Create torque models
                basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap(
                            bodyMap, torqueMap, bodiesToIntegrate );

                // Define propagator settings variables.
                SelectedAccelerationMap accelerationMap;
                std::vector< std::string > translationalBodiesToPropagate;
                std::vector< std::string > translationalCentralBodies;

                // Define acceleration model settings.
                std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfPhobos;
                if( dynamicsType == 0 )
                {
                    accelerationsOfPhobos[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
                }
                else if( dynamicsType == 1 )
                {
                    accelerationsOfPhobos[ "Mars" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
                }
                else if( dynamicsType == 2 )
                {
                    accelerationsOfPhobos[ "Mars" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 2, 2 ) );
                }


                if( dynamicsType > 0 )
                {
                    accelerationsOfPhobos[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
                    accelerationsOfPhobos[ "Deimos" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
                    accelerationsOfPhobos[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
                }
                accelerationMap[ "Phobos" ] = accelerationsOfPhobos;

                translationalBodiesToPropagate.push_back( "Phobos" );
                translationalCentralBodies.push_back( "Mars" );

                // Create acceleration models
                basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                            bodyMap, accelerationMap, translationalBodiesToPropagate, translationalCentralBodies );

                // Define integrator settings.
                double timeStep = 60.0 * std::pow( 2.0, integrationType );

                std::string propagationDataFilePrefix = propagationDataFilePrefixBase + "_step_" +
                        std::to_string( static_cast< int >( std::round( timeStep ) ) );

                std::shared_ptr< IntegratorSettings< > > integratorSettings =
                        std::make_shared< RungeKuttaVariableStepSizeSettings< > >
                        ( rungeKuttaVariableStepSize, initialEphemerisTime, timeStep,
                          RungeKuttaCoefficients::rungeKuttaFehlberg78, timeStep, timeStep, 1.0, 1.0 );

                // Define list of dependent variables to save.
                std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
                dependentVariablesList.push_back(
                            std::make_shared< SingleDependentVariableSaveSettings >( keplerian_state_dependent_variable, "Phobos", "Mars" ) );

                if( dynamicsType > 0 || propagationType == 2 )
                {
                    dependentVariablesList.push_back(
                                std::make_shared< SingleDependentVariableSaveSettings >(
                                    body_fixed_relative_spherical_position, "Phobos", "Mars" ) );
                    dependentVariablesList.push_back(
                                std::make_shared< SingleDependentVariableSaveSettings >(
                                    body_fixed_relative_spherical_position, "Mars", "Phobos" ) );
                    dependentVariablesList.push_back(
                                std::make_shared< SingleDependentVariableSaveSettings >(
                                    euler_angles_to_body_fixed_313, "Phobos" ) );
                }
                // Create object with list of dependent variables
                std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
                        std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );




                // Define propagator settings.
                std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
                        std::make_shared< RotationalStatePropagatorSettings< double > >
                        ( torqueModelMap, bodiesToIntegrate, initialRotationalState, std::make_shared< PropagationTimeTerminationSettings >(
                              finalEphemerisTime ), dependentVariablesToSave, printInterval );

                Eigen::VectorXd initialTranslationalState = //Eigen::VectorXd( 6 );
                        //initialTranslationalState <<-1989892.877155822, -9287092.947239853, 558710.5104264524, 1843.221947897575, -444.6826783220159, -917.140452833435;
                        propagators::getInitialStatesOfBodies(
                            translationalBodiesToPropagate, translationalCentralBodies, bodyMap, initialEphemerisTime );
                std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                        std::make_shared< TranslationalStatePropagatorSettings< double > >
                        ( translationalCentralBodies, accelerationModelMap, translationalBodiesToPropagate,
                          initialTranslationalState, finalEphemerisTime, cowell, dependentVariablesToSave, printInterval );

                std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > >  propagatorSettingsList;
                propagatorSettingsList.push_back( translationalPropagatorSettings );
                propagatorSettingsList.push_back( rotationalPropagatorSettings );

                std::shared_ptr< MultiTypePropagatorSettings< double > > coupledPropagatorSettings =
                        std::make_shared< MultiTypePropagatorSettings< double > >(
                            propagatorSettingsList,
                            std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime ),
                            dependentVariablesToSave, printInterval );

                std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings;
                if( propagationType == 0 )
                {
                    propagatorSettings = translationalPropagatorSettings;
                }
                else if( propagationType == 1 )
                {
                    propagatorSettings = rotationalPropagatorSettings;
                }
                else if( propagationType == 2 )
                {
                    propagatorSettings = coupledPropagatorSettings;
                }

                if( propagateDynamicsOnly )
                {

                    SingleArcDynamicsSimulator< > dynamicsSimulator =
                            SingleArcDynamicsSimulator< >(
                                bodyMap, integratorSettings, propagatorSettings );

                    input_output::writeDataMapToTextFile(
                                dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),
                                propagationDataFilePrefix + "StateOnly.dat", dataFolder );
                    input_output::writeDataMapToTextFile(
                                dynamicsSimulator.getDependentVariableHistory( ),
                                propagationDataFilePrefix + "StateOnlyDependent.dat", dataFolder );
                }

                if( performPerturbationAnalysis )
                {
                    double positionPerturbation = 0.1;
                    double velocityPerturbation = 1.0E-4;

                    double orientationPerturbation = 1.0E-5;
                    double rotationRatePerturbation = 1.0E-10;
                    double zRotationRatePerturbation = 1.0E-13;


                    SingleArcDynamicsSimulator< > dynamicsSimulator =
                            SingleArcDynamicsSimulator< >(
                                bodyMap, integratorSettings, propagatorSettings, false );


                    if( propagationType == 0 || propagationType == 2 )
                    {
                        Eigen::VectorXd originalInitialTranslationalState = initialTranslationalState;
                        for( int element = 0; element < 6; element++ )
                        {
                            double perturbationMagnitude;
                            if( element < 3 )
                            {
                                perturbationMagnitude = positionPerturbation;
                            }
                            else
                            {
                                perturbationMagnitude = velocityPerturbation;
                            }
                            ;
                            for( int perturbation = 1; perturbation <= 2; perturbation++ )
                            {
                                initialTranslationalState = originalInitialTranslationalState;
                                initialTranslationalState( element ) += perturbation * perturbationMagnitude;
                                translationalPropagatorSettings->resetInitialStates( initialTranslationalState );
                                coupledPropagatorSettings->resetInitialStates(
                                            createCombinedInitialState( coupledPropagatorSettings->propagatorSettingsMap_ ) );

                                dynamicsSimulator.integrateEquationsOfMotion( propagatorSettings->getInitialStates( ) );
                                input_output::writeDataMapToTextFile(
                                            dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),
                                            propagationDataFilePrefix + "PerturbedTransStateUp" + std::to_string( perturbation ) +
                                            "_Element" + std::to_string( element ) + ".dat", dataFolder );
                                input_output::writeDataMapToTextFile(
                                            dynamicsSimulator.getDependentVariableHistory( ),
                                            propagationDataFilePrefix + "PerturbedTransDependentUp" + std::to_string( perturbation )+
                                            "_Element" + std::to_string( element ) + ".dat", dataFolder );

                                initialTranslationalState = originalInitialTranslationalState;
                                initialTranslationalState( element ) -= perturbation * perturbationMagnitude;
                                translationalPropagatorSettings->resetInitialStates( initialTranslationalState );
                                coupledPropagatorSettings->resetInitialStates(
                                            createCombinedInitialState( coupledPropagatorSettings->propagatorSettingsMap_ ) );

                                dynamicsSimulator.integrateEquationsOfMotion( propagatorSettings->getInitialStates( ) );
                                input_output::writeDataMapToTextFile(
                                            dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),
                                            propagationDataFilePrefix + "PerturbedTransStateDown" + std::to_string( perturbation )+
                                            "_Element" + std::to_string( element ) + ".dat", dataFolder );
                                input_output::writeDataMapToTextFile(
                                            dynamicsSimulator.getDependentVariableHistory( ),
                                            propagationDataFilePrefix + "PerturbedTransDependentDown" + std::to_string( perturbation )+
                                            "_Element" + std::to_string( element ) + ".dat", dataFolder );


                                initialTranslationalState = originalInitialTranslationalState;
                                translationalPropagatorSettings->resetInitialStates( initialTranslationalState );
                                coupledPropagatorSettings->resetInitialStates(
                                            createCombinedInitialState( coupledPropagatorSettings->propagatorSettingsMap_ ) );

                            }
                        }
                    }

                    if( propagationType == 1 || propagationType == 2 )
                    {
                        Eigen::VectorXd originalInitialRotationalState = initialRotationalState;

                        for( int element = 0; element < 7; element++ )
                        {
                            double perturbationMagnitude;
                            if( element < 3 )
                            {
                                perturbationMagnitude = orientationPerturbation;
                            }
                            else if( element < 6)
                            {
                                perturbationMagnitude = rotationRatePerturbation;
                            }
                            else
                            {
                                perturbationMagnitude = zRotationRatePerturbation;
                            }

                            perturbationMagnitude = perturbationMagnitude;

                            for( int perturbation = 1; perturbation <= 2; perturbation++ )
                            {

                                initialRotationalState = originalInitialRotationalState;
                                initialRotationalState( element ) += perturbation * perturbationMagnitude;
                                rotationalPropagatorSettings->resetInitialStates( initialRotationalState );
                                coupledPropagatorSettings->resetInitialStates(
                                            createCombinedInitialState( coupledPropagatorSettings->propagatorSettingsMap_ ) );

                                dynamicsSimulator.integrateEquationsOfMotion( propagatorSettings->getInitialStates( ) );
                                input_output::writeDataMapToTextFile(
                                            dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),
                                            propagationDataFilePrefix + "PerturbedRotStateUp" + std::to_string( perturbation )+
                                            "_Element" + std::to_string( element ) + ".dat", dataFolder );
                                input_output::writeDataMapToTextFile(
                                            dynamicsSimulator.getDependentVariableHistory( ),
                                            propagationDataFilePrefix + "PerturbedRotDependentUp" + std::to_string( perturbation )+
                                            "_Element" + std::to_string( element ) + ".dat", dataFolder );

                                initialRotationalState = originalInitialRotationalState;
                                initialRotationalState( element ) -= perturbation * perturbationMagnitude;
                                rotationalPropagatorSettings->resetInitialStates( initialRotationalState );
                                coupledPropagatorSettings->resetInitialStates(
                                            createCombinedInitialState( coupledPropagatorSettings->propagatorSettingsMap_ ) );

                                dynamicsSimulator.integrateEquationsOfMotion( propagatorSettings->getInitialStates( ) );
                                input_output::writeDataMapToTextFile(
                                            dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),
                                            propagationDataFilePrefix + "PerturbedRotStateDown" + std::to_string( perturbation )+
                                            "_Element" + std::to_string( element ) + ".dat", dataFolder );
                                input_output::writeDataMapToTextFile(
                                            dynamicsSimulator.getDependentVariableHistory( ),
                                            propagationDataFilePrefix + "PerturbedRotDependentDown" + std::to_string( perturbation )+
                                            "_Element" + std::to_string( element ) + ".dat", dataFolder );

                                initialRotationalState = originalInitialRotationalState;
                                rotationalPropagatorSettings->resetInitialStates( initialRotationalState );
                                coupledPropagatorSettings->resetInitialStates(
                                            createCombinedInitialState( coupledPropagatorSettings->propagatorSettingsMap_ ) );


                            }
                        }
                    }
                }

                if( performParameterPerturbationAnalysis )
                {

                    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;

                    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >
                                              ( "Phobos", gravitational_parameter ) );
                    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >
                                              ( "Deimos", gravitational_parameter ) );
                    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >
                                              ( "Mars", gravitational_parameter ) );
                    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >(
                                                  "Phobos", ground_station_position, "Lander" ) );


                    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                                  2, 0, 2, 2, "Phobos", spherical_harmonics_cosine_coefficient_block ) );
                    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                                  2, 1, 2, 2, "Phobos", spherical_harmonics_sine_coefficient_block ) );
                    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                                  2, 0, 2, 2, "Mars", spherical_harmonics_cosine_coefficient_block ) );
                    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                                  2, 1, 2, 2, "Mars", spherical_harmonics_sine_coefficient_block ) );

                    // Create parameters
                    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                            createParametersToEstimate( parameterNames, bodyMap );

                    Eigen::VectorXd parameterPerturbations = Eigen::VectorXd( 16 );
                    parameterPerturbations( 0 ) = 1000.0;
                    parameterPerturbations( 1 ) = 1000.0;
                    parameterPerturbations( 2 ) = 1000.0;

                    parameterPerturbations( 3 ) = 1.0;
                    parameterPerturbations( 4 ) = 1.0;
                    parameterPerturbations( 5 ) = 1.0;

                    parameterPerturbations( 6 ) = 1.0E-6;
                    parameterPerturbations( 7 ) = 1.0E-3;
                    parameterPerturbations( 8 ) = 1.0E-6;

                    parameterPerturbations( 9 ) = 1.0E-3;
                    parameterPerturbations( 10 ) = 1.0E-3;

                    parameterPerturbations( 11 ) = 1.0E-10;
                    parameterPerturbations( 12 ) = 1.0E-10;
                    parameterPerturbations( 13 ) = 1.0E-10;

                    parameterPerturbations( 14 ) = 1.0E-10;
                    parameterPerturbations( 15 ) = 1.0E-10;

                    Eigen::VectorXd initialParameters = parametersToEstimate->getFullParameterValues< double >( );
                    Eigen::VectorXd perturbedParameters = initialParameters;


                    SingleArcDynamicsSimulator< > dynamicsSimulator =
                            SingleArcDynamicsSimulator< >(
                                bodyMap, integratorSettings, propagatorSettings, false );


                    for( int element = 6; element < 11; element++ )
                    {
                        double perturbationMagnitude = parameterPerturbations( element );

                        std::cout<<"Index: "<<element - 6<<" =================================================== "<<std::endl;

                        for( int perturbation = 1; perturbation <= 2; perturbation++ )
                        {
                            perturbedParameters = initialParameters;
                            perturbedParameters( element ) += perturbation * perturbationMagnitude;
                            std::cout<<"UP "<<perturbation<<std::endl;

                            parametersToEstimate->resetParameterValues( perturbedParameters );

                            dynamicsSimulator.integrateEquationsOfMotion( propagatorSettings->getInitialStates( ) );
                            input_output::writeDataMapToTextFile(
                                        dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),
                                        propagationDataFilePrefix + "PerturbedParameterUp" + std::to_string( perturbation ) +
                                        "_Element" + std::to_string( element ) + ".dat", dataFolder );
//                            input_output::writeDataMapToTextFile(
//                                        dynamicsSimulator.getDependentVariableHistory( ),
//                                        propagationDataFilePrefix + "PerturbedTransDependentUp" + std::to_string( perturbation )+
//                                        "_Element" + std::to_string( element ) + ".dat", dataFolder );

                            std::cout<<"Down "<<perturbation<<std::endl;
                            perturbedParameters = initialParameters;
                            perturbedParameters( element ) -= perturbation * perturbationMagnitude;
                            parametersToEstimate->resetParameterValues( perturbedParameters );

                            dynamicsSimulator.integrateEquationsOfMotion( propagatorSettings->getInitialStates( ) );
                            input_output::writeDataMapToTextFile(
                                        dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),
                                        propagationDataFilePrefix + "PerturbedParameterDown" + std::to_string( perturbation )+
                                        "_Element" + std::to_string( element ) + ".dat", dataFolder );
//                            input_output::writeDataMapToTextFile(
//                                        dynamicsSimulator.getDependentVariableHistory( ),
//                                        propagationDataFilePrefix + "PerturbedTransDependentDown" + std::to_string( perturbation )+
//                                        "_Element" + std::to_string( element ) + ".dat", dataFolder );
                            perturbedParameters = initialParameters;
                            parametersToEstimate->resetParameterValues( perturbedParameters );
                        }
                    }



                }

                if( findDampedInitialState && propagationType > 0 )
                {
                    std::vector< double > dissipationTimes;

                    double currentDissipationTime = 0.025 * 86400.0;
                    while( currentDissipationTime < finalEphemerisTime / 10.0 )
                    {
                        dissipationTimes.push_back( currentDissipationTime );
                        currentDissipationTime *= 2.0;
                    }

                    std::vector< std::pair< std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >,
                            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > > > propagatedStates;
                    std::vector< std::pair< std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >,
                            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > > > dependentVariables;

                    Eigen::VectorXd dampedInitialState = getZeroProperModeRotationalState(
                                bodyMap, integratorSettings, propagatorSettings,
                                2.0 * mathematical_constants::PI / ( 0.3190365374 * 86400.0 ), dissipationTimes,
                                propagatedStates, dependentVariables, true, true, propagationDataFilePrefix, dataFolder + "/DampedAnalysis/" );
                    input_output::writeMatrixToFile(
                                dampedInitialState, propagationDataFilePrefix + "DampedInitialState.dat", 16, dataFolder  );
                }

                if( propagateVariationalEquations )
                {
                    std::vector< LinkEnds > linkEndsList;
                    LinkEnds currentLinkEnds;
                    currentLinkEnds[ transmitter ] = std::make_pair( "Earth", "" );
                    currentLinkEnds[ receiver ] = std::make_pair( "Phobos", "Lander" );
                    linkEndsList.push_back( currentLinkEnds );

                    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
                    linkEndsPerObservable[ one_way_range ].push_back( linkEndsList[ 0 ] );

                    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
                    if( propagationType == 0 || propagationType == 2 )
                    {
                        parameterNames.push_back(
                                    std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                        "Phobos", initialTranslationalState, "Mars" ) );
                    }

                    if( propagationType == 1 || propagationType == 2 )
                    {
                        parameterNames.push_back(
                                    std::make_shared< InitialRotationalStateEstimatableParameterSettings< double > >(
                                        "Phobos", initialRotationalState ) );
                    }

                    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >
                                              ( "Phobos", gravitational_parameter ) );
                    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >
                                              ( "Deimos", gravitational_parameter ) );
                    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >
                                              ( "Mars", gravitational_parameter ) );
                    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >(
                                                  "Phobos", ground_station_position, "Lander" ) );


                    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                                  2, 0, 2, 2, "Phobos", spherical_harmonics_cosine_coefficient_block ) );
                    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                                  2, 1, 2, 2, "Phobos", spherical_harmonics_sine_coefficient_block ) );
                    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                                  2, 0, 2, 2, "Mars", spherical_harmonics_cosine_coefficient_block ) );
                    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                                  2, 1, 2, 2, "Mars", spherical_harmonics_sine_coefficient_block ) );

                    // Create parameters
                    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                            createParametersToEstimate( parameterNames, bodyMap );
                    printEstimatableParameterEntries( parametersToEstimate );
                    integratorSettings->saveFrequency_ = 4;

                    SingleArcVariationalEquationsSolver< > variationalEquationsSolver =
                            SingleArcVariationalEquationsSolver< >(
                                bodyMap, integratorSettings, propagatorSettings, parametersToEstimate,
                                true, nullptr, false );

                    input_output::writeDataMapToTextFile(
                                variationalEquationsSolver.getNumericalVariationalEquationsSolution( )[ 0 ],
                            propagationDataFilePrefix + "StateTransition.dat", dataFolder );
                    input_output::writeDataMapToTextFile(
                                variationalEquationsSolver.getNumericalVariationalEquationsSolution( )[ 1 ],
                            propagationDataFilePrefix + "Sensitivity.dat", dataFolder  );
                    input_output::writeDataMapToTextFile(
                                variationalEquationsSolver.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( ),
                                propagationDataFilePrefix + "StateFromVar.dat", dataFolder  );


                }
            }
        }
    }
}

