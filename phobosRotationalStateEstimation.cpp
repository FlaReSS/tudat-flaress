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
                                       -3600.0, 1.0 * 86400.0 + 3600.0, 120.0,
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
    bodyMap[ "Phobos" ]->setBodyInertiaTensor( phobosInertiaTensor,
                                               ( phobosInertiaTensor( 0, 0 ) + phobosInertiaTensor( 1, 1 ) + phobosInertiaTensor( 2, 2 ) ) / 3.0 );
    std::cout<<"Scaled MOI "<<bodyMap[ "Phobos" ]->getScaledMeanMomentOfInertia( )<<std::endl;

    bodyMap[ "Phobos" ]->setShapeModel(
                std::make_shared< SphericalBodyShapeModel >( 15.0E3 ) );
    std::cout<<"Scaled MOI "<<bodyMap[ "Phobos" ]->getScaledMeanMomentOfInertia( )<<std::endl;

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
    std::cout<<"Scaled MOI "<<bodyMap[ "Phobos" ]->getScaledMeanMomentOfInertia( )<<std::endl;

    bodyMap[ "Phobos" ]->setGravityFieldModel(
                std::make_shared< gravitation::SphericalHarmonicsGravityField >(
                    phobosGravitationalParameter, phobosReferenceRadius, phobosCosineGravityFieldCoefficients,
                    phobosSineGravityFieldCoefficients, "IAU_Phobos",
                    std::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodyMap.at( "Phobos" ), true ) ) );
    std::cout<<"Scaled MOI "<<bodyMap[ "Phobos" ]->getScaledMeanMomentOfInertia( )<<std::endl;

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
    bool findInitialState = false;
    bool propagateVariationalEquationsOnly = false;
    bool setPropagationDataFromFiles = true;
    std::string propagationDataFilePrefix = "phobosCoupledNominalWithMoi";


    std::string dataFolder = "/home/dominic/Documents/Articles/PhobosCoupledDynamics/Data/Estimation/";
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "meta_PhobosPCK.tm" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "NOE-4-2015-b.bsp" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de430.bsp" );

    // Define time range of test.
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 3.0 * 365 *86400.0;

    // Retrieve list of body objects.
    NamedBodyMap bodyMap = getTestBodyMap( 9376.0E3, 0, finalEphemerisTime );
    createGroundStation( bodyMap.at( "Phobos" ), "Lander", ( Eigen::Vector3d( ) << 0.1, 0.35, 0.0 ).finished( ), geodetic_position );
    setGlobalFrameBodyEphemerides( bodyMap, "Mars", "ECLIPJ2000" );

    // Set torques between bodies that are to be taken into account.
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Phobos" );

    // Define mean motion (equal to rotation rate).
    Eigen::VectorXd initialRotationalState =
            propagators::getInitialRotationalStateOfBody(
                "Phobos", "ECLIPJ2000", bodyMap, initialEphemerisTime );

    // Create torque models
    SelectedTorqueMap torqueMap;
    torqueMap[ "Phobos" ][ "Mars" ].push_back( std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );
    torqueMap[ "Phobos" ][ "Sun" ].push_back( std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );
    torqueMap[ "Phobos" ][ "Jupiter" ].push_back( std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );
    torqueMap[ "Phobos" ][ "Deimos" ].push_back( std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );

    // Create torque models
    basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap(
                bodyMap, torqueMap, bodiesToIntegrate );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > translationalBodiesToPropagate;
    std::vector< std::string > translationalCentralBodies;

    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfPhobos;
    //        accelerationsOfPhobos[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    //        accelerationsOfPhobos[ "Mars" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
    accelerationsOfPhobos[ "Mars" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 2, 2 ) );

    accelerationsOfPhobos[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfPhobos[ "Deimos" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfPhobos[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );

    accelerationMap[ "Phobos" ] = accelerationsOfPhobos;

    translationalBodiesToPropagate.push_back( "Phobos" );
    translationalCentralBodies.push_back( "Mars" );

    // Create acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, translationalBodiesToPropagate, translationalCentralBodies );

    // Define integrator settings.
    double timeStep = 60.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettings< > >
            ( rungeKuttaVariableStepSize, initialEphemerisTime, timeStep,
              RungeKuttaCoefficients::rungeKuttaFehlberg78, timeStep, timeStep, 1.0, 1.0 );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >( keplerian_state_dependent_variable, "Phobos", "Mars" ) );
    dependentVariablesList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_relative_spherical_position, "Phobos", "Mars" ) );
    dependentVariablesList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_relative_spherical_position, "Mars", "Phobos" ) );
    dependentVariablesList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    euler_angles_to_body_fixed_313, "Phobos" ) );
    // Create object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );


    // Define propagator settings.
    std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
            std::make_shared< RotationalStatePropagatorSettings< double > >
            ( torqueModelMap, bodiesToIntegrate, initialRotationalState, std::make_shared< PropagationTimeTerminationSettings >(
                  finalEphemerisTime ) );

    Eigen::VectorXd initialTranslationalState =
            propagators::getInitialStatesOfBodies(
                translationalBodiesToPropagate, translationalCentralBodies, bodyMap, initialEphemerisTime );
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( translationalCentralBodies, accelerationModelMap, translationalBodiesToPropagate,
              initialTranslationalState, finalEphemerisTime, cowell, dependentVariablesToSave, 365.25 / 10.0 * 86400.0 );


    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > >  propagatorSettingsList;
    propagatorSettingsList.push_back( translationalPropagatorSettings );
    propagatorSettingsList.push_back( rotationalPropagatorSettings );

    std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings =
            std::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsList,
                std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime ),
                dependentVariablesToSave, 365.25 / 10.0 * 86400.0 );

    std::cout<<"Initial state "<<std::setprecision( 16 )<<propagatorSettings->getInitialStates( ).transpose( )<<std::endl;

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
    else if( findInitialState )
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
                    propagatedStates, dependentVariables, true, true, propagationDataFilePrefix, dataFolder );

        input_output::writeMatrixToFile(
                    dampedInitialState, propagationDataFilePrefix + "DampedInitialState.dat", 16, dataFolder  );
    }
    else
    {
        std::vector< LinkEnds > linkEndsList;
        LinkEnds currentLinkEnds;
        currentLinkEnds[ transmitter ] = std::make_pair( "Earth", "" );
        currentLinkEnds[ receiver ] = std::make_pair( "Phobos", "Lander" );
        linkEndsList.push_back( currentLinkEnds );

        std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
        linkEndsPerObservable[ one_way_range ].push_back( linkEndsList[ 0 ] );

        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
        parameterNames.push_back(
                    std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                        "Phobos", initialTranslationalState, "Mars" ) );
        parameterNames.push_back(
                    std::make_shared< InitialRotationalStateEstimatableParameterSettings< double > >(
                        "Phobos", initialRotationalState ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >
                                  ( "Phobos", gravitational_parameter ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >
                                  ( "Deimos", gravitational_parameter ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >
                                  ( "Phobos", mean_moment_of_inertia ) );
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

        if( propagateVariationalEquationsOnly )
        {
            integratorSettings->initialTimeStep_ = integratorSettings->initialTimeStep_;
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
        else
        {

            std::cout<<"FULL SIMULATION "<<std::endl;

            observation_models::ObservationSettingsMap observationSettingsMap;
            for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
                 linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
            {
                ObservableType currentObservable = linkEndIterator->first;

                std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
                for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
                {
                    observationSettingsMap.insert( std::make_pair( currentLinkEndsList.at( i ),
                                                                   std::make_shared< ObservationSettings >( currentObservable ) ) );
                }
            }
            std::cout<<"Propagating"<<std::endl;


            OrbitDeterminationManager< double, double > orbitDeterminationManager =
                    OrbitDeterminationManager< double, double >(
                        bodyMap, parametersToEstimate, observationSettingsMap,
                        integratorSettings, propagatorSettings, !setPropagationDataFromFiles );

            if( setPropagationDataFromFiles )
            {
                std::map< double, Eigen::MatrixXd > stateTransitionHistory =
                        input_output::readMatrixHistoryFromFile< double, double >(
                            13, 13, dataFolder + propagationDataFilePrefix + "StateTransition.dat" );
                std::map< double, Eigen::MatrixXd > sensitivityHistory =
                        input_output::readMatrixHistoryFromFile< double, double >(
                            13, 16, dataFolder + propagationDataFilePrefix + "Sensitivity.dat" );
                std::map< double, Eigen::VectorXd > stateHistory =
                        input_output::readVectorHistoryFromFile< double, double >(
                            13, dataFolder + propagationDataFilePrefix + "StateFromVar.dat" );
                std::dynamic_pointer_cast< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                            std::dynamic_pointer_cast< SingleArcVariationalEquationsSolver< > >(
                                orbitDeterminationManager.getVariationalEquationsSolver( ) )->getStateTransitionMatrixInterface( ) )->
                        updateMatrixInterpolators(
                            interpolators::createOneDimensionalInterpolator(
                                stateTransitionHistory, std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) ),
                            interpolators::createOneDimensionalInterpolator(
                                sensitivityHistory, std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) ) );
                std::dynamic_pointer_cast< SingleArcVariationalEquationsSolver< > >(
                            orbitDeterminationManager.getVariationalEquationsSolver( ) )->getDynamicsSimulator()->manuallySetAndProcessRawNumericalEquationsOfMotionSolution(
                            stateHistory, std::map< double, Eigen::VectorXd >( ), true );

                propagationDataFilePrefix += "WithPrecomputed";

            }


            std::vector< double > observationTimes;
            double currentTime = initialEphemerisTime + 1800.0;
            double observationTimeStep = 60.0;
            while( currentTime < finalEphemerisTime - 1800.0 )
            {
                for( int i = 0; i < 60; i++ )
                {
                    observationTimes.push_back( currentTime );
                    currentTime += observationTimeStep;
                }
                currentTime += 86400.0;
            }

            std::map< ObservableType, std::map< LinkEnds, std::shared_ptr< ObservationSimulationTimeSettings< double > > > > measurementSimulationInput;
            for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
                 linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
            {
                ObservableType currentObservable = linkEndIterator->first;
                std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
                for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
                {
                    measurementSimulationInput[ currentObservable ][ currentLinkEndsList.at( i ) ] =
                            std::make_shared< TabulatedObservationSimulationTimeSettings< double > >(
                                receiver, observationTimes );
                }
            }

            // Create observation viability settings and calculators
            std::vector< std::shared_ptr< ObservationViabilitySettings > > observationViabilitySettings;

            observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                        body_avoidance_angle, std::make_pair( "Earth", "" ), "Sun",
                                                        5.0 * mathematical_constants::PI / 180.0 ) );
            observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                        body_avoidance_angle, std::make_pair( "Phobos", "" ), "Sun",
                                                        5.0 * mathematical_constants::PI / 180.0 ) );
            observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                        body_occultation, std::make_pair( "Earth", "" ), "Mars" ) );
            observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                        body_occultation, std::make_pair( "Phobos", "" ), "Mars" ) );
            PerObservableObservationViabilityCalculatorList viabilityCalculators = createObservationViabilityCalculators(
                        bodyMap, linkEndsPerObservable, observationViabilitySettings );

            typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
            typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > > SingleObservablePodInputType;
            typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

            // Simulate observations
            PodInputDataType observationsAndTimes = simulateObservations< double, double >(
                        measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ),
                        viabilityCalculators );

            // Perturb parameter estimate
            Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
                    parametersToEstimate->template getFullParameterValues< double >( );
            Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
            std::cout<<"Truth: "<<( truthParameters ).transpose( )<<std::endl;

            int parameterSize = initialParameterEstimate.rows( );

            // Define estimation input
            std::shared_ptr< PodInput< double, double  > > podInput =
                    std::make_shared< PodInput< double, double > >(
                        observationsAndTimes, parameterSize,
                        Eigen::MatrixXd::Zero( parameterSize, parameterSize ) );
            podInput->defineEstimationSettings( false, false, true, true, true, true );

            // Perform estimation
            std::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                        podInput, std::make_shared< EstimationConvergenceChecker >( 0 ) );

            std::cout<<( podOutput->parameterEstimate_ - truthParameters ).transpose( )<<std::endl;

            input_output::writeMatrixToFile( podOutput->normalizedInformationMatrix_,
                                             propagationDataFilePrefix + "InformationMatrix.dat", 16, dataFolder );
            input_output::writeMatrixToFile( podOutput->getCorrelationMatrix( ),
                                             propagationDataFilePrefix + "Correlations.dat", 16, dataFolder );
            input_output::writeMatrixToFile( podOutput->inverseNormalizedCovarianceMatrix_,
                                             propagationDataFilePrefix + "InverseNormalizedCovariance.dat", 16, dataFolder );
            input_output::writeMatrixToFile( podOutput->getFormalErrorVector( ),
                                             propagationDataFilePrefix + "FormalEstimationError.dat", 16, dataFolder );
            input_output::writeMatrixToFile( truthParameters,
                                             propagationDataFilePrefix + "TruthParameters.dat", 16, dataFolder );
            //            input_output::writeMatrixToFile( podOutput->residuals_,
            //                                             propagationDataFilePrefix + "Residuals.dat", 16, dataFolder );
            input_output::writeMatrixToFile( podOutput->informationMatrixTransformationDiagonal_,
                                             propagationDataFilePrefix + "ParameterNormalization.dat", 16, dataFolder );
            //            input_output::writeMatrixToFile( podOutput->getResidualHistoryMatrix( ),
            //                                             propagationDataFilePrefix + "ResidualHistory.dat", 16, dataFolder );
            //            input_output::writeMatrixToFile( podOutput->getParameterHistoryMatrix( ),
            //                                             propagationDataFilePrefix + "ParameterHistory.dat", 16, dataFolder );
            input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                                 getConcatenatedTimeVector( podInput->getObservationsAndTimes( ) ) ),
                                             propagationDataFilePrefix + "ObservationTimes.dat", 16, dataFolder );
            input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ),
                                             propagationDataFilePrefix + "Observations.dat", 16, dataFolder );
        }
    }
}


