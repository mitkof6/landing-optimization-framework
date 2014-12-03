/* -------------------------------------------------------------------------- *
 *                     OpenSim:  explore_control.cpp                          *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2014 Stanford University and the Authors                *
 * Author(s): Carmichael Ong, Ajay Seth                                       *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

//==============================================================================
// This file generates a main() to explore the control parameters (desired length
// length gain, speed gain) of muscle stretch reflex controllers to decrese
// ankle inversion.
// This will build an executable that can take the following arguments:
//
// -m <OptimizedModel.osim>: model file to be loaded
// -csfile <OptimizedControllerSet.xml>: initial guess of controller set parameters for optimizer
// -kfile <OptimizedStates.sto>: trajectories of desired kinematics
// -ti <0.0> default time from which to start simulating
// -tf <0.1> default simulation end time
//
// Usage: exploreControl.exe -m optimizedDropModel.osim 
//          -csfile optimizedControllerSet.xml -kfile optimizedStates.sto
//			
//
//
// The -opt flag is used when an optimization should be run. If the -opt flag is
// not present, then a single forward simulation with analyses will be performed.
//
// Parallelization of the CMA optimizer is done by using MPI. If ENABLE_MPI was
// checked in CMake, then you can call on the executable to use more threads using:
//
// mpiexec -n <numthreads> dropLandingSim_nodal.exe -opt -lambda 50 ...
//==============================================================================

#include <ctime>  // clock(), clock_t, CLOCKS_PER_SEC
#include <OpenSim\OpenSim.h>
#include <DelayedPathReflexController.h>
#include "CMAOptimizer.h"

#ifdef USE_MPI
#include <mpi.h>
#endif 

#define DIETTAG 2

using namespace OpenSim;
using namespace SimTK;
using namespace std;

// Global variables
int evalCount = 0;
double bestSoFar = SimTK::Infinity;

double integratorTolerance = 1.0e-4;
const static double initialTime = 0.0, finalTime = 0.1;
const static double reflexDelay = 0.06; // 60 miliseconds
const static double platformAngle = 30; // degrees
const static std::string resultDir = "";
// variables related to objective function

// variables related to parameters
double maxDesiredFiberLength = 1.5;
double minDesiredFiberLength = 0.5;
double maxGain = 2;
double minGain = 0;

ofstream fwdLog("fwdLog.txt", ofstream::out);



//=============================================================================
// Uitility Methods
//=============================================================================

#include "time.h"
#include <string>

// string formatting helper utility
template <typename T>
std::string to_string(T const& value) {
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}

// crop path and .osim extension from filename
std::string cropFileName(std::string completeFileName) {

    //unsigned foundSlash = completeFileName.find_last_of("/\\");
    //std::string fileName = completeFileName.substr(foundSlash);
    unsigned foundPeriod = completeFileName.find_last_of(".");
    std::string croppedName = completeFileName.substr(0, foundPeriod);
    //sstr << value;
    return croppedName;
}

std::string get_date(void)
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer, 80, "%Y-%m-%d %X", timeinfo);
    puts(buffer);

    return std::string(buffer);
}


void simulateDropHeightAndCocontraction(Model* model, double height, double slope, double activation, string dir = "")
{
    // get the component of gravity in the y direction
    double g = model->getGravity()[1];

    // use gravity and height to compute the velocity after falling a distance equal to the input height.

    double landingVelocity = g*sqrt(-2 * height / g);

    Coordinate& ty = model->updCoordinateSet().get("pelvis_ty");
    ty.set_locked(false);
    ty.setDefaultSpeedValue(landingVelocity);
    /*ty.set_locked(false);
    ty.setDefaultValue(-1*height);
    ty.set_locked(true);*/

    Coordinate& rx = model->updCoordinateSet().get("platform_rx");
    rx.set_locked(false);
    rx.setDefaultValue(Pi*slope / 180.0);
    rx.set_locked(true);


    // set all invertor excitations to the predetermined value
    PrescribedController* invControls = dynamic_cast<PrescribedController*> (&model->updControllerSet().get("inverter_controls_r"));
    invControls->set_isDisabled(false);
    for (int i = 0; i<invControls->getActuatorSet().getSize(); i++)
    {

        OpenSim::Constant currentAct(activation);
        invControls->prescribeControlForActuator(i, currentAct.clone());

        OpenSim::ActivationFiberLengthMuscle* muscle = dynamic_cast<OpenSim::ActivationFiberLengthMuscle*> (&invControls->updActuators().get(i));
        if (muscle)
        {
            muscle->setDefaultActivation(activation);
        }
    }

    // set all evertor excitations to the predetermined value
    PrescribedController* evControls = dynamic_cast<PrescribedController*> (&model->updControllerSet().get("everter_controls_r"));
    evControls->set_isDisabled(false);
    for (int i = 0; i<evControls->getActuatorSet().getSize(); i++)
    {
        OpenSim::Constant currentAct(activation);
        evControls->prescribeControlForActuator(i, currentAct.clone());

        OpenSim::ActivationFiberLengthMuscle* muscle = dynamic_cast<OpenSim::ActivationFiberLengthMuscle*> (&evControls->updActuators().get(i));
        if (muscle)
        {
            muscle->setDefaultActivation(activation);
        }
    }

    // turn off other everter and inverter controllers
    model->updControllerSet().get("AnkleEverterReflexes").set_isDisabled(true);
    model->updControllerSet().get("AnkleInverterReflexes").set_isDisabled(true);
    model->updControllerSet().get("AnkleEverterDelayedReflexes").set_isDisabled(true);
    model->updControllerSet().get("AnkleInverterDelayedReflexes").set_isDisabled(true);


    // Create a force reporter
    ForceReporter* reporter = new ForceReporter(model);
    model->addAnalysis(reporter);

    SimTK::State& osim_state = model->initSystem();
    //ty.setValue(osim_state, -1*height);
    ty.setSpeedValue(osim_state, landingVelocity);
    model->getMultibodySystem().realize(osim_state, Stage::Velocity);
    model->equilibrateMuscles(osim_state);

    // Create the integrator for integrating system dynamics
    SimTK::RungeKuttaMersonIntegrator integrator(model->getMultibodySystem());
    integrator.setAccuracy(1.0e-6);
    integrator.setMaximumStepSize(1);
    integrator.setMinimumStepSize(1e-8);

    // Create the manager managing the forward integration and its outputs
    Manager manager(*model, integrator);
    // Integrate from initial time to final time
    double ti = 0;
    double tf = 0.3; //300 ms after landing

    manager.setInitialTime(ti);
    manager.setFinalTime(tf);
    cout << "\nIntegrating from " << ti << " to " << tf << endl;
    manager.integrate(osim_state);

    //////////////////////////////
    // SAVE THE RESULTS TO FILE //
    //////////////////////////////
    // Save the model states from forward integration
    Storage statesDegrees(manager.getStateStorage());

    char trial_name[100];
    sprintf(trial_name, "%sheight_%.1f_incline_%.1f_activation_%.1f", dir, height, slope, activation);
    statesDegrees.print(string(trial_name) + string("_states.sto"));

    // Save the forces
    reporter->getForceStorage().print(string(trial_name) + string("_forces.sto"));

    // Save the controls
    model->printControlStorage(trial_name + string("_controls.sto"));

    // save the model 
    model->print(string(trial_name) + string("_model.osim"));
}
//______________________________________________________________________________

void ExploreCoactivation(Model* model, State& initial_state, double ti, double tf, double slope, double activation, string dir = "")
{



    Coordinate& rx = model->updCoordinateSet().get("platform_rx");
    rx.set_locked(false);
    rx.setDefaultValue(Pi*slope / 180.0);
    rx.set_locked(true);

    ControllerSet& controllers = model->updControllerSet();

    // set all invertor excitations to the predetermined value
    PrescribedController* invControls = dynamic_cast<PrescribedController*> (&model->updControllerSet().get("inverter_controls_r"));
    invControls->set_isDisabled(false);
    for (int i = 0; i<invControls->getActuatorSet().getSize(); i++)
    {

        OpenSim::Constant currentAct(activation);
        invControls->prescribeControlForActuator(i, currentAct.clone());

        OpenSim::ActivationFiberLengthMuscle* muscle = dynamic_cast<OpenSim::ActivationFiberLengthMuscle*> (&invControls->updActuators().get(i));
        if (muscle)
        {
            muscle->setDefaultActivation(activation);
        }
    }

    // set all evertor excitations to the predetermined value
    PrescribedController* evControls = dynamic_cast<PrescribedController*> (&model->updControllerSet().get("everter_controls_r"));
    evControls->set_isDisabled(false);
    for (int i = 0; i<evControls->getActuatorSet().getSize(); i++)
    {
        OpenSim::Constant currentAct(activation);
        evControls->prescribeControlForActuator(i, currentAct.clone());

        OpenSim::ActivationFiberLengthMuscle* muscle = dynamic_cast<OpenSim::ActivationFiberLengthMuscle*> (&evControls->updActuators().get(i));
        if (muscle)
        {
            muscle->setDefaultActivation(activation);
        }
    }

    // turn off other everter and inverter controllers
    //model->updControllerSet().get("AnkleEverterDelayedReflexes").set_isDisabled(true);
    //model->updControllerSet().get("AnkleInverterDelayedReflexes").set_isDisabled(true);

    // Create a force reporter
    ForceReporter* reporter = new ForceReporter(model);
    model->addAnalysis(reporter);

    SimTK::State& osim_state = model->initSystem();
    osim_state.setQ(initial_state.getQ());
    osim_state.setU(initial_state.getU());

    rx.set_locked(false);
    rx.setValue(osim_state, Pi*slope / 180.0);
    rx.set_locked(true);
    model->equilibrateMuscles(osim_state);

    // Create the integrator for integrating system dynamics
    SimTK::RungeKuttaMersonIntegrator integrator(model->getMultibodySystem());
    integrator.setAccuracy(integratorTolerance);
    integrator.setMaximumStepSize(1);
    integrator.setMinimumStepSize(1e-8);

    // Create the manager managing the forward integration and its outputs
    Manager manager(*model, integrator);
    // Integrate from initial time to final time

    manager.setInitialTime(ti);
    manager.setFinalTime(tf);
    cout << "\nIntegrating from " << ti << " to " << tf << endl;
    manager.integrate(osim_state);

    //////////////////////////////
    // SAVE THE RESULTS TO FILE //
    //////////////////////////////
    // Save the model states from forward integration
    Storage statesDegrees(manager.getStateStorage());

    char trial_name[100];
    sprintf(trial_name, "%sincline_%.1f_activation_%.1f", dir, slope, activation);

    statesDegrees.print(string(trial_name) + string("_states.sto"));

    // Save the forces
    reporter->getForceStorage().print(string(trial_name) + string("_forces.sto"));

    // Save the controls
    model->printControlStorage(trial_name + string("_controls.sto"));

    // save the model 
    model->print(string(trial_name) + string("_model.osim"));

}
//______________________________________________________________________________

void ExploreReflexes(Model* model, State& initial_state, double ti, double tf, double slope, double gain, double delay, string dir = "")
{

    

    Coordinate& rx = model->updCoordinateSet().get("platform_rx");
    rx.set_locked(false);
    rx.setDefaultValue(Pi*slope / 180.0);
    rx.set_locked(true);

    ControllerSet& controllers = model->updControllerSet();

    for (int i = 0; i < controllers.getSize(); i++)
    {

        DelayedPathReflexController* reflex = dynamic_cast<DelayedPathReflexController*>(&controllers.get(i));
        if (reflex) {
            reflex->set_isDisabled(false);
            reflex->set_gain(gain);
            reflex->set_delay(delay);
        }
    }

        // Create a force reporter
    ForceReporter* reporter = new ForceReporter(model);
    model->addAnalysis(reporter);
    
    SimTK::State& osim_state = model->initSystem();
    osim_state.setQ(initial_state.getQ());
    osim_state.setU(initial_state.getU());
    
    rx.set_locked(false);
    rx.setValue(osim_state, Pi*slope / 180.0);
    rx.set_locked(true);
    model->equilibrateMuscles(osim_state);

    // Create the integrator for integrating system dynamics
    SimTK::RungeKuttaMersonIntegrator integrator(model->getMultibodySystem());
    integrator.setAccuracy(integratorTolerance);
    integrator.setMaximumStepSize(1);
    integrator.setMinimumStepSize(1e-8);

    // Create the manager managing the forward integration and its outputs
    Manager manager(*model, integrator);
    // Integrate from initial time to final time

    manager.setInitialTime(ti);
    manager.setFinalTime(tf);
    cout << "\nIntegrating from " << ti << " to " << tf << endl;
    manager.integrate(osim_state);

    //////////////////////////////
    // SAVE THE RESULTS TO FILE //
    //////////////////////////////
    // Save the model states from forward integration
    Storage statesDegrees(manager.getStateStorage());

    char trial_name[100];
    sprintf(trial_name, "%sincline_%.1f_gain_%.1f_delay_%.2f", dir, slope, gain, delay);

    statesDegrees.print(string(trial_name) + string("_states.sto"));

    // Save the forces
    reporter->getForceStorage().print(string(trial_name) + string("_forces.sto"));

    // Save the controls
    model->printControlStorage(trial_name + string("_controls.sto"));

    // save the model 
    model->print(string(trial_name) + string("_model.osim"));

}
//______________________________________________________________________________

int main(int argc, char* argv[])
{

    try {
        std::clock_t startTime = std::clock();	

        string modelFile = "bestModel.osim";
        string ctrlFile = "reflexControllers.xml";
        string kinematicsFile = "DesiredKinematics.sto";
        double ti = initialTime;
        double tf = finalTime; 

        // Options
        for (int i = 0; i < argc; i++) {
        
            if (strcmp(argv[i], "-m") == 0) {
                modelFile = argv[i+1];
            }

            if (strcmp(argv[i], "-csfile") == 0) {
                ctrlFile = argv[i+1];
            }

            if (strcmp(argv[i], "-kfile") == 0) {
                kinematicsFile = argv[i + 1];
            }

            if (strcmp(argv[i], "-ti") == 0) {
                ti = atof(argv[i + 1]);
            }

            if (strcmp(argv[i], "-tf") == 0) {
                tf = atof(argv[i + 1]);
            }

        }

        

        OpenSim::Model osimModel(modelFile);
        Vec3 platform_origin(0.0, 0.0, 0.29);
        Vec3 twist_y(0.0, 0.20943951023931953, 0.0);
        osimModel.updJointSet().get("ground_platform").setLocationInParent(platform_origin);
        //osimModel.updJointSet().get("ground_platform").setOrientationInParent(twist_y);

        osimModel.updForceSet().get("foot_floor_l").set_isDisabled(true);
        CoordinateSet& coords = osimModel.updCoordinateSet();
        //coords.get("platform_rx").set_locked(false);
        //coords.get("platform_rx").setDefaultValue(platformAngle*Pi/180);
        /*
        coords.get("platform_rx").setDefaultValue(0.0);
        coords.get("platform_rx").set_locked(true);
        coords.get("platform_ry").set_locked(false);
        coords.get("platform_ry").setDefaultValue(0.0);
        coords.get("platform_ry").set_locked(true);
        coords.get("platform_rz").set_locked(false);
        coords.get("platform_rz").setDefaultValue(0.0);
        coords.get("platform_rz").set_locked(true);
        coords.get("platform_ty").set_locked(false);
        coords.get("platform_ty").setDefaultValue(0.0);
        coords.get("platform_ty").set_locked(true);
        */


        Array<double> actLevels;

        actLevels.append(0.0);
        actLevels.append(0.1);
        actLevels.append(0.2);
        actLevels.append(0.3);
        actLevels.append(0.4);
        actLevels.append(0.5);
        actLevels.append(0.6);
        actLevels.append(0.7);
        actLevels.append(0.8);
        actLevels.append(0.9);
        actLevels.append(1.0);


        Array<double> reflexGains;

        reflexGains.append(0.0);
        reflexGains.append(0.5);
        reflexGains.append(1.0);
        reflexGains.append(2.0);
        reflexGains.append(5.0);
        reflexGains.append(10.0);
        //reflexGains.append(100.0);
        //reflexGains.append(1000.0);

        

        ControllerSet inputControllerSet(osimModel,ctrlFile);
        cout << "Found " << inputControllerSet.getSize() << " controllers." << endl;
        
        for (int i = 0; i < inputControllerSet.getSize(); i++)
        {
            
            DelayedPathReflexController* reflex = dynamic_cast<DelayedPathReflexController*>(&inputControllerSet.get(i));
            if (reflex) {
                for (int i = 0; i<reflex->getActuatorSet().getSize(); i++)
                {
                    OpenSim::ActivationFiberLengthMuscle* muscle = dynamic_cast<OpenSim::ActivationFiberLengthMuscle*> (&reflex->updActuators().get(i));
                    if (muscle)
                    {
                        muscle->setDefaultActivation(0.0);
                    }
                }
                reflex->set_isDisabled(true);
                osimModel.addController(reflex->clone());
            }
        }
        
        for (int i = 0; i < inputControllerSet.getSize(); i++)
        {

            PrescribedController* prescribed = dynamic_cast<PrescribedController*>(&inputControllerSet.get(i));
            if (prescribed) {
                prescribed->set_isDisabled(true);
                osimModel.addController(prescribed->clone());
            }
        }
        ControllerSet& controllers = osimModel.updControllerSet();
        cout << "Model has " << controllers.getSize() << " controllers." << endl;
        int numControllers, numParameters;
        numControllers = controllers.getSize();
        //optLog << "numControllers = " << numControllers << endl;
        
        
        State& osimState = osimModel.initSystem();

        /* Read desired kinematics*/
        Storage k_desired = Storage(kinematicsFile);
        
        Storage *q_desired, *u_desired;
        osimModel.getSimbodyEngine().formCompleteStorages(osimState, k_desired, q_desired, u_desired);
        cout << "Joint velocities estimated." << endl;
        if (q_desired->isInDegrees())
        {
            cout << "Converting coordinates from degrees to radians." << endl;
            osimModel.getSimbodyEngine().convertDegreesToRadians(*q_desired);
        }
        cout << "All desired coordinates are in radians." << endl;

        if (u_desired->isInDegrees())
        {
            cout << "Converting velocities from degrees to radians." << endl;
            osimModel.getSimbodyEngine().convertDegreesToRadians(*u_desired);
        }
        cout << "All desired velocities are in radians." << endl;
        Vector qi(osimState.getNQ(), 0.0), ui(osimState.getNU(),0.0);
        // get and set initial position
        q_desired->getDataAtTime(ti, osimState.getNQ(), qi);
        //cout << "Setting q's." << endl;
        //qi.dump();

        osimState.setQ(qi);
        // get and set initial velocity
        u_desired->getDataAtTime(ti, osimState.getNU(), ui);
        //cout << "Setting u's" << endl;
        //ui.dump();
        osimState.setU(ui);
        cout << "Initial position and velocity set." << endl;
        State initial_state = osimState;
        
        for (int j = 0; j<actLevels.getSize(); j++)
        {
            cout << "\tCo-contraction level " << actLevels[j]
                << " (" << j + 1 << "/" << actLevels.getSize() << "):" << endl;

            Model* modelCopy = osimModel.clone();
            ExploreCoactivation(modelCopy, initial_state, ti, tf, platformAngle, actLevels[j], resultDir);
        }
        

        for (int j = 0; j<reflexGains.getSize(); j++)
        {
            cout << "\tReflex Gain " << reflexGains[j]
                << " (" << j + 1 << "/" << reflexGains.getSize() << "):" << endl;
            Model* modelCopy = osimModel.clone();
            ExploreReflexes(modelCopy, initial_state, ti, tf, platformAngle, reflexGains[j], reflexDelay,resultDir);
            //simulateDropHeightAndReflexes(dropModel, dropHeights[i], platformAngle, reflexGains[j], reflexDelay);
        }

        clock_t endTime = std::clock();
        fwdLog << "computeTime: " << endTime - startTime << "ms" << endl;
    }
    catch (const std::exception& ex)
    {
        std::cout << ex.what() << std::endl;
        return 1;
    }

    catch (...)
    {
        std::cout << "UNRECOGNIZED EXCEPTION" << std::endl;
        return 1;
    }

    std::cout << "\nSimulations Completed\n";
    return 0;
    // End of main() routine.
    
}

