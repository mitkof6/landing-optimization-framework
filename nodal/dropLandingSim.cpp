/* -------------------------------------------------------------------------- *
 *                     OpenSim:  dropLandingSim.cpp                           *
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
// This file generates a main() to optimize the control parameters (desired length
// length gain, speed gain) of muscle stretch controllers to match a given motion.
// This will build an executable that can take the following arguments:
//
// -m <InputDropModel.osim>: model file to be loaded
// -csfile <InitialControllerSet.xml>: initial guess of controller set parameters for optimizer
// -kfile <DesiredKinematics.sto>: trajectories of desired kinematics
// -ti <0.0> default time from which to start simulating
// -tf <0.1> default simulation end time
// -opt: flag to run optimization
// -resume: flag to run optimization but resume using resumecmaes.dat file
// -lambda <int>: CMA parameter, set number of samples per outer iteration
// -sigma <double>: CMA parameter, set initial constant of covariance matrix
// -maxIters <int>: optimizer parameter, set max outer iterations
//
// Usage: dropLandingSim_nodal.exe -m InputDropModel.osim 
//          -csfile InitialControllerSet.xml -kfile DesiredKinematics.sto
//			-ti 0.0 -tf 0.1 -opt -lambda 50 -sigma 0.1 -maxIters 2000
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
#include <MuscleFiberStretchController.h>
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

// variables related to objective function

// variables related to parameters
double maxDesiredFiberLength = 1.5;
double minDesiredFiberLength = 0.5;
double maxGain = 2;
double minGain = 0;
// Optimizer settings and flags
int lambda = 30;
double sigma = 0.1;
int maxIterations = 1;
bool doOpt = false;
bool doResume = false;

//Vec3 toeInFootPosition = Vec3(0.14, 0.0, 0.0);
ofstream optLog("optLog.txt", ofstream::out);
ofstream fwdLog("fwdLog.txt", ofstream::out);
std::clock_t optimizerStartTime;


//=============================================================================
// Uitility Methods
//=============================================================================

#include "time.h"
#include <string>

#define MAX_DATE 19

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

double gainToParameter(double gainValue)
{
	return (gainValue - minGain) / (maxGain - minGain);
}

double parameterToGain(double parameterValue)
{
	return parameterValue*(maxGain - minGain) + minGain;
}

double desiredFiberLengthToparameter(double lengthValue)
{
	return (lengthValue - minDesiredFiberLength) / (maxDesiredFiberLength - minDesiredFiberLength);
}

double parameterToDesiredFiberLength(double parameterValue)
{
	return parameterValue*(maxDesiredFiberLength - minDesiredFiberLength) + minDesiredFiberLength;
}

// need to convert strech controller's 3 control parameters into normalized optimizer parameters

Vector getNormalizedControllerParameters(const OpenSim::ControllerSet& cset)
{
	Vector parameters = Vector(3*cset.getSize());
	int k = 0;
	for (int i = 0; i < cset.getSize(); i++)
	{
		OpenSim::MuscleFiberStretchController* stretch = dynamic_cast<OpenSim::MuscleFiberStretchController*>(&cset.get(i));
		if (stretch)
		{
			parameters(k) = desiredFiberLengthToparameter(stretch->get_normalized_rest_length());
			k++;
			parameters(k) = gainToParameter(stretch->get_gain_length());
			k++;
			parameters(k) = gainToParameter(stretch->get_gain_velocity());
			k++;
		}
	}
	return parameters;
}

void pushNormalizedParametersToControllers(Vector inParameters, OpenSim::ControllerSet& cset)
{
	int k = 0;
	for (int i = 0; i < cset.getSize(); i++)
	{
		OpenSim::MuscleFiberStretchController* stretch = dynamic_cast<OpenSim::MuscleFiberStretchController*>(&cset.get(i));
		if (stretch)
		{
			stretch->set_normalized_rest_length(parameterToDesiredFiberLength(inParameters(k)));
			k++;
			stretch->set_gain_length(parameterToGain(inParameters(k)));
			k++;
			stretch->set_gain_velocity(parameterToGain(inParameters(k)));
			k++;
		}
	}
}

double evalRMSerror(const Storage& q_out, const Storage& ref_q);

//=============================================================================
// OPTIMIZER SYSTEM: LandingOptimizationSystem
// Defines a constructor and objective function for the optimization
//=============================================================================
class LandingOptimizationSystem : public OptimizerSystem {
   public:

	   /* Constructor class. Parameters passed are accessed in the objectiveFunc() class. */
	   LandingOptimizationSystem(int numParameters, State& s, OpenSim::Model& aModel, OpenSim::Manager& manager, OpenSim::Storage& q_desired): 
             OptimizerSystem(numParameters),
			 si(s),
		     osimModel(aModel),
			 manager(manager),
             ref_q(q_desired)
	   {
		   // Get a copy of the state to reuse
           cout << "Instantiating LandingOptimization System." << endl;
		   _sCopy = s;           
		   
	   }
   

	   
    void evalSimulation( const Vector &newParameters, double &f ) const {
		std::clock_t startTime = std::clock();

        // reuse the initial state;
		_sCopy = si;
		State& s = _sCopy;

		//update the controllers
		ControllerSet& controllers = osimModel.updControllerSet();
		pushNormalizedParametersToControllers(newParameters, controllers);
		
        //newParameters.dump();
        //osimModel.print("tested model.osim");

		osimModel.equilibrateMuscles(s);

        // must purge state storage from previous runs of the simulation.
        manager.getStateStorage().purge();
		manager.integrate(s);

		Storage& q_out = manager.getStateStorage();
        //q_out.print("last simulated states.sto");

		double rms_error = evalRMSerror(q_out, ref_q);
        //cout << "rms error = " << rms_error << endl;
        
		f = rms_error;

   }

	int objectiveFunc(const Vector &newControls, bool new_coefficients, Real& f ) const {
		evalSimulation(newControls, f);
		
		return 0;  
	}

private:
	OpenSim::Model& osimModel;
	OpenSim::Manager& manager;
	// the initial state that all simulations must begin with
	const State& si;
	// keep a working copy of the state so each objective func eval does not have to make a copy
	mutable State _sCopy;
	// the reference kinematics (desired kinematics) for comparison
	const Storage& ref_q;
 };


//______________________________________________________________________________
/**
 * Define an optimization problem that finds a set of stretch controller parameters 
 * to minimize the error between model kinematis and desired kinematics
 */
int main(int argc, char* argv[])
{
#ifdef USE_MPI
		int myrank; 
		int rc = MPI_Init(&argc, &argv);
		if (rc != MPI_SUCCESS) {
			printf ("Error starting MPI program. Terminating.\n");
			MPI_Abort(MPI_COMM_WORLD, rc);
		}

		int numtasks;
		numtasks = 0;
		myrank = 0;
		MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
		MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
		printf ("Number of processes= %d My rank= %d\n", numtasks,myrank);
#endif
//	try {
		std::clock_t startTime = std::clock();	

		string modelFile = "optimizedModel.osim";
		string ctrlFile = "InitialControllerSet.xml";
		string kinematicsFile = "DesiredKinematics.sto";
		double ti = 0.0;
		double tf = 0.1; // 100 milliseconds

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

			if (strcmp(argv[i], "-resume") == 0) {
				doResume = true;
			}
        
			if (strcmp(argv[i], "-maxIters") == 0) {
				maxIterations = atoi(argv[i+1]);
			}

			if (strcmp(argv[i], "-lambda") == 0) {
				lambda = atoi(argv[i+1]);
			}

			if (strcmp(argv[i], "-sigma") == 0) {
				sigma = atof(argv[i+1]);
			}

			if (strcmp(argv[i], "-opt") == 0) {
				doOpt = true; 
			}

			if (strcmp(argv[i], "-ti") == 0) {
				ti = atof(argv[i + 1]);
			}

			if (strcmp(argv[i], "-tf") == 0) {
				tf = atof(argv[i + 1]);
			}

		}

        

		OpenSim::Model osimModel(modelFile);
        if (!doOpt) { osimModel.setUseVisualizer(true); }
		// For safety, empty the controllerSet
		osimModel.updControllerSet().clearAndDestroy();

		ControllerSet inputControllerSet(osimModel,ctrlFile);
        cout << "Found " << inputControllerSet.getSize() << " controllers." << endl;
        
		for (int i = 0; i < inputControllerSet.getSize(); i++)
		{
            
			MuscleFiberStretchController* stretch = dynamic_cast<MuscleFiberStretchController*>(&inputControllerSet.get(i));
			if (stretch) { 
                osimModel.addController(stretch->clone());
            }
		}
        

		ControllerSet& controllers = osimModel.updControllerSet();
        cout << "Model has " << controllers.getSize() << " controllers." << endl;
		// The number of parameters is (numController*3) because each controller has a desired length, length gain, and speed gain
		int numControllers, numParameters;
		numControllers = controllers.getSize();
		optLog << "numControllers = " << numControllers << endl;
		numParameters = numControllers*3;
		optLog << "numParameters = " << numParameters << endl;

		/* Define initial values for controllerParameters. */
		
		// initialize controller parameters
		Vector controllerParameters = getNormalizedControllerParameters(controllers);
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
        osimModel.equilibrateMuscles(osimState);
		// Create the integrator for the simulation.
		RungeKuttaMersonIntegrator integrator(osimModel.getMultibodySystem());
		//SemiExplicitEuler2Integrator integrator(osimModel.getMultibodySystem());
		integrator.setAccuracy(integratorTolerance);

		// Create a manager to run the simulation. Can change manager options to save run time and memory or print more information
        
        Manager manager(osimModel, integrator);
		manager.setWriteToStorage(true);
		manager.setPerformAnalyses(false);

		// Integrate from initial time to final time and integrate
		manager.setInitialTime(ti);
		manager.setFinalTime(tf);	

        cout << "Manager Configured:" << endl;
        cout << "start time = " << ti << endl;
        cout << "end time = " << tf << endl;


        int failed = 0;
			
		if (doOpt) 
		{
            cout << "Creating optimization system." << endl;
            // Create the OptimizationSystem. Initialize the objective function value "f".
            LandingOptimizationSystem sys(numParameters, osimState, osimModel, manager, *q_desired);

            cout << "Landing Optimization System Instantiated" << endl;
            Real f = NaN;

            // Set lower and upper bounds.  We've handled autonormalization of parameters in the optimizer
            Vector lower_bounds(numParameters, 0.0);
            Vector upper_bounds(numParameters, 1.0);
            sys.setParameterLimits(lower_bounds, upper_bounds);
            optLog << "lower_bounds = " << lower_bounds << endl;
            optLog << "upper_bounds = " << upper_bounds << endl << endl;

            cout << "Optimization System Configured" << endl;
                      

			CMAOptimizer opt(sys);
			optLog << "using CMA Optimizer" << endl;
			//SimTK::Optimizer opt(sys, InteriorPoint); optLog << "using IPOPT" << endl;
			//SimTK::Optimizer opt(sys, LBFGSB); optLog << "using LBFGSB" << endl;

			optLog << "maxIterations = " << maxIterations << endl;
			optLog << "lambda = " << lambda << endl;
			optLog << "sigma = " << sigma << endl;
			optLog << "integratorTolerance = " << integratorTolerance << endl;


			//opt.setConvergenceTolerance(integratorTolerance);
			//opt.useNumericalGradient(true, 10*integratorTolerance);

			opt.setMaxIterations(maxIterations);


			if (doResume) {opt.setAdvancedBoolOption("resume", true);}
			opt.setAdvancedIntOption("lambda", lambda);
			opt.setAdvancedRealOption("sigma", sigma);
#ifdef USE_MPI
			opt.setAdvancedBoolOption("enableMPI", true);
			if (numtasks == 1) 
				opt.setAdvancedBoolOption("enableMPI", false);
#else
			opt.setAdvancedBoolOption("enableMPI", false);
#endif
			
			try{
				
				f = opt.optimize(controllerParameters);
			}
			catch (const std::exception& ex){
				std::cout << ex.what() << std::endl;
				failed = 1;
			}
	
			cout << "Elapsed time = " << (std::clock()-startTime)/CLOCKS_PER_SEC << "s" << endl;

			ControllerSet& optimizedControllers = osimModel.updControllerSet();
			pushNormalizedParametersToControllers(controllerParameters, optimizedControllers);

			osimModel.print("optimizedModel.osim");
			optimizedControllers.print("optimizedControllers.xml");
			

            string date = get_date();
            replace(date.begin(), date.end(), ':', '-');
            cout << "Completion date and time: " << date << endl;

            osimModel.print(string("optimizedModel ") + date + string(".osim"));
            optimizedControllers.print(string("optimizedControllers ") + date + string(".xml"));
			cout << "\nMinimum Objective Function Value = " << f << endl;

		} else {

			// Run a single integration either specific control file specified without -opt tag

            //initialState.setY(osimState.getY());
            //osimModel.equilibrateMuscles(initialState);
			// Adding analyses
			ForceReporter* forceReporter = new ForceReporter(&osimModel);
			forceReporter->setName("ForceReporter");
			forceReporter->includeConstraintForces(true);
			osimModel.addAnalysis(forceReporter);

			/*
			PointKinematics* pointKinematics_toe = new PointKinematics(&osimModel);
			pointKinematics_toe->setName("PointKinematics_TOE");
			pointKinematics_toe->setBodyPoint("foot", toeInFootPosition);
			osimModel.addAnalysis(pointKinematics_toe);

			PointKinematics* pointKinematics_heel = new PointKinematics(&osimModel);
			pointKinematics_heel->setName("PointKinematics_HEEL");
			pointKinematics_heel->setBodyPoint("foot", Vec3(0.0));
			osimModel.addAnalysis(pointKinematics_heel);
			*/
			startTime = std::clock();

			manager.setWriteToStorage(true);
			manager.setPerformAnalyses(true);

			// Integrate from initial time to final time and integrate
			manager.setInitialTime(ti);
			manager.setFinalTime(tf);
            cout << "Saving current model" << endl;

            osimModel.print("model_to_integrate.osim");
            //osimModel.setUseVisualizer(true);
            

            cout << "Integrating from " << ti << " to " << tf << " seconds." << endl;
			manager.integrate(osimState);

			Storage states(manager.getStateStorage());
			states.print("states.sto");

			Storage forceStorage = forceReporter->getForceStorage();
			forceStorage.print("forceReporter.sto");

			clock_t endTime = std::clock();
			fwdLog << "computeTime: " << endTime - startTime << "ms" << endl;

			double rms_error = evalRMSerror(states, *q_desired);
			fwdLog << "kinematic RMS error = " << rms_error;

            int pause;
            cout << "Enter any letter key to close...";
            cin >> pause;
		}

#ifdef USE_MPI
	for (int rank = 1; rank < numtasks; ++rank) {
		MPI_Send(0, 0, MPI_INT, rank, DIETTAG, MPI_COMM_WORLD);
	}
	MPI_Finalize();
#endif
	
	// End of main() routine.
	return failed;
}

double evalRMSerror(const Storage& q_out, const Storage& ref_q)
{

	Array<std::string> refCoordinates = ref_q.getColumnLabels();
    
	int refSize = refCoordinates.getSize();

    Array<std::string> outCoordinates = q_out.getColumnLabels();
    int outSize = q_out.getColumnLabels().getSize();

	Array<double> time_stamps;
	q_out.getTimeColumn(time_stamps);
	
	Array<int> out_indices;
	for (int i = 0; i < refSize; i++)
	{
		out_indices.append(q_out.getColumnLabels().findIndex(refCoordinates.get(i)));
	}
	double rms_error = 0.0;
	double diff = 0.0;
	double time = 0.0;
	Vector refdata = Vector(refSize);
	Vector outdata = Vector(outSize);
    cout << time_stamps.getSize() << " time steps" << endl;
    
	for (int i = 0; i < time_stamps.getSize(); i++)
	{
		time = time_stamps.get(i);
		ref_q.getDataAtTime(time, refSize, refdata);
		q_out.getDataAtTime(time, outSize, outdata);
		for (int j = 1; j < refSize; j++)
		{
			string coord = refCoordinates.get(j);
			diff = refdata(j - 1) - outdata(out_indices.get(j - 1));
			rms_error += diff*diff;
		}
	}
    cout << "rms error = " << rms_error << endl;
	return rms_error;
}
