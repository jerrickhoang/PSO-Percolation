

// for random numbers
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Random;
import java.text.SimpleDateFormat;



public class PSO {


	// ****************  MISCELLANEOUS	  ******************

	// for random numbers
	public static boolean useSeedForRand;
	public static int seed;   
	public static Random rand;


	// ****************  PSO   ******************

	public static int totalNumRuns;
	public static int currentRunNum;

	// iterations are counted in terms of function evaluations (FEs)
	// to facilitate runs where we are counting those instead of iterations
	public static int totalNumIters;
	public static int currentIterNum;
	public static int numItersPerOutputInterval;
	public static int numInitialFEsIgnored;

	// are we counting iterations or function evaluations?
	public static boolean useIterations;

	public static int totalNumFEs;
	public static int currentFENum;
	public static int numFEsPerOutputInterval;
	public static int numFEIntervals;

	public static int numIntervalsForDataOutput;

	// this is global static because it is needed in Particle.java to set
	// the iteration number of the currentSolution and the personal best solution
//	public static int currentNumParticles;


	// shape of the neighborhood
	// FNN neighborhoods must be distinguished from standard neighborhoods so we know
	// how to do updates in the Particle class
	// in Solé & Miramontes paper, it is Moore
	public static enum Topology {
		GBEST, RING, vonNEUMANN, MOORE, FNN_GBEST, FNN_RING, FNN_vonNEUMANN, FNN_MOORE,
	}
	public static int[] numRowsVonNeumannAndMooreList;
	public static int[] numColsVonNeumannAndMooreList;


	// include self?
	// same issue for any kind of neighborhood;
	// in Solé & Miramontes paper, it is to include self
	public static enum SelfModel {
		INCLUDE_SELF, NOT_INCLUDE_SELF
	}

	// who has influence in the neighborhood?
	// this is *not* an issue for the operation of the FNN, it is only an 
	// an issue wrt how the neighborhood constructed by the FNN is used;
	// it is not specified in Solé & Miramontes paper because
	// the are not using the FNN for PSO
	public static enum InfluenceModel {
		NEIGH_BEST, FIPS
	}

	// is the neighborhood a lattice or a torus?
	// the default is:
	//		1) a lattice for the functioning of the FNN
	//		2) a torus for the PSO 
	public static enum BoundaryModel {
		LATTICE, TORUS
	}

	// this is an issue for both the PSO and the FNN,
	// the assumption is that it will always be the same for both
	public static enum ActivityModel {
		ALL_NEURONS, ONLY_ACTIVE_NEURONS, SIMILAR_ACTIVITY_STATUS_NEURONS
	}

	
	
	// these PSO characteristics will almost always be torus and only-active-neurons,
	// so none of these has a for-loop for testing multiple possibilities
	public static BoundaryModel currentPSOBoundaryModel; 
	public static ActivityModel currentPSOActivityModel; 
	
	
	// these FNN characteristics will almost always be set to Moore, includeSelf, lattice, only-active-neurons.
	// so none of these has a for-loop for testing multiple possibilities
	public static Topology currentFNNTopology; 
	public static SelfModel currentFNNSelfModel; 
	public static BoundaryModel currentFNNBoundaryModel; 
	public static ActivityModel currentFNNActivityModel;


	
	
	// the usual PSO parameters
	public static double neighborhoodTheta;
	public static double personalTheta;
	public static double theta;
	public static double constrictionFactor;


	// data files
	// this one stores all the parameters of a run
	public static boolean summaryDataFileNeeded; 
	// this one stores min, max, and quartile info
	public static boolean intervalDataFileNeeded;
	// this one stores final function values for Mann-Whitney calculation
	public static boolean allFinalFuncValsDataFileNeeded;
	// this one stores final function values the plots of final function values by number of particles
	public static boolean finalFuncValsByNumParticlesDataFileNeeded;


	public static boolean useIndividualGains = false;
	public static double lowGain;
	public static double highGain;
	public static double gainIncrement;
	public static double gainDecrement;
	
	public static double totalNumGainIncreases = 0;
	public static double totalNumGainDecreases = 0;
	public static double sumGainValues = 0;
	
	public static boolean useNewSpontActMechanism = false;
	
	
	// FNN
//	public static int numIterationsDiscarded;
//	public static int numIterationsData;
//	public static int firstIterationData;

	

	public static void main(String[] args) {

		try {

			// create date string for screen output and output file names
			SimpleDateFormat dateformatter = new SimpleDateFormat("yyyy-MM-dd--hh-mm-ss-a");
			Calendar date = Calendar.getInstance();
			String dateString = dateformatter.format(date.getTime());

			// to output window
			System.out.println("RUNNING CODE ON " + dateString + "\n");			


			// create a PrintWriter that writes to System.out so I can use methods that
			// write to a file to write to the output window; gets closed at the end of
			// this (main) method
			PrintWriter outputWindow = new PrintWriter(System.out);

			String folder = "results/";
			
			summaryDataFileNeeded = true;    
			// if summary data file needed, name it, create it, and print the date to it;
			// this will contain summary data for *all* the runs that are done
			PrintWriter summaryDataFile = null;
			if (summaryDataFileNeeded) {
				String fileDataDescription = folder + "SUMMARY-DATA-";
				summaryDataFile = new PrintWriter(new FileWriter(fileDataDescription + "-" + dateString));
				summaryDataFile.println(dateString);
				summaryDataFile.println();
			}

			allFinalFuncValsDataFileNeeded = true;
			// if summary data file needed, name it, create it, and print the date to it;
			// this will contain summary data for *all* the runs that are done
			PrintWriter allFinalFuncValsDataFile = null;
			if (allFinalFuncValsDataFileNeeded) {
				String fileDataDescription = folder + "FINAL-FUNC-VALS-";
				allFinalFuncValsDataFile = new PrintWriter(new FileWriter(fileDataDescription + "-" + dateString));
				allFinalFuncValsDataFile.println(dateString);
				allFinalFuncValsDataFile.println();
			}

			// NOTE:   
			// NOTE: if we are only testing one level of #particles, this is automatically set to false
			//        (right after the section where we select the #particles levels to be tested)
			finalFuncValsByNumParticlesDataFileNeeded = true;
			// if files with final function values for each number of particles are needed, 
			// they will be named and created later because there will be a separate file 
			// for each run done, but declare a PrintWriter for the data files now;
			// create arrays to:
			// 1) save the names of the output files that need to go 
			//    into the gnuplot script created at the end of the run, and 
			// 2) save the final function values
			PrintWriter finalFuncValByNumParticlesDataFile = null;
			String[][][][][][][][] finalFuncValsByNumParticlesOutputFileNames = null;
			double[][][][][][][][][] finalFuncValues = null;
			if (finalFuncValsByNumParticlesDataFileNeeded) {
				finalFuncValsByNumParticlesOutputFileNames = new String[8][8][2][2][8][10][4][4];
				finalFuncValues = new double[8][8][8][2][2][8][10][4][4];
			}
			
			
			intervalDataFileNeeded = true;
			// if files with interval data are needed, they will be named and created later 
			// because there will be a separate file for each run done, but declare a
			// PrintWriter for the data files now;
			// create an array to save the names of the output files that need to go 
			// into the gnuplot script created at the end of the run
			PrintWriter intervalDataFile = null;
			String[][][][][][][][][] iterationOutputFileNames = null;
			if (intervalDataFileNeeded) {
				iterationOutputFileNames = new String[8][8][8][2][2][8][10][4][4];
			}

			
			boolean gnuplotScriptsByGainsNeeded = true;
			
			
			// this has never gotten used!!
			useSeedForRand = false;
			seed = 8778;   
			rand = useSeedForRand? new Random(seed): new Random();


			// standard PSO parameters
			neighborhoodTheta = 2.05;		
			personalTheta = 2.05;
			theta = neighborhoodTheta + personalTheta;
			constrictionFactor = 2.0 / (theta - 2.0 + Math.sqrt(theta*theta - 4.0*theta));


			// how many runs?
			totalNumRuns = 50;  


			// PSO (not FNN) TOPOLOGY
			Topology[] topologiesList = { Topology.GBEST, Topology.RING, Topology.vonNEUMANN, Topology.MOORE, 
					Topology.FNN_GBEST, Topology.FNN_RING, Topology.FNN_vonNEUMANN, Topology.FNN_MOORE };
			boolean doGBEST = false;
			boolean doRING = false;
			boolean doVonNEUMANN = true;
			boolean doMOORE =  true;
			boolean doFNN_GBEST = false;
			boolean doFNN_RING = false;
			boolean doFNN_VonNEUMANN = true;
			boolean doFNN_MOORE =  true;
			boolean[] doTopology = { doGBEST, doRING, doVonNEUMANN, doMOORE, doFNN_GBEST, doFNN_RING, doFNN_VonNEUMANN, doFNN_MOORE };

			// SELF MODEL
			SelfModel[] selfModelsList = { SelfModel.INCLUDE_SELF, SelfModel.NOT_INCLUDE_SELF };
			boolean doINCLUDE_SELF = true;
			boolean doNOT_INCLUDE_SELF = true;
			boolean[] doSelfModel = { doINCLUDE_SELF, doNOT_INCLUDE_SELF };

			// INFLUENCE MODEL
			InfluenceModel[] influenceModelsList = { InfluenceModel.NEIGH_BEST, InfluenceModel.FIPS };
			boolean doNEIGH_BEST = true;
			boolean doFIPS = true;
			boolean[] doInfluenceModel = { doNEIGH_BEST, doFIPS };


			// for the PSO
			currentPSOBoundaryModel = PSO.BoundaryModel.TORUS; 
			currentPSOActivityModel = PSO.ActivityModel.ALL_NEURONS; 
			
			// for the FNN
			currentFNNTopology = PSO.Topology.FNN_MOORE; 
			currentFNNSelfModel = PSO.SelfModel.INCLUDE_SELF; 
			currentFNNBoundaryModel = PSO.BoundaryModel.LATTICE; 
			currentFNNActivityModel = PSO.ActivityModel.ALL_NEURONS; 

			
			// from when these had for-loops to test multiple possibilities
			// BOUNDARY
//			PSO.BoundaryModel[] boundariesList = { BoundaryModel.LATTICE, BoundaryModel.TORUS };
//			boolean doLattice = true;
//			boolean doTorus = false;
//			boolean[] doBoundaryModel = { doLattice, doTorus };

			// ACTIVITY
//			PSO.ActivityModel[] activitiesList = { ActivityModel.ALL_NEURONS, ActivityModel.ONLY_ACTIVE_NEURONS };
//			boolean doAllNeurons = false;
//			boolean doOnlyActiveNeurons = true;
//			boolean[] doActivityModel = { doAllNeurons, doOnlyActiveNeurons };

			
			//  FNN parameters *************************  BEGIN

			// LATTICE SIZE
			int[] latticeSizesList = { 4, 5, 6, 7, 8, 9, 10, 12 }; 
			boolean doLatticeSize1 = false;   
			boolean doLatticeSize2 = false;   
			boolean doLatticeSize3 = false;   
			boolean doLatticeSize4 = false;    
			boolean doLatticeSize5 = false;    
			boolean doLatticeSize6 = true;    
			boolean doLatticeSize7 = false;    
			boolean doLatticeSize8 = false;    
			boolean[] doLatticeSize = { doLatticeSize1, doLatticeSize2, doLatticeSize3, doLatticeSize4, doLatticeSize5, doLatticeSize6, doLatticeSize7, doLatticeSize8 };

			// NOTE: density is not explicitly set; it is numParticles / (latticeSize * latticeSize)

			// GAIN
//			double[] gainsList = { 0.08, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 }; 
//			double[] gainsList = { 0.2, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 }; 
//			double[] gainsList = { 0.14, 0.20 }; 
//			double[] gainsList = { 0.25 }; 
//			double[] gainsList = { 0.206 }; 
//			double[] gainsList = { 0.2 }; 
//			double[] gainsList = { 100.0 }; 
			double[] gainsList = { 0.0 }; 
//			double[] gainsList = { 0.0001, 0.01, 0.14, 0.16, 0.20, 0.21, 0.23 }; 

			boolean doGain1 = true;   
			boolean doGain2 = false;   
			boolean doGain3 = false;   
			boolean doGain4 = false;    
			boolean doGain5 = false;    
			boolean doGain6 = false;    
			boolean doGain7 = false;    
			boolean doGain8 = false;    
			boolean doGain9 = false;    
			boolean doGain10 = false;    
			boolean[] doGain = { doGain1, doGain2, doGain3, doGain4, doGain5, doGain6, doGain7, doGain8, doGain9, doGain10 };

			// ***************************************************
			// ***************************************************
			// ***************************************************
			useIndividualGains = false;  
			
			lowGain = 0.1;
			highGain = 0.5;
			
			gainIncrement = 0.001;
			gainDecrement = 0.0001;
			
			totalNumGainIncreases = 0;
			totalNumGainDecreases = 0;
			sumGainValues = 0;

			// ***************************************************
			// ***************************************************
			// ***************************************************

			// ###################################################
			// ###################################################
			// ###################################################

			useNewSpontActMechanism = false;
			
			FluidNN.INITIAL_ACTIVATION_LOW_LEVEL = 0.0;
			FluidNN.INITIAL_ACTIVATION_HIGH_LEVEL = 1.0;
			FluidNN.INITIAL_ACTIVATION_RANGE = FluidNN.INITIAL_ACTIVATION_HIGH_LEVEL - FluidNN.INITIAL_ACTIVATION_LOW_LEVEL;
						
			// ###################################################
			// ###################################################
			// ###################################################

			
			// SPONTANEOUS ACTIVATION LEVEL
			double[] spontActLevelsList = { 1e-6, 0.15, 0.2, 0.25 }; 
			boolean doSpontActLevel1 = true;   
			boolean doSpontActLevel2 = false;   
			boolean doSpontActLevel3 = false;   
			boolean doSpontActLevel4 = false;    
			boolean[] doSpontActLevel = { doSpontActLevel1, doSpontActLevel2, doSpontActLevel3, doSpontActLevel4 };

			// SPONTANEOUS ACTIVATION PROBABILITY
			double[] spontActProbsList = { 1e-4, 1e-3, 1e-2, 0.025 }; 
			boolean doSpontActProb1 = true;   
			boolean doSpontActProb2 = false;   
			boolean doSpontActProb3 = false;   
			boolean doSpontActProb4 = false;    
			boolean[] doSpontActProb = { doSpontActProb1, doSpontActProb2, doSpontActProb3, doSpontActProb4 };

			// thresholds
			// always 0.0 in S&M tests
			double sumNeighborActivationsThreshold = 0.0;
			// always 1e-16 in S&M tests
			double activationThreshold = 1e-16; 
			
			//  FNN parameters *************************  END


			// FUNCTION
			boolean doSchwefel = false;
			boolean doRastrigin = true;  
			boolean doAckley = true;     
			boolean doGriewank =  true;   
			boolean doPenalFunc1 = true;   
			boolean doPenalFunc2 = true;   
			boolean doSphere = true;
			boolean doRosenbrock = true;    
			boolean[] doFunction = { doSchwefel, doRastrigin, doAckley, doGriewank, doPenalFunc1, doPenalFunc2, doSphere, doRosenbrock };

			// NOTE: the location of the optimum is shifted by a random amount that is a fraction of the 
			// distance from the optimum location to the edge of the search space.  in the case that the 
			// optimum location is not the same distance from both edges, it is a fraction of the shorter 
			// distance


			// DIMENSIONS
			int[] numDimsList = { 3, 10, 30, 50 };
			boolean doNumDims1 = false;  
			boolean doNumDims2 = false;   
			boolean doNumDims3 = true;  
			boolean doNumDims4 = false;  
			boolean[] doNumDims = { doNumDims1, doNumDims2, doNumDims3, doNumDims4 };


			// NUM PARTICLES
			int[] numParticlesList = { 4, 10, 20, 30, 40, 50, 60, 70 }; 
			boolean doNumParticles1 = false;   
			boolean doNumParticles2 = false;   
			boolean doNumParticles3 = false;   
			boolean doNumParticles4 = false;    
			boolean doNumParticles5 = true;    
			boolean doNumParticles6 = false;    
			boolean doNumParticles7 = false;    
			boolean doNumParticles8 = false;    
			boolean[] doNumParticles = { doNumParticles1, doNumParticles2, doNumParticles3, doNumParticles4, doNumParticles5, doNumParticles6, doNumParticles7, doNumParticles8 };

			// ********************************************************************************************************************************
			// ********************************************************************************************************************************
			// ********************************************************************************************************************************
			// if we are only testing one level of #particles, don't bother collecting data to be plotted against the number of particles
			int numParticleLevelsTested = 0;
			for (int i = 0 ; i < doNumParticles.length ; ++i) {
				if (doNumParticles[i])
					++numParticleLevelsTested;
			}
			if (numParticleLevelsTested == 1) {
				finalFuncValsByNumParticlesDataFileNeeded = false;
			}
			// ********************************************************************************************************************************
			// ********************************************************************************************************************************
			// ********************************************************************************************************************************


			// we need to specify the number of rows and columns for the von Neumann
			// and Moore neighborhoods for each number of particles
			numRowsVonNeumannAndMooreList = new int[numParticlesList.length];
			numColsVonNeumannAndMooreList = new int[numParticlesList.length];
			
			numRowsVonNeumannAndMooreList[0] = 2;
			numColsVonNeumannAndMooreList[0] = 2;

			numRowsVonNeumannAndMooreList[1] = 2;
			numColsVonNeumannAndMooreList[1] = 5;

			numRowsVonNeumannAndMooreList[2] = 4;
			numColsVonNeumannAndMooreList[2] = 5;

			numRowsVonNeumannAndMooreList[3] = 5;
			numColsVonNeumannAndMooreList[3] = 6;

			numRowsVonNeumannAndMooreList[4] = 5;
			numColsVonNeumannAndMooreList[4] = 8;

			numRowsVonNeumannAndMooreList[5] = 5;
			numColsVonNeumannAndMooreList[5] = 10;

			numRowsVonNeumannAndMooreList[6] = 6;
			numColsVonNeumannAndMooreList[6] = 10;

			numRowsVonNeumannAndMooreList[7] = 7;
			numColsVonNeumannAndMooreList[7] = 10;


			
			// NUM ITERATIONS AND NUM FUNCTION EVALUATIONS
			// this is were we start dealing with the fact that we might be counting iterations,
			// or we might be counting function evaluations; we need to have the ability to specify both
			useIterations = true;
			if (useIterations) {
				numItersPerOutputInterval = 100;
			}
			else {
				numFEsPerOutputInterval = 200;
			}

			// NUM ITERATIONS
			int[] numIterationsList = { 10000, 30000, 50000, 100000, 200000 }; 
			boolean doNumIterations1 = true;       //  10,000
			boolean doNumIterations2 = false;      //  30,000
			boolean doNumIterations3 = false;      //  50,000
			boolean doNumIterations4 = false;      // 100,000
			boolean doNumIterations5 = false;      // 200,000
			boolean[] doNumIterations = { doNumIterations1, doNumIterations2, doNumIterations3, doNumIterations4, doNumIterations5 };

			// NUM FUNCTION EVALUATIONS
			int[] numFEsList = { 10000, 30000, 50000, 100000, 200000 }; 
			boolean doNumFEs1 = false;      //  10,000
			boolean doNumFEs2 = false;      //  30,000
			boolean doNumFEs3 = false;      //  50,000
			boolean doNumFEs4 = false;      // 100,000
			boolean doNumFEs5 = false;      // 200,000
			boolean[] doNumFEs = { doNumFEs1, doNumFEs2, doNumFEs3, doNumFEs4, doNumFEs5 };



			// Here's where we start all the nested for-loops to run all the parameter combinations specified above


			// TOPOLOGY
			// ========
			for (int topoIndex = 0 ; topoIndex < topologiesList.length ; ++topoIndex) {

				if (!doTopology[topoIndex])
					continue;

				Topology currentPSOTopology = topologiesList[topoIndex]; 
				//			System.out.println("currentTopology = " + currentTopology);



				// INCLUDING SELF?
				// ===============
				for (int selfIndex = 0 ; selfIndex < selfModelsList.length ; ++selfIndex) {

					if (!doSelfModel[selfIndex])
						continue;

					SelfModel currentPSOSelfModel =  selfModelsList[selfIndex]; 
					//				System.out.println("                  currentSelfModel = " + currentSelfModel);



					// INFLUENCE MODEL
					// ===============
					for (int influenceIndex = 0 ; influenceIndex < influenceModelsList.length ; ++influenceIndex) {

						if (!doInfluenceModel[influenceIndex])
							continue;

						InfluenceModel currentPSOInfluenceModel = influenceModelsList[influenceIndex]; 
						//					System.out.println("                        currentInfluenceModel = " + currentInfluenceModel);



						// BOUNDARY MODEL
						// ==============
//						for (int boundIndex = 0 ; boundIndex < boundariesList.length ; ++boundIndex) {
//
//							if ((isFNNTopology(currentTopology) && !doBoundaryModel[boundIndex]) || (isStaticTopology(currentTopology) && boundIndex > 0))
//								continue;
//
//							FNN_BoundaryModel currentFNNBoundaryModel = boundariesList[boundIndex]; 
//							//						System.out.println("   currentFNNBoundaryModel = " + currentFNNBoundaryModel);



							// ACTIVITY MODEL
							// ==============
//							for (int activityIndex = 0 ; activityIndex < activitiesList.length ; ++activityIndex) {
//
//								if ((isFNNTopology(currentTopology) && !doActivityModel[activityIndex]) || (isStaticTopology(currentTopology) && activityIndex > 0))
//									continue;
//
//								ActivityModel currentFNNActivityModel = activitiesList[activityIndex]; 
//								//							System.out.println("      currentFNNActivityModel = " + currentFNNActivityModel);



								// LATTICE SIZE
								// ============
								for (int latticeIndex = 0 ; latticeIndex < latticeSizesList.length ; ++latticeIndex) {

									if ((isFNNTopology(currentPSOTopology) && !doLatticeSize[latticeIndex]) || (isStaticTopology(currentPSOTopology) && latticeIndex > 0))
										continue;

									int currentLatticeSideSize = latticeSizesList[latticeIndex]; 
									//								System.out.println("                     latticeSize = " + latticeSize);



									// GAIN
									// ========
									for (int gainIndex = 0 ; gainIndex < gainsList.length ; ++gainIndex) {

										if ((isFNNTopology(currentPSOTopology) && !doGain[gainIndex]) || (isStaticTopology(currentPSOTopology) && gainIndex > 0))
											continue;

										double currentGain = gainsList[gainIndex]; 
										//									System.out.println("         currentGain = " + currentGain);



										// SPONTANEOUS ACTIVATION LEVEL
										// ============================
										for (int spontLevelIndex = 0 ; spontLevelIndex < spontActLevelsList.length ; ++spontLevelIndex) {

											if ((isFNNTopology(currentPSOTopology) && !doSpontActLevel[spontLevelIndex]) || (isStaticTopology(currentPSOTopology) && spontLevelIndex > 0))
												continue;

											double currentSpontActLevel = spontActLevelsList[spontLevelIndex]; 
											//										System.out.println("            spontaneousActivationLevel = " + spontaneousActivationLevel);



											// SPONTANEOUS ACTIVATION PROBABILITY
											// ==================================
											for (int spontProbIndex = 0 ; spontProbIndex < spontActProbsList.length ; ++spontProbIndex) {

												if ((isFNNTopology(currentPSOTopology) && !doSpontActProb[spontProbIndex]) || (isStaticTopology(currentPSOTopology) && spontProbIndex > 0))
													continue;

												double currentSpontActProb = spontActProbsList[spontProbIndex]; 
												//											System.out.println("               spontaneousActivationProbability = " + spontaneousActivationProbability);



												// TEST FUNCTIONS
												// ==============
												for (int currentFunctionNum = 0 ; currentFunctionNum < doFunction.length ; ++currentFunctionNum) {

													if (!doFunction[currentFunctionNum])
														continue;

													//												System.out.println("                           currentFunctionNum = " + currentFunctionNum);


													// # DIMENSIONS
													// ============
													for (int numDimensionsIndex = 0 ; numDimensionsIndex < doNumDims.length ; ++numDimensionsIndex ) {

														if (!doNumDims[numDimensionsIndex])
															continue;

														int currentNumDimensions = numDimsList[numDimensionsIndex];

														
														
														
														
														String topologyAndInfluence = getPSOTopologyAndInfluence(currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel,
																currentPSOBoundaryModel, currentPSOActivityModel);
														
														String FNNparameters = "";
														if (isFNNTopology(currentPSOTopology)) {
															FNNparameters = getFNNDataString(currentLatticeSideSize, currentGain, currentSpontActLevel, currentSpontActProb);
														}
														String fileName = TestFunctions.getShortFunctionName(currentFunctionNum) + 
																		"-d" + currentNumDimensions + 
																		"-it" + totalNumIters + 
																		"-" + topologyAndInfluence + 
																		FNNparameters + 
																		".txt";

														if (finalFuncValsByNumParticlesDataFileNeeded) {
															finalFuncValByNumParticlesDataFile = new PrintWriter(new FileWriter(folder + "final-func-val-data-" + fileName));
															finalFuncValByNumParticlesDataFile.println("# " + dateString);
															finalFuncValByNumParticlesDataFile.println("# ");
															outputTestParameters(finalFuncValByNumParticlesDataFile, currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel, 
																	currentPSOBoundaryModel, currentPSOActivityModel, 
																	currentFunctionNum, currentNumDimensions, 999, totalNumRuns, totalNumFEs, totalNumIters,	
																	currentLatticeSideSize, currentGain, sumNeighborActivationsThreshold, activationThreshold, 
																	currentSpontActLevel, currentSpontActProb);		

															finalFuncValsByNumParticlesOutputFileNames[currentFunctionNum]
															                                           [topoIndex]
															                                            [selfIndex]
															                                             [influenceIndex]
															                                              [latticeIndex]
															                                               [gainIndex]
															                                                [spontLevelIndex]
															                                                 [spontProbIndex] = "final-func-val-data-" + fileName;
														}

														
														// # PARTICLES
														// ===========
														for (int numParticlesIndex = 0 ; numParticlesIndex < doNumParticles.length ; ++numParticlesIndex ) {

															if (!doNumParticles[numParticlesIndex])
																continue;

															int currentNumParticles = numParticlesList[numParticlesIndex];
															// ??????
															// currentNumParticles is a global so this info can be accessed in Particle.java to set
															// the iteration number of the currentSolution and the personal best solution
															currentNumParticles = numParticlesList[numParticlesIndex];

															


															// # ITERATIONS
															// ============
															int numIterationsIndex;
															int numFEsIndex;
															for (numIterationsIndex = 0, numFEsIndex = 0 ; numIterationsIndex < doNumIterations.length && numFEsIndex < doNumFEs.length ; ++numIterationsIndex, ++numFEsIndex ) {
																// the second condition = if we're using iterations we don't need to do a run for every level of FEs                            ??????????????????????????

																// if using iterations
																if (useIterations && !doNumIterations[numIterationsIndex])
																	continue;

																// if using function evaluations
																if (!useIterations && !doNumFEs[numFEsIndex])
																	continue;


																currentFENum = 0;
																currentIterNum = 0;

																// this puts everything into number of FEs terms so we can just use that for iteration control
																// and data output control
																if (useIterations) {
																	totalNumIters = numIterationsList[numIterationsIndex];

																	// to get that many iterations (not counting the numParticles FEs used to evaluate the
																	// particles when the swarm is created, i.e. iteration 0), we need an additional numParticles FEs
																	totalNumFEs = (totalNumIters * currentNumParticles) + currentNumParticles;

																	// numFEsPerIntervalForOutput needs to be calculated based on numIterationsPerIntervalForOutput
																	numFEsPerOutputInterval = currentNumParticles * numItersPerOutputInterval;
																	// we need this because the first numParticles FEs are in iteration 0 and shouldn't be counted 
																	// for *iterations* reporting
																	numInitialFEsIgnored = currentNumParticles;

																	numIntervalsForDataOutput = totalNumIters / numItersPerOutputInterval;

																	//								System.out.println("totalNumIterations = " + totalNumIterations);
																	//								System.out.println("numIterationsPerIntervalForOutput = " + numIterationsPerIntervalForOutput);
																	//								System.out.println("numParticles = " + numParticles);
																	//								System.out.println("totalNumFunctionEvaluations = " + totalNumFunctionEvaluations);
																	//								System.out.println("numFEsPerIntervalForOutput = " + numFEsPerIntervalForOutput);
																	//								System.out.println("numInitialFEsIgnored = " + numInitialFEsIgnored);
																	//								System.out.println("numIntervalsForDataOutput = " + numIntervalsForDataOutput);
																}

																else {
																	totalNumFEs = numFEsList[numFEsIndex];

																	// numFEsPerIntervalForOutput does not need to be calculated since it is assigned
																	// a value above where we indicate whether to use iterations or FEs

																	// no initial FEs are ignored because, unlike iterations, where we don't start counting iterations
																	// until after the swarm has been created and numParticles FEs have been used in the process,
																	// here we start counting FEs right away
																	numInitialFEsIgnored = 0;
																	numIntervalsForDataOutput = totalNumFEs / numFEsPerOutputInterval;
																	// would not expect the following to ever be an issue since we are setting both 
																	// totalNumFunctionEvaluations and numFEsPerIntervalForOutput, and it is unlikely
																	// that we would set them so that totalNumFunctionEvaluations is not a multiple of 
																	// numFEsPerIntervalForOutput, but if that is the case we would need an extra interval
																	// for data output to take care of the "extra" FEs
																	if (totalNumFEs % numFEsPerOutputInterval != 0)
																		++numIntervalsForDataOutput;
																}



																// FNN
																//										numIterations = 11000;
																//										numIterationsDiscarded = 1000;
																//										numIterationsData = numIterations - numIterationsDiscarded;
																//										firstIterationData = numIterationsDiscarded + 1;



																// create the DataOutput objects in iterationData once and for all and
																// just copy the values/errors in the globalBests returned by update into these 
																// DataOutput objects
																//
																// NOTE: we need to add 1 to numIntervalsForDataOutput because: 
																//       whether we are doing output by # of iterations or output by # of FEs, 
																//       the output with index 0,
																//       which is right after the swarm is created if counting iterations, 
																//           or before the first FE if counting FEs 
																//       is not wanted
																// perhaps a better way to look at it is that output interval indices are *not* 
																// zero-based, so we need that extra (useless) item in the index 0 spot
																DataOutput[][] intervalData = new DataOutput[numIntervalsForDataOutput+1][totalNumRuns];
																for (int i = 0 ; i < intervalData.length ; ++i) {
																	for (int r = 0 ; r < totalNumRuns ; ++r) {
																		intervalData[i][r] = new DataOutput();
																	}
																}
																// can create the IntervalSummaryData objects once and for all because we will just be
																// recalculating their contents after every series of runs
																IntervalSummaryData[] intervalSummaryData = new IntervalSummaryData[numIntervalsForDataOutput+1];
																for (int i = 0 ; i < intervalSummaryData.length ; ++i) {
																	intervalSummaryData[i] = new IntervalSummaryData();
																}
																


																// this will put EACH CASE:
																//	1) topology
																//  2) self model
																//  3) influence model
																//  4) function
																//  5) # dimensions 
																//  6) # particles
																//  7) # iterations
																// in its own file
																if (intervalDataFileNeeded) {
																	
																	topologyAndInfluence = getPSOTopologyAndInfluence(currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel,
																			currentPSOBoundaryModel, currentPSOActivityModel);
																	
																	FNNparameters = "";
																	if (isFNNTopology(currentPSOTopology)) {
																		FNNparameters = getFNNDataString(currentLatticeSideSize, currentGain, currentSpontActLevel, currentSpontActProb);
																	}
																	
																	fileName = TestFunctions.getShortFunctionName(currentFunctionNum) + 
																					"-d" + currentNumDimensions + 
																					"-p" + currentNumParticles + 
																					"-it" + totalNumIters + 
																					"-" + topologyAndInfluence + 
																					FNNparameters + 
																					".txt";
																	
																	intervalDataFile = new PrintWriter(new FileWriter(folder + "iter-data-" + fileName));
																	intervalDataFile.println("# " + dateString);
																	intervalDataFile.println("# ");
																	outputTestParameters(intervalDataFile, currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel, 
																			currentPSOBoundaryModel, currentPSOActivityModel, 
																			currentFunctionNum, currentNumDimensions, currentNumParticles, totalNumRuns, totalNumFEs, totalNumIters,	
																			currentLatticeSideSize, currentGain, sumNeighborActivationsThreshold, activationThreshold, 
																			currentSpontActLevel, currentSpontActProb);		

																	iterationOutputFileNames[currentFunctionNum]
																				                [numParticlesIndex]
																				                 [topoIndex]
																				                  [selfIndex]
																				                   [influenceIndex]
																				                    [latticeIndex]
																				                     [gainIndex]
																				                      [spontLevelIndex]
																				                       [spontProbIndex] = "iter-data-" + fileName;

																}



																if (allFinalFuncValsDataFileNeeded) {
																	allFinalFuncValsDataFile.println();
																	allFinalFuncValsDataFile.println();
																	allFinalFuncValsDataFile.println();
																	outputTestParameters(allFinalFuncValsDataFile, currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel, 
																			currentPSOBoundaryModel, currentPSOActivityModel, 
																			currentFunctionNum, currentNumDimensions, currentNumParticles, totalNumRuns, totalNumFEs, totalNumIters,	
																			currentLatticeSideSize, currentGain, sumNeighborActivationsThreshold, activationThreshold, 
																			currentSpontActLevel, currentSpontActProb);		
																	allFinalFuncValsDataFile.println();
																}
																
																


																// haven't been too interested in this yet....
																long startTimeAllRuns = System.currentTimeMillis();  


																// RUNS
																// ====
																for (currentRunNum = 0 ; currentRunNum < totalNumRuns ; ++currentRunNum) {


																	// for each run generate a random shift of the location of the optimum in that function's search space
																	double shiftVectorAmount = TestFunctions.SHIFT_RANGE[currentFunctionNum] * rand.nextDouble();
																	if (rand.nextDouble() < 0.5) 
																		shiftVectorAmount *= -1.0;

																	// finally, we can set the function parameters because we know the shift for that run
																	TestFunctions.setShiftVector(currentFunctionNum, new DoubleVector(currentNumDimensions, shiftVectorAmount));

																	// initialize numFunctionEvaluations here so that FEs that occur in creating the swarm will be counted
																	currentFENum = 0;

																	// create the swarm
																	// numParticlesIndex is sent because it needs to be used in the creation of the vonNeumann
																	// and Moore neighborhoods to access the arrays above that indicate the number of rows/cols
																	//
																	// need to send all FNN parameters (even though they will not be used for other topologies)
																	// because I am going to want to test FNN-PSO over a range of parameter values that is set 
																	// in PSO.java
																	Swarm swarm = new Swarm (currentNumParticles, numParticlesIndex, currentFunctionNum, currentNumDimensions, 
																			currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel,
																			currentLatticeSideSize, currentGain, sumNeighborActivationsThreshold, activationThreshold, 
																			currentSpontActLevel, currentSpontActProb, intervalData);



																	// ITERATIONS
																	// ==========
																	// - we are going solely by the number of FEs used;
																	// - we calculate the number of FEs we need to get the number of iterations we want, and stop when we have used those up;
																	// - so we are only keeping track of iterations for the purpose of things like deciding when to restructure a 
																	//		dynamic topology (which is not part of the Basic-PSO code, but will be added in subsequent proejects);
																	// - NOTE: currentFunctionEvaluationNum is initialized to 0 above before the swarm is created because some FEs
																	//         will be used when the particles are created and evaluated, and it is incremented in TestFunctions 
																	//         when the function is evaluated in the evalWithError method

																	for ( ; currentFENum < totalNumFEs ; ) {

																		// one round of asynchronous updates is an iteration;
																		// currentIterationNum starts at 0, so needs to be incremented *before* the call to asynchronousUpdate
																		++currentIterNum;

//																		swarm.update(currentFunctionNum, currentTopology, currentSelfModel, currentInfluenceModel, 
//																				currentBoundaryModel, currentActivityModel, intervalData);
																		swarm.update(currentFunctionNum, currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel, intervalData);

																	}  // END ITERATIONS FOR-LOOP 


																	double finalFuncValue = swarm.globalBest.getFunctionValue();
																	finalFuncValues[currentFunctionNum]
																	                [numParticlesIndex]
																	                 [topoIndex]
																	                  [selfIndex]
																	                   [influenceIndex]
																	                    [latticeIndex]
																	                     [gainIndex]
																	                      [spontLevelIndex]
																	                       [spontProbIndex] = finalFuncValue;
																	if (allFinalFuncValsDataFileNeeded) {
																		allFinalFuncValsDataFile.println(finalFuncValue);
																	}
																	
																}  // END RUNS FOR-LOOP 


																// TIMING STATS
																long endTimeAllRuns = System.currentTimeMillis();
																double secondsPerRun = ((endTimeAllRuns - startTimeAllRuns) / 1000.0) / totalNumRuns;

																
																// CALCULATIONS FOR OUTPUT
																calculateSummaryDataOverRuns(intervalData, intervalSummaryData, currentFunctionNum, currentNumDimensions);

																
																// SUMMARY DATA OUTPUT TO OUTPUT WINDOW
																outputTestParameters(outputWindow, currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel, 
																		currentPSOBoundaryModel, currentPSOActivityModel, 
																		currentFunctionNum, currentNumDimensions, currentNumParticles, totalNumRuns, totalNumFEs, totalNumIters,	
																		currentLatticeSideSize, currentGain, sumNeighborActivationsThreshold, activationThreshold, 
																		currentSpontActLevel, currentSpontActProb);		

																if (isFNNTopology(currentPSOTopology) && useIndividualGains) {
																	outputWindow.println("--------------------------------------------------------");
																	outputWindow.println("totalNumGainIncreases = " + totalNumGainIncreases);
																	outputWindow.println("totalNumGainDecreases = " + totalNumGainDecreases);
																	outputWindow.println("sumGainValues = " + sumGainValues);
																	double totalNumberGains = totalNumGainIncreases + totalNumGainDecreases;
																	outputWindow.println("average gain = " + (sumGainValues / totalNumberGains));
																	outputWindow.println("--------------------------------------------------------");
																	outputWindow.println("#");
																}

																outputSummaryData(outputWindow, intervalSummaryData, secondsPerRun);								

																
																// SUMMARY DATA OUTPUT TO FILE
																// NOTE:  the file is not closed here because the summary data for *all* cases is put in this file, 
																//        so it's not closed until *all* loops are done
																if (summaryDataFileNeeded) {
																	outputTestParameters(summaryDataFile, currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel, 
																			currentPSOBoundaryModel, currentPSOActivityModel, 
																			currentFunctionNum, currentNumDimensions, currentNumParticles, totalNumRuns, totalNumFEs, totalNumIters,	
																			currentLatticeSideSize, currentGain, sumNeighborActivationsThreshold, activationThreshold, 
																			currentSpontActLevel, currentSpontActProb);		

																	if (isFNNTopology(currentPSOTopology) && useIndividualGains) {
																		summaryDataFile.println("--------------------------------------------------------");
																		summaryDataFile.println("totalNumGainIncreases = " + totalNumGainIncreases);
																		summaryDataFile.println("totalNumGainDecreases = " + totalNumGainDecreases);
																		summaryDataFile.println("sumGainValues = " + sumGainValues);
																		double totalNumberGains = totalNumGainIncreases + totalNumGainDecreases;
																		summaryDataFile.println("average gain = " + (sumGainValues / totalNumberGains));
																		summaryDataFile.println("--------------------------------------------------------");
																		summaryDataFile.println("#");
																	}

																	outputSummaryData(summaryDataFile, intervalSummaryData, secondsPerRun);								

																}

																
																if (allFinalFuncValsDataFileNeeded) {
																	allFinalFuncValsDataFile.println();
																	outputSummaryData(allFinalFuncValsDataFile, intervalSummaryData, secondsPerRun);								
																}

																
																if (finalFuncValsByNumParticlesDataFileNeeded) {
																	outputFinalFuncValue(finalFuncValByNumParticlesDataFile, intervalSummaryData, currentNumParticles);
																}
																
																if (useIndividualGains) {
																	totalNumGainIncreases = 0;
																	totalNumGainDecreases = 0;
																	sumGainValues = 0;
																}

																
																// INTERVAL DATA OUTPUT TO FILE
																if (intervalDataFileNeeded) {
																	outputIterationData(intervalDataFile, intervalSummaryData);
																	// close the file here because we create a separate file for every case (function, #dimensions, #particles, etc.)
																	intervalDataFile.close();
																}


															}  // END NUM ITERATIONS


														}  // END NUM PARTICLES FOR-LOOP 

														
														if (finalFuncValsByNumParticlesDataFileNeeded) {
															finalFuncValByNumParticlesDataFile.close();
														}
														


													}  // END NUM DIMENSIONS FOR-LOOP 		


												} // END FUNCTION NUM FOR-LOOP 


											}  // END SPONTANEOUS ACTIVATION PROBABILITY FOR-LOOP


										}  // END SPONTANEOUS ACTIVATION LEVEL FOR-LOOP


									}  // END GAIN FOR-LOOP  


								} // END LATTICE SIZE FOR-LOOP     


//							} // END ACTIVITY FOR-LOOP      


//						} // END BOUNDARY FOR-LOOP     


					} // END INFLUENCE MODEL FOR-LOOP       


				} // INCLUDE SELF FOR-LOOP        


			} // END TOPOLOGY FOR-LOOP 


			// close this file after *all* cases are done running (since 
			// output from all cases goes to this file)
			if (summaryDataFileNeeded) {
				summaryDataFile.close();   
			}
			

			if (allFinalFuncValsDataFileNeeded) {
				allFinalFuncValsDataFile.close();
			}


			// close the gnuplot script file after *all* cases are done running (since 
			// the script has a "plot" line for all topos for each other combination of params)
			if (intervalDataFileNeeded) {
				
				PrintWriter gnuplotScriptFile = new PrintWriter(new FileWriter(folder + "0-gnuplot-script-intervals-" + dateString + ".txt"));
				gnuplotScriptFile.println("#  " + dateString + "\n\n");

				gnuplotScriptFile.println("set logscale y \n\n");

				for (int currentFunctionNum = 0 ; currentFunctionNum < doFunction.length ; ++currentFunctionNum) {
					if (!doFunction[currentFunctionNum])
						continue;

					for (int numPartsIndex = 0 ; numPartsIndex < numParticlesList.length ; ++numPartsIndex) {
						if (!doNumParticles[numPartsIndex])
							continue;

						for (int latticeIndex = 0 ; latticeIndex < latticeSizesList.length ; ++latticeIndex) {
							if (!doLatticeSize[latticeIndex])
								continue;

							for (int gainIndex = 0 ; gainIndex < gainsList.length ; ++gainIndex) {
								if (!doGain[gainIndex])
									continue;

								for (int spontLevelIndex = 0 ; spontLevelIndex < spontActLevelsList.length ; ++spontLevelIndex) {
									if (!doSpontActLevel[spontLevelIndex])
										continue;

									for (int spontProbIndex = 0 ; spontProbIndex < spontActProbsList.length ; ++spontProbIndex) {
										if (!doSpontActProb[spontProbIndex])
											continue;

										
										// ********************************************
										//   ALL of the cases
										
										gnuplotScriptFile.println("# " + TestFunctions.getFunctionName(currentFunctionNum));
										gnuplotScriptFile.println("# num particles = " + numParticlesList[numPartsIndex]);
										gnuplotScriptFile.println("# lattice size = " + latticeSizesList[latticeIndex]);
										gnuplotScriptFile.println("# gain = " + gainsList[gainIndex]);
										gnuplotScriptFile.println("# spont act level = " + spontActLevelsList[spontLevelIndex]);
										gnuplotScriptFile.println("# spont prob level = " + spontActProbsList[spontProbIndex] );
										gnuplotScriptFile.println();
										gnuplotScriptFile.println("# ALL CASES");
										gnuplotScriptFile.println();
										gnuplotScriptFile.print("#plot [][] ");

										boolean firstFile = true;
										
										// print out s-pso cases
										for (int topoIndex = 0 ; topoIndex < topologiesList.length ; ++topoIndex) {
											if (!doTopology[topoIndex] || topoIndex >= 4) 
												continue;

											for (int selfIndex = 0 ; selfIndex < selfModelsList.length ; ++selfIndex) {
												if (!doSelfModel[selfIndex])
													continue;

												for (int influenceIndex = 0 ; influenceIndex < influenceModelsList.length ; ++influenceIndex) {
													if (!doInfluenceModel[influenceIndex])
														continue;

													String thisFileName =  iterationOutputFileNames[currentFunctionNum]
													                                       [numPartsIndex]
													                                        [topoIndex]
													                                         [selfIndex]
													                                          [influenceIndex]
													                                           [0]
													                                            [0]
													                                             [0]
													                                              [0];
													if (thisFileName != null) {
														if (firstFile) {
															gnuplotScriptFile.print("\"" + thisFileName + "\" using 1:2 with linesp");
															firstFile = false;
														}
														else {
															gnuplotScriptFile.print(", \"" + thisFileName + "\" using 1:2 with linesp");															
														}
													}
												}
											}
										}

										// print out fnn-pso cases
										for (int topoIndex = 0 ; topoIndex < topologiesList.length ; ++topoIndex) {
											if (!doTopology[topoIndex] || topoIndex < 4) 
												continue;

											for (int selfIndex = 0 ; selfIndex < selfModelsList.length ; ++selfIndex) {
												if (!doSelfModel[selfIndex])
													continue;

												for (int influenceIndex = 0 ; influenceIndex < influenceModelsList.length ; ++influenceIndex) {
													if (!doInfluenceModel[influenceIndex])
														continue;


													String thisFileName =  iterationOutputFileNames[currentFunctionNum]
													                                       [numPartsIndex]
													                                        [topoIndex]
													                                         [selfIndex]
													                                          [influenceIndex]
													                                           [latticeIndex]
													                                            [gainIndex]
													                                             [spontLevelIndex]
													                                              [spontProbIndex];
													if (thisFileName != null) {
														if (firstFile) {
															gnuplotScriptFile.print("\"" + thisFileName + "\" using 1:2 with linesp");
															firstFile = false;
														}
														else {
															gnuplotScriptFile.print(", \"" + thisFileName + "\" using 1:2 with linesp");															
														}
													}
												}
											}
										}
										// close the "plot"
										gnuplotScriptFile.print("\n");

										//   ALL of the cases
										// ********************************************

										
										// skip 2 lines
										gnuplotScriptFile.print("\n\n");

										
										// ++++++++++++++++++++++++++++++++++++++++++++
										//   just yesSelf cases
										
										gnuplotScriptFile.println("# JUST yesSelf CASES");
										gnuplotScriptFile.println();
										gnuplotScriptFile.print("#plot [][] ");

										firstFile = true;
										
										// print out s-pso cases
										for (int topoIndex = 0 ; topoIndex < topologiesList.length ; ++topoIndex) {
											if (!doTopology[topoIndex] || topoIndex >= 4) 
												continue;

											for (int selfIndex = 0 ; selfIndex < selfModelsList.length ; ++selfIndex) {
												if (!doSelfModel[selfIndex])
													continue;

												for (int influenceIndex = 0 ; influenceIndex < influenceModelsList.length ; ++influenceIndex) {
													if (!doInfluenceModel[influenceIndex])
														continue;

													String thisFileName =  iterationOutputFileNames[currentFunctionNum]
													                                       [numPartsIndex]
													                                        [topoIndex]
													                                         [selfIndex]
													                                          [influenceIndex]
													                                           [0]
													                                            [0]
													                                             [0]
													                                              [0];
													if (thisFileName != null) {
														if (thisFileName.contains("yesSELF")) {
															if (firstFile) {
																gnuplotScriptFile.print("\"" + thisFileName + "\" using 1:2 with linesp");
																firstFile = false;
															}
															else {
																gnuplotScriptFile.print(", \"" + thisFileName + "\" using 1:2 with linesp");															
															}
														}
													}
												}
											}
										}

										// print out fnn-pso cases
										for (int topoIndex = 0 ; topoIndex < topologiesList.length ; ++topoIndex) {
											if (!doTopology[topoIndex] || topoIndex < 4) 
												continue;

											for (int selfIndex = 0 ; selfIndex < selfModelsList.length ; ++selfIndex) {
												if (!doSelfModel[selfIndex])
													continue;

												for (int influenceIndex = 0 ; influenceIndex < influenceModelsList.length ; ++influenceIndex) {
													if (!doInfluenceModel[influenceIndex])
														continue;


													String thisFileName =  iterationOutputFileNames[currentFunctionNum]
													                                       [numPartsIndex]
													                                        [topoIndex]
													                                         [selfIndex]
													                                          [influenceIndex]
													                                           [latticeIndex]
													                                            [gainIndex]
													                                             [spontLevelIndex]
													                                              [spontProbIndex];
													if (thisFileName != null) {
														if (thisFileName.contains("yesSELF")) {
															if (firstFile) {
																gnuplotScriptFile.print("\"" + thisFileName + "\" using 1:2 with linesp");
																firstFile = false;
															}
															else {
																gnuplotScriptFile.print(", \"" + thisFileName + "\" using 1:2 with linesp");															
															}
														}
													}
												}
											}
										}
										// close the "plot"
										gnuplotScriptFile.print("\n");
										

										gnuplotScriptFile.print("\n\n\n\n");
										gnuplotScriptFile.println("####################################################################################");
										gnuplotScriptFile.print("\n\n\n\n");

										//   just yesSelf cases
										// ++++++++++++++++++++++++++++++++++++++++++++

										
										
										
										
										
									}
								}
							}
						}
					}
				}

				gnuplotScriptFile.close();
			}

			
			// ===============================================================================================================================================
			// ===============================================================================================================================================
			// ===============================================================================================================================================
			// ===============================================================================================================================================
			// ===============================================================================================================================================
			// ===============================================================================================================================================
			
			
			
			if (finalFuncValsByNumParticlesDataFileNeeded) {
				
				PrintWriter gnuplotScriptFile = new PrintWriter(new FileWriter(folder + "0-gnuplot-script-final-func-values-" + dateString + ".txt"));
				gnuplotScriptFile.println("#  " + dateString + "\n\n");

				gnuplotScriptFile.println("set logscale y \n\n");

				for (int currentFunctionNum = 0 ; currentFunctionNum < doFunction.length ; ++currentFunctionNum) {
					if (!doFunction[currentFunctionNum])
						continue;

//					for (int numPartsIndex = 0 ; numPartsIndex < numParticlesList.length ; ++numPartsIndex) {
//						if (!doNumParticles[numPartsIndex])
//							continue;

						for (int latticeIndex = 0 ; latticeIndex < latticeSizesList.length ; ++latticeIndex) {
							if (!doLatticeSize[latticeIndex])
								continue;

							for (int gainIndex = 0 ; gainIndex < gainsList.length ; ++gainIndex) {
								if (!doGain[gainIndex])
									continue;

								for (int spontLevelIndex = 0 ; spontLevelIndex < spontActLevelsList.length ; ++spontLevelIndex) {
									if (!doSpontActLevel[spontLevelIndex])
										continue;

									for (int spontProbIndex = 0 ; spontProbIndex < spontActProbsList.length ; ++spontProbIndex) {
										if (!doSpontActProb[spontProbIndex])
											continue;

										
										
										// ********************************************
										//   ALL of the cases
										
										gnuplotScriptFile.println("# " + TestFunctions.getFunctionName(currentFunctionNum));
//										gnuplotScriptFile.println("# num particles = " + numParticlesList[numPartsIndex]);
										gnuplotScriptFile.println("# lattice size = " + latticeSizesList[latticeIndex]);
										gnuplotScriptFile.println("# gain = " + gainsList[gainIndex]);
										gnuplotScriptFile.println("# spont act level = " + spontActLevelsList[spontLevelIndex]);
										gnuplotScriptFile.println("# spont prob level = " + spontActProbsList[spontProbIndex] );
										gnuplotScriptFile.println();
										gnuplotScriptFile.println("# ALL CASES");
										gnuplotScriptFile.println();
										gnuplotScriptFile.print("#plot [][] ");

										boolean firstFile = true;
										
										// print out s-pso cases
//										for (int numPartsIndex = 0 ; numPartsIndex < numParticlesList.length ; ++numPartsIndex) {
//											if (!doNumParticles[numPartsIndex])
//												continue;

										for (int topoIndex = 0 ; topoIndex < topologiesList.length ; ++topoIndex) {
											if (!doTopology[topoIndex] || topoIndex >= 4) 
												continue;

											for (int selfIndex = 0 ; selfIndex < selfModelsList.length ; ++selfIndex) {
												if (!doSelfModel[selfIndex])
													continue;

												for (int influenceIndex = 0 ; influenceIndex < influenceModelsList.length ; ++influenceIndex) {
													if (!doInfluenceModel[influenceIndex])
														continue;

													String thisFileName =  finalFuncValsByNumParticlesOutputFileNames[currentFunctionNum]
//													                                       [numPartsIndex]
													                                        [topoIndex]
													                                         [selfIndex]
													                                          [influenceIndex]
													                                           [0]
													                                            [0]
													                                             [0]
													                                              [0];
													if (thisFileName != null) {
														if (firstFile) {
															gnuplotScriptFile.print("\"" + thisFileName + "\" using 1:2 with linesp");
															firstFile = false;
														}
														else {
															gnuplotScriptFile.print(", \"" + thisFileName + "\" using 1:2 with linesp");															
														}
													}
												}
											}
										}
//										}

										// print out fnn-pso cases
//										for (int numPartsIndex = 0 ; numPartsIndex < numParticlesList.length ; ++numPartsIndex) {
//											if (!doNumParticles[numPartsIndex])
//												continue;

										for (int topoIndex = 0 ; topoIndex < topologiesList.length ; ++topoIndex) {
											if (!doTopology[topoIndex] || topoIndex < 4) 
												continue;

											for (int selfIndex = 0 ; selfIndex < selfModelsList.length ; ++selfIndex) {
												if (!doSelfModel[selfIndex])
													continue;

												for (int influenceIndex = 0 ; influenceIndex < influenceModelsList.length ; ++influenceIndex) {
													if (!doInfluenceModel[influenceIndex])
														continue;


													String thisFileName =  finalFuncValsByNumParticlesOutputFileNames[currentFunctionNum]
//													                                       [numPartsIndex]
													                                        [topoIndex]
													                                         [selfIndex]
													                                          [influenceIndex]
													                                           [latticeIndex]
													                                            [gainIndex]
													                                             [spontLevelIndex]
													                                              [spontProbIndex];
													if (thisFileName != null) {
														if (firstFile) {
															gnuplotScriptFile.print("\"" + thisFileName + "\" using 1:2 with linesp");
															firstFile = false;
														}
														else {
															gnuplotScriptFile.print(", \"" + thisFileName + "\" using 1:2 with linesp");															
														}
													}
												}
											}
										}
										// close the "plot"
										gnuplotScriptFile.print("\n");

//										}

										//   ALL of the cases
										// ********************************************

										
										// skip 2 lines
										gnuplotScriptFile.print("\n\n");

										
										
										// ++++++++++++++++++++++++++++++++++++++++++++
										//   just s-pso cases


										gnuplotScriptFile.println("# JUST s-pso CASES");
										gnuplotScriptFile.println();
										gnuplotScriptFile.print("#plot [][] ");


										firstFile = true;
										
										// print out s-pso cases
//										for (int numPartsIndex = 0 ; numPartsIndex < numParticlesList.length ; ++numPartsIndex) {
//											if (!doNumParticles[numPartsIndex])
//												continue;

										for (int topoIndex = 0 ; topoIndex < topologiesList.length ; ++topoIndex) {
											if (!doTopology[topoIndex] || topoIndex >= 4) 
												continue;

											for (int selfIndex = 0 ; selfIndex < selfModelsList.length ; ++selfIndex) {
												if (!doSelfModel[selfIndex])
													continue;

												for (int influenceIndex = 0 ; influenceIndex < influenceModelsList.length ; ++influenceIndex) {
													if (!doInfluenceModel[influenceIndex])
														continue;

													String thisFileName =  finalFuncValsByNumParticlesOutputFileNames[currentFunctionNum]
//													                                       [numPartsIndex]
													                                        [topoIndex]
													                                         [selfIndex]
													                                          [influenceIndex]
													                                           [0]
													                                            [0]
													                                             [0]
													                                              [0];
													if (thisFileName != null) {
														if (firstFile) {
															gnuplotScriptFile.print("\"" + thisFileName + "\" using 1:2 with linesp");
															firstFile = false;
														}
														else {
															gnuplotScriptFile.print(", \"" + thisFileName + "\" using 1:2 with linesp");															
														}
													}
												}
											}
										}
										// close the "plot"
										gnuplotScriptFile.print("\n");

//										}

										//   just s-pso cases
										// ++++++++++++++++++++++++++++++++++++++++++++

										
										// skip 2 lines
										gnuplotScriptFile.print("\n\n");

										
										// ++++++++++++++++++++++++++++++++++++++++++++
										//   just fnn-pso cases

										gnuplotScriptFile.println("# JUST fnn-pso CASES");
										gnuplotScriptFile.println();
										gnuplotScriptFile.print("#plot [][] ");

										firstFile = true;
										
										// print out fnn-pso cases
//										for (int numPartsIndex = 0 ; numPartsIndex < numParticlesList.length ; ++numPartsIndex) {
//											if (!doNumParticles[numPartsIndex])
//												continue;

										for (int topoIndex = 0 ; topoIndex < topologiesList.length ; ++topoIndex) {
											if (!doTopology[topoIndex] || topoIndex < 4) 
												continue;

											for (int selfIndex = 0 ; selfIndex < selfModelsList.length ; ++selfIndex) {
												if (!doSelfModel[selfIndex])
													continue;

												for (int influenceIndex = 0 ; influenceIndex < influenceModelsList.length ; ++influenceIndex) {
													if (!doInfluenceModel[influenceIndex])
														continue;


													String thisFileName =  finalFuncValsByNumParticlesOutputFileNames[currentFunctionNum]
//													                                       [numPartsIndex]
													                                        [topoIndex]
													                                         [selfIndex]
													                                          [influenceIndex]
													                                           [latticeIndex]
													                                            [gainIndex]
													                                             [spontLevelIndex]
													                                              [spontProbIndex];
													if (thisFileName != null) {
														if (firstFile) {
															gnuplotScriptFile.print("\"" + thisFileName + "\" using 1:2 with linesp");
															firstFile = false;
														}
														else {
															gnuplotScriptFile.print(", \"" + thisFileName + "\" using 1:2 with linesp");															
														}
													}
												}
											}
										}
//										}
										// close the "plot"
										gnuplotScriptFile.print("\n");
										

										gnuplotScriptFile.print("\n\n\n\n");
										gnuplotScriptFile.println("####################################################################################");
										gnuplotScriptFile.print("\n\n\n\n");

										//   just fnn-pso cases
										// ++++++++++++++++++++++++++++++++++++++++++++


										
										
										
									}
								}
							}
						}
					}
//				}

				gnuplotScriptFile.close();
				
			}
			
			
			// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			//	BY GAINS
			
			if (gnuplotScriptsByGainsNeeded) {

				PrintWriter gnuplotScriptFile = new PrintWriter(new FileWriter(folder + "0-gnuplot-script-gains-" + dateString + ".txt"));
				gnuplotScriptFile.println("#  " + dateString + "\n\n");

				gnuplotScriptFile.println("set logscale y \n\n");

				for (int currentFunctionNum = 0 ; currentFunctionNum < doFunction.length ; ++currentFunctionNum) {
					if (!doFunction[currentFunctionNum])
						continue;

					for (int numPartsIndex = 0 ; numPartsIndex < numParticlesList.length ; ++numPartsIndex) {
						if (!doNumParticles[numPartsIndex])
							continue;

						for (int latticeIndex = 0 ; latticeIndex < latticeSizesList.length ; ++latticeIndex) {
							if (!doLatticeSize[latticeIndex])
								continue;

							for (int spontLevelIndex = 0 ; spontLevelIndex < spontActLevelsList.length ; ++spontLevelIndex) {
								if (!doSpontActLevel[spontLevelIndex])
									continue;

								for (int spontProbIndex = 0 ; spontProbIndex < spontActProbsList.length ; ++spontProbIndex) {
									if (!doSpontActProb[spontProbIndex])
										continue;
												
									for (int topoIndex = 0 ; topoIndex < topologiesList.length ; ++topoIndex) {

										if (!doTopology[topoIndex] || topoIndex < 4) 
											continue;
											
										for (int selfIndex = 0 ; selfIndex < selfModelsList.length ; ++selfIndex) {
											if (!doSelfModel[selfIndex])
												continue;

											for (int influenceIndex = 0 ; influenceIndex < influenceModelsList.length ; ++influenceIndex) {
												if (!doInfluenceModel[influenceIndex])
													continue;
												
												String topo = 
													getPSOTopologyAndInfluence(topologiesList[topoIndex], selfModelsList[selfIndex], influenceModelsList[influenceIndex], 
															currentPSOBoundaryModel, currentPSOActivityModel);
				
												gnuplotScriptFile.println("# " + topo);
												gnuplotScriptFile.println();
												gnuplotScriptFile.println("# " + TestFunctions.getFunctionName(currentFunctionNum));
												gnuplotScriptFile.println("# num particles = " + numParticlesList[numPartsIndex]);
												gnuplotScriptFile.println("# lattice size = " + latticeSizesList[latticeIndex]);
												gnuplotScriptFile.println("# spont act level = " + spontActLevelsList[spontLevelIndex]);
												gnuplotScriptFile.println("# spont prob level = " + spontActProbsList[spontProbIndex] );
												gnuplotScriptFile.println();
												gnuplotScriptFile.print("#plot [][] ");

												boolean firstFile = true;

												for (int gainIndex = 0 ; gainIndex < gainsList.length ; ++gainIndex) {

													if (!doGain[gainIndex])
														continue;


													String thisFileName =  iterationOutputFileNames[currentFunctionNum]
													                                                [numPartsIndex]
													                                                 [topoIndex]
													                                                  [selfIndex]
													                                                   [influenceIndex]
													                                                    [latticeIndex]
													                                                     [gainIndex]
													                                                      [spontLevelIndex]
													                                                       [spontProbIndex];
													if (thisFileName != null) {
														if (firstFile) {
															gnuplotScriptFile.print("\"" + thisFileName + "\" using 1:2 with linesp");
															firstFile = false;
														}
														else {
															gnuplotScriptFile.print(", \"" + thisFileName + "\" using 1:2 with linesp");															
														}
													}
												}
												
												// close the "plot"
												gnuplotScriptFile.print("\n\n\n\n");
												
											}
										}

									}


									gnuplotScriptFile.print("\n\n\n\n");
									gnuplotScriptFile.println("####################################################################################");
									gnuplotScriptFile.print("\n\n\n\n");

								}
							}
						}
					}
				}

				gnuplotScriptFile.close();
			}

			
			//	BY GAINS
			// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			
			
			
			
			// outputWindow is System.out
			outputWindow.println("DONE");	
			outputWindow.close();

			
		}  // try

		catch (Exception e)
		{
		    System.err.println("Caught exception: " + e.getMessage());
		}


	}


	public static boolean isStaticTopology (Topology currentTopology) {  

		return currentTopology == Topology.GBEST ||
		currentTopology == Topology.RING ||
		currentTopology == Topology.vonNEUMANN ||
		currentTopology == Topology.MOORE;

	}


	public static boolean isFNNTopology (Topology currentTopology) {  

		return currentTopology == Topology.FNN_GBEST ||
		currentTopology == Topology.FNN_RING ||
		currentTopology == Topology.FNN_vonNEUMANN ||
		currentTopology == Topology.FNN_MOORE;

	}


	public static String getTopology(int topologyIndex) {
		
		if (topologyIndex == 0){
			return "GBEST";
		}
		else if (topologyIndex == 1){
			return "RING";
		}
		else if (topologyIndex == 2){
			return "VON NEUMANN";
		}
		else if (topologyIndex == 3){
			return "MOORE";
		}
		else if (topologyIndex == 4){
			return "FNN GBEST";
		}
		else if (topologyIndex == 5){
			return "FNN RING";
		}
		else if (topologyIndex == 6){
			return "FNN VON NEUMANN";
		}
		else if (topologyIndex == 7){
			return "FNN MOORE";
		}

		return "UNKNOWN TOPOLOGY";
	}
	
	
	public static String getPSOTopologyAndInfluence(Topology currentPSOTopology, SelfModel currentPSOSelfModel, InfluenceModel currentPSOInfluenceModel,
			BoundaryModel currentPSOBoundaryModel, ActivityModel currentPSOActivityModel) {

		String topologyAndInfluence = "";

		if (currentPSOTopology == Topology.GBEST) {
			topologyAndInfluence += "GB";
		}
		else if (currentPSOTopology == Topology.RING){
			topologyAndInfluence += "RI";
		}
		else if (currentPSOTopology == Topology.vonNEUMANN) {
			topologyAndInfluence += "vN";
		}
		else if (currentPSOTopology == Topology.MOORE){
			topologyAndInfluence += "MO";
		}
		else if (currentPSOTopology == Topology.FNN_GBEST) {
			topologyAndInfluence += "FNN_GB";
		}
		else if (currentPSOTopology == Topology.FNN_RING){
			topologyAndInfluence += "FNN_RI";
		}
		else if (currentPSOTopology == Topology.FNN_vonNEUMANN) {
			topologyAndInfluence += "FNN_vN";
		}
		else if (currentPSOTopology == Topology.FNN_MOORE){
			topologyAndInfluence += "FNN_MO";
		}
		else {
			topologyAndInfluence += "UNKNOWN_PSO_TOPOLOGY";
		}
		

		if (currentPSOSelfModel == SelfModel.INCLUDE_SELF) {
			topologyAndInfluence += "-yesSELF";
		}
		else if (currentPSOSelfModel == SelfModel.NOT_INCLUDE_SELF) {
			topologyAndInfluence += "-noSELF";
		}
		else {
			topologyAndInfluence += "UNKNOWN_PSO_SELF_MODEL";
		}

		
		if (currentPSOInfluenceModel == InfluenceModel.NEIGH_BEST) {
			topologyAndInfluence += "-nBEST";
		}
		else if (currentPSOInfluenceModel == InfluenceModel.FIPS){
			topologyAndInfluence += "-FIPS";
		}
		else {
			topologyAndInfluence += "UNKNOWN_PSO_INFLUENCE_MODEL";
		}


		if (isFNNTopology(currentPSOTopology)) {
			if (currentPSOBoundaryModel == BoundaryModel.LATTICE) {
				topologyAndInfluence += "-LAT";
			}
			else if (currentPSOBoundaryModel == BoundaryModel.TORUS) {
				topologyAndInfluence += "-TOR";
			}
			else {
				topologyAndInfluence += "UNKNOWN_PSO_BOUNDARY_MODEL";
			}
			

			if (currentPSOActivityModel == ActivityModel.ALL_NEURONS) {
				topologyAndInfluence += "-ALLN";
			}
			else if (currentPSOActivityModel == ActivityModel.ONLY_ACTIVE_NEURONS) {
				topologyAndInfluence += "-ACTN";
			}
			else if (currentPSOActivityModel == ActivityModel.SIMILAR_ACTIVITY_STATUS_NEURONS) {
				topologyAndInfluence += "-SIMACTN";
			}
			else {
				topologyAndInfluence += "UNKNOWN_PSO_ACTIVITY_MODEL";
			}
			
		}
		
		return topologyAndInfluence;

	}

	
	
	public static String getFNNParameters() {

		String topologyAndInfluence = "";

		
		if (currentFNNTopology == Topology.FNN_vonNEUMANN) {
			topologyAndInfluence += "vN";
		}
		else if (currentFNNTopology == Topology.FNN_MOORE){
			topologyAndInfluence += "MO";
		}
		else {
			topologyAndInfluence += "UNKNOWN_FNN_TOPOLOGY";
		}
		
		
		if (currentFNNSelfModel == SelfModel.INCLUDE_SELF) {
			topologyAndInfluence += "-yesSELF";
		}
		else if (currentFNNSelfModel == SelfModel.NOT_INCLUDE_SELF) {
			topologyAndInfluence += "-noSELF";
		}
		else {
			topologyAndInfluence += "UNKNOWN_FNN_SELF_MODEL";
		}

		
		if (currentFNNBoundaryModel == BoundaryModel.LATTICE) {
			topologyAndInfluence += "-LAT";
		}
		else if (currentFNNBoundaryModel == BoundaryModel.TORUS) {
			topologyAndInfluence += "-TOR";
		}
		else {
			topologyAndInfluence += "UNKNOWN_FNN_BOUNDARY_MODEL";
		}

		
		if (currentFNNActivityModel == ActivityModel.ALL_NEURONS) {
			topologyAndInfluence += "-ALLN";
		}
		else if (currentFNNActivityModel == ActivityModel.ONLY_ACTIVE_NEURONS) {
			topologyAndInfluence += "-ACTN";
		}
		else {
			topologyAndInfluence += "UNKNOWN_FNN_ACTIVITY_MODEL";
		}
		
		return topologyAndInfluence;

	}

	

	public static String getFNNDataString(int latticeSize, double gain, double spontActLevel, double spontActProb) {

		String FNNdata = "";
		
		FNNdata += "-lat" + latticeSize;
		FNNdata += "-g" + gain;
		FNNdata += "-sal" + spontActLevel;
		FNNdata += "-sap" + spontActProb;
		
		return FNNdata;
	}
	
	
	public static void outputTestParameters(PrintWriter resultsOutputFile, Topology currentPSOTopology, SelfModel currentPSOSelfModel, InfluenceModel currentPSOInfluenceModel, 
			BoundaryModel currentPSOBoundaryModel, ActivityModel currentPSOActivityModel,
			int currentFunctionNum, int numDimensions, int numParticles, int numRuns, int totalNumFunctionEvaluations, int numIterations,
			int latticeSideSizeFNN, double gain, double sumNeighborActivationsThreshold, double activationThreshold, 
			double spontaneousActivationLevel, double spontaneousActivationProbability) {


		resultsOutputFile.println("# " + TestFunctions.getFunctionName(currentFunctionNum) + ", location of optimum shifted randomly for each run ");

		resultsOutputFile.println("# "+ getPSOTopologyAndInfluence(currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel, currentPSOBoundaryModel, currentPSOActivityModel));
		
		if (isFNNTopology(currentPSOTopology)) {

			resultsOutputFile.println("# FNN DYNAMICS: " + getFNNParameters());
			
			resultsOutputFile.print("# lattice = "+ latticeSideSizeFNN + "x" + latticeSideSizeFNN + ", ");
			resultsOutputFile.printf("implied density = %4.2f, ", (double) numParticles / (latticeSideSizeFNN * latticeSideSizeFNN));
			resultsOutputFile.print("gain = "+ gain + ", ");
			resultsOutputFile.print("sum neigh acts thresh = "+ sumNeighborActivationsThreshold + ", ");
			resultsOutputFile.print("act thresh = "+ activationThreshold + ", ");
			resultsOutputFile.print("spont act level = "+ spontaneousActivationLevel + ", ");
			resultsOutputFile.println("spont act prob = "+ spontaneousActivationProbability);

			resultsOutputFile.println("# initial activation range = [ " + FluidNN.INITIAL_ACTIVATION_LOW_LEVEL + ", " + FluidNN.INITIAL_ACTIVATION_HIGH_LEVEL + " ]");	

			if (useNewSpontActMechanism) {
				resultsOutputFile.println("# new spontaneous activation mechanism");
			}
			else {
				resultsOutputFile.println("# old spontaneous activation mechanism");
			}

		}

		
		resultsOutputFile.println("# " + numDimensions + " dimensions");

		resultsOutputFile.println("# " + numParticles + " particles");

		resultsOutputFile.print("# ");
		if (useIterations)
			resultsOutputFile.println(numIterations + " iterations");
		else
			resultsOutputFile.println(totalNumFunctionEvaluations + " function evaluations");

		resultsOutputFile.println("# neighborhood theta = " + neighborhoodTheta + ", personal theta = " + personalTheta +
				", overall theta = " + theta + ", constriction = " + constrictionFactor);			

		resultsOutputFile.print("# ");
		if (useSeedForRand)
			resultsOutputFile.println("random number generator seed = " + seed);			
		else
			resultsOutputFile.println("random number generator seed = no seed");

		resultsOutputFile.println("# " + numRuns + " runs");

		resultsOutputFile.println("# ");

		resultsOutputFile.flush();

	}




	// print out the function value summary data for the last iteration
	private static void outputSummaryData(PrintWriter summaryDataFile, IntervalSummaryData intervalSummaryData[], double secondsPerRun) {

		summaryDataFile.printf("%s %.4e \n", "Average function value over runs = ", intervalSummaryData[intervalSummaryData.length-1].getAverageFunctionValue());
		summaryDataFile.printf("%s %.4e \n", "Standard deviation over runs     = ", intervalSummaryData[intervalSummaryData.length-1].getStdDevFunctionValue());
		summaryDataFile.printf("%s %.4e \n", "Minimum function value over runs = ", intervalSummaryData[intervalSummaryData.length-1].getMinimumFunctionValue());
		summaryDataFile.printf("%s %.4e \n", "Median function value over runs  = ", intervalSummaryData[intervalSummaryData.length-1].getMedianFunctionValue());
		summaryDataFile.printf("%s %.4e \n", "Maximum function value over runs = ", intervalSummaryData[intervalSummaryData.length-1].getMaximumFunctionValue());
		//	summaryDataFile.printf("%s %4.3e \n", "RMSE over runs                   = ", iterationSummaryData[iterationSummaryData.length-1].getRootMeanSqrErrFunctionValue());
		summaryDataFile.println();		
		summaryDataFile.println("Time per run: " + secondsPerRun + " seconds");		
		summaryDataFile.println("--------------------------------------------------------------------------------------------------------------------------------------------------------------------");
		summaryDataFile.println();

		summaryDataFile.flush();

	}

	
	private static void outputFinalFuncValue(PrintWriter finalFuncValDataFile, IntervalSummaryData intervalSummaryData[], int numParticles) {

		finalFuncValDataFile.printf("%d  %.4e \n", numParticles, intervalSummaryData[intervalSummaryData.length-1].getAverageFunctionValue());
//		summaryDataFile.printf("%s %.4e \n", "Average function value over runs = ", intervalSummaryData[intervalSummaryData.length-1].getAverageFunctionValue());
//		summaryDataFile.printf("%s %.4e \n", "Standard deviation over runs     = ", intervalSummaryData[intervalSummaryData.length-1].getStdDevFunctionValue());
//		summaryDataFile.printf("%s %.4e \n", "Minimum function value over runs = ", intervalSummaryData[intervalSummaryData.length-1].getMinimumFunctionValue());
//		summaryDataFile.printf("%s %.4e \n", "Median function value over runs  = ", intervalSummaryData[intervalSummaryData.length-1].getMedianFunctionValue());
//		summaryDataFile.printf("%s %.4e \n", "Maximum function value over runs = ", intervalSummaryData[intervalSummaryData.length-1].getMaximumFunctionValue());
		//	summaryDataFile.printf("%s %4.3e \n", "RMSE over runs                   = ", iterationSummaryData[iterationSummaryData.length-1].getRootMeanSqrErrFunctionValue());
//		summaryDataFile.println();		
//		summaryDataFile.println("Time per run: " + secondsPerRun + " seconds");		
//		summaryDataFile.println("--------------------------------------------------------------------------------------------------------------------------------------------------------------------");
//		summaryDataFile.println();
//
//		summaryDataFile.flush();

	}



	// print out the function and error summary data for each iteration at which it was collected
	private static void outputIterationData(PrintWriter intervalSummaryDataFile, IntervalSummaryData intervalSummaryData[]) {


		intervalSummaryDataFile.println("#  val = function value,   err = absolute value error");

		if (useIterations)
			intervalSummaryDataFile.print("# iter      ");
		else
			intervalSummaryDataFile.print("# #FEs      ");

		intervalSummaryDataFile.print("mean-val    std-dev-val     min-val        q1-val      median-val      q3-val       max-val");
		intervalSummaryDataFile.println("       mean-err      min-err        q1-err      median-err      q3-err       max-err");


		intervalSummaryDataFile.println("#-------------------------------------------------------------------------------------------" + 
		"-------------------------------------------------------------------------------------------------");

		// start at interval i = 1 because we have set it up so that interval numbers are *not* zero-based
		// NOTE:  the length of the intervalSummaryData array is (numIntervalsForDataOutput+1), so
		// the last intervalNum, i.e. intervalSummaryData.length - 1, is actually numIntervalsForDataOutput
		for (int i = 1 ; i < intervalSummaryData.length ; ++i) {	

			if (useIterations) 
				intervalSummaryDataFile.printf("%6d     ", i * numItersPerOutputInterval);
			else
				intervalSummaryDataFile.printf("%6d     ", i * numFEsPerOutputInterval);

			// mean function value
			intervalSummaryDataFile.printf("%.4e    ", 
					intervalSummaryData[i].getAverageFunctionValue());
			// std dev function value
			intervalSummaryDataFile.printf("%.4e    ", 
					intervalSummaryData[i].getStdDevFunctionValue());
			// min function value
			intervalSummaryDataFile.printf("%.4e    ", 
					intervalSummaryData[i].getMinimumFunctionValue());
			// q1 function value
			intervalSummaryDataFile.printf("%.4e    ", 
					intervalSummaryData[i].getFirstQuartileFunctionValue());
			// q2 function value
			intervalSummaryDataFile.printf("%.4e    ", 
					intervalSummaryData[i].getMedianFunctionValue());
			// q3 function value
			intervalSummaryDataFile.printf("%.4e    ", 
					intervalSummaryData[i].getThirdQuartileFunctionValue());
			// max function value
			intervalSummaryDataFile.printf("%.4e    ", 
					intervalSummaryData[i].getMaximumFunctionValue());


			// mean abs val error
			intervalSummaryDataFile.printf("%.4e    ", 
					intervalSummaryData[i].getAverageAbsValError());
			// min abs val error
			intervalSummaryDataFile.printf("%.4e    ", 
					intervalSummaryData[i].getMinimumAbsValError());
			// q1 abs val error
			intervalSummaryDataFile.printf("%.4e    ", 
					intervalSummaryData[i].getFirstQuartileAbsValError());
			// q2 abs val error
			intervalSummaryDataFile.printf("%.4e    ", 
					intervalSummaryData[i].getMedianAbsValError());
			// q3 abs val error
			intervalSummaryDataFile.printf("%.4e    ", 
					intervalSummaryData[i].getThirdQuartileAbsValError());
			// max abs val error
			intervalSummaryDataFile.printf("%.4e    ", 
					intervalSummaryData[i].getMaximumAbsValError());


			intervalSummaryDataFile.println(); 

		}

		intervalSummaryDataFile.println("#-------------------------------------------------------------------------------------------" + 
		"-------------------------------------------------------------------------------------------------");

		if (useIterations)
			intervalSummaryDataFile.print("# iter      ");
		else
			intervalSummaryDataFile.print("# #FEs      ");

		intervalSummaryDataFile.print("mean-val    std-dev-val     min-val        q1-val      median-val      q3-val       max-val");
		intervalSummaryDataFile.println("       mean-err      min-err        q1-err      median-err      q3-err       max-err");

		intervalSummaryDataFile.flush();

	}

	
	



	// this has to be done before any summary statistics are printed to the screen or sent to a file
	// iterationData contains data for every run for every iteration
	// 
	private static void calculateSummaryDataOverRuns(DataOutput[][] intervalData, IntervalSummaryData[] intervalSummaryData, 
			int functionNum, int numDimensions) {

		for (int i = 0 ; i < intervalSummaryData.length ; ++i) {
			// we need to send iterationData[i], the data over all runs for iteration i, to the 
			// summarizeData method of intervalSummaryData[i], which is the IntervalSummaryData object 
			// for iteration i.  the method uses that data to calculate statistics (average, min, max, etc.) for that 
			// iteration over all the runs and set those values in the IntervalSummaryData object
			intervalSummaryData[i].summarizeData(intervalData[i], functionNum, numDimensions);
		}

	}




}





//for (int currentFunctionNum = 0 ; currentFunctionNum < doFunction.length ; ++currentFunctionNum) {
//	if (!doFunction[currentFunctionNum])
//		continue;
//
//	for (int numPartsIndex = 0 ; numPartsIndex < numParticlesList.length ; ++numPartsIndex) {
//		if (!doNumParticles[numPartsIndex])
//			continue;
//
//		
//		
//		
//		for (int latticeIndex = 0 ; latticeIndex < latticeSizesList.length ; ++latticeIndex) {
////			if (!doLatticeSize[latticeIndex])
////				continue;
//
//			for (int gainIndex = 0 ; gainIndex < gainsList.length ; ++gainIndex) {
////				if (!doGain[gainIndex])
////					continue;
//
//				for (int spontLevelIndex = 0 ; spontLevelIndex < spontActLevelsList.length ; ++spontLevelIndex) {
////					if (!doSpontActLevel[spontLevelIndex])
////						continue;
//
//					for (int spontProbIndex = 0 ; spontProbIndex < spontActProbsList.length ; ++spontProbIndex) {
////						if (!doSpontActProb[spontProbIndex])
////							continue;
//						
//						boolean firstTime = true;
//						
//						for (int topoIndex = 0 ; topoIndex < topologiesList.length ; ++topoIndex) {
//							for (int selfIndex = 0 ; selfIndex < selfModelsList.length ; ++selfIndex) {
//								for (int influenceIndex = 0 ; influenceIndex < influenceModelsList.length ; ++influenceIndex) {
//
//									if ((!doLatticeSize[latticeIndex] || !doGain[gainIndex] || !doSpontActLevel[spontLevelIndex] || !doSpontActProb[spontProbIndex] ||
//											topoIndex < 4 ||  !doToplogy[topoIndex])
//											&& 
//											(latticeIndex != 0 || gainIndex != 0 || spontLevelIndex != 0 || spontProbIndex != 0 || 
//													topoIndex >= 4 || 
//													!doToplogy[topoIndex])) {
//										continue;
//									}
//									
//									if (firstTime) {
//										gnuplotScriptFile.println("# " + TestFunctions.getFunctionName(currentFunctionNum));
//										gnuplotScriptFile.println("# num particles = " + numParticlesList[numPartsIndex]);
//										gnuplotScriptFile.println("# lattice size = " + latticeSizesList[latticeIndex]);
//										gnuplotScriptFile.println("# gain = " + gainsList[gainIndex]);
//										gnuplotScriptFile.println("# spont act level = " + spontActLevelsList[spontLevelIndex]);
//										gnuplotScriptFile.println("# spont prob level = " + spontActProbsList[spontProbIndex] );
//										gnuplotScriptFile.print("#plot [][] ");
//										firstTime = false;
//									}
//
//									String thisFileName =  outputFileNames[currentFunctionNum]
//									                                       [numPartsIndex]
//									                                        [topoIndex]
//									                                         [selfIndex]
//									                                          [influenceIndex]
//									                                           [latticeIndex]
//									                                            [gainIndex]
//									                                             [spontLevelIndex]
//									                                              [spontProbIndex];
//									if (thisFileName != null) {
//										gnuplotScriptFile.print("\"" + thisFileName + "\" using 1:2 with linesp, ");
//									}
//								}
//							}
//						}
//						if (!firstTime)
//							gnuplotScriptFile.print("\n\n\n");
//
//					}
//				}
//			}
//		}
//	}
//}




//
//for (int currentFunctionNum = 0 ; currentFunctionNum < doFunction.length ; ++currentFunctionNum) {
//	if (!doFunction[currentFunctionNum])
//		continue;
//
//	for (int numPartsIndex = 0 ; numPartsIndex < numParticlesList.length ; ++numPartsIndex) {
//		if (!doNumParticles[numPartsIndex])
//			continue;
//
//		for (int latticeIndex = 0 ; latticeIndex < latticeSizesList.length ; ++latticeIndex) {
//			if (!doLatticeSize[latticeIndex])
//				continue;
//
//			for (int gainIndex = 0 ; gainIndex < gainsList.length ; ++gainIndex) {
//				if (!doGain[gainIndex])
//					continue;
//
//				for (int spontLevelIndex = 0 ; spontLevelIndex < spontActLevelsList.length ; ++spontLevelIndex) {
//					if (!doSpontActLevel[spontLevelIndex])
//						continue;
//
//					for (int spontProbIndex = 0 ; spontProbIndex < spontActProbsList.length ; ++spontProbIndex) {
//						if (!doSpontActProb[spontProbIndex])
//							continue;
//
//						gnuplotScriptFile.println("# " + TestFunctions.getFunctionName(currentFunctionNum));
//						gnuplotScriptFile.println("# num particles = " + numParticlesList[numPartsIndex]);
//						gnuplotScriptFile.println("# lattice size = " + latticeSizesList[latticeIndex]);
//						gnuplotScriptFile.println("# gain = " + gainsList[gainIndex]);
//						gnuplotScriptFile.println("# spont act level = " + spontActLevelsList[spontLevelIndex]);
//						gnuplotScriptFile.println("# spont prob level = " + spontActProbsList[spontProbIndex] );
//						gnuplotScriptFile.print("#plot [][] ");
//
//						
//						for (int topoIndex = 0 ; topoIndex < topologiesList.length ; ++topoIndex) {
//							for (int selfIndex = 0 ; selfIndex < selfModelsList.length ; ++selfIndex) {
//								for (int influenceIndex = 0 ; influenceIndex < influenceModelsList.length ; ++influenceIndex) {
//
//									String thisFileName =  outputFileNames[currentFunctionNum]
//									                                       [numPartsIndex]
//									                                        [topoIndex]
//									                                         [selfIndex]
//									                                          [influenceIndex]
//									                                           [latticeIndex]
//									                                            [gainIndex]
//									                                             [spontLevelIndex]
//									                                              [spontProbIndex];
//									if (thisFileName != null) {
//										gnuplotScriptFile.print("\"" + thisFileName + "\" using 1:2 with linesp, ");
//									}
//								}
//							}
//						}
//						gnuplotScriptFile.print("\n\n\n");
//
//					}
//				}
//			}
//		}
//	}
//}
//
//
//





	
//	
//	
//	
//}
//gnuplotScriptFile.println("\"" + outputFileNames[0][numOutputFilesFunc0 - 1] + "\" using 1:2 with linesp");
//gnuplotScriptFile.print("\n\n\n");
//
//
//// make a plot for each function
//if (numOutputFilesFunc0 > 0) {
//	gnuplotScriptFile.print("#plot [][] ");
//	for (int caseNum = 0 ; caseNum < numOutputFilesFunc0 - 1 ; ++caseNum) {
//		gnuplotScriptFile.print("\"" + outputFileNames[0][caseNum] + "\" using 1:2 with linesp, ");
//	}
//	gnuplotScriptFile.println("\"" + outputFileNames[0][numOutputFilesFunc0 - 1] + "\" using 1:2 with linesp");
//	gnuplotScriptFile.print("\n\n\n");
//}
//if (numOutputFilesFunc1 > 0) {
//	gnuplotScriptFile.print("#plot [][] ");
//	for (int caseNum = 0 ; caseNum < numOutputFilesFunc1 - 1 ; ++caseNum) {
//		gnuplotScriptFile.print("\"" + outputFileNames[1][caseNum] + "\" using 1:2 with linesp, ");
//	}
//	gnuplotScriptFile.println("\"" + outputFileNames[1][numOutputFilesFunc1 - 1] + "\" using 1:2 with linesp");
//	gnuplotScriptFile.print("\n\n\n");
//}
//if (numOutputFilesFunc2 > 0) {
//	gnuplotScriptFile.print("#plot [][] ");
//	for (int caseNum = 0 ; caseNum < numOutputFilesFunc2 - 1 ; ++caseNum) {
//		gnuplotScriptFile.print("\"" + outputFileNames[2][caseNum] + "\" using 1:2 with linesp, ");
//	}
//	gnuplotScriptFile.println("\"" + outputFileNames[2][numOutputFilesFunc2 - 1] + "\" using 1:2 with linesp");
//	gnuplotScriptFile.print("\n\n\n");
//}
//if (numOutputFilesFunc3 > 0) {
//	gnuplotScriptFile.print("#plot [][] ");
//	for (int caseNum = 0 ; caseNum < numOutputFilesFunc3 - 1 ; ++caseNum) {
//		gnuplotScriptFile.print("\"" + outputFileNames[3][caseNum] + "\" using 1:2 with linesp, ");
//	}
//	gnuplotScriptFile.println("\"" + outputFileNames[3][numOutputFilesFunc3 - 1] + "\" using 1:2 with linesp");
//	gnuplotScriptFile.print("\n\n\n");
//}
//if (numOutputFilesFunc4 > 0) {
//	gnuplotScriptFile.print("#plot [][] ");
//	for (int caseNum = 0 ; caseNum < numOutputFilesFunc4 - 1 ; ++caseNum) {
//		gnuplotScriptFile.print("\"" + outputFileNames[4][caseNum] + "\" using 1:2 with linesp, ");
//	}
//	gnuplotScriptFile.println("\"" + outputFileNames[4][numOutputFilesFunc4 - 1] + "\" using 1:2 with linesp");
//	gnuplotScriptFile.print("\n\n\n");
//}
//if (numOutputFilesFunc5 > 0) {
//	gnuplotScriptFile.print("#plot [][] ");
//	for (int caseNum = 0 ; caseNum < numOutputFilesFunc5 - 1 ; ++caseNum) {
//		gnuplotScriptFile.print("\"" + outputFileNames[5][caseNum] + "\" using 1:2 with linesp, ");
//	}
//	gnuplotScriptFile.println("\"" + outputFileNames[5][numOutputFilesFunc5 - 1] + "\" using 1:2 with linesp");
//	gnuplotScriptFile.print("\n\n\n");
//}
//if (numOutputFilesFunc6 > 0) {
//	gnuplotScriptFile.print("#plot [][] ");
//	for (int caseNum = 0 ; caseNum < numOutputFilesFunc6 - 1 ; ++caseNum) {
//		gnuplotScriptFile.print("\"" + outputFileNames[6][caseNum] + "\" using 1:2 with linesp, ");
//	}
//	gnuplotScriptFile.println("\"" + outputFileNames[6][numOutputFilesFunc6 - 1] + "\" using 1:2 with linesp");
//	gnuplotScriptFile.print("\n\n\n");
//}
//if (numOutputFilesFunc7 > 0) {
//	gnuplotScriptFile.print("#plot [][] ");
//	for (int caseNum = 0 ; caseNum < numOutputFilesFunc7 - 1 ; ++caseNum) {
//		gnuplotScriptFile.print("\"" + outputFileNames[7][caseNum] + "\" using 1:2 with linesp, ");
//	}
//	gnuplotScriptFile.println("\"" + outputFileNames[7][numOutputFilesFunc7 - 1] + "\" using 1:2 with linesp");
//	gnuplotScriptFile.print("\n\n\n");
//}


// don't bother making a plot that includes *all* cases
//gnuplotScriptFile.print("plot [][] ");
//for (int i = 0 ; i < numOutputFiles - 1 ; ++i) {
//	gnuplotScriptFile.print("\"" + outputFileNames[i] + "\" using 1:2 with lines, ");
//}
//gnuplotScriptFile.println("\"" + outputFileNames[numOutputFiles-1] + "\" using 1:2 with lines");

// don't bother making a plot for each cases
//gnuplotScriptFile.print("\n\n\n");
//for (int i = 0 ; i < numOutputFiles ; ++i) {
//	gnuplotScriptFile.println("#plot [][] \"" + outputFileNames[i] + "\" using 1:2 with lines\n");
//}

//System.out.println ("HERE 1");
//for (int numPartsIndex = 0 ; numPartsIndex < numParticlesList.length ; ++numPartsIndex) {
//	if (doNumParticles[numPartsIndex]) {
//		gnuplotScriptFile.println("# num particles = " + numParticlesList[numPartsIndex]);
//		gnuplotScriptFile.print("#plot [][] ");
//		for (int currentFunctionNum = 0 ; currentFunctionNum < doFunction.length ; ++currentFunctionNum) {
//			for (int topoIndex = 0 ; topoIndex < topologiesList.length ; ++topoIndex) {
//				for (int selfIndex = 0 ; selfIndex < selfModelsList.length ; ++selfIndex) {
//					for (int influenceIndex = 0 ; influenceIndex < influenceModelsList.length ; ++influenceIndex) {
//						for (int latticeIndex = 0 ; latticeIndex < latticeSizesList.length ; ++latticeIndex) {
//							for (int gainIndex = 0 ; gainIndex < gainsList.length ; ++gainIndex) {
//								for (int spontLevelIndex = 0 ; spontLevelIndex < spontActLevelsList.length ; ++spontLevelIndex) {
//									for (int spontProbIndex = 0 ; spontProbIndex < spontActProbsList.length ; ++spontProbIndex) {
//										String thisFileName =  outputFileNames[currentFunctionNum]
//										                                       [numPartsIndex]
//										                                        [topoIndex]
//										                                         [selfIndex]
//										                                          [influenceIndex]
//										                                           [latticeIndex]
//										                                            [gainIndex]
//										                                             [spontLevelIndex]
//										                                              [spontProbIndex];
//										//													System.out.println ("BEFORE IF thisFileName = " + thisFileName);
//
//										if (thisFileName != null) {
//											//														System.out.println ("INSIDE IF thisFileName = " + thisFileName);
//											gnuplotScriptFile.print("\"" + thisFileName + "\" using 1:2 with linesp, ");
//										}
//									}
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//		gnuplotScriptFile.print("\n\n\n");
//	}
//}
//
//for (int currentFunctionNum = 0 ; currentFunctionNum < doFunction.length ; ++currentFunctionNum) {
//	if (doFunction[currentFunctionNum]) {
//		gnuplotScriptFile.println("# function = " + TestFunctions.getFunctionName(currentFunctionNum));
//		gnuplotScriptFile.print("#plot [][] ");
//		for (int numPartsIndex = 0 ; numPartsIndex < numParticlesList.length ; ++numPartsIndex) {
//			for (int topoIndex = 0 ; topoIndex < topologiesList.length ; ++topoIndex) {
//				for (int selfIndex = 0 ; selfIndex < selfModelsList.length ; ++selfIndex) {
//					for (int influenceIndex = 0 ; influenceIndex < influenceModelsList.length ; ++influenceIndex) {
//						for (int latticeIndex = 0 ; latticeIndex < latticeSizesList.length ; ++latticeIndex) {
//							for (int gainIndex = 0 ; gainIndex < gainsList.length ; ++gainIndex) {
//								for (int spontLevelIndex = 0 ; spontLevelIndex < spontActLevelsList.length ; ++spontLevelIndex) {
//									for (int spontProbIndex = 0 ; spontProbIndex < spontActProbsList.length ; ++spontProbIndex) {
//										String thisFileName =  outputFileNames[currentFunctionNum]
//										                                       [numPartsIndex]
//										                                        [topoIndex]
//										                                         [selfIndex]
//										                                          [influenceIndex]
//										                                           [latticeIndex]
//										                                            [gainIndex]
//										                                             [spontLevelIndex]
//										                                              [spontProbIndex];
//										//													System.out.println ("BEFORE IF thisFileName = " + thisFileName);
//
//										if (thisFileName != null) {
//											//														System.out.println ("INSIDE IF thisFileName = " + thisFileName);
//											gnuplotScriptFile.print("\"" + thisFileName + "\" using 1:2 with linesp, ");
//										}
//									}
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//		gnuplotScriptFile.print("\n\n\n");
//	}
//}
//
//
//for (int topoIndex = 0 ; topoIndex < topologiesList.length ; ++topoIndex) {
//	if (doToplogy[topoIndex]) {
//		gnuplotScriptFile.println("# topology = " + getTopology(topoIndex));
//		gnuplotScriptFile.print("#plot [][] ");
//		for (int currentFunctionNum = 0 ; currentFunctionNum < doFunction.length ; ++currentFunctionNum) {
//		for (int numPartsIndex = 0 ; numPartsIndex < numParticlesList.length ; ++numPartsIndex) {
//			for (int selfIndex = 0 ; selfIndex < selfModelsList.length ; ++selfIndex) {
//					for (int influenceIndex = 0 ; influenceIndex < influenceModelsList.length ; ++influenceIndex) {
//						for (int latticeIndex = 0 ; latticeIndex < latticeSizesList.length ; ++latticeIndex) {
//							for (int gainIndex = 0 ; gainIndex < gainsList.length ; ++gainIndex) {
//								for (int spontLevelIndex = 0 ; spontLevelIndex < spontActLevelsList.length ; ++spontLevelIndex) {
//									for (int spontProbIndex = 0 ; spontProbIndex < spontActProbsList.length ; ++spontProbIndex) {
//										String thisFileName =  outputFileNames[currentFunctionNum]
//										                                       [numPartsIndex]
//										                                        [topoIndex]
//										                                         [selfIndex]
//										                                          [influenceIndex]
//										                                           [latticeIndex]
//										                                            [gainIndex]
//										                                             [spontLevelIndex]
//										                                              [spontProbIndex];
//										//													System.out.println ("BEFORE IF thisFileName = " + thisFileName);
//
//										if (thisFileName != null) {
//											//														System.out.println ("INSIDE IF thisFileName = " + thisFileName);
//											gnuplotScriptFile.print("\"" + thisFileName + "\" using 1:2 with linesp, ");
//										}
//									}
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//		gnuplotScriptFile.print("\n\n\n");
//	}
//}
//
//for (int topoIndex = 0 ; topoIndex < topologiesList.length ; ++topoIndex) {
//	for (int currentFunctionNum = 0 ; currentFunctionNum < doFunction.length ; ++currentFunctionNum) {
//		for (int selfIndex = 0 ; selfIndex < selfModelsList.length ; ++selfIndex) {
//			for (int influenceIndex = 0 ; influenceIndex < influenceModelsList.length ; ++influenceIndex) {
//				for (int latticeIndex = 0 ; latticeIndex < latticeSizesList.length ; ++latticeIndex) {
//					for (int gainIndex = 0 ; gainIndex < gainsList.length ; ++gainIndex) {
//						for (int spontLevelIndex = 0 ; spontLevelIndex < spontActLevelsList.length ; ++spontLevelIndex) {
//							for (int spontProbIndex = 0 ; spontProbIndex < spontActProbsList.length ; ++spontProbIndex) {
//								for (int numPartsIndex = 0 ; numPartsIndex < numParticlesList.length ; ++numPartsIndex) {
//									gnuplotScriptFile.println("# by num particles");
//									gnuplotScriptFile.print("#plot [][] ");
//									if (doNumParticles[numPartsIndex]) {
//										String thisFileName =  outputFileNames[currentFunctionNum]
//										                                       [numPartsIndex]
//										                                        [topoIndex]
//										                                         [selfIndex]
//										                                          [influenceIndex]
//										                                           [latticeIndex]
//										                                            [gainIndex]
//										                                             [spontLevelIndex]
//										                                              [spontProbIndex];
//										//													System.out.println ("BEFORE IF thisFileName = " + thisFileName);
//
//										if (thisFileName != null) {
//											//														System.out.println ("INSIDE IF thisFileName = " + thisFileName);
//											gnuplotScriptFile.print("\"" + thisFileName + "\" using 1:2 with linesp, ");
//										}
//									}
//									gnuplotScriptFile.print("\n\n\n");
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//}



//private static void showSummaryData(RunIterationStatistics intervalSummaryDataOverRuns[], double secondsPerRun) {
//
//		System.out.printf("%s %.4e \n", "Average function value over runs = ", intervalSummaryDataOverRuns[intervalSummaryDataOverRuns.length-1].getAverageFunctionValue());
//		System.out.printf("%s %.4e \n", "Standard deviation over runs     = ", intervalSummaryDataOverRuns[intervalSummaryDataOverRuns.length-1].getStdDevFunctionValue());
//		System.out.printf("%s %.4e \n", "Minimum function value over runs = ", intervalSummaryDataOverRuns[intervalSummaryDataOverRuns.length-1].getMinimumFunctionValue());
//		System.out.printf("%s %.4e \n", "Median function value over runs  = ", intervalSummaryDataOverRuns[intervalSummaryDataOverRuns.length-1].getMedianFunctionValue());
//		System.out.printf("%s %.4e \n", "Maximum function value over runs = ", intervalSummaryDataOverRuns[intervalSummaryDataOverRuns.length-1].getMaximumFunctionValue());
//		//	System.out.printf("%s %4.3e \n", "RMSE over runs                   = ", iterationSummaryDataOverRuns[totalNumIterations].getRootMeanSqrErrFunctionValue());
//		System.out.println();		
//		System.out.println("Time per run: " + secondsPerRun + " seconds");		
//		System.out.println("--------------------------------------------------------------------------------------------------------------------------------------------------------------------");
//		System.out.println();
//
//
//	}




//public static void showTypeOfTest(Topology currentTopology, SelfModel currentSelfModel, InfluenceModel currentInfluenceModel, 
//int currentFunctionNum, int numDimensions, int numParticles, int numRuns, int totalNumFunctionEvaluations, int numIterations) {
//
//
//
//
//System.out.println();
//
//System.out.print(TestFunctions.getFunctionName(currentFunctionNum));
//System.out.print(", ");
//
//System.out.print("location of optimum shifted randomly for each run ");
//System.out.println();
//
//
//System.out.print("# ");
//
//System.out.println(getTopologyAndInfluence(currentTopology, currentSelfModel, currentInfluenceModel));
//
////if (currentTopology == Topology.GBEST) {
////System.out.print("GBEST");
////}
////else if (currentTopology == Topology.RING){
////System.out.print("RING");
////}
////else if (currentTopology == Topology.vonNEUMANN) {
////System.out.print("vonNEUMANN");
////}
////else if (currentTopology == Topology.MOORE){
////System.out.print("MOORE");
////}
////
////if (currentSelfModel == SelfModel.INCLUDE_SELF) {
////System.out.print(", including SELF");
////}
////else {
////System.out.print(", NOT including SELF");
////}
////
////if (currentInfluenceModel == InfluenceModel.NEIGH_BEST) {
////System.out.print(", using NEIGHBORHOOD BEST");
////}
////else if (currentInfluenceModel == InfluenceModel.FIPS){
////System.out.print(", using FIPS");
////}
////
////System.out.println();
//
//
//System.out.print("# ");
//
//System.out.print(numDimensions + " dimensions");
//
//System.out.println();
//System.out.print("# ");
//
//System.out.print(numParticles + " particles");
//
//System.out.println();
//System.out.print("# ");
//
//if (useIterations)
//System.out.print(numIterations + " iterations");
//else
//System.out.print(totalNumFunctionEvaluations + " function evaluations");
//
//System.out.println();
//System.out.print("# ");
//
//System.out.print("global theta = " + neighborhoodTheta + ", personal theta = " + personalTheta +
//	", overall theta = " + theta + ", constriction = " + constrictionFactor);			
//
//System.out.println();
//System.out.print("# ");
//
//if (useSeedForRand)
//System.out.print("random number generator seed = " + seed);			
//else
//System.out.print("random number generator seed = no seed");
//
//System.out.println();
//System.out.print("# ");
//
//System.out.println(numRuns + " runs");
//
//System.out.println();
//
//}
