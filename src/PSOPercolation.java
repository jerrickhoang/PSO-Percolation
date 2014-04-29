import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Random;

public class PSOPercolation {
	
	// ****************  MISCELLANEOUS	  ******************

	// for random numbers
	public static boolean useSeedForRand;
	public static int seed;   
	public static Random rand = new Random();

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

	public static enum Topology {
		GBEST, RING, vonNEUMANN, MOORE, FNN_GBEST, FNN_RING, FNN_vonNEUMANN, FNN_MOORE,
	}
	public static int[] numRowsVonNeumannAndMooreList;
	public static int[] numColsVonNeumannAndMooreList;

	public static enum SelfModel {
		INCLUDE_SELF, NOT_INCLUDE_SELF
	}

	public static enum InfluenceModel {
		NEIGH_BEST, FIPS
	}

	public static enum BoundaryModel {
		LATTICE, TORUS
	}

	public static enum ActivityModel {
		ALL_NEURONS, ONLY_ACTIVE_NEURONS, SIMILAR_ACTIVITY_STATUS_NEURONS
	}
	
	public static BoundaryModel currentPSOBoundaryModel; 
	public static ActivityModel currentPSOActivityModel; 
	public static double neighborhoodTheta;
	public static double personalTheta;
	public static double theta;
	public static double constrictionFactor;

	public static boolean summaryDataFileNeeded = true; 
	public static boolean intervalDataFileNeeded;
	public static boolean allFinalFuncValsDataFileNeeded = true;
	public static boolean finalFuncValsByNumParticlesDataFileNeeded;


	
	public static void main(String[] args) {

		finalFuncValsByNumParticlesDataFileNeeded = true;
		intervalDataFileNeeded = true;

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

		
		PrintWriter outputWindow = new PrintWriter(System.out);

		double[][][][][][][][][] finalFuncValues = null;
		if (finalFuncValsByNumParticlesDataFileNeeded) {
			finalFuncValues = new double[8][8][8][2][2][8][10][4][4];
		}
		// standard PSO parameters
		neighborhoodTheta = 2.05;		
		personalTheta = 2.05;
		theta = neighborhoodTheta + personalTheta;
		constrictionFactor = 2.0 / (theta - 2.0 + Math.sqrt(theta*theta - 4.0*theta));

		// how many runs?
		totalNumRuns = 50;  
		SimpleDateFormat dateformatter = new SimpleDateFormat("yyyy-MM-dd--hh-mm-ss-a");
		Calendar date = Calendar.getInstance();
		String dateString = dateformatter.format(date.getTime());

		System.out.println("RUNNING CODE ON " + dateString + "\n");
		String folder = "results/";
		PrintWriter allFinalFuncValsDataFile = createrWriterToAllFinalFuncValsDataFile(folder, dateString);
		PrintWriter summaryDataFile = createWriterToSummaryFile(folder, dateString);
		
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
		currentPSOBoundaryModel = BoundaryModel.TORUS; 
		currentPSOActivityModel = ActivityModel.ALL_NEURONS; 

		Topology currentPSOTopology = topologiesList[7];
		SelfModel currentPSOSelfModel =  selfModelsList[0]; 
		InfluenceModel currentPSOInfluenceModel = influenceModelsList[0]; 
		int currentNumDimensions = 10;
		int currentLatticeSideSize = 10;
		int currentNumParticles = 50;
		int numParticlesIndex = 5;
		double probability = 0.5;
		
		long startTimeAllRuns = System.currentTimeMillis();  
		
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
		int[] numIterationsList = { 1000, 30000, 50000, 100000, 200000 }; 
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

		//Probability 
		double[] probabilityList = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
		
		

		// # ITERATIONS
		// ============
		int numIterationsIndex;
		int numFEsIndex;
		for (int probabilityIndex = 0; probabilityIndex < probabilityList.length; probabilityIndex ++) {
			probability = probabilityList[probabilityIndex];
			for (int currentFunctionNum = 0 ; currentFunctionNum < doFunction.length ; ++currentFunctionNum) {
				if (!doFunction[currentFunctionNum])
					continue;
				for (int numDimensionsIndex = 0 ; numDimensionsIndex < doNumDims.length ; ++numDimensionsIndex ) {
					if (!doNumDims[numDimensionsIndex])
						continue;
					currentNumDimensions = numDimsList[numDimensionsIndex];
					for (numIterationsIndex = 0, numFEsIndex = 0 ; numIterationsIndex < doNumIterations.length && numFEsIndex < doNumFEs.length ; ++numIterationsIndex, ++numFEsIndex ) {
			
						if (useIterations && !doNumIterations[numIterationsIndex])
							continue;
			
						if (!useIterations && !doNumFEs[numFEsIndex])
							continue;
			
			
						currentFENum = 0;
						currentIterNum = 0;
			
						if (useIterations) {
							totalNumIters = numIterationsList[numIterationsIndex];
			
							totalNumFEs = (totalNumIters * currentNumParticles) + currentNumParticles;
			
							numFEsPerOutputInterval = currentNumParticles * numItersPerOutputInterval;
			
							numInitialFEsIgnored = currentNumParticles;
			
							numIntervalsForDataOutput = totalNumIters / numItersPerOutputInterval;
			
						}
			
						else {
							totalNumFEs = numFEsList[numFEsIndex];
			
							numInitialFEsIgnored = 0;
							numIntervalsForDataOutput = totalNumFEs / numFEsPerOutputInterval;
			
							if (totalNumFEs % numFEsPerOutputInterval != 0)
								++numIntervalsForDataOutput;
						}
					
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
						
						for (currentRunNum = 0 ; currentRunNum < totalNumRuns ; ++currentRunNum) {
				
							runExperiment(finalFuncValues, allFinalFuncValsDataFile,
									currentPSOTopology, currentPSOSelfModel,
									currentPSOInfluenceModel, currentFunctionNum,
									currentNumDimensions, currentLatticeSideSize,
									currentNumParticles, numParticlesIndex, probability,
									intervalData);
							
						}  // END RUNS FOR-LOOP 
						
				
						long endTimeAllRuns = System.currentTimeMillis();
						double secondsPerRun = ((endTimeAllRuns - startTimeAllRuns) / 1000.0) / totalNumRuns;
						
						writeToOutputFiles(outputWindow, allFinalFuncValsDataFile,
								summaryDataFile, currentPSOTopology, currentPSOSelfModel,
								currentPSOInfluenceModel, currentFunctionNum,
								currentNumDimensions, currentLatticeSideSize,
								currentNumParticles, intervalData, intervalSummaryData,
								secondsPerRun);
					}
				}
			}
		}
	}

	private static void runExperiment(double[][][][][][][][][] finalFuncValues,
			PrintWriter allFinalFuncValsDataFile, Topology currentPSOTopology,
			SelfModel currentPSOSelfModel,
			InfluenceModel currentPSOInfluenceModel, int currentFunctionNum,
			int currentNumDimensions, int currentLatticeSideSize,
			int currentNumParticles, int numParticlesIndex, double probability,
			DataOutput[][] intervalData) {
		double shiftVectorAmount = TestFunctions.SHIFT_RANGE[currentFunctionNum] * rand.nextDouble();
		if (rand.nextDouble() < 0.5) 
			shiftVectorAmount *= -1.0;
		TestFunctions.setShiftVector(currentFunctionNum, new DoubleVector(currentNumDimensions, shiftVectorAmount));
		currentFENum = 0;

		Swarm swarm = new Swarm (currentNumParticles, numParticlesIndex, currentFunctionNum, currentNumDimensions, 
				currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel,
				currentLatticeSideSize, intervalData, probability);

		for ( ; currentFENum < totalNumFEs ; ) {
			++currentIterNum;
			swarm.update(currentFunctionNum, currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel, intervalData);

		}  // END ITERATIONS FOR-LOOP 


		double finalFuncValue = swarm.globalBest.getFunctionValue();
		finalFuncValues[currentFunctionNum]
		                [numParticlesIndex]
		                 [0]
		                  [0]
		                   [0]
		                    [0]
		                     [0]
		                      [0]
		                       [0] = finalFuncValue;
		if (allFinalFuncValsDataFileNeeded) {
			allFinalFuncValsDataFile.println(finalFuncValue);
		}
	}

	private static void writeToOutputFiles(PrintWriter outputWindow,
			PrintWriter allFinalFuncValsDataFile, PrintWriter summaryDataFile,
			Topology currentPSOTopology, SelfModel currentPSOSelfModel,
			InfluenceModel currentPSOInfluenceModel, int currentFunctionNum,
			int currentNumDimensions, int currentLatticeSideSize,
			int currentNumParticles, DataOutput[][] intervalData,
			IntervalSummaryData[] intervalSummaryData, double secondsPerRun) {
		calculateSummaryDataOverRuns(intervalData, intervalSummaryData, currentFunctionNum, currentNumDimensions);

		outputTestParameters(outputWindow, currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel, 
				currentPSOBoundaryModel, currentPSOActivityModel, 
				currentFunctionNum, currentNumDimensions, currentNumParticles, totalNumRuns, totalNumFEs, totalNumIters,	
				currentLatticeSideSize);		

		outputSummaryData(outputWindow, intervalSummaryData, secondsPerRun);								

		if (summaryDataFileNeeded) {
			outputTestParameters(summaryDataFile, currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel, 
					currentPSOBoundaryModel, currentPSOActivityModel, 
					currentFunctionNum, currentNumDimensions, currentNumParticles, totalNumRuns, totalNumFEs, totalNumIters,	
					currentLatticeSideSize);		

			outputSummaryData(summaryDataFile, intervalSummaryData, secondsPerRun);								

		}

		
		if (allFinalFuncValsDataFileNeeded) {
			allFinalFuncValsDataFile.println();
			outputSummaryData(allFinalFuncValsDataFile, intervalSummaryData, secondsPerRun);								
		}

		
//			if (finalFuncValsByNumParticlesDataFileNeeded) {
//				outputFinalFuncValue(finalFuncValByNumParticlesDataFile, intervalSummaryData, currentNumParticles);
//			}
		
		
//			if (intervalDataFileNeeded) {
//				outputIterationData(intervalDataFile, intervalSummaryData);
//				intervalDataFile.close();
//			}
	}
	
	public static void outputTestParameters(PrintWriter resultsOutputFile, Topology currentPSOTopology, SelfModel currentPSOSelfModel, InfluenceModel currentPSOInfluenceModel, 
			BoundaryModel currentPSOBoundaryModel, ActivityModel currentPSOActivityModel,
			int currentFunctionNum, int numDimensions, int numParticles, int numRuns, int totalNumFunctionEvaluations, int numIterations,
			int latticeSideSizeFNN) {


		resultsOutputFile.println("# " + TestFunctions.getFunctionName(currentFunctionNum) + ", location of optimum shifted randomly for each run ");

		resultsOutputFile.println("# "+ getPSOTopologyAndInfluence(currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel, currentPSOBoundaryModel, currentPSOActivityModel));

		resultsOutputFile.print("# lattice = "+ latticeSideSizeFNN + "x" + latticeSideSizeFNN + ", ");
		resultsOutputFile.printf("implied density = %4.2f, ", (double) numParticles / (latticeSideSizeFNN * latticeSideSizeFNN));


		resultsOutputFile.println("# initial activation range = [ " + FluidNN.INITIAL_ACTIVATION_LOW_LEVEL + ", " + FluidNN.INITIAL_ACTIVATION_HIGH_LEVEL + " ]");	

		
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

	public static boolean isFNNTopology (Topology currentTopology) {  

		return currentTopology == Topology.FNN_GBEST ||
		currentTopology == Topology.FNN_RING ||
		currentTopology == Topology.FNN_vonNEUMANN ||
		currentTopology == Topology.FNN_MOORE;

	}
	
	public static PrintWriter createWriterToSummaryFile(String folder, String dateString) {
		PrintWriter summaryDataFile = null;
		if (summaryDataFileNeeded) {
			String fileDataDescription = folder + "SUMMARY-DATA-";
			try {
				summaryDataFile = new PrintWriter(new FileWriter(fileDataDescription + "-" + dateString));
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			summaryDataFile.println(dateString);
			summaryDataFile.println();
		}
		return summaryDataFile;
	}
	
	public static PrintWriter createrWriterToAllFinalFuncValsDataFile(String folder, String dateString) {
		PrintWriter allFinalFuncValsDataFile = null;
		if (allFinalFuncValsDataFileNeeded) {
			String fileDataDescription = folder + "FINAL-FUNC-VALS-";
			try {
				allFinalFuncValsDataFile = new PrintWriter(new FileWriter(fileDataDescription + "-" + dateString));
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			allFinalFuncValsDataFile.println(dateString);
			allFinalFuncValsDataFile.println();
		}
		return allFinalFuncValsDataFile;
	}

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

	private static void calculateSummaryDataOverRuns(DataOutput[][] intervalData, IntervalSummaryData[] intervalSummaryData, 
			int functionNum, int numDimensions) {

		for (int i = 0 ; i < intervalSummaryData.length ; ++i) {
			intervalSummaryData[i].summarizeData(intervalData[i], functionNum, numDimensions);
		}

	}

	private static void outputFinalFuncValue(PrintWriter finalFuncValDataFile, IntervalSummaryData intervalSummaryData[], int numParticles) {

		finalFuncValDataFile.printf("%d  %.4e \n", numParticles, intervalSummaryData[intervalSummaryData.length-1].getAverageFunctionValue());

	}

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

	
	

}
