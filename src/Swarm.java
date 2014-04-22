



public class Swarm {

	public static Particle[] particles;  
	public static FluidNN fnn;

	public static Solution globalBest;
	private static Neighborhood[] allNeighs;



	public Swarm (int numParticles, 
			int numParticlesIndex, 
			int functionNum, 
			int numDimensions, 	
			PSO.Topology currentPSOTopology, 
			PSO.SelfModel currentPSOSelfModel, 
			PSO.InfluenceModel currentPSOInfluenceModel,
			int latticeSideSizeFNN, 
			double gain, 		
			double sumNeighActThresh, 
			double activationThreshold, 	
			double spontActLevel, 
			double spontActProb,
			DataOutput[][] intervalRunData) {

		// create arrays to hold particle and neighborhoods
		particles = new Particle[numParticles];   
		allNeighs = new Neighborhood[numParticles];		

		// array needed to get back function evaluation results from Particle constructor
		double[] initParticleData = new double[2];    

		// create the particles

		particles[0] = new Particle(functionNum, numDimensions, 0, initParticleData);

		// first one is the current global best
		int globalBestParticleNum = 0;
		double globalBestValue = initParticleData[TestFunctions.VAL_INDEX];
		double globalBestError = initParticleData[TestFunctions.ERR_INDEX];

		// getPosition does not return a copy, but that's okay because Solution constructor makes a copy of the DoubleVector sent in
		// NOTE: false means it is not an approximated value
		globalBest = new Solution(particles[0].getPosition(), globalBestValue, globalBestError, 0, globalBestParticleNum, false);


		for (int particleID = 1 ; particleID < particles.length ; particleID++) {

			particles[particleID] = new Particle(functionNum, numDimensions, particleID, initParticleData);

			// evaluate
			double particleValue = initParticleData[TestFunctions.VAL_INDEX];
			double particleError = initParticleData[TestFunctions.ERR_INDEX];

			if (particleValue < globalBestValue) {  
				globalBest.copyFromPosition(particles[particleID].getPosition());
				globalBest.setFunctionValue(particleValue);
				globalBest.setError(particleError);
				globalBest.setIterationCreated(0);
				globalBest.setParticleID(particleID);
				globalBest.setApproximated(false);
			}

			// data is collected based on number of function evaluations; if we are counting iterations, things are set up in PSO.java
			// so that the number of FEs in each interval is equal to the number of iterations we want in each interval
			// NOTE: PSO.numInitialFEsIgnored is > 0 only if we are counting iterations, since the FEs when the particles are
			// created should not count towards the FEs used to count iterations
			if ((PSO.currentFENum - PSO.numInitialFEsIgnored) % PSO.numFEsPerOutputInterval == 0) {
				int intervalIndex = (PSO.currentFENum - PSO.numInitialFEsIgnored) / PSO.numFEsPerOutputInterval;				
				intervalRunData[intervalIndex][PSO.currentRunNum].copyDataFromSolution(globalBest);
				//				System.out.println("CONST: putting item at " + (PSO.currentFunctionEvaluationNum - PSO.numInitialFEsIgnored) + " FEs (not including FEs to create swarm) in index " + intervalRunDataIndex + " of intervalRunData");
			}

		}


		// if FNN topology, create the Fluid Neural Net and connect the swarm to it
		if (PSO.isFNNTopology(currentPSOTopology)) {

			fnn = new FluidNN(latticeSideSizeFNN, latticeSideSizeFNN, numParticles, gain, 
					sumNeighActThresh, activationThreshold, spontActLevel, spontActProb);

			connectSwarmToFNN();
			//		printParticlesAndNeurons();
			//		fnn.printNeuronsAndParticles();
			//		checkParticlesAndNeurons();
		}


		// create neighborhoods, but only for static topologies;
		// if topology is dynamic, we will need to reconstruct it every iteration             *********************************
		if (PSO.isStaticTopology(currentPSOTopology)) {

			// these two calls create, for each particle:
			//		1) its Neighborhood, and 
			//		2) a list of Neighborhoods that contain it 


			// create the neighborhoods
			// this is where the heavy lifting gets done, i.e. creating the actual list of particles for each neighborhood
			createNeighborhoods(currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel, numParticlesIndex);


			//		System.out.println();
			//		System.out.println();
			//		showAllNeighborhoods();
			//		System.out.println();
			//		System.out.println();

			// each particle has a list of the neighborhoods it is in
			createNeighLists();
		}

	}

	
	
	
	
	

	public void update (int currentFunctionNum, PSO.Topology currentPSOTopology, PSO.SelfModel currentPSOSelfModel, PSO.InfluenceModel currentPSOInfluenceModel,
			DataOutput[][] intervalData) {

		Neighborhood[] particleNeighs = null;
		if (PSO.isFNNTopology(currentPSOTopology)) {
			// updateNeurons returns the particle Neighborhoods that were calculated when the activation levels were updated;
			// send that to updateParticles so they can be used there instead of recalculating the neighborhoods
			particleNeighs = updateNeurons(PSO.currentFNNTopology, PSO.currentFNNSelfModel, PSO.currentFNNBoundaryModel, PSO.currentFNNActivityModel);
		}

		// although updateNeurons returns an array of the Neighborhoods calculated when neuron activations were updated, we no  
		// longer send it to updateParticles because the topology/selfModel/influenceModel of the PSO can be different from that 
		// of the FNN, so the neighborhoods will not be the same, i.e. the PSO will have to calculate its own neighborhoods
		// in the getNeighBestPositionFNN and getFIPSAccelerationFNN methods in the Neighborhood class
		// (code that does send array of the Neighborhoods calculated when neuron activations were updated
		//    to updateParticles is in the two commented-out method below:
		//   1) updateSendNeighborhoods
		//   2) updateParticlesSendNeighborhoods
		updateParticles(currentFunctionNum, currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel, intervalData);

	}
	
	
	
	public void updateParticles (int currentFunctionNum, PSO.Topology currentPSOTopology, PSO.SelfModel currentPSOSelfModel, PSO.InfluenceModel currentPSOInfluenceModel, 
			DataOutput[][] intervalData) {

		for (int particleID = 0 ; particleID < particles.length ; particleID++) {

			Solution newSolution = null;			
			newSolution = particles[particleID].update(currentFunctionNum, currentPSOTopology, currentPSOSelfModel, currentPSOInfluenceModel);  

			
			double newParticleValue = newSolution.getFunctionValue();      

			if (newParticleValue < globalBest.getFunctionValue()) {

				globalBest.copyFromPosition(newSolution.getPosition());
				globalBest.setFunctionValue(newParticleValue);
				globalBest.setError(newSolution.getError());
				globalBest.setIterationCreated(PSO.currentIterNum);
				globalBest.setParticleID(particleID);
				globalBest.setApproximated(false);

			} // if new global best

			// data is collected based on number of function evaluations; if we are counting iterations, things are set up in PSO.java
			// so that the number of FEs in each interval is equal to the number of iterations we want in each interval
			// NOTE: PSO.numInitialFEsIgnored is > 0 only if we are counting iterations, since the FEs when the particles are
			// created should not count towards the FEs used to count iterations
			if ((PSO.currentFENum - PSO.numInitialFEsIgnored) % PSO.numFEsPerOutputInterval == 0) {
				int intervalIndex = (PSO.currentFENum - PSO.numInitialFEsIgnored) / PSO.numFEsPerOutputInterval;
				intervalData[intervalIndex][PSO.currentRunNum].copyDataFromSolution(globalBest);
				//				System.out.println("putting item at " + (PSO.currentFunctionEvaluationNum - PSO.numInitialFEsIgnored) + " FEs (not including FEs to create swarm) in index " + intervalRunDataIndex + " of intervalRunData");
			}

		}  // end for-loop particles

	}   


	
	// return the particle Neighborhoods that were calculated when the activation levels were updated
	public Neighborhood[] updateNeurons (PSO.Topology currentFNNTopology, PSO.SelfModel currentFNNSelfModel, 
			PSO.BoundaryModel currentFNNBoundaryModel, PSO.ActivityModel currentFNNActivityModel) {
		
		Neighborhood[] particleNeighs = fnn.moveAndUpdateNeurons(currentFNNTopology, currentFNNSelfModel, currentFNNBoundaryModel, currentFNNActivityModel);

		return particleNeighs;
	}




//	public void updateSendNeighborhoods (int currentFunctionNum, PSO.Topology currentTopology, PSO.SelfModel currentSelfModel, PSO.InfluenceModel currentInfluenceModel,
//			PSO.BoundaryModel currentFNNBoundaryModel, PSO.ActivityModel currentFNNActivityModel, DataOutput[][] intervalData) {
//
//		Neighborhood[] particleNeighs = null;
//		if (PSO.isFNNTopology(currentTopology)) {
//			// updateNeurons returns the particle Neighborhoods that were calculated when the activation levels were updated;
//			// send that to updateParticles so they can be used there instead of recalculating the neighborhoods
//			particleNeighs = updateNeurons(currentTopology, currentSelfModel, currentFNNBoundaryModel, currentFNNActivityModel);
//		}
//
//		updateParticles(currentFunctionNum, currentTopology, currentSelfModel, currentInfluenceModel, 
//				currentFNNBoundaryModel, currentFNNActivityModel, particleNeighs, intervalData);
//
//	}
//
//
//	public void updateParticlesSendNeighborhoods (int currentFunctionNum, PSO.Topology currentTopology, PSO.SelfModel currentSelfModel, PSO.InfluenceModel currentInfluenceModel, 
//			PSO.BoundaryModel currentFNNBoundaryModel, PSO.ActivityModel currentFNNActivityModel, Neighborhood[] particleNeighs,
//			DataOutput[][] intervalData) {
//
//		for (int particleID = 0 ; particleID < particles.length ; particleID++) {
//
//			Solution newSolution = null;
//			// if it's a static topology, we won't be constructing the array of neighborhoods, one for each 
//			// particle and Particle.udpate will not try to use such an array, so just send null
//			if (PSO.isStaticTopology(currentTopology)) {
//				newSolution = particles[particleID].update(currentFunctionNum, currentTopology, currentSelfModel, currentInfluenceModel,
//						currentFNNBoundaryModel, currentFNNActivityModel, null);  
//			}
//			// otherwise, send the Neighborhood for the Particle being updated
//			else if (PSO.isFNNTopology(currentTopology)) {
//				newSolution = particles[particleID].update(currentFunctionNum, currentTopology, currentSelfModel, currentInfluenceModel,
//						currentFNNBoundaryModel, currentFNNActivityModel, particleNeighs[particleID]);  
//			}
//
//			double newParticleValue = newSolution.getFunctionValue();      
//
//			if (newParticleValue < globalBest.getFunctionValue()) {
//
//				globalBest.copyFromPosition(newSolution.getPosition());
//				globalBest.setFunctionValue(newParticleValue);
//				globalBest.setError(newSolution.getError());
//				globalBest.setIterationCreated(PSO.currentIterNum);
//				globalBest.setParticleID(particleID);
//				globalBest.setApproximated(false);
//
//			} // if new global best
//
//			// data is collected based on number of function evaluations; if we are counting iterations, things are set up in PSO.java
//			// so that the number of FEs in each interval is equal to the number of iterations we want in each interval
//			// NOTE: PSO.numInitialFEsIgnored is > 0 only if we are counting iterations, since the FEs when the particles are
//			// created should not count towards the FEs used to count iterations
//			if ((PSO.currentFENum - PSO.numInitialFEsIgnored) % PSO.numFEsPerOutputInterval == 0) {
//				int intervalIndex = (PSO.currentFENum - PSO.numInitialFEsIgnored) / PSO.numFEsPerOutputInterval;
//				intervalData[intervalIndex][PSO.currentRunNum].copyDataFromSolution(globalBest);
//				//				System.out.println("putting item at " + (PSO.currentFunctionEvaluationNum - PSO.numInitialFEsIgnored) + " FEs (not including FEs to create swarm) in index " + intervalRunDataIndex + " of intervalRunData");
//			}
//
//		}  // end for-loop particles
//
//	}   


	
	

	public void connectSwarmToFNN () {

		for (int particleID = 0 ; particleID < particles.length ; particleID++) {

			Particle p = particles[particleID];
			Neuron n = fnn.getNeuron(particleID);

			p.setNeuron(n);
			n.setParticle(p);
		}

	}



	public void createNeighborhoods(PSO.Topology currentTopology, PSO.SelfModel currentSelfModel, PSO.InfluenceModel currentInfluenceModel, int numParticlesIndex) {

		// this is bad for GBEST, because it's going to create a distinct Neighborhood
		// object for every particle, in spite of the fact that there is only one neighborhood...
		// but go with it for now

		for (int particleID = 0 ; particleID < particles.length ; particleID++) {

			// create the actual neighborhood with list of references to Particles in neighborhood
			Neighborhood neigh = new Neighborhood(particles, numParticlesIndex, particleID, currentTopology, currentSelfModel, currentInfluenceModel);	
			allNeighs[particleID] = neigh;
			particles[particleID].setNeighborhood(neigh);
		}

	}


	// for each particle, create an array of neighborhoods that it is in
	public void createNeighLists() {

		// do it for every particle
		for (int particleID = 0 ; particleID < particles.length ; particleID++) {
			Particle thisParticle = particles[particleID];

			// make an array of the neighborhoods the particle is in by going through all the 
			// neighborhoods and checking if this particle is in them; make the array big enough 
			// to hold *all* neighborhoods
			// NOTE: this way of doing it takes care of *any* neighborhood type (e.g. asymmetric 
			//       neighborhoods), since we are actually going through each neighborhood and
			//       checking who is in there *after* the neighborhoods have been created
			int numNeighsContainingParticle = 0;
			Neighborhood[] neighsContainingParticle = new Neighborhood[allNeighs.length];

			for (int neighID = 0 ; neighID < allNeighs.length ; neighID++) {	
				Neighborhood thisNeigh = allNeighs[neighID];
				if (thisNeigh.containsParticle(thisParticle)) {
					neighsContainingParticle[numNeighsContainingParticle++] = thisNeigh;
				}
			}

			// now create the right size array neighsContainingParticle in the particle
			// NOTE: this is bad for the gbest neighborhood, since every particle will have the
			//   same list of all particles (unless the neighborhood does not include self)
			thisParticle.initializeNeighsList(numNeighsContainingParticle);			
			// go through the first numNeighsContainingParticle elements in the array and
			// put the neighborhoods contained therein in the neighsContainingParticle array
			// of the particle
			for (int n = 0 ; n < numNeighsContainingParticle ; ++n) {
				thisParticle.addNeigh(neighsContainingParticle[n]);
			}
		}
	}



	public void showAllNeighborhoods () {

		for (int n = 0 ; n < allNeighs.length ; ++n) {
			showNeighborhood(allNeighs[n]);
		}

	}


	public static void showNeighborhood (Neighborhood neigh) {

		Particle[] particlesInNeighborhood = neigh.getNeighParticles();

		System.out.println("neighborhoodID = " + neigh.getNeighID() + " includes particles:");
		for (int i = 0 ; i < particlesInNeighborhood.length ; i++) {
			System.out.println("particleID = " + particlesInNeighborhood[i].getParticleID());
		}

	}


	public void printParticles () {

		for(int p = 0 ; p < particles.length ; ++p) {
			System.out.println("particle " + p + " = " + particles[p]);
		}

	}



	public void printParticlesAndNeurons() {

		for(int p = 0 ; p < particles.length ; ++p) {
			System.out.println("particle " + p + " = " + particles[p] + " and has neuron " + particles[p].getNeuron());
		}

	}


	public void checkParticlesAndNeurons() {

		for(int p = 0 ; p < particles.length ; ++p) {
			System.out.println("particle " + p + " = " + particles[p] + 
					"   and has neuron " + particles[p].getNeuron() + 
					" which has particle \n              " + particles[p].getNeuron().getParticle() +
					" which has neuron " + particles[p].getNeuron().getParticle().getNeuron());
			System.out.println();
		}

	}



}




