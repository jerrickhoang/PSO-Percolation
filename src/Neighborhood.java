


public class Neighborhood {


	private int particleID;       // the particle whose neighborhood this is
	private int neighID;
	private Particle[] neighParticles;  
	private int nextNeighIndex;   // only used when adding neighborhoods; keeps track of index of next one



	// this is used when we are creating neighborhoods in FluidNN;
	// it is called from FluidNN.getSumActivations because I am calculating the neighborhood of neurons
	// for a given neuron in order to sum the activation levels in the neighborhood;
	// while we have the neighborhood of Neurons, we'll create the neighborhood of Particles using
	// this constructor and the addNeighbor method in this class;
	// this constructor just creates an array big enough to hold all the neighbors -- 
	// neighbors are added using the addNeighbor;
	// I'm tagging the neighborhood with the particleID, but that's not really necessary, since
	// these neighborhoods are being put into an array in FluidNN where the array indices 
	// correspond to particle IDS
	public Neighborhood(int size, int particleID) {
		this.particleID = particleID;
		this.neighID = particleID;
		nextNeighIndex = 0;
		neighParticles = new Particle[size];
	}

	
	// after the particles are created in the OptSwarm constructor, the constructor calls the createNeighborhoods method, 
	// which calls this constructor for each particle to create the Neighborhood in that particle ("neigh")
	public Neighborhood (Particle[] particles, int numParticlesIndex, int particleID, 
			PSO.Topology currentTopology, PSO.SelfModel currentSelfModel, PSO.InfluenceModel currentInfluenceModel) {


		this.particleID = particleID;
		this.neighID = particleID;
		nextNeighIndex = 0;

		if (currentTopology == PSO.Topology.GBEST) {

			if (currentSelfModel == PSO.SelfModel.INCLUDE_SELF) {
				neighParticles = new Particle[particles.length];
				for (int partID = 0 ; partID < particles.length ; partID++) {
					neighParticles[partID] = particles[partID];
				}
			}

			else {
				int nextParticleIndex = 0;  // need this to keep track when self is not included and that partID is skipped
				neighParticles = new Particle[particles.length-1];
				for (int partID = 0 ; partID < particles.length ; partID++) {
					if (partID != particleID)
						neighParticles[nextParticleIndex++] = particles[partID];
				}

			}


		}

		else if (currentTopology == PSO.Topology.RING) {

			// the non-self particles are always the first two in the list; 
			// if self is included, make the array size 3 and put the self at the end
			if (currentSelfModel == PSO.SelfModel.INCLUDE_SELF) {
				neighParticles = new Particle[3];
				neighParticles[2] = particles[particleID];
			}
			else {
				neighParticles = new Particle[2];
			}

			int leftIndex = particleID == 0? particles.length - 1: particleID - 1;				
			int rightIndex = particleID == particles.length - 1? 0: particleID + 1;

			neighParticles[0] = particles[leftIndex];
			neighParticles[1] = particles[rightIndex];

		}

		else if (currentTopology == PSO.Topology.vonNEUMANN) {

			// the non-self particles are always the first four in the list; 
			// if self is included, make the array size 5 and put the self at the end
			if (currentSelfModel == PSO.SelfModel.INCLUDE_SELF) {
				neighParticles = new Particle[5];
				neighParticles[4] = particles[particleID];
			}
			else {
				neighParticles = new Particle[4];
			}

			// get the dimensions of the torus from PSO, where they are set by hand
			// for each possible # of particles
			int numRowsVonNeumann = PSO.numRowsVonNeumannAndMooreList[numParticlesIndex];
			int numColsVonNeumann = PSO.numColsVonNeumannAndMooreList[numParticlesIndex];

			int row = particleID  / numColsVonNeumann;
			int col = particleID  % numColsVonNeumann;

			int northRow = row - 1 < 0? numRowsVonNeumann - 1: row - 1;
			int northCol = col;
			int particleNum = (northRow * numColsVonNeumann) + northCol;
			neighParticles[0] = particles[particleNum];

			int eastRow = row;
			int eastCol = (col + 1) % numColsVonNeumann;
			particleNum = (eastRow * numColsVonNeumann) + eastCol;
			neighParticles[1] = particles[particleNum];

			int southRow = (row + 1) % numRowsVonNeumann;
			int southCol = col;
			particleNum = (southRow * numColsVonNeumann) + southCol;
			neighParticles[2] = particles[particleNum];

			int westRow = row;
			int westCol = col - 1 < 0? numColsVonNeumann - 1: col - 1;
			particleNum = (westRow * numColsVonNeumann) + westCol;
			neighParticles[3] = particles[particleNum];			

		}

		else if (currentTopology == PSO.Topology.MOORE) {

			// the non-self particles are always the first eight in the list; 
			// if self is included, make the array size 9 and put the self at the end
			if (currentSelfModel == PSO.SelfModel.INCLUDE_SELF) {
				neighParticles = new Particle[9];
				neighParticles[8] = particles[particleID];
			}
			else {
				neighParticles = new Particle[8];
			}

			// get the dimensions of the torus from PSO, where they are set by hand
			// for each possible # of particles
			int numRowsMoore = PSO.numRowsVonNeumannAndMooreList[numParticlesIndex];
			int numColsMoore = PSO.numColsVonNeumannAndMooreList[numParticlesIndex];

			int row = particleID  / numColsMoore;
			int col = particleID  % numColsMoore;

			int nextParticleIndex = 0;
			for (int rDelta = -1 ; rDelta <= 1 ; ++rDelta) {
				for (int cDelta = -1 ; cDelta <= 1 ; ++cDelta) {

					// don't do this for the particle itself
					if (rDelta != 0 || cDelta != 0) {
						int neighRow = row + rDelta;
						int neighCol = col + cDelta;

						if (neighRow < 0)
							neighRow = numRowsMoore - 1;
						else if (neighRow == numRowsMoore) {
							neighRow = 0;
						}

						if (neighCol < 0)
							neighCol = numColsMoore - 1;
						else if (neighCol == numColsMoore) {
							neighCol = 0;
						}

						int particleNum = (neighRow * numColsMoore) + neighCol;
						neighParticles[nextParticleIndex++] = particles[particleNum];
					}
				}
			}
		}
	}


	
	// **********************************************************************************************************************************
	// these three methods do *not* need to be static because they will be called only when the neighborhood is static (in the 
	// neighborhood sense, i.e. not changing each iteration), in which case there will be a Neighborhood object
	// through which this can be called
	public DoubleVector getVectorToNeighBestPosition(Particle particle, PSO.Topology currentTopology, PSO.SelfModel currentSelfModel) {

		DoubleVector neighBestPosition = getNeighBestPosition(currentTopology, currentSelfModel);
		DoubleVector vectorToNeighBestPosition = DoubleVector.sub(neighBestPosition, particle.getPosition());

		return vectorToNeighBestPosition;

	}

	
	public DoubleVector getNeighBestPosition(PSO.Topology currentTopology, PSO.SelfModel currentSelfModel) {

		// if it's the standard gbest including self, it's faster to just return the
		// true global best that we're keeping track of in Swarm.java
		if (currentTopology == PSO.Topology.GBEST && currentSelfModel == PSO.SelfModel.INCLUDE_SELF) {
			return Swarm.globalBest.getPosition();
		}

		// whether the self is included was dealt with when the neighborhood was created;
		// if the self is not supposed to be in the neighborhood, it's not (see constructor for details)
		// so just go through the list of particles and find the best one
		Particle bestParticle = null;
		double bestPBestFuncVal = Double.MAX_VALUE;

		int numPartsInNeigh = neighParticles.length;
		for (int p = 0 ; p < numPartsInNeigh ; ++p) {
			Particle nextParticle = neighParticles[p];
			double nextPBestFuncVal = nextParticle.getPersonalBest().getFunctionValue();
			if (nextPBestFuncVal < bestPBestFuncVal) {
				bestParticle = nextParticle;
				bestPBestFuncVal = nextPBestFuncVal;
			}
		}

		return bestParticle.getPersonalBest().getPosition();
	}

	
	public DoubleVector getFIPSAcceleration (Particle particle) {

		DoubleVector position = particle.getPosition();
		DoubleVector acceleration = new DoubleVector(position.size(), 0.0);

		int numPartsInNeigh = neighParticles.length;
		double componentTheta = PSO.theta / numPartsInNeigh;

		// whether the self is included was dealt with when the neighborhood was created;
		// if the self is not supposed to be in the neighborhood, it's not (see constructor for details)
		// so just go through the list of particles and do the standard FIPS calculation
		for (int p = 0 ; p < numPartsInNeigh ; ++p) {
			Particle nextParticle = neighParticles[p];
			DoubleVector vectorToNextPBest = DoubleVector.sub(nextParticle.getPersonalBest().getPosition(), position);
			vectorToNextPBest.multRandomScalar(0.0, componentTheta);
			acceleration.addVector(vectorToNextPBest);
		}

		return acceleration;
	}

	// **********************************************************************************************************************************

	
	
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// THESE VERSIONS CONSTRUCT THE NEIGHBORHOOD
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// these three methods *do* need to be static because they will be called when the neighborhood is changing each iteration
	// (using a fluid neural network), in which case there will not be a Neighborhood object through which a method
	// can be called
	public static DoubleVector getVectorToNeighBestPositionFNN(Particle particle, PSO.Topology currentPSOTopology, PSO.SelfModel currentPSOSelfModel) {

		DoubleVector neighBestPosition = getNeighBestPositionFNN(particle, currentPSOTopology, currentPSOSelfModel);
		DoubleVector vectorToNeighBestPosition = DoubleVector.sub(neighBestPosition, particle.getPosition());

		return vectorToNeighBestPosition;

	}

	
	public static DoubleVector getNeighBestPositionFNN(Particle particle, PSO.Topology currentPSOTopology, PSO.SelfModel currentPSOSelfModel) {

		// topology and selfModel are sent as parameters because there are loops for them in PSO.java to allow us to try different possibilities
		// and we need to send which one is currently being tested
		// the boundary model and activity model are anticipated to be torus and only-active-neurons so much of the time that these values are
		// assigned to static variables in PSO.java.  they can be changed, but there is no loop for either of them that would allow us to test
		// different possibilities, so we just access the static variables here
		Neuron[] neighNeurons = Swarm.fnn.getNeighborhood(particle.getNeuron(), currentPSOTopology, currentPSOSelfModel, PSO.currentPSOBoundaryModel, PSO.currentPSOActivityModel);
		int numNeuronsInNeigh = neighNeurons.length;
		
		
		Particle bestParticle = null;
		double bestPBestFuncVal = Double.MAX_VALUE;

		for (int n = 0 ; n < numNeuronsInNeigh ; ++n) {
			Particle nextParticle = neighNeurons[n].getParticle();
			
			// the is actually determined in the call to getNeighborhood; if the self is not to be included, the 
			// neuron corresponding to this particle is not added into the neighborhood that is returned; this 
			// code is just double-checking
			if (currentPSOSelfModel == PSO.SelfModel.NOT_INCLUDE_SELF && nextParticle.getParticleID() == particle.getParticleID())
				continue;
			
			double nextPBestFuncVal = nextParticle.getPersonalBest().getFunctionValue();
			if (nextPBestFuncVal < bestPBestFuncVal) {
				bestParticle = nextParticle;
				bestPBestFuncVal = nextPBestFuncVal;
			}
		}

		// the particle itself may be the only particle in the neighborhood, and, if INCLUDE_SELF
		// is the case, it will be the "bestParticle" and its pbest position will be returned, BUT
		// if the particle itself is the only particle in the neighborhood, and NOT_INCLUDE_SELF
		// is the case, then there are *no* particles in the neighborhood so bestParticle is still null, so send back
		// position of the particle itself (not its pbest position); this means that in getVectorToNeighBestPositionPM, the
		// vectorToNeighBestPosition returned will be the particle's position minus its position,
		// i.e. the zero vector, which is correct since, if the neighborhood is empty, the neighborhood best 
		// component added to the acceleration being computed in the update method in Particle.java 
		// should be zero
		if (bestParticle == null)
			return particle.getPosition();
		
		return bestParticle.getPersonalBest().getPosition();

	}

	public static DoubleVector getFIPSAccelerationFNN (Particle particle, PSO.Topology currentPSOTopology, PSO.SelfModel currentPSOSelfModel) {

		// topology and selfModel are sent as parameters because there are loops for them in PSO.java to allow us to try different possibilities
		// and we need to send which one is currently being tested
		// the boundary model and activity model are anticipated to be torus and only-active-neurons so much of the time that these values are
		// assigned to static variables in PSO.java.  they can be changed, but there is no loop for either of them that would allow us to test
		// different possibilities, so we just access the static variables here
		Neuron[] neighNeurons = Swarm.fnn.getNeighborhood(particle.getNeuron(), currentPSOTopology, currentPSOSelfModel, PSO.currentPSOBoundaryModel, PSO.currentPSOActivityModel);
		int numNeuronsInNeigh = neighNeurons.length;
		double componentTheta = PSO.theta / numNeuronsInNeigh;

		DoubleVector position = particle.getPosition();
		DoubleVector acceleration = new DoubleVector(position.size(), 0.0);
		
		for (int n = 0 ; n < numNeuronsInNeigh ; ++n) {
			Particle nextParticle = neighNeurons[n].getParticle();
			
			// the is actually determined in the call to getNeighborhood; if the self is not to be included, the 
			// neuron corresponding to this particle is not added into the neighborhood that is returned; this 
			// code is just double-checking
			if (currentPSOSelfModel == PSO.SelfModel.NOT_INCLUDE_SELF && nextParticle.getParticleID() == particle.getParticleID())
				continue;
			
			DoubleVector vectorToNextPBest = DoubleVector.sub(nextParticle.getPersonalBest().getPosition(), position);
			vectorToNextPBest.multRandomScalar(0.0, componentTheta);
			acceleration.addVector(vectorToNextPBest);
		}

		return acceleration;
	}

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	
	
//	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	// THESE VERSIONS DO *NOT* CONSTRUCT THE NEIGHBORHOOD; INSTEAD THEY USE THE NEIGHBORHOOD PASSED IN THROUGH particleNeigh,
//	// WHICH IS THE NEIGHBORHOOD FOR THIS PARTICLE THAT WAS CONSTRUCTED WHEN ITS ACTIVATION STATUS WAS UPDATED
//	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	// these three methods *do* need to be static because it will be called when the neighborhood is changing each iteration
//	// (using a fluid neural network), in which case there will not be a Neighborhood object through which a method
//	// can be called
//	public static DoubleVector getVectorToNeighBestPositionFNN(Particle particle, PSO.SelfModel currentSelfModel, Neighborhood particleNeigh) {
//
//		DoubleVector neighBestPosition = getNeighBestPositionFNN(particle, currentSelfModel, particleNeigh);
//		DoubleVector vectorToNeighBestPosition = DoubleVector.sub(neighBestPosition, particle.getPosition());
//
//		return vectorToNeighBestPosition;
//
//	}
//
//	
//	public static DoubleVector getNeighBestPositionFNN(Particle particle, PSO.SelfModel currentSelfModel, Neighborhood particleNeigh) {
//
//		Particle[] theseNeighParticles = particleNeigh.getNeighParticles();
//		int numParticlesInNeigh = theseNeighParticles.length;
//		
//		
//		Particle bestParticle = null;
//		double bestPBestFuncVal = Double.MAX_VALUE;
//
//		for (int p = 0 ; p < numParticlesInNeigh ; ++p) {
//			Particle nextParticle = theseNeighParticles[p];
//			
//			// the is actually determined in the call to getNeighborhood; if the self is not to be included, the 
//			// neuron corresponding to this particle is not added into the neighborhood that is returned; this 
//			// code is just double-checking
//			if (currentSelfModel == PSO.SelfModel.NOT_INCLUDE_SELF && nextParticle.getParticleID() == particle.getParticleID())
//				continue;
//			
//			double nextPBestFuncVal = nextParticle.getPersonalBest().getFunctionValue();
//			if (nextPBestFuncVal < bestPBestFuncVal) {
//				bestParticle = nextParticle;
//				bestPBestFuncVal = nextPBestFuncVal;
//			}
//		}
//
//		// the particle itself may be the only particle in the neighborhood, and, if INCLUDE_SELF
//		// is the case, it will be the "bestParticle" and its pbest position will be returned, BUT
//		// if the particle itself is the only particle in the neighborhood, and NOT_INCLUDE_SELF
//		// is the case, then there are *no* particles in the neighborhood so bestParticle is still null, so send back
//		// position of the particle itself (not its pbest position); this means that in getVectorToNeighBestPositionPM, the
//		// vectorToNeighBestPosition returned will be the particle's position minus its position,
//		// i.e. the zero vector, which is correct since, if the neighborhood is empty, the neighborhood best 
//		// component added to the acceleration being computed in the update method in Particle.java 
//		// should be zero
//		if (bestParticle == null)
//			return particle.getPosition();
//		
//		return bestParticle.getPersonalBest().getPosition();
//
//	}
//
//	
//	public static DoubleVector getFIPSAccelerationFNN (Particle particle, PSO.SelfModel currentSelfModel, Neighborhood particleNeigh) {
//
//		Particle[] theseNeighParticles = particleNeigh.getNeighParticles();
//		int numParticlesInNeigh = theseNeighParticles.length;
//		double componentTheta = PSO.theta / numParticlesInNeigh;
//
//		DoubleVector position = particle.getPosition();
//		DoubleVector acceleration = new DoubleVector(position.size(), 0.0);
//		
//		for (int p = 0 ; p < numParticlesInNeigh ; ++p) {
//			Particle nextParticle = theseNeighParticles[p];
//			
//			// the is actually determined in the call to getNeighborhood; if the self is not to be included, the 
//			// neuron corresponding to this particle is not added into the neighborhood that is returned; this 
//			// code is now redundant (1/13/14), but leaving it in for now
//			if (currentSelfModel == PSO.SelfModel.NOT_INCLUDE_SELF && nextParticle.getParticleID() == particle.getParticleID())
//				continue;
//			
//			DoubleVector vectorToNextPBest = DoubleVector.sub(nextParticle.getPersonalBest().getPosition(), position);
//			vectorToNextPBest.multRandomScalar(0.0, componentTheta);
//			acceleration.addVector(vectorToNextPBest);
//		}
//
//		return acceleration;
//	}
//
//	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	
		
	
	public boolean containsParticle (Particle particle) {

		for (int p = 0 ; p < neighParticles.length ; ++p) {
			if (neighParticles[p] == particle)	{
				return true;
			}
		}

		return false;
	}



	public void addNeighbor(Particle particle) {
		neighParticles[nextNeighIndex++] = particle;
	}



	
	public Particle[] getNeighParticles() {

		return neighParticles;

	}

	
	public int sizeNeigh() {
		return neighParticles.length;
	}


	public int getNeighID () {
		return neighID;
	}

	


}
