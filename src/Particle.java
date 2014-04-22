



public class Particle {

	private int particleID;

	private DoubleVector position;  

	private DoubleVector velocity;
	//	private DoubleVector nextVelocity;

	private Solution currSolution;
	private Solution personalBest;

	// the neighborhood that this particle is in
	// (which is also included in neighsContainingParticle (next variable)
	private Neighborhood staticNeigh;
	// a list of neighborhoods that contain this particle
	private Neighborhood[] neighsContainingParticle;
	// used in addNeigh when adding neighbors into the list;
	// keeps track of where the next one should go
	private int nextNeighIndex = 0;

	private Neuron neuron;

	private boolean pBestImproved = false;
	

	public Particle(int functionNum, int numDimensions, int particleID, double[] sendBackResults) {

		this.particleID = particleID;

		position = new DoubleVector(numDimensions);
		for(int i = 0 ; i < position.size() ; ++i) {
			position.set(i, TestFunctions.INIT_MIN_VALS[functionNum] + (TestFunctions.INIT_RANGES[functionNum] * PSO.rand.nextDouble()));
		}

		double[] results = TestFunctions.evalWithError(position, functionNum);
		// need to get the function value and error back to the Swarm constructor,
		// so we can determine the initial global best
		sendBackResults[TestFunctions.VAL_INDEX] =  results[TestFunctions.VAL_INDEX];
		sendBackResults[TestFunctions.ERR_INDEX] =  results[TestFunctions.ERR_INDEX];

		// can send position itself because the Solution constructor makes a copy of the position DoubleVector sent in
		currSolution = new Solution(position, results[TestFunctions.VAL_INDEX], results[TestFunctions.ERR_INDEX], 0, particleID, false);
		personalBest = currSolution.getCopy();

		// start with small random velocity
		velocity = new DoubleVector(numDimensions);
		double minSpeed = 0.0;
		double speedRange = 0.0;
		// don't let the initial speed be greater than a small amount
		if (TestFunctions.UNIVERSAL_SPEED_RANGE < TestFunctions.SPEED_RANGES[functionNum]) {
			minSpeed = TestFunctions.UNIVERSAL_MIN_INIT_SPEED;
			speedRange = TestFunctions.UNIVERSAL_SPEED_RANGE;
		}
		else {
			minSpeed = TestFunctions.SPEED_MIN_VALS[functionNum];
			speedRange = TestFunctions.SPEED_RANGES[functionNum];			
		}
		for(int i = 0 ; i < velocity.size() ; ++i) {
			velocity.set(i, minSpeed + (speedRange * PSO.rand.nextDouble()));
		}


	}


	// create the neighborhood list of the specified size
	public void initializeNeighsList (int numNeighs) {
		neighsContainingParticle = new Neighborhood[numNeighs];
	}


	// used in createNeighLists in Swarm.java when creating, for each particle, a list
	// of Neighborhoods that it is a member of
	public void addNeigh (Neighborhood neigh) {
		neighsContainingParticle[nextNeighIndex++] = neigh;
	}



	// the velocity and position update
	public Solution update(int functionNum, PSO.Topology currentPSOTopology, PSO.SelfModel currentPSOSelfModel, PSO.InfluenceModel currentPSOInfluenceModel) {


		// acceleration starts at 0.0
		DoubleVector acceleration = new DoubleVector(position.size(), 0.0);

		// static topology
		if (PSO.isStaticTopology(currentPSOTopology)) {

			// the heavy lifting is done in the Neighborhood class;
			// much of it has already been done when the neighborhoods were created once an for all
			// before the iterations started
			if (currentPSOInfluenceModel == PSO.InfluenceModel.NEIGH_BEST) {

				// neighborhood best component
				DoubleVector neighborhoodComponent = staticNeigh.getVectorToNeighBestPosition(this, currentPSOTopology, currentPSOSelfModel);
				neighborhoodComponent.multRandomScalar(0.0, PSO.neighborhoodTheta);
				acceleration.addVector(neighborhoodComponent);

				// personal component
				DoubleVector personalComponent = DoubleVector.sub(personalBest.getPosition(), position);
				personalComponent.multRandomScalar(0.0, PSO.personalTheta);
				acceleration.addVector(personalComponent);
			}

			// NOTE:  FIPS does not have an explicit personal component in addition to factoring in the 
			//        particle's pbest (along with all the other pbests in the neighborhood) 
			else if (currentPSOInfluenceModel == PSO.InfluenceModel.FIPS) {

				acceleration.addVector(staticNeigh.getFIPSAcceleration(this));	

			}
		}

		// dynamic topology
		
		else if (PSO.isFNNTopology(currentPSOTopology)) {

			// the heavy lifting is done in the Neighborhood class; none of it was done before the
			// iterations started -- all of it is done dynamically on each iteration
			if (currentPSOInfluenceModel == PSO.InfluenceModel.NEIGH_BEST) {

				// neighborhood best component
				DoubleVector neighborhoodComponent = Neighborhood.getVectorToNeighBestPositionFNN(this, currentPSOTopology, currentPSOSelfModel);
				
				neighborhoodComponent.multRandomScalar(0.0, PSO.neighborhoodTheta);
				acceleration.addVector(neighborhoodComponent);
				
				// personal component
				DoubleVector personalComponent = DoubleVector.sub(personalBest.getPosition(), position);
				personalComponent.multRandomScalar(0.0, PSO.personalTheta);
				acceleration.addVector(personalComponent);
			}

			// NOTE:  FIPS does not have an explicit personal component in addition to factoring in the 
			//        particle's pbest (along with all the other pbests in the neighborhood) 
			else if (currentPSOInfluenceModel == PSO.InfluenceModel.FIPS) {

				acceleration.addVector(Neighborhood.getFIPSAccelerationFNN(this, currentPSOTopology, currentPSOSelfModel));	
				
			}

		}
		
		else {
			System.out.println("error: topology unspecifed in Particle.update");
			System.exit(-1);
		}

		
		
		// update the velocity and apply the constriction factor
		velocity.addVector(acceleration);
		velocity.multScalar(PSO.constrictionFactor);



		// bound velocity
		for (int i = 0 ; i < velocity.size() ; ++i) {
			if (velocity.get(i) < TestFunctions.SPEED_MIN_VALS[functionNum])
				velocity.set(i, TestFunctions.SPEED_MIN_VALS[functionNum]);
			else if (velocity.get(i) > TestFunctions.SPEED_MAX_VALS[functionNum])
				velocity.set(i, TestFunctions.SPEED_MAX_VALS[functionNum]);
		}


		// move the particle 
		position.addVector(velocity); 


		// evaluate the new position and set currentSolution
		double[] results = TestFunctions.evalWithError(position, functionNum);
		double newPositionValue = results[TestFunctions.VAL_INDEX];
		double newPositionError = results[TestFunctions.ERR_INDEX];

		currSolution.copyFromPosition(position);
		currSolution.setFunctionValue(newPositionValue);
		currSolution.setError(newPositionError);
		currSolution.setIterationCreated(PSO.currentIterNum);

		// update the personal best, if necessary
		pBestImproved = false;
		if (newPositionValue < personalBest.getFunctionValue()) {

			personalBest.copyFromPosition(position);
			personalBest.setFunctionValue(newPositionValue);
			personalBest.setError(newPositionError);
			personalBest.setIterationCreated(PSO.currentIterNum);
			
			pBestImproved = true;
			
			// pbest improved, so decrease gain if it's an fnn topology and we're using individual gains
			if (PSO.isFNNTopology(currentPSOTopology) && PSO.useIndividualGains) {
				
				//			this.getNeuron().setGain(PSO.lowGain);
				++PSO.totalNumGainDecreases;
				double reducedGain = this.getNeuron().getGain() - PSO.gainDecrement;
				if (reducedGain < 0.0) {
					reducedGain = 0.0;
				}
				this.getNeuron().setGain(reducedGain);
				
			}
			
		}
		// pbest did not improve, so increase gain if it's an fnn topology and we're using individual gains
		else if (PSO.isFNNTopology(currentPSOTopology) && PSO.useIndividualGains) {
			pBestImproved = false;   // not actually needed, since pBestImproved is set to 
								     // false before we check for improvement 
			
			//			this.getNeuron().setGain(PSO.highGain);
			++PSO.totalNumGainIncreases;
			double increasedGain = this.getNeuron().getGain() + PSO.gainIncrement;
			if (increasedGain > 1.0) {
				increasedGain = 1.0;
			}
			this.getNeuron().setGain(increasedGain);

		}

		
		if (PSO.isFNNTopology(currentPSOTopology)) {
			PSO.sumGainValues += this.getNeuron().getGain();
		}
		
		
		return currSolution;
	}		

	
	

	public boolean pBestImproved () {
		return pBestImproved;
	}
	

	public int getParticleID () {
		return particleID;
	}


	public DoubleVector getPosition() {
		return position;
	}



	public Solution getPersonalBest() {
		return personalBest;
	}



	public void setPersonalBest(Solution personalBest) {
		this.personalBest = personalBest;
	}



	public Neuron getNeuron() {
		return neuron;
	}



	public void setNeuron(Neuron neuron) {
		this.neuron = neuron;
	}




	public Neighborhood getNeighborhood() {
		return staticNeigh;
	}



	public void setNeighborhood(Neighborhood neigh) {
		this.staticNeigh = neigh;
	}




	public void printPosition() {
		position.print();
	}


	public void printlnPosition() {
		position.println();
	}




}


