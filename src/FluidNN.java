import java.util.Random;

/**
 * Implements a Fluid Neural Network as described in:
 * 		Sole and Miramontes
 * 		"Information at the edge of chaos in fluid neural networks"
 * 		Physica D 80 (1995) 171-180
 * 
 * Another relevant paper:
 * 		Delgado and Sole
 * 		"Mean-field theory of fluid neural networks"
 * 		Physcial Review E 57(2) (1998) 2204-11
 */

/**
 * @author Stephen Majercik
 * 12/03/13
 *
 */
public class FluidNN {

	// a list of the neurons
	private Neuron[] neuronList;
	// the nodes/automata of an FNN live on a lattice/grid
	private Neuron[][] grid;
	// size of grid
	private int numRows;
	private int numCols;
	// number of neurons
	private int numNeurons;

	// this is the "coupling matrix" in Sole & Miramontes:
	//
	// J_ij is an arbitrary function of S_i and S_j, the activation levels
	// of neurons i and j, and is a multiplier in the activation formula:
	// in the paper, however, J is just a 2x2 matrix of constants:
	//
	//    lambda_11  lambda 12
	//    lambda_21  lambda 22
	//
	// where:
	//   J_ij = lambda_11 if both i and j are active 
	//   J_ij = lambda_12 if i is active and j is inactive 
	//   J_ij = lambda_21 if i is inactive and j is active 
	//   J_ij = lambda_22 if both i and j are inactive 
	//
	// in the paper, in all their experiments, all lambdas were 1.0
	//
	// NOTE: J is 2x2 only because we are defining only two states
	// for a neuron; in general, J will be k x k, where k is the number
	// of states a neuron can be in
	//
	// for now, we will follow the paper
	private double[][] J = { { 1.0, 1.0 }, { 1.0, 1.0 } };


	// NOTE: all of the values below are defaults and will usually be over-ridden (with
	//       the exception of the constants specifying the range of initial
	//       activation levels)

	// this is factor that dials the activation level up and down;
	// it is applied to the activation sum before the "squashing
	//  function" (tanh, in this case) is applied
	private double gain = 0.0;

	// new neurons have an activation level between 0.0 and 1.0
	public static double INITIAL_ACTIVATION_LOW_LEVEL = 0.0;
	public static double INITIAL_ACTIVATION_HIGH_LEVEL = 1.0;
	public static double INITIAL_ACTIVATION_RANGE = 
		INITIAL_ACTIVATION_HIGH_LEVEL- INITIAL_ACTIVATION_LOW_LEVEL;

	// there appear to be two thresholds in the Sole & Miramontes paper:
	//   1) one that is used in the sum of neighbor activations calculation
	//		(but is set to 0.0 in the paper)
	private double sumNeighborActivationsThreshold = 0.0;

	// a neuron can become active in one of two ways:
	//   1) its activation level exceeds a threshold:
	//      (in Sole & Miramontes this is 1e-16)
	private double activationThreshold = 1e-16;
	//   2) it is inactive and becomes active spontaneously
	//		at a specified level with a specified probability
	//      In the Sole & Miramontes paper, the level is given 
	//		as 0.1; the probability is not given, but in the
	//		Delgado and Sole paper cited above, they experiment
	//		with values in the range of 0.00003 to 0.3, with
	// 		0.1 appearing to be the "default"
	private double spontaneousActivationLevel = 0.1;
	private double spontaneousActivationProbability = 0.1;   



	// create a FluidNN of a given size (dimensions), but with no neurons
	public FluidNN (int numRows, int numCols) {

		this.numRows = numRows;
		this.numCols = numCols;
		grid = new Neuron[numRows][numCols];
		this.numNeurons = 0;

	}

	public FluidNN (int numRows, int numCols, int numNeurons, double gain,
			double sumNeighborActivationsThreshold, double activationThreshold, 
			double spontaneousActivationLevel, double spontaneousActivationProbability) {

		this.numRows = numRows;
		this.numCols = numCols;
		neuronList = new Neuron[numNeurons];
		grid = new Neuron[numRows][numCols];
		this.numNeurons = numNeurons;
		
		
		// ***************************************************************************
		// ***************************************************************************
		// ***************************************************************************
		// ***************************************************************************
		randomPopulate(numNeurons, gain);
		
//		randomPopulateConnected(numNeurons, PSO.currentFNNBoundaryModel, gain);
		// ***************************************************************************
		// ***************************************************************************
		// ***************************************************************************
		// ***************************************************************************

		
		this.gain = gain;
		this.sumNeighborActivationsThreshold = sumNeighborActivationsThreshold;
		this.activationThreshold = activationThreshold;
		this.spontaneousActivationLevel = spontaneousActivationLevel;
		this.spontaneousActivationProbability = spontaneousActivationProbability;

	}
	
	
	// create a FluidNN of a given size with a given number of neurons (randomly placed)
	// and with specified values for the parameters
	public FluidNN (int numRows, int numCols, int numNeurons, double p) {

		this.numRows = numRows;
		this.numCols = numCols;
		neuronList = new Neuron[numNeurons];
		grid = new Neuron[numRows][numCols];
		this.numNeurons = numNeurons;
		
		randomSequentialPopulate(numNeurons, p);
	}

	public void randomSequentialPopulate(int numNeurons, double probability) {
		if (numNeurons > numRows * numCols) {
			System.out.println("Grid is too small for " + numNeurons + " neurons");
			System.exit(0);
		}
		Random r = new Random();
		int count = 0;
		int n = 0;
		Neuron tempNeuron;
		do {
//			n++;
//			if (n > 1000) {
//				System.out.println(numRows + " " + numCols + " " + probability + " " + numNeurons);
//				System.out.println("Initialization took too long");
//				System.exit(0);
//			}
			neuronList = new Neuron[numNeurons];
			count = 0;
			for (int i = 0; i < numRows; i ++) {
				for (int j = 0; j < numCols; j++) {
					grid[i][j] = null;
					if (r.nextDouble() < probability) {
						double initialActivationLevel = INITIAL_ACTIVATION_LOW_LEVEL + (r.nextDouble() * INITIAL_ACTIVATION_RANGE);
						tempNeuron = new Neuron(count, i, j, initialActivationLevel, 0.0, this);
						grid[i][j] = tempNeuron;
						if (count < numNeurons) neuronList[count] = tempNeuron;
						count ++;
					}
				}
			}
		} while(count != numNeurons);
	}


	// put a specified number of neurons in the net at random locations; 
	// assumes that the grid is empty and that number of neurons is less 
	// than the number of cells in the grid
	public void randomPopulate(int numNeurons, double currentGain) {

		for (int neuronID = 0 ; neuronID < numNeurons ; ++neuronID) {
			Neuron newNeuron = addRandomNeuron(neuronID, currentGain);
			neuronList[neuronID] = newNeuron;
			// null if grid is full; should not happen
			if (newNeuron == null) { 
				System.out.println("error: grid too small in randomPopulate");
				System.exit(0);
			}
		}
	}

	
	// uses randomPopulate to get neurons onto the grid, but then iterates,
	// looking for isolated neurons and relocating them so that they have
	// a neighbor
	public void randomPopulateConnected(int numNeurons, PSO.BoundaryModel currentFNNBoundaryModel, double currentGain) {

		randomPopulate(numNeurons, currentGain);
		
		Neuron isolatedNeuron = getIsolatedNeuron(currentFNNBoundaryModel);
		while (isolatedNeuron != null) {
			relocateIsolatedNeuron(isolatedNeuron, currentFNNBoundaryModel);
			isolatedNeuron = getIsolatedNeuron(currentFNNBoundaryModel);
		}
		
	}


	
	// returns a neuron that is isolated, if there is one
	public Neuron getIsolatedNeuron (PSO.BoundaryModel currentFNNBoundaryModel) {
		
		
		for (int i = 0 ; i < neuronList.length ; ++i) {
			Neuron n = neuronList[i];
			if (!hasNeighbors(n, currentFNNBoundaryModel)) {
				return n;			}
		}

		return null;
		
	}

	

	// relocate the neuron sent in to a place where it has neighbors
	public void relocateIsolatedNeuron (Neuron isolatedNeuron, PSO.BoundaryModel currentFNNBoundaryModel) {
		
		int oldRow = isolatedNeuron.getRow();
		int oldCol = isolatedNeuron.getCol();
		
		// get random row and col
		int newRow = PSO.rand.nextInt(numRows);
		int newCol = PSO.rand.nextInt(numCols);

		// kludgey way of finding an empty location with neighbors
		while (grid[newRow][newCol] != null || !hasNeighbors(newRow, newCol, currentFNNBoundaryModel)) {
			newRow = PSO.rand.nextInt(numRows);
			newCol = PSO.rand.nextInt(numCols);
		}

		// reset the neuron's row and column
		isolatedNeuron.setRow(newRow);
		isolatedNeuron.setCol(newCol);

		// move the neuron to its new location
		grid[newRow][newCol] = isolatedNeuron;
		// be sure to set the old location to null
		grid[oldRow][oldCol] = null;
		
	}


	
	// find a random empty location in the grid (if one exists; if not, 
	// return null), put a neuron there, and return that neuron
	public Neuron addRandomNeuron(int neuronID, double currentGain) {

		if (gridFull()) {
			return null;
		}

		// get random row and col
		int r = PSO.rand.nextInt(numRows);
		int c = PSO.rand.nextInt(numCols);

		// kludgey way of finding an empty location
		while (grid[r][c] != null) {
			r = PSO.rand.nextInt(numRows);
			c = PSO.rand.nextInt(numCols);
		}

		// random initial activation level
		double initialActivationLevel = INITIAL_ACTIVATION_LOW_LEVEL + (PSO.rand.nextDouble() * INITIAL_ACTIVATION_RANGE);
		grid[r][c] = new Neuron(neuronID, r, c, initialActivationLevel, currentGain, this);

		return grid[r][c];

	}


	
	// move all the neurons and update their activation levels;
	// return the particle Neighborhoods assembled in updateActivationLevels
	public Neighborhood[] moveAndUpdateNeurons(PSO.Topology currentFNNTopology, PSO.SelfModel currentFNNSelfModel, 
			PSO.BoundaryModel currentFNNBoundaryModel, PSO.ActivityModel currentFNNActivityModel) {

		// I am going to move them first, since I want to use the neighborhoods that are calculated 
		// during activation updates for the PSO particle updates;
		//   1) tests indicate that it doesn't matter (see 1/23-28/14 entry 0-NOTES-CODE-AND-TESTS in FluidNeuralNetworks)
		//   2) this may not be the order in which Solé & Miramontes do it, but it shouldn't matter
		// moving will ALWAYS used Moore, includeSelf, and only active neurons; allow for the possibility that it might use a torus, instead of a lattice
		moveAllMoore(currentFNNBoundaryModel);

		// update activations
		Neighborhood[] particleNeighs = updateActivationLevels(currentFNNTopology, currentFNNSelfModel, currentFNNBoundaryModel, currentFNNActivityModel);

		return particleNeighs;

	}



	// update activation levels according to the Sole & Miramontes paper;
	// assuming that new activation levels must be calculated entirely from 
	// current activation level, i.e. new levels are *not* used as they
	// become available
	// since we will be computing the neighborhoods and will need those in Particle.udpate, put
	// them in an array and return those
	public Neighborhood[] updateActivationLevels(PSO.Topology currentFNNTopology, PSO.SelfModel currentFNNSelfModel, 
			PSO.BoundaryModel currentFNNBoundaryModel, PSO.ActivityModel currentFNNActivityModel) {

		// create an array of the neighborhoods of each particle to send back for eventual use in
		// Particle.update; although Java uses pass-by value, we are changing the *contents* of the
		// array. so we're okay
		Neighborhood[] particleNeighs = new Neighborhood[neuronList.length];
		
		double[] newActivationLevels = new double[neuronList.length];

		for (int n = 0 ; n < neuronList.length ; ++n) {
			double sumAct = getSumActivations(neuronList[n], currentFNNTopology, currentFNNSelfModel, currentFNNBoundaryModel, currentFNNActivityModel, particleNeighs);  
			
			if (PSO.useIndividualGains) {
				newActivationLevels[n] = Math.tanh(neuronList[n].getGain() * (sumAct - sumNeighborActivationsThreshold));
			}
			else {
				newActivationLevels[n] = Math.tanh(gain * (sumAct - sumNeighborActivationsThreshold));
			}
		}

		// using the raw activation levels just computed, update the neurons' activity status, which 
		// includes both its activation level and whether it is actually "active"
		for (int n = 0 ; n < neuronList.length ; ++n) {
			if (PSO.useNewSpontActMechanism) {
				neuronList[n].updateActivationStatusNewSpontActMechanism(newActivationLevels[n], currentFNNBoundaryModel);
			}
			else {
				neuronList[n].updateActivationStatus(newActivationLevels[n]);
			}
		}

		return particleNeighs;

	}




	// sum activations according to Sole & Miramontes paper; not clear whether
	// the neuron itself is included, but the Delgado & Sole paper seems to
	// indicate that it is, so we do
	public double getSumActivations(Neuron neuron, PSO.Topology currentFNNTopology, PSO.SelfModel currentFNNSelfModel, 
			PSO.BoundaryModel currentFNNBoundaryModel, PSO.ActivityModel currentFNNActivityModel, Neighborhood[] particleNeighs) {

		// NOTE: in the Sole & Miramontes paper, it appears that the activation
		//       level of the neuron itself is also multiplied by a J factor;
		//       it will always be lambda_11 ("both" active) or lambda_22 
		//       ("both" inactive)

		// this is a list of the Neurons that are neighbors
		
		Neuron[] neighbors = getNeighborhood(neuron, currentFNNTopology, currentFNNSelfModel, currentFNNBoundaryModel, currentFNNActivityModel);

		// sum the activations, but ALSO
		// put together the neighborhood of corresponding particles and send that back up the line,
		// so Swarm.update can send it to Swarm.updateParticles, which sends it to Particle.update,
		// which send it to the methods in Neighborhood that need this to do their calculations;
		// previously those methods in Neighborhood recalculated the neighborhood -- a waste, since
		// it has already been done;
		// send particleID to the constructor since we want the 
		// the Neighborhood to have the ID of the particle (though not strictly necessary, since
		// the Neighborhoods are being put into an array of Neighborhoods in updateActivationLevels,
		// and the array indices correspond to the particle IDS
		int particleID = neuron.getParticle().getParticleID();
		Neighborhood thisNeigh = new Neighborhood(neighbors.length, particleID);
		
		double sumActivations = 0.0;
		for (int n = 0 ; n < neighbors.length ; ++n) {
			thisNeigh.addNeighbor(neighbors[n].getParticle());
			sumActivations += getJValueStandard(neuron, neighbors[n]) * neighbors[n].getActivationLevel();
		}
		
		// put the Neighborhood for this particle into the array
		particleNeighs[particleID] = thisNeigh;
		
		return sumActivations;

	}


	
	// returns a list of neurons in the neighborhood
	public Neuron[] getNeighborhood(Neuron neuron, PSO.Topology topology, PSO.SelfModel selfModel, 
			PSO.BoundaryModel boundaryModel, PSO.ActivityModel activityModel) {

		Neuron[] neighborhood = null;

		if (topology == PSO.Topology.FNN_GBEST) {		
			neighborhood = getGbestNeighbors(neuron, selfModel, boundaryModel, activityModel);			
		}

		else if (topology == PSO.Topology.FNN_RING) {
			neighborhood = getRingNeighbors(neuron, selfModel, boundaryModel, activityModel);				
		}

		else if (topology == PSO.Topology.FNN_vonNEUMANN) {
			neighborhood = getVonNeumannNeighbors(neuron, selfModel, boundaryModel, activityModel);				
		}

		else if (topology == PSO.Topology.FNN_MOORE) {
			neighborhood = getMooreNeighbors(neuron, selfModel, boundaryModel, activityModel);			
		}
		else {
			System.out.println("error:  unknown topology in FluidNN.getNeigborhood"); 
			System.exit(-1);
		}

		return neighborhood;
	}
	
	public Neuron[] getNeighborhood(Neuron neuron, PSOPercolation.Topology topology, PSOPercolation.SelfModel selfModel, 
			PSOPercolation.BoundaryModel boundaryModel, PSOPercolation.ActivityModel activityModel) {

		Neuron[] neighborhood = null;

		if (topology == PSOPercolation.Topology.FNN_MOORE) {
			neighborhood = getMooreNeighbors(neuron, selfModel, boundaryModel, activityModel);			
		}
		else {
			System.out.println("error:  unknown topology in FluidNN.getNeigborhood"); 
			System.exit(-1);
		}

		return neighborhood;
	}
	


	
	// NOTE: we don't actually need the boundary model, since the neighborhood is the entire grid, but include it for 
	// the sake of uniformity with respect to the other getXXXNeighborhood methods;
	// returns a list of neurons in the neighborhood
	public Neuron[] getGbestNeighbors (Neuron neuron, PSO.SelfModel selfModel, PSO.BoundaryModel boundaryModel, 
			PSO.ActivityModel activityModel) {

		// create an array that can hold the maximum number of neurons in the neighborhood, including the neuron itself;
		// in the GBEST case, it could be all the neurons
		Neuron[] possibleNeighbors = new Neuron[neuronList.length];
		int numNeigh = 0;

		// don't need to call getNeighbor because boundary model is irrelevant (since the neighborhood is
		// everything) and we can take care of self model and activity model right here
		for (int n = 0 ; n < neuronList.length ; ++n) {
			
			// if we're not including self and the potential neighbor is the neuron itself, skip it
			if (selfModel == PSO.SelfModel.NOT_INCLUDE_SELF && neuronList[n] == neuron) {
				continue;
			}
			
			// if we're only including active neurons and the potential neighbor is not active, skip it
			if (activityModel == PSO.ActivityModel.ONLY_ACTIVE_NEURONS && !neuronList[n].active()) {
				continue;
			}
			
			// if we're only including neurons that have the same activation status, and the
			// potential neighbor doesn't, skip it
			if (activityModel == PSO.ActivityModel.SIMILAR_ACTIVITY_STATUS_NEURONS) {
				boolean sameActivationStatus =  (neuron.active() && neuronList[n].active()) || 
												(!neuron.active() && !neuronList[n].active());
				if (!sameActivationStatus) {	
					continue;
				}
			}
			
			possibleNeighbors[numNeigh++] = neuronList[n];
		}

		// if all the possible neurons are in the neighborhood, can just return the array 
		// they are already in
		if (numNeigh == neuronList.length) {   
			return possibleNeighbors;
		}

		// otherwise transfer the members of the neighborhood to an array that's exactly large 
		// enough to hold them
		Neuron[] actualNeighbors = new Neuron[numNeigh];
		for (int i = 0 ; i < numNeigh ; ++i) {
			actualNeighbors[i] = possibleNeighbors[i];
		}

		return actualNeighbors;

	}



	// this isn't really RING; it's treating each row of the grid/lattice as a ring, which is not quite right,
	// but it's reasonable (although not an exact analogue to standard RING);
	// in any case, RING in the context of FNN neighborhoods does not make a whole lot of sense;
	// returns a list of neurons in the neighborhood
	public Neuron[] getRingNeighbors (Neuron neuron, PSO.SelfModel selfModel, PSO.BoundaryModel boundaryModel, 
			PSO.ActivityModel activityModel) {

		// create an array that can hold the maximum number of neurons in the neighborhood, including the neuron itself
		Neuron[] possibleNeighbors = new Neuron[3];
		int numNeigh = 0;

		// center, i.e. the particle itself
		// NOTE: we could potentially check here if the currentSelfModel is INCLUDE_SELF, but we would also have to 
		// check the current activity model; isNeighbor does both, so call isNeighbor 
		int rDelta = 0;
		int cDelta = 0;
		Neuron possibleNeighbor = getNeighbor(neuron, rDelta, cDelta, selfModel, boundaryModel, activityModel);
		if (possibleNeighbor != null)
			possibleNeighbors[numNeigh++] = possibleNeighbor;

		// right
		rDelta = 0;
		cDelta = 1;
		possibleNeighbor = getNeighbor(neuron, rDelta, cDelta, selfModel, boundaryModel, activityModel);
		if (possibleNeighbor != null)
			possibleNeighbors[numNeigh++] = possibleNeighbor;

		// left
		rDelta = 0;
		cDelta = -1;
		possibleNeighbor = getNeighbor(neuron, rDelta, cDelta, selfModel, boundaryModel, activityModel);
		if (possibleNeighbor != null)
			possibleNeighbors[numNeigh++] = possibleNeighbor;

		// if all the possible neurons are in the neighborhood, can just return the array 
		// they are already in
		if (numNeigh == 3) {   // 3 includes the neuron whose neighbors we are determining
			return possibleNeighbors;
		}

		// otherwise transfer the members of the neighborhood to an array that's exactly large 
		// enough to hold them
		Neuron[] actualNeighbors = new Neuron[numNeigh];
		for (int i = 0 ; i < numNeigh ; ++i) {
			actualNeighbors[i] = possibleNeighbors[i];
		}

		return actualNeighbors;

	}



	// returns a list of neurons in the neighborhood
	public Neuron[] getVonNeumannNeighbors (Neuron neuron, PSO.SelfModel selfModel, PSO.BoundaryModel boundaryModel, 
			PSO.ActivityModel activityModel) {


		// create an array that can hold the maximum number of neurons in the neighborhood, including the neuron itself
		Neuron[] possibleNeighbors = new Neuron[5];
		int numNeigh = 0;

		// center, i.e. the particle itself
		// NOTE: we could potentially check here if the currentSelfModel is INCLUDE_SELF, but we would also have to 
		// check the current activity model; isNeighbor does both, so call isNeighbor 
		int rDelta = 0;
		int cDelta = 0;
		Neuron possibleNeighbor = getNeighbor(neuron, rDelta, cDelta, selfModel, boundaryModel, activityModel);
		if (possibleNeighbor != null)
			possibleNeighbors[numNeigh++] = possibleNeighbor;

		// north
		rDelta = -1;
		cDelta = 0;
		possibleNeighbor = getNeighbor(neuron, rDelta, cDelta, selfModel, boundaryModel, activityModel);
		if (possibleNeighbor != null)
			possibleNeighbors[numNeigh++] = possibleNeighbor;

		// east
		rDelta = 0;
		cDelta = 1;
		possibleNeighbor = getNeighbor(neuron, rDelta, cDelta, selfModel, boundaryModel, activityModel);
		if (possibleNeighbor != null)
			possibleNeighbors[numNeigh++] = possibleNeighbor;

		// south
		rDelta = 1;
		cDelta = 0;
		possibleNeighbor = getNeighbor(neuron, rDelta, cDelta, selfModel, boundaryModel, activityModel);
		if (possibleNeighbor != null)
			possibleNeighbors[numNeigh++] = possibleNeighbor;

		// west
		rDelta = 0;
		cDelta = -1;
		possibleNeighbor = getNeighbor(neuron, rDelta, cDelta, selfModel, boundaryModel, activityModel);
		if (possibleNeighbor != null)
			possibleNeighbors[numNeigh++] = possibleNeighbor;

		// if all the possible neurons are in the neighborhood, can just return the array 
		// they are already in
		if (numNeigh == 5) {   // 5 includes the neuron whose neighbors we are determining
			return possibleNeighbors;
		}

		// otherwise transfer the members of the neighborhood to an array that's exactly large 
		// enough to hold them
		Neuron[] actualNeighbors = new Neuron[numNeigh];
		for (int i = 0 ; i < numNeigh ; ++i) {
			actualNeighbors[i] = possibleNeighbors[i];
		}

		return actualNeighbors;


	}


	// returns a list of neurons in the neighborhood
	public Neuron[] getMooreNeighbors (Neuron neuron, PSO.SelfModel selfModel, PSO.BoundaryModel boundaryModel, 
			PSO.ActivityModel activityModel) {

		// create an array that can hold the maximum number of neurons in the neighborhood, including the neuron itself
		Neuron[] possibleNeighbors = new Neuron[9];
		int numNeigh = 0;

		// this will also consider the neuron itself;
		// NOTE: we could potentially check here if the currentSelfModel is INCLUDE_SELF, but we would also have to 
		// check the current activity model; isNeighbor does both, so call isNeighbor 
		for (int rDelta = -1 ; rDelta <= 1 ; ++rDelta) {
			for (int cDelta = -1 ; cDelta <= 1 ; ++cDelta) {
				Neuron possibleNeighbor = getNeighbor(neuron, rDelta, cDelta, selfModel, boundaryModel, activityModel);
				if (possibleNeighbor != null)
					possibleNeighbors[numNeigh++] = possibleNeighbor;
			}
		}

		// if all the possible neurons are in the neighborhood, can just return the array 
		// they are already in
		if (numNeigh == 9) {   // 9 includes the neuron whose neighbors we are determining
			return possibleNeighbors;
		}

		// otherwise transfer the members of the neighborhood to an array that's exactly large 
		// enough to hold them
		Neuron[] actualNeighbors = new Neuron[numNeigh];
		for (int i = 0 ; i < numNeigh ; ++i) {
			actualNeighbors[i] = possibleNeighbors[i];
		}

		return actualNeighbors;

	}

	// returns a list of neurons in the neighborhood
	public Neuron[] getMooreNeighbors (Neuron neuron, PSOPercolation.SelfModel selfModel, PSOPercolation.BoundaryModel boundaryModel, 
			PSOPercolation.ActivityModel activityModel) {

		// create an array that can hold the maximum number of neurons in the neighborhood, including the neuron itself
		Neuron[] possibleNeighbors = new Neuron[9];
		int numNeigh = 0;

		// this will also consider the neuron itself;
		// NOTE: we could potentially check here if the currentSelfModel is INCLUDE_SELF, but we would also have to 
		// check the current activity model; isNeighbor does both, so call isNeighbor 
		for (int rDelta = -1 ; rDelta <= 1 ; ++rDelta) {
			for (int cDelta = -1 ; cDelta <= 1 ; ++cDelta) {
				Neuron possibleNeighbor = getNeighbor(neuron, rDelta, cDelta, selfModel, boundaryModel, activityModel);
				if (possibleNeighbor != null)
					possibleNeighbors[numNeigh++] = possibleNeighbor;
			}
		}

		// if all the possible neurons are in the neighborhood, can just return the array 
		// they are already in
		if (numNeigh == 9) {   // 9 includes the neuron whose neighbors we are determining
			return possibleNeighbors;
		}

		// otherwise transfer the members of the neighborhood to an array that's exactly large 
		// enough to hold them
		Neuron[] actualNeighbors = new Neuron[numNeigh];
		for (int i = 0 ; i < numNeigh ; ++i) {
			actualNeighbors[i] = possibleNeighbors[i];
		}

		return actualNeighbors;

	}

	
	// checks to see if there is a neuron at (row + rDelta, col + cDelta),wrapping around if
	// the boundary model is a torus, and, if there is, returns that neuron if it is consistent
	// with the self model and the activity model
	public Neuron getNeighbor (Neuron neuron, int rDelta, int cDelta, PSO.SelfModel selfModel, 
			PSO.BoundaryModel boundaryModel, PSO.ActivityModel activityModel) {

		int row = neuron.getRow();
		int col = neuron.getCol();

		Neuron neighbor = null;

		if (selfModel == PSO.SelfModel.INCLUDE_SELF) {

			if (boundaryModel == PSO.BoundaryModel.TORUS) {

				int neighRow = rowWrap(row + rDelta);
				int neighCol = colWrap(col + cDelta);

				if (activityModel == PSO.ActivityModel.ONLY_ACTIVE_NEURONS) {
					if ((grid[neighRow][neighCol] != null) && grid[neighRow][neighCol].active())
						neighbor = grid[neighRow][neighCol];
				}

				else if (activityModel == PSO.ActivityModel.ALL_NEURONS) {				
					if (grid[neighRow][neighCol] != null)
						neighbor = grid[neighRow][neighCol];
				}
				
				else if (activityModel == PSO.ActivityModel.SIMILAR_ACTIVITY_STATUS_NEURONS) {
					if (grid[neighRow][neighCol] != null) {
						boolean sameActivityStatus = (neuron.active() && grid[neighRow][neighCol].active()) || 
						  							(!neuron.active() && !grid[neighRow][neighCol].active());
						if (sameActivityStatus)
							neighbor = grid[neighRow][neighCol];
					}
				}

			}				

			else if (boundaryModel == PSO.BoundaryModel.LATTICE) {

				int neighRow = row + rDelta;
				int neighCol = col + cDelta;

				if (activityModel == PSO.ActivityModel.ONLY_ACTIVE_NEURONS) {
					if (legalCell(neighRow, neighCol) && (grid[neighRow][neighCol] != null) && grid[neighRow][neighCol].active())
						neighbor = grid[neighRow][neighCol];
				}

				else if (activityModel == PSO.ActivityModel.ALL_NEURONS) {				
					if (legalCell(neighRow, neighCol) && (grid[neighRow][neighCol] != null))
						neighbor = grid[neighRow][neighCol];
				}
				
				else if (activityModel == PSO.ActivityModel.SIMILAR_ACTIVITY_STATUS_NEURONS) {
					if (legalCell(neighRow, neighCol) && grid[neighRow][neighCol] != null) {
						boolean sameActivityStatus = (neuron.active() && grid[neighRow][neighCol].active()) || 
						  							(!neuron.active() && !grid[neighRow][neighCol].active());
						if (sameActivityStatus)
							neighbor = grid[neighRow][neighCol];
					}
				}

			}

		}

		else if (selfModel == PSO.SelfModel.NOT_INCLUDE_SELF) {

			if (boundaryModel == PSO.BoundaryModel.TORUS) {

				int neighRow = rowWrap(row + rDelta);
				int neighCol = colWrap(col + cDelta);

				if (activityModel == PSO.ActivityModel.ONLY_ACTIVE_NEURONS) {
					if ((neighRow != row || neighCol != col) && (grid[neighRow][neighCol] != null) && grid[neighRow][neighCol].active())
						neighbor = grid[neighRow][neighCol];
				}

				else if (activityModel == PSO.ActivityModel.ALL_NEURONS) {				
					if ((neighRow != row || neighCol != col) && (grid[neighRow][neighCol] != null))
						neighbor = grid[neighRow][neighCol];
				}
				
				else if (activityModel == PSO.ActivityModel.SIMILAR_ACTIVITY_STATUS_NEURONS) {
					if ((neighRow != row || neighCol != col) && grid[neighRow][neighCol] != null) {
						boolean sameActivityStatus = (neuron.active() && grid[neighRow][neighCol].active()) || 
						  							(!neuron.active() && !grid[neighRow][neighCol].active());
						if (sameActivityStatus)
							neighbor = grid[neighRow][neighCol];
					}
				}

			}
			else if (boundaryModel == PSO.BoundaryModel.LATTICE) {

				int neighRow = row + rDelta;
				int neighCol = col + cDelta;

				if (activityModel == PSO.ActivityModel.ONLY_ACTIVE_NEURONS) {
					if ((neighRow != row || neighCol != col) && legalCell(neighRow, neighCol) && (grid[neighRow][neighCol] != null) && grid[neighRow][neighCol].active())
						neighbor = grid[neighRow][neighCol];
				}

				else if (activityModel == PSO.ActivityModel.ALL_NEURONS) {				
					if ((neighRow != row || neighCol != col) && legalCell(neighRow, neighCol) && (grid[neighRow][neighCol] != null))
						neighbor = grid[neighRow][neighCol];
				}
				
				else if (activityModel == PSO.ActivityModel.SIMILAR_ACTIVITY_STATUS_NEURONS) {
					if ((neighRow != row || neighCol != col) && legalCell(neighRow, neighCol) && (grid[neighRow][neighCol] != null)) {
						boolean sameActivityStatus = (neuron.active() && grid[neighRow][neighCol].active()) || 
						  							(!neuron.active() && !grid[neighRow][neighCol].active());
						if (sameActivityStatus)
							neighbor = grid[neighRow][neighCol];
					}
				}

			}

		}

		// could be null
		return neighbor;

	}

	// checks to see if there is a neuron at (row + rDelta, col + cDelta),wrapping around if
	// the boundary model is a torus, and, if there is, returns that neuron if it is consistent
	// with the self model and the activity model
	public Neuron getNeighbor (Neuron neuron, int rDelta, int cDelta, PSOPercolation.SelfModel selfModel, 
			PSOPercolation.BoundaryModel boundaryModel, PSOPercolation.ActivityModel activityModel) {

		int row = neuron.getRow();
		int col = neuron.getCol();

		Neuron neighbor = null;

		if (selfModel == PSOPercolation.SelfModel.INCLUDE_SELF) {
				

			if (boundaryModel == PSOPercolation.BoundaryModel.LATTICE) {

				int neighRow = row + rDelta;
				int neighCol = col + cDelta;

				if (activityModel == PSOPercolation.ActivityModel.ONLY_ACTIVE_NEURONS) {
					if (legalCell(neighRow, neighCol) && (grid[neighRow][neighCol] != null) && grid[neighRow][neighCol].active())
						neighbor = grid[neighRow][neighCol];
				}

				else if (activityModel == PSOPercolation.ActivityModel.ALL_NEURONS) {				
					if (legalCell(neighRow, neighCol) && (grid[neighRow][neighCol] != null))
						neighbor = grid[neighRow][neighCol];
				}
				
				else if (activityModel == PSOPercolation.ActivityModel.SIMILAR_ACTIVITY_STATUS_NEURONS) {
					if (legalCell(neighRow, neighCol) && grid[neighRow][neighCol] != null) {
						boolean sameActivityStatus = (neuron.active() && grid[neighRow][neighCol].active()) || 
						  							(!neuron.active() && !grid[neighRow][neighCol].active());
						if (sameActivityStatus)
							neighbor = grid[neighRow][neighCol];
					}
				}
			}

		}

		// could be null
		return neighbor;

	}



	
	// n1's activation is being calculated; n2 is the neighbor
	public double getJValueStandard(Neuron n1, Neuron n2) {

		double JValue = 0.0;

		boolean isActive1 = n1.active();
		boolean isActive2 = n2.active();

		if (isActive1 && isActive2) {
			JValue = J[0][0];
		}

		else if (isActive1 && !isActive2) {
			JValue = J[0][1];
		}

		else if (!isActive1 && isActive2) {
			JValue = J[1][0];
		}

		else if (!isActive1 && !isActive2) {
			JValue = J[1][1];
		}

		return JValue;

	}

	
	// n1's activation is being calculated; n2 is the neighbor
	public double getJValuePSOpBestRelativeValues(Neuron n1, Neuron n2) {

		double JValue = 0.0;

		Particle p1 = n1.getParticle();
		Particle p2 = n2.getParticle();
		
		// if n1's particle has a better pbest than its neighbor's particle,
		// increase n1's activation sum by a greater amount
		if (p1.getPersonalBest().getFunctionValue() < p2.getPersonalBest().getFunctionValue()) {
			JValue = 1.5;
		}

		// if n1's neighbor's particle has a better pbest than n1's particle,
		// increase n1's activation level by a lesser amount
		else if (p1.getPersonalBest().getFunctionValue() > p2.getPersonalBest().getFunctionValue()) {
			JValue = 0.5;
		}

		// if both n1's particle and its neighbor's particle have the same pbest value,
		// increase n1's activation level by the activation of its neighbor
		else {
			JValue = 1.0;
		}

		return JValue;

	}

	// n1's activation is being calculated; n2 is the neighbor
	public double getJValuePSOpBestIproved(Neuron n1, Neuron n2) {

		double JValue = 0.0;

		Particle p1 = n1.getParticle();
		Particle p2 = n2.getParticle();
		
		// if n1's particle has a better pbest than its neighbor's particle,
		// increase n1's activation sum by a greater amount
		if (p1.pBestImproved()) {
			JValue = 2.0;
		}
		else {
			JValue = 1.0;			
		}
		
		return JValue;

	}

	

	// move the neurons:
	// not in parallel, instead sweep through the grid, moving them sequentially;
	// this means that a neuron that could have moved before the sweep started,
	// might not be able to move when its turn comes
	// NOTE: movement will always be in a Moore neighborhood at this point;
	//   movement in a vonNeumann neighborhood is also a possibility, but what
	//   I'm doing depends strongly on the dynamics/behavior of the FNN model
	//   in the Solé & Miramontes paper, which uses a Moore neighborhood
	public void moveAllMoore(PSO.BoundaryModel currentFNNBoundaryModel) {

		for (int i = 0 ; i < neuronList.length ; ++i) {
			Neuron n = neuronList[i];
			if (n.active()) {
				moveMoore(n, currentFNNBoundaryModel);
			}
		}


		//		int numMoved = 0;
		//		for (int i = 0 ; i < neuronList.length ; ++i) {
		//			Neuron n = neuronList[i];
		//			if (n.active()) {
		//				if (moveMoore(n, currentBoundaryModel))
		//					++numMoved;
		//			}
		//		}

	}



	// move the neuron at [row,col] to a random location in the Moore neighborborhood;
	// if there is an empty neighbor cell, a move is *always* made;
	// this is in contrast to picking a neighbor cell randomly and then moving there
	// if it is not occupied (and if it is occupied, no movement takes place)
	// returns true if move was made
	public boolean moveMoore(Neuron neuron, PSO.BoundaryModel currentFNNBoundaryModel) {

		int row = neuron.getRow();
		int col = neuron.getCol();

		if (noMovePossible(row, col, currentFNNBoundaryModel)) {
			//			System.out.println("no move possible");
			return false;
		}

		// get random changes in row and column (-1, 0, or +1 in both cases)
		int newRow = 0;
		int newCol = 0;


		if (currentFNNBoundaryModel == PSO.BoundaryModel.LATTICE) {
			int rowChange = PSO.rand.nextInt(3) - 1;
			int colChange = PSO.rand.nextInt(3) - 1;
			newRow = row + rowChange;
			newCol = col + colChange;
			// must:
			//   1) not be the same location it is in already
			//   2) not an illegal location (no wrap-around)
			//   3) not already occupied
			// we know a move is possible, since we checked that at the beginning
			while ((rowChange == 0 && colChange == 0) || 
					!legalCell(newRow, newCol) || 
					(legalCell(newRow, newCol) && grid[newRow][newCol] != null)) {
				rowChange = PSO.rand.nextInt(3) - 1;
				colChange = PSO.rand.nextInt(3) - 1;
				newRow = row + rowChange;
				newCol = col + colChange;
			}
		}

		else if (currentFNNBoundaryModel == PSO.BoundaryModel.TORUS) {
			int rowChange = PSO.rand.nextInt(3) - 1;
			int colChange = PSO.rand.nextInt(3) - 1;
			newRow = rowWrap(row + rowChange);
			newCol = colWrap(col + colChange);
			// must:
			//   1) not be the same location it is in already
			//   2) not already occupied
			// we know a move is possible, since we checked that at the beginning
			while ((rowChange == 0 && colChange == 0) || grid[newRow][newCol] != null) {
				rowChange = PSO.rand.nextInt(3) - 1;
				colChange = PSO.rand.nextInt(3) - 1;
				newRow = rowWrap(row + rowChange);
				newCol = colWrap(col + colChange);

			}

		}


		// In an earlier paper (Miramontes, Solé, and Goodwin, "Collective behavior of random-activated mobile cellular automata," 1993)
		//  a neuron is allowed 6 attempts to move, i.e. pick a neighboring space and move to it if unoccupied; the latter paper (1995) I am using
		//  as my implementation guide, does not seem to have this limit, but it is unclear how it works
		//						boolean foundSpace = false;
		//						int newRow = 0;
		//						int newCol = 0;
		//						for (int i = 1 ; i <= 6 && !foundSpace; ++i) {
		//							int rowChange = FNN.rand.nextInt(3) - 1;
		//							int colChange = FNN.rand.nextInt(3) - 1;
		//							newRow = row + rowChange;
		//							newCol = col + colChange;
		//							// must:
		//							//   1) not be the same location it is already in
		//							//   2) not an illegal location (no wrap-around)
		//							//   3) not already occupied
		//							while ((rowChange == 0 && colChange == 0) || 
		//									!legalCell(newRow, newCol)) {
		//								rowChange = FNN.rand.nextInt(3) - 1;
		//								colChange = FNN.rand.nextInt(3) - 1;
		//								newRow = row + rowChange;
		//								newCol = col + colChange;
		//							}
		//				
		//							if (grid[newRow][newCol] == null) {
		//								foundSpace = true;
		//							}
		//						}
		//				
		//						if (!foundSpace)
		//							return false;

		
		// reset the neuron's row and column
		neuron.setRow(newRow);
		neuron.setCol(newCol);
		// mark it as having moved, so we do not move it again if we encounter it
		// in its new position on our sweep through the grid
		//		neuron.setMoved(true);

		// move the neuron to its new location
		grid[newRow][newCol] = neuron;
		// be sure to set the old location to null
		grid[row][col] = null;

		return true;
		//		System.out.println("node moved to [" + newRow + "][" + newCol + "]");

	}



	// check whether given values for row and column are legal
	public boolean legalCell(int row, int col) {

		return row >= 0 && row < numRows && col >= 0 && col < numCols;

	}




	// if one of the 8 Moore neighbor cells is empty, return true; otherwise false.
	// assume that automata cannot wrap around when moving
	public boolean noMovePossible(int row, int col, PSO.BoundaryModel boundaryCondition) {

		int newRow;
		int newCol;

		if (boundaryCondition == PSO.BoundaryModel.LATTICE) {
			for (int rDelta = -1 ; rDelta <= 1 ; ++rDelta) {
				for (int cDelta = -1 ; cDelta <= 1 ; ++cDelta) {
					newRow = row + rDelta;
					newCol = col + cDelta;
					// ignore the neuron at [row][col], since that's the neuron trying to move
					if ((newRow != row || newCol != col) && legalCell(newRow, newCol) && grid[newRow][newCol] == null)
						return false;
				}
			}
		}

		else if (boundaryCondition == PSO.BoundaryModel.TORUS) {
			for (int rDelta = -1 ; rDelta <= 1 ; ++rDelta) {
				for (int cDelta = -1 ; cDelta <= 1 ; ++cDelta) {
					newRow = rowWrap(row + rDelta);
					newCol = colWrap(col + cDelta);
					// ignore the neuron at [row][col], since that's the neuron trying to move
					// since we're using a torus, the new row and col will always be a legal cell
					if ((newRow != row || newCol != col) && grid[newRow][newCol] == null)
						return false;
				}
			}
		}

		return true;

	}


	
	// does neuron sent in have neighbors?
	public boolean hasNeighbors(Neuron neuron, PSO.BoundaryModel currentFNNBoundaryModel) {
		
		int row = neuron.getRow();
		int col = neuron.getCol();
		
		int newRow;
		int newCol;

		if (currentFNNBoundaryModel == PSO.BoundaryModel.LATTICE) {
			for (int rDelta = -1 ; rDelta <= 1 ; ++rDelta) {
				for (int cDelta = -1 ; cDelta <= 1 ; ++cDelta) {
					newRow = row + rDelta;
					newCol = col + cDelta;
					// ignore the neuron at [row][col]
					if ((newRow != row || newCol != col) && legalCell(newRow, newCol) && grid[newRow][newCol] != null)
						return true;
				}
			}
		}

		else if (currentFNNBoundaryModel == PSO.BoundaryModel.TORUS) {
			for (int rDelta = -1 ; rDelta <= 1 ; ++rDelta) {
				for (int cDelta = -1 ; cDelta <= 1 ; ++cDelta) {
					newRow = rowWrap(row + rDelta);
					newCol = colWrap(col + cDelta);
					// ignore the neuron at [row][col], since that's the neuron trying to move
					// since we're using a torus, the new row and col will always be a legal cell
					if ((newRow != row || newCol != col) && grid[newRow][newCol] != null)
						return true;
				}
			}
		}

		return false;

	}

	// does row and col location sent in have neighbors?
	public boolean hasNeighbors(int row, int col, PSO.BoundaryModel currentFNNBoundaryModel) {
		
		int newRow;
		int newCol;

		if (currentFNNBoundaryModel == PSO.BoundaryModel.LATTICE) {
			for (int rDelta = -1 ; rDelta <= 1 ; ++rDelta) {
				for (int cDelta = -1 ; cDelta <= 1 ; ++cDelta) {
					newRow = row + rDelta;
					newCol = col + cDelta;
					// ignore the neuron at [row][col]
					if ((newRow != row || newCol != col) && legalCell(newRow, newCol) && grid[newRow][newCol] != null)
						return true;
				}
			}
		}

		else if (currentFNNBoundaryModel == PSO.BoundaryModel.TORUS) {
			for (int rDelta = -1 ; rDelta <= 1 ; ++rDelta) {
				for (int cDelta = -1 ; cDelta <= 1 ; ++cDelta) {
					newRow = rowWrap(row + rDelta);
					newCol = colWrap(col + cDelta);
					// ignore the neuron at [row][col], since that's the neuron trying to move
					// since we're using a torus, the new row and col will always be a legal cell
					if ((newRow != row || newCol != col) && grid[newRow][newCol] != null)
						return true;
				}
			}
		}

		return false;

	}

	
	
//	// grid full?
	public boolean gridFull() {

		for (int r = 0 ; r < grid.length ; ++r) {
			for (int c = 0 ; c < grid[r].length ; ++c) {
				if (grid[r][c] == null)
					return false;
			}
		}

		return true;

	}


	// get the number of active neurons; note that we are not checking whether
	// the activation level is > 0.0, because a neuron could have a non-zero
	// activation level, but still be inactive (I think -- this is a question
	// I have for Solé)
//	public int numActiveNeurons() {
//
//		int numActiveNeurons = 0;
//		for (int i = 0 ; i < neuronList.length ; ++i) {
//			if (neuronList[i].active())
//				++numActiveNeurons;
//		}
//
//		return numActiveNeurons;
//
//	}


	public Neuron getRandomNeuron () {

		return neuronList[PSO.rand.nextInt(neuronList.length)];

	}



	// row wrap-around
	public int rowWrap (int row) {
		if (row < 0)
			return numRows - 1;

		if (row == numRows) {
			return 0;
		}

		return row;
	}


	// column wrap-around
	public int colWrap (int col) {
		if (col < 0)
			return numCols - 1;

		if (col == numCols) {
			return 0;
		}

		return col;
	}



	// print activation levels
	public void printActivationLevels() {

		for (int r = 0 ; r < grid.length ; ++r) {
			for (int c = 0 ; c < grid[r].length ; ++c) {
				if (grid[r][c] != null) {
					System.out.printf("%8.4f   ", grid[r][c].getActivationLevel());
				}
				else {
					System.out.printf(" _________ ");					
				}
			}
			System.out.println();
		}
		System.out.println();
		System.out.println();
		System.out.println();

	}


	// print active status
	public void printActiveStatus() {

		for (int r = 0 ; r < grid.length ; ++r) {
			for (int c = 0 ; c < grid[r].length ; ++c) {
				if (grid[r][c] != null) {
					if (grid[r][c].active()) {
						System.out.printf("1  ");
					}
					else {
						System.out.printf("0  ");
					}
				}
				else {
					System.out.printf("-  ");					
				}
			}
			System.out.println();
		}
		System.out.println();
		System.out.println();
		System.out.println();

	}


	// print active status
	public void printActiveStatusWithIDs() {

		for (int r = 0 ; r < grid.length ; ++r) {
			for (int c = 0 ; c < grid[r].length ; ++c) {
				Neuron n = grid[r][c];
				if (n != null) {
					if (n.active()) {
						System.out.printf(" %2d ", n.getID());
					}
					else {
						System.out.printf(" %2d ", 0);
					}
				}
				else {
					System.out.printf(" -- ");
				}
			}
			System.out.println();
		}
		System.out.println();
		System.out.println();
		System.out.println();

	}


	public void printNeurons () {

		for(int n = 0 ; n < neuronList.length ; ++n) {
			System.out.println("neuron " + n + " = " + neuronList[n]);
		}

	}


	public void printNeuronsAndParticles () {

		for(int n = 0 ; n < neuronList.length ; ++n) {
			System.out.println("neuron " + n + " = " + neuronList[n] + " and has particle " + neuronList[n].getParticle());
		}

	}



	// getters and setters

	public Neuron getNeuron (int neuronID) {
		return neuronList[neuronID];
	}


	public double getGain() {
		return gain;
	}

	public void setGain(double gain) {
		this.gain = gain;
	}


	public double getSumNeighborActivationsThreshold() {
		return sumNeighborActivationsThreshold;
	}

	public void setSumNeighborActivationsThreshold(double sumNeighborActivationsThreshold) {
		this.sumNeighborActivationsThreshold = sumNeighborActivationsThreshold;
	}


	public double getActivationThreshold() {
		return activationThreshold;
	}

	public void setActivationThreshold(double activationThreshold) {
		this.activationThreshold = activationThreshold;
	}


	public double getSpontaneousActivationLevel() {
		return spontaneousActivationLevel;
	}

	public void setSpontaneousActivationLevel(double spontaneousActivationLevel) {
		this.spontaneousActivationLevel = spontaneousActivationLevel;
	}


	public double getSpontaneousActivationProbability() {
		return spontaneousActivationProbability;
	}

	public void setSpontaneousActivationProbability(double spontaneousActivationProbability) {
		this.spontaneousActivationProbability = spontaneousActivationProbability;
	}





}




