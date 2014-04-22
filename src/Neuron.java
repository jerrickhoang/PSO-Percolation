/**
 * The objects that are in a Fluid Neural Network
 */

/**
 * @author Stephen Majercik
 * 12/03/13
 *
 */


public class Neuron {

	// neuron has ID for tracking purposes
	private int ID;

	// an node has an activation level, but it is "active" only if the activation 
	// level exceeds a threshold  or if it becomes active "spontaneously" (see the 
	// setActivationLevel method) 
	private double activationLevel;
	private boolean active;

	// nodes live on a grid, so each one has a row and column
	private int row;
	private int col;



//	private double[] activationLevelHistory = new double[PSO.numIterationsData];
//	private int[] activeInactiveHistory = new int[PSO.numIterationsData];
	private int iteration;

	private Particle particle;
	
	// need this to be able to check if the neuron has neighbors
	private FluidNN fnn;
	

	private double gain;
	
	
	// basic constructor
	public Neuron() {

		iteration = 0;

		activationLevel = 0.0;
//		activationLevelHistory[0] = 0.0;
//		activeInactiveHistory[0] = 0;
		active = false;

		row = 0;
		col = 0;

	}


	// create at given location
	public Neuron(int ID, int row, int col) {

		this.ID = ID;
		iteration = 0;

		activationLevel = 0.0;
		//		activationLevelHistory[0] = 0.0;
		//		activeInactiveHistory[0] = 0;
		active = false;

		this.row = row;
		this.col = col;

	}


	// create at given location with given activation level;
	// mark it as active right away if it exceeds the activation
	// threshold
	public Neuron(int ID, int row, int col, double activationLevel, double currentGain, FluidNN fnn) {

		this.ID = ID;
		iteration = 0;

		this.fnn = fnn;

		this.activationLevel = activationLevel;
		//		activationLevelHistory[0] = activationLevel;
		if (activationLevel > fnn.getActivationThreshold()) {
			active = true;
			//			activeInactiveHistory[0] = 1;
		}
		else {
			active = false;
			//			activeInactiveHistory[0] = 0;			
		}		

		this.row = row;
		this.col = col;

		
		gain = currentGain;
		
	}


	// update the activationLevel:
	// set the new activation level
	// if it exceeds the activation level threshold, it becomes active
	//    else, it can become active spontaneously 
	//
	// All of the values that govern this are in FluidNN.java:
	//   	private static double activationThreshold = 1e-16;
	//      private static double spontaneousActivationLevel = 0.1;
	//      private static double spontaneousActivationProbability = 0.1; 
	public void updateActivationStatus(double activationLevel) {


		// set this here, not in the "if" below because if the activationLevel does not
		// exceed the threshold *and* it does not spontaneously activate, it should still
		// have this activation level (used by the nodes that it is a neighbor of)
		this.activationLevel = activationLevel;

		// does it become active? Note that activity does not carry over between iterations
		active = false;

		// exceeds threshold
		if (activationLevel > fnn.getActivationThreshold()) {
			this.activationLevel = activationLevel;
			active = true;
		}
		// else may become active spontaneously
		else if (PSO.rand.nextDouble() < fnn.getSpontaneousActivationProbability()) {
			this.activationLevel = fnn.getSpontaneousActivationLevel();
			active = true;
		}	

		// for simulation of Fernandes work
//		// ************************************************************************************************************************************
//		active = true;
//		// ************************************************************************************************************************************

		// for case where neurons never move after initial placement
//		// ************************************************************************************************************************************
//		active = false;
//		// ************************************************************************************************************************************

	}

	
	// "NEW" SPONTANEOUS ACTIVATION MECHANISM: a particle can spontaneously activate ONLY when it is isolated
	// update the activationLevel:
	// set the new activation level
	// if it exceeds the activation level threshold, it becomes active
	//    else, IF IT IS ISOLATED it can become active spontaneously 
	//
	// All of the values that govern this are in FluidNN.java:
	//   	private static double activationThreshold = 1e-16;
	//      private static double spontaneousActivationLevel = 0.1;
	//      private static double spontaneousActivationProbability = 0.1; 
	public void updateActivationStatusNewSpontActMechanism (double activationLevel, PSO.BoundaryModel currentFNNBoundaryModel) {

		// even if it doesn't become active, its level should be updated
		this.activationLevel = activationLevel;
		
		// does it become active? Note that activity does not carry over between iterations
		active = false;
		
		// does current activation level suffice to  activate it?
		if (activationLevel > fnn.getActivationThreshold()) {
			this.activationLevel = activationLevel;
			active = true;
		}
		
		// if current activation level not high enough to activate it and it's isolated, can spontaneously activate
		else if (!fnn.hasNeighbors(this, currentFNNBoundaryModel)) {
			if (PSO.rand.nextDouble() < fnn.getSpontaneousActivationProbability()) {
				this.activationLevel = fnn.getSpontaneousActivationLevel();
				active = true;
			}	
		}

	}


	// returns probability that neuron n was in a particular state, given the history
//	public static double probabilityOfState(Neuron n, int state) {
//
//		int[] activeInactiveHistory = n.getActiveInactiveHistory();
//
//		double numItersInState = 0.0;
//		for (int i = 0 ; i < activeInactiveHistory.length ; ++i) {
//			if (activeInactiveHistory[i] == state)
//				++numItersInState;
//		}
//		return numItersInState / activeInactiveHistory.length;
//
//
//	}
//
//
//	// returns probability that neurons n1 and n2 were in the specified states, given the history
//	public static double probabilityOfJointStates(Neuron n1, int state1, Neuron n2, int state2) {
//
//		int[] activeInactiveHistory1 = n1.getActiveInactiveHistory();
//		int[] activeInactiveHistory2 = n2.getActiveInactiveHistory();
//
//		double numItersInJointStates = 0.0;
//		for (int i = 0 ; i < activeInactiveHistory1.length ; ++i) {
//			if (activeInactiveHistory1[i] == state1 && activeInactiveHistory2[i] == state2)
//				++numItersInJointStates;
//		}
//		return numItersInJointStates / activeInactiveHistory1.length;		
//
//	}
//
//	public int[] getActiveInactiveHistory() {
//		return activeInactiveHistory;
//	}

	
	// getters and setters

	public double getGain () {
		return gain;
	}

	
	public void setGain (double gain) {
		this.gain = gain;
	}

	
	public int getID() {
		return ID;
	}

	public void setID(int ID) {
		this.ID = ID;
	}

	public double getActivationLevel() {
		return activationLevel;
	}

	public void setActivationLevel(double activationLevel) {
		this.activationLevel = activationLevel;
	}


	public int getRow() {
		return row;
	}

	public void setRow(int row) {
		this.row = row;
	}


	public int getCol() {
		return col;
	}

	public void setCol(int col) {
		this.col = col;
	}


	public boolean active() {
		return active;
	}

	public void setActive(boolean active) {
		this.active = active;
	}


	public Particle getParticle() {
		return particle;
	}


	public void setParticle(Particle particle) {
		this.particle = particle;
	}




}
