
public class Solution {

	
	private DoubleVector position;
	private double value;
	private double error;
	private int iterationCreated;
	private int particleID;
	private boolean approximated = false;
	
	
	
	public Solution () {
		this.position = null;
		this.value = 0.0;
		this.error = 0.0;
		this.iterationCreated = 0;
		this.particleID = 0;
		this.approximated = false;
	}

	public Solution (DoubleVector position, double value, double error, int iterationNumCreated, int particleID, boolean approximated) {
		this.position = position.getCopy();
		this.value = value;
		this.error = error;
		this.iterationCreated = iterationNumCreated;
		this.particleID = particleID;			
		this.approximated = approximated;
	}
	
	
	public Solution getCopy() {
		return  new Solution(position, value, error, iterationCreated, particleID, approximated);
	}

	
	public void copyFrom(Solution s) {
		this.position = s.getPositionCopy();
		this.value = s.getFunctionValue();
		this.error = s.getError();
		this.iterationCreated = s.getIterationCreated();
		this.particleID = s.getParticleID();
		this.approximated = s.getApproximated();
	}
	
	
	public double distance (Solution s) {
		return position.distance(s.getPosition());
		
	}
	
//	public void averageIn(Solution sol, int functionNum) {
//		position.averageIn(sol.getPosition());
//		double[] results = TestFunctions.evalWithError(position, functionNum);
//		value = results[TestFunctions.RESULTS_VAL_INDEX];
//		error = results[TestFunctions.RESULTS_ERR_INDEX];
////		value = Double.MAX_VALUE;
////		error = Double.MAX_VALUE;
////		iterationNumCreated = PSO.iterationNum;
//		particleID = -1;
//	}
	
	public int getParticleID() {
		return particleID;
	}

	public void setParticleID(int particleID) {
		this.particleID = particleID;
	}


	// return the actual position
	public DoubleVector getPosition() {
		return position;
	}

	// return a copy of the position
	public DoubleVector getPositionCopy() {
		return position.getCopy();
	}

	// copy from a given position to this position
	public void copyFromPosition(DoubleVector inPosition) {
		this.position.copyFrom(inPosition);
	}

	
	public double getFunctionValue() {
		return value;
	}

	public void setFunctionValue(double value) {
		this.value = value;
	}

	
	public double getError() {
		return error;
	}

	public double getAbsValError() {
		return Math.abs(error);
	}

	public void setError(double error) {
		this.error = error;
	}

	public int getIterationCreated() {
		return iterationCreated;
	}
	
	public void setIterationCreated(int iterationCreated) {
		this.iterationCreated = iterationCreated;
	}

	
	public boolean getApproximated() {
		return approximated;
	}

	public void setApproximated(boolean approximated) {
		this.approximated = approximated;
	}

//	
//	public double calcSolutionValue (int functionNum) {
//		double[] results = TestFunctions.evalWithError(position , functionNum);
//		return results[TestFunctions.RESULTS_VAL_INDEX];
//	}
//
//	public double calcSolutionError (int functionNum) {
//		double[] results = TestFunctions.evalWithError(position , functionNum);
//		return results[TestFunctions.RESULTS_ERR_INDEX];
//	}


	
	
	
	public void print() {
		position.print();
		System.out.printf("%s%2d%s%.8e%s%.8e%s%d", "  particleID = ", particleID, "  val = ", value,"  err = ", error, "  iter = ", iterationCreated);
	}

	public void println() {
		print();
		System.out.println();
	}



//	public void printExp() {
//		position.printExp();
//		System.out.printf("%s%2d%s%.8e%s%.8e%s%d", "  particleID = ", particleID, "  val = ", value,"  err = ", error, "  iter = ", iterationCreated);
//	}
//
//	public void printlnExp() {
//		printExp();
//		System.out.println();
//	}






}



//
//public class Solution {
//
//	
//	private DoubleVector position;
//	private double value;
//	private double error;
//	private int particleID;
//
//	
//	public Solution () {
//		this.position = null;
//		this.value = 0.0;
//		this.error = 0.0;
//		this.particleID = 0;
//	}
//
//	public Solution (DoubleVector position, double value, double error, int particleID) {
//		this.position = position.getCopy();
//		this.value = value;
//		this.error = error;
//		this.particleID = particleID;	
//	}
//	
//	public Solution getCopy() {
//		return new Solution(position.getCopy(), value, error, particleID);
//	}
//
//	
//	public void copyFrom(Solution s) {
//		position = s.getCopyPositionVector();
//		value = s.getValue();
//		error = s.getError();
//		particleID = s.getParticleID();
//	}
//	
//	
//	public int getParticleID() {
//		return particleID;
//	}
//
//	public void setParticleID(int particleID) {
//		this.particleID = particleID;
//	}
//
//
//	public DoubleVector getPositionVector() {
//		return position;
//	}
//
//	public DoubleVector getCopyPositionVector() {
//		return position.getCopy();
//	}
//
//	public void copyToPositionVector(DoubleVector inPosition) {
//		this.position.copyFrom(inPosition);
//	}
//
//	
//	public double getValue() {
//		return value;
//	}
//
//	public void setValue(double value) {
//		this.value = value;
//	}
//
//	
//	public double getError() {
//		return error;
//	}
//
//	public double getAbsValError() {
//		return Math.abs(error);
//	}
//
//	public void setError(double error) {
//		this.error = error;
//	}
//
//
//	public double evalSolutionValue (int functionNum) {
//		double[] results = TestFunctions.evalWithError(position , functionNum);
//		return results[TestFunctions.RESULTS_VAL_INDEX];
//	}
//
//	public double evalSolutionError (int functionNum) {
//		double[] results = TestFunctions.evalWithError(position , functionNum);
//		return results[TestFunctions.RESULTS_ERR_INDEX];
//	}
//
//
//	public void println() {
//		position.print();
//		System.out.printf("%s%2d%s%.8e%s%.8e", "  particleID = ", particleID, "  val = ", value,"  err = ", error);
//		System.out.println();
//	}
//
//
//	public void print() {
//		position.print();
//		System.out.printf("%s%2d%s%.8e%s%.8e", "  particleID = ", particleID, "  val = ", value,"  err = ", error);
//	}
//
//
//
//
//}
