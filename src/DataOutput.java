
public class DataOutput {

	// holds the function value and error for a single run
	
	private double functionValue;
	private double error;
	private double absValError;

	
	public DataOutput () {
		this.functionValue = 0.0;
		this.error = 0.0;
		this.absValError = 0.0;
	}

	public DataOutput (double value, double error, double absValError) {
		this.functionValue = value;
		this.error = error;
		this.absValError = absValError;
	}
	
	public DataOutput getCopy() {
		return new DataOutput(functionValue, error, absValError);
	}

	
	public void copyFrom (DataOutput fromDataOutput) {
		this.functionValue = fromDataOutput.getFunctionValue();
		this.error = fromDataOutput.getError();
		this.absValError = fromDataOutput.getAbsValError();
	}

		
	public void copyDataFromSolution(Solution s) {
		functionValue = s.getFunctionValue();
		error = s.getError();
		absValError = Math.abs(error);

	}


	
	public double getFunctionValue() {
		return functionValue;
	}

	public void setFunctionValue(double functionValue) {
		this.functionValue = functionValue;
	}

	
	public double getError() {
		return error;
	}

	public void setError(double error) {
		this.error = error;
	}

	
	public double getAbsValError() {
		return absValError;
	}

	public void setAbsValError(double absValError) {
		this.absValError = absValError;
	}

	

}
