
public class TestFunctions {


	// 2-element arrays are frequently used to hold the function value and error
	// after a solution has been evaluated; these are the indices of those
	public static final int VAL_INDEX = 0;
	public static final int ERR_INDEX = 1;


	public static final int NUMBER_OPT_FUNCTIONS = 8;

	// the initialization and search space ranges for all but the two Penalized Functions are taken from:
	//  J.J.Liang, A.K.Qin, P.N.Suganthan, and S.Baskar
	// "Comprehensive learning particle swarm optimizer for global optimization of multimodal functions"
	// IEEE Transactions on Evolutionary Computation, v. 10, no. 3, pp. 281-295, 2006.
	//
	// the initialization and search space ranges for the two Penalized Functions are taken from:
	//  M.E.H.Pederson
	// "Good parameters for particle swarm optimization"
	// Hvass Labs, Tech. Report No. HL1001, 2010.

	// -------------------------------------------------------------------------------------------------------------------------------

	//  Schwefel Problem 2.26
	// 	minimum is 0.0, which occurs at (-420.9687,...,-420.9687)
	public static final int SCHWEFEL_FUNCTION_NUM = 0;

	private static final double SCHWEFEL_INIT_MIN_VAL = -500.0;
	private static final double SCHWEFEL_INIT_MAX_VAL = 500.0;
	private static final double SCHWEFEL_INIT_RANGE = SCHWEFEL_INIT_MAX_VAL - SCHWEFEL_INIT_MIN_VAL;

	private static final double SCHWEFEL_SEARCH_SPACE_MIN_VAL = -500.0;
	private static final double SCHWEFEL_SEARCH_SPACE_MAX_VAL = 500.0;
	private static final double SCHWEFEL_SEARCH_SPACE_RANGE =SCHWEFEL_SEARCH_SPACE_MAX_VAL - SCHWEFEL_SEARCH_SPACE_MIN_VAL;

	private static final double SCHWEFEL_SPEED_MIN_VAL = SCHWEFEL_SEARCH_SPACE_MIN_VAL;
	private static final double SCHWEFEL_SPEED_MAX_VAL = SCHWEFEL_SEARCH_SPACE_MAX_VAL;
	private static final double SCHWEFEL_SPEED_RANGE = SCHWEFEL_SPEED_MAX_VAL - SCHWEFEL_SPEED_MIN_VAL;

	private static double SCHWEFEL_OPT_VALUE = 0.0;
	private static DoubleVector schwefelShiftVector;

	private static final double SCHWEFEL_OPT_COORD = -420.9687;
	public static final double SCHWEFEL_SHIFT_RANGE = -(SCHWEFEL_SEARCH_SPACE_MIN_VAL - SCHWEFEL_OPT_COORD);

	// -------------------------------------------------------------------------------------------------------------------------------

	//  Rastrigin Function
	// 	minimum is 0.0, which occurs at (0.0,...,0.0)
	public static final int RASTRIGIN_FUNCTION_NUM = 1;

	private static final double RASTRIGIN_INIT_MIN_VAL = -5.12;
	private static final double RASTRIGIN_INIT_MAX_VAL = 2.0;
	private static final double RASTRIGIN_INIT_RANGE = RASTRIGIN_INIT_MAX_VAL - RASTRIGIN_INIT_MIN_VAL;

	private static final double RASTRIGIN_SEARCH_SPACE_MIN_VAL = -5.12;
	private static final double RASTRIGIN_SEARCH_SPACE_MAX_VAL = 5.12;
	private static final double RASTRIGIN_SEARCH_SPACE_RANGE = RASTRIGIN_SEARCH_SPACE_MAX_VAL - RASTRIGIN_SEARCH_SPACE_MIN_VAL;

	private static final double RASTRIGIN_SPEED_MIN_VAL = RASTRIGIN_SEARCH_SPACE_MIN_VAL;
	private static final double RASTRIGIN_SPEED_MAX_VAL = RASTRIGIN_SEARCH_SPACE_MAX_VAL;
	private static final double RASTRIGIN_SPEED_RANGE = RASTRIGIN_SPEED_MAX_VAL - RASTRIGIN_SPEED_MIN_VAL;

	private static final double RASTRIGIN_OPT_VALUE = 0.0;
	private static DoubleVector rastriginShiftVector;

	public static final double RASTRIGIN_OPT_COORD = 0.0; 
	public static final double RASTRIGIN_SHIFT_RANGE = RASTRIGIN_SEARCH_SPACE_MAX_VAL - RASTRIGIN_OPT_COORD;

	// -------------------------------------------------------------------------------------------------------------------------------

	//  Ackley Function
	// 	minimum is 0.0, which occurs at (0.0,...,0.0)
	public static final int ACKLEY_FUNCTION_NUM = 2;

	// have also seen (10.0, 20.0) for initialization and (-32.0, 32.0) for search range
	private static final double ACKLEY_INIT_MIN_VAL = -32.768; 
	private static final double ACKLEY_INIT_MAX_VAL = 16.0;   
	//	private static final double ACKLEY_INIT_MIN_VAL = 15.0; 
	//	private static final double ACKLEY_INIT_MAX_VAL = 30.0;   
	private static final double ACKLEY_INIT_RANGE = ACKLEY_INIT_MAX_VAL - ACKLEY_INIT_MIN_VAL;

	private static final double ACKLEY_SEARCH_SPACE_MIN_VAL = -32.768;   
	private static final double ACKLEY_SEARCH_SPACE_MAX_VAL = 32.768;   
	//	private static final double ACKLEY_SEARCH_SPACE_MIN_VAL = -30.0;   
	//	private static final double ACKLEY_SEARCH_SPACE_MAX_VAL = 30.0;   
	private static final double ACKLEY_SEARCH_SPACE_RANGE = ACKLEY_SEARCH_SPACE_MAX_VAL - ACKLEY_SEARCH_SPACE_MIN_VAL;

	private static final double ACKLEY_SPEED_MIN_VAL = ACKLEY_SEARCH_SPACE_MIN_VAL;
	private static final double ACKLEY_SPEED_MAX_VAL = ACKLEY_SEARCH_SPACE_MAX_VAL;
	private static final double ACKLEY_SPEED_RANGE = ACKLEY_SPEED_MAX_VAL - ACKLEY_SPEED_MIN_VAL;

	private static final double ACKLEY_OPT_VALUE = 0.0;
	private static DoubleVector ackleyShiftVector;

	public static final double ACKLEY_OPT_COORD = 0.0; 
	public static final double ACKLEY_SHIFT_RANGE = ACKLEY_SEARCH_SPACE_MAX_VAL - ACKLEY_OPT_COORD;

	// -------------------------------------------------------------------------------------------------------------------------------

	// Griewank function
	// 	minimum is 0.0, which occurs at (0.0,...,0.0)
	public static final int GRIEWANK_FUNCTION_NUM = 3;

	// have also seen (300.0, 600.0) for initialization and (-600.0, 600.0) for search range
	private static final double GRIEWANK_INIT_MIN_VAL = -600.0;
	private static final double GRIEWANK_INIT_MAX_VAL = 200.0;
	private static final double GRIEWANK_INIT_RANGE = GRIEWANK_INIT_MAX_VAL - GRIEWANK_INIT_MIN_VAL;

	private static final double GRIEWANK_SEARCH_SPACE_MIN_VAL = -600.0;
	private static final double GRIEWANK_SEARCH_SPACE_MAX_VAL = 600.0;
	private static final double GRIEWANK_SEARCH_SPACE_RANGE = GRIEWANK_SEARCH_SPACE_MAX_VAL -  GRIEWANK_SEARCH_SPACE_MIN_VAL;

	private static final double GRIEWANK_SPEED_MIN_VAL = GRIEWANK_SEARCH_SPACE_MIN_VAL;
	private static final double GRIEWANK_SPEED_MAX_VAL = GRIEWANK_SEARCH_SPACE_MAX_VAL;
	private static final double GRIEWANK_SPEED_RANGE = GRIEWANK_SPEED_MAX_VAL - GRIEWANK_SPEED_MIN_VAL;

	private static final double GRIEWANK_OPT_VALUE = 0.0;
	private static DoubleVector griewankShiftVector;

	public static final double GRIEWANK_OPT_COORD = 0.0; 
	public static final double GRIEWANK_SHIFT_RANGE = GRIEWANK_SEARCH_SPACE_MAX_VAL - GRIEWANK_OPT_COORD;

	// -------------------------------------------------------------------------------------------------------------------------------

	// Penalized Function 1
	// 	minimum is 0.0, which occurs at (1.0,...,1.0)   // I think this is wrong even though it's what everyone says...
	// 	minimum is 0.0, which occurs at (-1.0,...,-1.0)
	public static final int PENALIZED_FUNCTION_1_NUM = 4;

	private static final double PENALIZED_FUNCTION_1_INIT_MIN_VAL = 5.0;
	private static final double PENALIZED_FUNCTION_1_INIT_MAX_VAL = 50.0;
	private static final double PENALIZED_FUNCTION_1_INIT_RANGE = PENALIZED_FUNCTION_1_INIT_MAX_VAL - PENALIZED_FUNCTION_1_INIT_MIN_VAL;

	private static final double PENALIZED_FUNCTION_1_SEARCH_SPACE_MIN_VAL = -50.0;
	private static final double PENALIZED_FUNCTION_1_SEARCH_SPACE_MAX_VAL = 50.0;
	private static final double PENALIZED_FUNCTION_1_SEARCH_SPACE_RANGE = PENALIZED_FUNCTION_1_SEARCH_SPACE_MAX_VAL -  PENALIZED_FUNCTION_1_SEARCH_SPACE_MIN_VAL;

	private static final double PENALIZED_FUNCTION_1_SPEED_MIN_VAL = PENALIZED_FUNCTION_1_SEARCH_SPACE_MIN_VAL;
	private static final double PENALIZED_FUNCTION_1_SPEED_MAX_VAL = PENALIZED_FUNCTION_1_SEARCH_SPACE_MAX_VAL;
	private static final double PENALIZED_FUNCTION_1_SPEED_RANGE = PENALIZED_FUNCTION_1_SPEED_MAX_VAL - PENALIZED_FUNCTION_1_SPEED_MIN_VAL;

	private static final double PENALIZED_FUNCTION_1_OPT_VALUE = 0.0;
	private static DoubleVector penalizedFunction1ShiftVector;

	public static final double PENALIZED_FUNCTION_1_OPT_COORD = -1.0; 
	public static final double PENALIZED_FUNCTION_1_SHIFT_RANGE = -(PENALIZED_FUNCTION_1_SEARCH_SPACE_MIN_VAL - PENALIZED_FUNCTION_1_OPT_COORD);

	// -------------------------------------------------------------------------------------------------------------------------------

	// Penalized Function 2
	// 	minimum is 0.0, which occurs at (1.0,...,1.0)
	public static final int PENALIZED_FUNCTION_2_NUM = 5;

	private static final double PENALIZED_FUNCTION_2_INIT_MIN_VAL = 5.0;
	private static final double PENALIZED_FUNCTION_2_INIT_MAX_VAL = 50.0;
	private static final double PENALIZED_FUNCTION_2_INIT_RANGE = PENALIZED_FUNCTION_2_INIT_MAX_VAL - PENALIZED_FUNCTION_2_INIT_MIN_VAL;

	private static final double PENALIZED_FUNCTION_2_SEARCH_SPACE_MIN_VAL = -50.0;
	private static final double PENALIZED_FUNCTION_2_SEARCH_SPACE_MAX_VAL = 50.0;
	private static final double PENALIZED_FUNCTION_2_SEARCH_SPACE_RANGE = PENALIZED_FUNCTION_2_SEARCH_SPACE_MAX_VAL -  PENALIZED_FUNCTION_2_SEARCH_SPACE_MIN_VAL;

	private static final double PENALIZED_FUNCTION_2_SPEED_MIN_VAL = PENALIZED_FUNCTION_2_SEARCH_SPACE_MIN_VAL;
	private static final double PENALIZED_FUNCTION_2_SPEED_MAX_VAL = PENALIZED_FUNCTION_2_SEARCH_SPACE_MAX_VAL;
	private static final double PENALIZED_FUNCTION_2_SPEED_RANGE = PENALIZED_FUNCTION_2_SPEED_MAX_VAL - PENALIZED_FUNCTION_2_SPEED_MIN_VAL;

	private static final double PENALIZED_FUNCTION_2_OPT_VALUE = 0.0;
	private static DoubleVector penalizedFunction2ShiftVector;

	public static final double PENALIZED_FUNCTION_2_OPT_COORD = 1.0; 
	public static final double PENALIZED_FUNCTION_2_SHIFT_RANGE = PENALIZED_FUNCTION_2_SEARCH_SPACE_MAX_VAL - PENALIZED_FUNCTION_2_OPT_COORD;

	// -------------------------------------------------------------------------------------------------------------------------------

	// Sphere function
	// 	minimum is 0.0, which occurs at (0.0,...,0.0)
	public static final int SPHERE_FUNCTION_NUM = 6;

	private static final double SPHERE_INIT_MIN_VAL = -100.0;
	private static final double SPHERE_INIT_MAX_VAL = 50.0;
	//	private static final double SPHERE_INIT_MIN_VAL = 50.0;      // from Pedersen paper: used when testing function shift of 25
	//	private static final double SPHERE_INIT_MAX_VAL = 100.0;     // from Pedersen paper: used when testing function shift of 25
	private static final double SPHERE_INIT_RANGE = SPHERE_INIT_MAX_VAL - SPHERE_INIT_MIN_VAL;

	private static final double SPHERE_SEARCH_SPACE_MIN_VAL = -100.0;
	private static final double SPHERE_SEARCH_SPACE_MAX_VAL = 100.0;
	//	private static final double SPHERE_SEARCH_SPACE_MIN_VAL = -100.0;
	//	private static final double SPHERE_SEARCH_SPACE_MAX_VAL = 100.0;
	private static final double SPHERE_SEARCH_SPACE_RANGE = SPHERE_SEARCH_SPACE_MAX_VAL - SPHERE_SEARCH_SPACE_MIN_VAL;

	private static final double SPHERE_SPEED_MIN_VAL = SPHERE_SEARCH_SPACE_MIN_VAL;
	private static final double SPHERE_SPEED_MAX_VAL = SPHERE_SEARCH_SPACE_MAX_VAL;
	private static final double SPHERE_SPEED_RANGE = SPHERE_SPEED_MAX_VAL - SPHERE_SPEED_MIN_VAL;

	private static double SPHERE_OPT_VALUE = 0.0;
	private static DoubleVector sphereShiftVector;

	public static final double SPHERE_OPT_COORD = 0.0; 
	public static final double SPHERE_SHIFT_RANGE = SPHERE_SEARCH_SPACE_MAX_VAL - SPHERE_OPT_COORD;

	// -------------------------------------------------------------------------------------------------------------------------------

	//  Rosenbrock Function
	// 	minimum is 0.0, which occurs at (1.0,...,1.0)
	public static final int ROSENBROCK_FUNCTION_NUM = 7;

	private static final double ROSENBROCK_INIT_MIN_VAL = -2.048;
	private static final double ROSENBROCK_INIT_MAX_VAL = 2.048;
	//	private static final double ROSENBROCK_INIT_MIN_VAL = 15.0;        // from Pedersen paper: used when testing function shift of 25
	//	private static final double ROSENBROCK_INIT_MAX_VAL = 30.0;        // from Pedersen paper: used when testing function shift of 25
	private static final double ROSENBROCK_INIT_RANGE = ROSENBROCK_INIT_MAX_VAL - ROSENBROCK_INIT_MIN_VAL;

	private static final double ROSENBROCK_SEARCH_SPACE_MIN_VAL = -2.048;
	private static final double ROSENBROCK_SEARCH_SPACE_MAX_VAL = 2.048;
	//	private static final double ROSENBROCK_SEARCH_SPACE_MIN_VAL = -100.0;     // from Pedersen paper: used when testing function shift of 25
	//	private static final double ROSENBROCK_SEARCH_SPACE_MAX_VAL = 100.0;      // from Pedersen paper: used when testing function shift of 25
	private static final double ROSENBROCK_SEARCH_SPACE_RANGE =ROSENBROCK_SEARCH_SPACE_MAX_VAL - ROSENBROCK_SEARCH_SPACE_MIN_VAL;

	private static final double ROSENBROCK_SPEED_MIN_VAL = ROSENBROCK_SEARCH_SPACE_MIN_VAL;
	private static final double ROSENBROCK_SPEED_MAX_VAL = ROSENBROCK_SEARCH_SPACE_MAX_VAL;
	private static final double ROSENBROCK_SPEED_RANGE = ROSENBROCK_SPEED_MAX_VAL - ROSENBROCK_SPEED_MIN_VAL;

	private static final double ROSENBROCK_OPT_VALUE = 0.0;
	private static DoubleVector rosenbrockShiftVector;

	public static final double ROSENBROCK_OPT_COORD = 1.0; 
	public static final double ROSENBROCK_SHIFT_RANGE = ROSENBROCK_SEARCH_SPACE_MAX_VAL - ROSENBROCK_OPT_COORD;


	// -------------------------------------------------------------------------------------------------------------------------------



	public static final double[] INIT_MIN_VALS  = 
	{ SCHWEFEL_INIT_MIN_VAL, RASTRIGIN_INIT_MIN_VAL, ACKLEY_INIT_MIN_VAL, GRIEWANK_INIT_MIN_VAL, 
		PENALIZED_FUNCTION_1_INIT_MIN_VAL, PENALIZED_FUNCTION_2_INIT_MIN_VAL, SPHERE_INIT_MIN_VAL, ROSENBROCK_INIT_MIN_VAL };
	public static final double[] INIT_RANGES  = 
	{ SCHWEFEL_INIT_RANGE, RASTRIGIN_INIT_RANGE, ACKLEY_INIT_RANGE, GRIEWANK_INIT_RANGE, 
		PENALIZED_FUNCTION_1_INIT_RANGE, PENALIZED_FUNCTION_2_INIT_RANGE, SPHERE_INIT_RANGE, ROSENBROCK_INIT_RANGE };


	public static final double[] SEARCH_SPACE_MAX_VALS  = 
	{ SCHWEFEL_SEARCH_SPACE_MAX_VAL, RASTRIGIN_SEARCH_SPACE_MAX_VAL, ACKLEY_SEARCH_SPACE_MAX_VAL, GRIEWANK_SEARCH_SPACE_MAX_VAL, 
		PENALIZED_FUNCTION_1_SEARCH_SPACE_MAX_VAL, PENALIZED_FUNCTION_2_SEARCH_SPACE_MAX_VAL, SPHERE_SEARCH_SPACE_MAX_VAL, ROSENBROCK_SEARCH_SPACE_MAX_VAL };

	public static final double[] SEARCH_SPACE_RANGES  = 
	{ SCHWEFEL_SEARCH_SPACE_RANGE, RASTRIGIN_SEARCH_SPACE_RANGE, ACKLEY_SEARCH_SPACE_RANGE, GRIEWANK_SEARCH_SPACE_RANGE, 
		PENALIZED_FUNCTION_1_SEARCH_SPACE_RANGE, PENALIZED_FUNCTION_2_SEARCH_SPACE_RANGE, SPHERE_SEARCH_SPACE_RANGE, ROSENBROCK_SEARCH_SPACE_RANGE };

	
	public static final double UNIVERSAL_MIN_INIT_SPEED = -2.0;
	public static final double UNIVERSAL_SPEED_RANGE = 4.0;

	public static final double[] SPEED_MIN_VALS  = 
	{ SCHWEFEL_SPEED_MIN_VAL, RASTRIGIN_SPEED_MIN_VAL, ACKLEY_SPEED_MIN_VAL, GRIEWANK_SPEED_MIN_VAL, 
		PENALIZED_FUNCTION_1_SPEED_MIN_VAL, PENALIZED_FUNCTION_2_SPEED_MIN_VAL, SPHERE_SPEED_MIN_VAL, ROSENBROCK_SPEED_MIN_VAL };

	public static final double[] SPEED_MAX_VALS  = 
	{ SCHWEFEL_SPEED_MAX_VAL, RASTRIGIN_SPEED_MAX_VAL, ACKLEY_SPEED_MAX_VAL, GRIEWANK_SPEED_MAX_VAL, 
		PENALIZED_FUNCTION_1_SPEED_MAX_VAL, PENALIZED_FUNCTION_2_SPEED_MAX_VAL, SPHERE_SPEED_MAX_VAL, ROSENBROCK_SPEED_MAX_VAL };

	public static final double[] SPEED_RANGES  = 
	{ SCHWEFEL_SPEED_RANGE, RASTRIGIN_SPEED_RANGE, ACKLEY_SPEED_RANGE, GRIEWANK_SPEED_RANGE, 
		PENALIZED_FUNCTION_1_SPEED_RANGE, PENALIZED_FUNCTION_2_SPEED_RANGE, SPHERE_SPEED_RANGE, ROSENBROCK_SPEED_RANGE };

	
	public static final double[] OPT_COORD  = 
	{ SCHWEFEL_OPT_COORD, RASTRIGIN_OPT_COORD, ACKLEY_OPT_COORD, GRIEWANK_OPT_COORD, 
		PENALIZED_FUNCTION_1_OPT_COORD, PENALIZED_FUNCTION_2_OPT_COORD, SPHERE_OPT_COORD, ROSENBROCK_OPT_COORD };

	
	public static final double[] SHIFT_RANGE  = 
	{ SCHWEFEL_SHIFT_RANGE, RASTRIGIN_SHIFT_RANGE, ACKLEY_SHIFT_RANGE, GRIEWANK_SHIFT_RANGE, 
		PENALIZED_FUNCTION_1_SHIFT_RANGE, PENALIZED_FUNCTION_2_SHIFT_RANGE, SPHERE_SHIFT_RANGE, ROSENBROCK_SHIFT_RANGE };

	
	public static final double[] OPT_VALUES  = 
	{ SCHWEFEL_OPT_VALUE, RASTRIGIN_OPT_VALUE, ACKLEY_OPT_VALUE, GRIEWANK_OPT_VALUE, 
		PENALIZED_FUNCTION_1_OPT_VALUE, PENALIZED_FUNCTION_2_OPT_VALUE, SPHERE_OPT_VALUE, ROSENBROCK_OPT_VALUE };




	// -------------------------------------------------------------------------------------------------------------------------------



	public static void setShiftVector (int functionNum, DoubleVector shiftVector) {

		if (functionNum == TestFunctions.SCHWEFEL_FUNCTION_NUM) {
			schwefelShiftVector = shiftVector.getCopy();
		}
		else if (functionNum == TestFunctions.ROSENBROCK_FUNCTION_NUM) {
			rosenbrockShiftVector = shiftVector.getCopy();
		}	
		else if (functionNum == TestFunctions.RASTRIGIN_FUNCTION_NUM) {
			rastriginShiftVector = shiftVector.getCopy();
		}	
		else if (functionNum == TestFunctions.ACKLEY_FUNCTION_NUM) {
			ackleyShiftVector = shiftVector.getCopy();
		}	
		else if (functionNum == TestFunctions.GRIEWANK_FUNCTION_NUM) {
			griewankShiftVector = shiftVector.getCopy();
		}	
		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_1_NUM) {
			penalizedFunction1ShiftVector = shiftVector.getCopy();
		}	
		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_2_NUM) {
			penalizedFunction2ShiftVector = shiftVector.getCopy();
		}
		else if (functionNum == TestFunctions.SPHERE_FUNCTION_NUM) {
			sphereShiftVector = shiftVector.getCopy();
		}

	}



	public static double getOptValue(int functionNum) {
		return OPT_VALUES[functionNum];
	}




	public static double[] evalWithError(DoubleVector v, int functionNum) {

		++PSO.currentFENum;
		//		System.out.println("                in TF  numFunctionEvaluations = " + PSO.numFunctionEvaluations);

		// gets assigned an array with two values:
		//	- the function value
		//	- the error
		double[] retArray = null;  

		if (functionNum == TestFunctions.SCHWEFEL_FUNCTION_NUM) {
			retArray = TestFunctions.evalSchwefelWithError(v);
		}
		else if (functionNum == TestFunctions.ROSENBROCK_FUNCTION_NUM) {
			retArray = TestFunctions.evalRosenbrockWithError(v);
		}	
		else if (functionNum == TestFunctions.RASTRIGIN_FUNCTION_NUM) {
			retArray = TestFunctions.evalRastriginWithError(v);
		}	
		else if (functionNum == TestFunctions.ACKLEY_FUNCTION_NUM) {
			retArray = TestFunctions.evalAckleyWithError(v);
		}	
		else if (functionNum == TestFunctions.GRIEWANK_FUNCTION_NUM) {
			retArray = TestFunctions.evalGriewankWithError(v);
		}	
		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_1_NUM) {
			retArray = TestFunctions.evalPenalizedFunction1WithError(v);
		}	
		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_2_NUM) {
			retArray = TestFunctions.evalPenalizedFunction2WithError(v);
		}
		else if (functionNum == TestFunctions.SPHERE_FUNCTION_NUM) {
			retArray = TestFunctions.evalSphereWithError(v);
		}	

		return retArray;

	}

	// does not increment PSO.numFunctionEvaluations 
	public static double[] evalWithErrorDoNotCountFE (DoubleVector v, int functionNum) {


		double[] retArray = null;  // new double[2];

		if (functionNum == TestFunctions.SCHWEFEL_FUNCTION_NUM) {
			retArray = TestFunctions.evalSchwefelWithError(v);
		}
		else if (functionNum == TestFunctions.ROSENBROCK_FUNCTION_NUM) {
			retArray = TestFunctions.evalRosenbrockWithError(v);
		}	
		else if (functionNum == TestFunctions.RASTRIGIN_FUNCTION_NUM) {
			retArray = TestFunctions.evalRastriginWithError(v);
		}	
		else if (functionNum == TestFunctions.ACKLEY_FUNCTION_NUM) {
			retArray = TestFunctions.evalAckleyWithError(v);
		}	
		else if (functionNum == TestFunctions.GRIEWANK_FUNCTION_NUM) {
			retArray = TestFunctions.evalGriewankWithError(v);
		}	
		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_1_NUM) {
			retArray = TestFunctions.evalPenalizedFunction1WithError(v);
		}	
		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_2_NUM) {
			retArray = TestFunctions.evalPenalizedFunction2WithError(v);
		}
		else if (functionNum == TestFunctions.SPHERE_FUNCTION_NUM) {
			retArray = TestFunctions.evalSphereWithError(v);
		}	

		return retArray;

	}


	

	//  Schwefel Problem 2.26
	// 	minimum is 0.0, which occurs at (-420.9687,...,-420.9687)

	public static double[] evalSchwefelWithError (DoubleVector v) {

		double value = evalSchwefel(v);
		double error = value - SCHWEFEL_OPT_VALUE;  // schwefelShiftedOptValue;
		double[] retArray = new double[2];
		retArray[VAL_INDEX] = value;
		retArray[ERR_INDEX] = error;
		return retArray;

	}


	public static double evalSchwefel (DoubleVector y) {

		DoubleVector v = DoubleVector.sub(y, schwefelShiftVector);

		double sum = 0;
		for(int i = 0 ; i < v.size() ; ++i) {
			double xi = v.get(i);
			sum += xi * Math.sin(Math.sqrt(Math.abs(xi)));
		}


		return 418.9829 * v.size() + sum;  // + schwefelBiasValue;

	}






	//  Rosenbrock Function
	// 	minimum is 0.0, which occurs at (1.0,...,1.0)

	public static double[] evalRosenbrockWithError (DoubleVector v) {

		double value = evalRosenbrock(v);
		double error = value - ROSENBROCK_OPT_VALUE;  //rosenbrockShiftedOptValue;
		double[] retArray = new double[2];
		retArray[VAL_INDEX] = value;
		retArray[ERR_INDEX] = error;
		return retArray;

	}

	public static double evalRosenbrock (DoubleVector y) {

		DoubleVector v = DoubleVector.sub(y, rosenbrockShiftVector);

		double retVal = 0;
		// NOTE:  LAST DIMENSION SHOULD *NOT* BE INCLUDED IN CALCULATION
		for(int i = 0 ; i < v.size() - 1 ; ++i) {
			double xi = v.get(i);
			double xiPlusOne = v.get(i+1);
			retVal += 100.0 * Math.pow(xiPlusOne - xi*xi, 2.0) + Math.pow(xi-1.0, 2.0);
		}

		return retVal; // + rosenbrockBiasValue;

	}






	//  Rastrigin Function
	// 	minimum is 0.0, which occurs at (0.0,...,0.0)

	public static double[] evalRastriginWithError (DoubleVector v) {

		double value = evalRastrigin(v);
		double error = value - RASTRIGIN_OPT_VALUE;   // rastriginShiftedOptValue;
		double[] retArray = new double[2];
		retArray[VAL_INDEX] = value;
		retArray[ERR_INDEX] = error;
		return retArray;

	}

	public static double evalRastrigin (DoubleVector y) {

		DoubleVector v = DoubleVector.sub(y, rastriginShiftVector);

		double retVal = 0;
		for(int i = 0 ; i < v.size() ; ++i) {
			double xi = v.get(i);
			retVal += xi*xi - 10.0*Math.cos(2.0*Math.PI*xi) + 10.0;
		}

		return retVal;  // + rastriginBiasValue;

	}


	



	//  Ackley Function
	// 	minimum is 0.0, which occurs at (0.0,...,0.0)

	public static double[] evalAckleyWithError (DoubleVector v) {

		double value = evalAckley(v);
		double error = value - ACKLEY_OPT_VALUE;  // ackleyShiftedOptValue;
		double[] retArray = new double[2];
		retArray[VAL_INDEX] = value;
		retArray[ERR_INDEX] = error;
		return retArray;

	}

	public static double evalAckley (DoubleVector y) {

		DoubleVector v = DoubleVector.sub(y, ackleyShiftVector);

		double firstSum = 0.0;
		double secondSum = 0.0;
		for(int i = 0 ; i < v.size() ; ++i) {
			double xi = v.get(i);
			firstSum += xi * xi;
			secondSum += Math.cos(2.0*Math.PI*xi);
		}

		return -20.0 * Math.exp(-0.2 * Math.sqrt(firstSum/v.size())) - 
		Math.exp(secondSum/v.size()) + 20.0 + Math.E;   // + ackleyBiasValue;

	}	


	




	// Griewank function
	// 	minimum is 0.0, which occurs at (0.0,...,0.0)

	public static double[] evalGriewankWithError (DoubleVector v) {

		double value = evalGriewank(v);
		double error = value - GRIEWANK_OPT_VALUE;  // griewankShiftedOptValue;
		double[] retArray = new double[2];
		retArray[VAL_INDEX] = value;
		retArray[ERR_INDEX] = error;
		return retArray;

	}

	public static double evalGriewank (DoubleVector y) {

		DoubleVector v = DoubleVector.sub(y, griewankShiftVector);

		double sumSquares = 0.0;
		double productCos = 1.0;
		for(int i = 0 ; i < v.size() ; ++i) {
			double xi = v.get(i);
			sumSquares += xi * xi;
			productCos *= Math.cos(xi/Math.sqrt(i+1));
		}

		return sumSquares/4000.0 - productCos + 1.0;   // + griewankBiasValue;

	}	

	

	



	// Penalized Function 1
	// 	minimum is 0.0, which occurs at (1.0,...,1.0)

	public static double[] evalPenalizedFunction1WithError (DoubleVector v) {

		double value = evalPenalizedFunction1(v);
		double error = value - PENALIZED_FUNCTION_1_OPT_VALUE;   // penalizedFunction1ShiftedOptValue;
		double[] retArray = new double[2];
		retArray[VAL_INDEX] = value;
		retArray[ERR_INDEX] = error;
		return retArray;

	}



	public static double evalPenalizedFunction1 (DoubleVector y) {

		DoubleVector x = DoubleVector.sub(y, penalizedFunction1ShiftVector);

		int xSize = x.size();

		double firstSum = 0.0;
		for(int i = 0 ; i < xSize - 1 ; ++i) {
			double yiAti = yi(x.get(i));
			double sinTerm = Math.sin(Math.PI*yi(x.get(i+1)));
			firstSum += (yiAti-1.0) * (yiAti-1.0) * (1.0 + 10.0*sinTerm*sinTerm);
		}

		double uSum = 0.0;
		for(int i = 0 ; i < xSize ; ++i) {
			uSum += u(x.get(i), 10.0, 100.0, 4.0);
		}

		double firstSinTerm = Math.sin(Math.PI*yi(x.get(0)));
		double yiAtnMinus1 = yi(x.get(xSize-1));

		// some papers have PI/30.0 (first coefficient), but more papers have PI/n  where n = number of dimensions
		return Math.PI/xSize * (10.0*firstSinTerm*firstSinTerm + firstSum + (yiAtnMinus1-1.0)*(yiAtnMinus1-1.0)) + uSum;   // + penalizedFunction1BiasValue;

	}	


	// needed in Penalized Function 1
	private static double yi(double xi) {
		return 1.0 + (xi + 1.0) / 4.0;
	}



	

	// Penalized Function 2
	// 	minimum is 0.0, which occurs at (1.0,...,1.0)

	public static double[] evalPenalizedFunction2WithError (DoubleVector v) {

		double value = evalPenalizedFunction2(v);
		double error = value - PENALIZED_FUNCTION_2_OPT_VALUE;  // penalizedFunction2ShiftedOptValue;
		double[] retArray = new double[2];
		retArray[VAL_INDEX] = value;
		retArray[ERR_INDEX] = error;
		return retArray;

	}



	public static double evalPenalizedFunction2 (DoubleVector y) {

		DoubleVector x = DoubleVector.sub(y, penalizedFunction2ShiftVector);

		int xSize = x.size();

		double firstSum = 0.0;
		for(int i = 0 ; i < xSize - 1 ; ++i) {
			double xi = x.get(i);
			double sinTerm = Math.sin(3.0*Math.PI*x.get(i+1));
			firstSum += (xi-1.0)*(xi-1.0) * (1.0 + sinTerm*sinTerm);
		}

		double uSum = 0.0;
		for(int i = 0 ; i < xSize ; ++i) {
			uSum += u(x.get(i), 5.0, 100.0, 4.0);
		}

		double firstSinTerm = Math.sin(3.0*Math.PI*x.get(0));
		double lastXi = x.get(xSize-1);
		double lastSinTerm = Math.sin(2.0*Math.PI*lastXi);
		return 0.1 * (firstSinTerm*firstSinTerm + firstSum + (lastXi-1.0)*(lastXi-1.0) * (1.0 + lastSinTerm*lastSinTerm)) + uSum ;  // + penalizedFunction2BiasValue;

	}	

	
	

	// needed in Penalized Function 1 and Penalized Function 2
	private static double u(double xi, double a, double k, double m) {

		double retVal = 0;
		if (xi > a) {
			retVal = k * Math.pow(xi-a, m);
		}
		else if (xi <= a && xi >= -a) {
			retVal = 0.0;
		}
		else {  // xi < -a
			retVal = k * Math.pow(-xi-a, m);
		}

		return retVal;

	}



	// Sphere function
	// 	minimum is 0.0, which occurs at (0.0,...,0.0)

	public static double[] evalSphereWithError (DoubleVector v) {

		double value = evalSphere(v);
		double error = value - SPHERE_OPT_VALUE;   // sphereShiftedOptValue;
		double[] retArray = new double[2];
		retArray[VAL_INDEX] = value;
		retArray[ERR_INDEX] = error;
		return retArray;

	}

	public static double evalSphere (DoubleVector y) {

		DoubleVector v = DoubleVector.sub(y, sphereShiftVector);

		double sumSquares = 0.0;
		for(int i = 0 ; i < v.size() ; ++i) {
			double xi = v.get(i);
			sumSquares += xi * xi;
		}

		return sumSquares;   //  + sphereBiasValue;

	}	


	
	


	public static String getFunctionName(int functionNum) {

		String retName = "";

		if (functionNum == TestFunctions.SCHWEFEL_FUNCTION_NUM) {
			retName = "SCHWEFEL_FUNCTION";
		}
		else if (functionNum == TestFunctions.ROSENBROCK_FUNCTION_NUM) {
			retName = "ROSENBROCK_FUNCTION";
		}	
		else if (functionNum == TestFunctions.RASTRIGIN_FUNCTION_NUM) {
			retName = "RASTRIGIN_FUNCTION";
		}	
		else if (functionNum == TestFunctions.ACKLEY_FUNCTION_NUM) {
			retName = "ACKLEY_FUNCTION";
		}	
		else if (functionNum == TestFunctions.GRIEWANK_FUNCTION_NUM) {
			retName = "GRIEWANK_FUNCTION";
		}	
		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_1_NUM) {
			retName = "PENALIZED_FUNCTION_1";
		}	
		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_2_NUM) {
			retName = "PENALIZED_FUNCTION_2";
		}
		else if (functionNum == TestFunctions.SPHERE_FUNCTION_NUM) {
			retName = "SPHERE_FUNCTION";
		}	

		return retName;

	}



	
	


	public static String getShortFunctionName(int functionNum) {

		String retName = "";

		if (functionNum == TestFunctions.SCHWEFEL_FUNCTION_NUM) {
			retName = "SCHWEF";
		}
		else if (functionNum == TestFunctions.ROSENBROCK_FUNCTION_NUM) {
			retName = "ROSE";
		}	
		else if (functionNum == TestFunctions.RASTRIGIN_FUNCTION_NUM) {
			retName = "RAST";
		}	
		else if (functionNum == TestFunctions.ACKLEY_FUNCTION_NUM) {
			retName = "ACK";
		}	
		else if (functionNum == TestFunctions.GRIEWANK_FUNCTION_NUM) {
			retName = "GRIE";
		}	
		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_1_NUM) {
			retName = "PEN1";
		}	
		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_2_NUM) {
			retName = "PEN2";
		}
		else if (functionNum == TestFunctions.SPHERE_FUNCTION_NUM) {
			retName = "SPHR";
		}	

		return retName;

	}



	
	
	public static void printlnFunctionName(int functionNum) {

		printFunctionName(functionNum);
		System.out.println();

	}

	public static void printFunctionName(int functionNum) {

		if (functionNum == TestFunctions.SCHWEFEL_FUNCTION_NUM) {
			System.out.print("SCHWEFEL_FUNCTION");
		}
		else if (functionNum == TestFunctions.ROSENBROCK_FUNCTION_NUM) {
			System.out.print("ROSENBROCK_FUNCTION");
		}	
		else if (functionNum == TestFunctions.RASTRIGIN_FUNCTION_NUM) {
			System.out.print("RASTRIGIN_FUNCTION");
		}	
		else if (functionNum == TestFunctions.ACKLEY_FUNCTION_NUM) {
			System.out.print("ACKLEY_FUNCTION");
		}	
		else if (functionNum == TestFunctions.GRIEWANK_FUNCTION_NUM) {
			System.out.print("GRIEWANK_FUNCTION");
		}	
		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_1_NUM) {
			System.out.print("PENALIZED_FUNCTION_1");
		}	
		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_2_NUM) {
			System.out.print("PENALIZED_FUNCTION_2");
		}
		else if (functionNum == TestFunctions.SPHERE_FUNCTION_NUM) {
			System.out.print("SPHERE_FUNCTION");
		}	

	}


	//	public static double eval (DoubleVector v) {
	//	
	//	double x = v.get(0);
	//	double y = v.get(1);
	//	return x * x + y * y; 
	//	
	//}

	//public static double eval (DoubleVector v) {
	//	
	//	double x = v.get(0);
	//	double y = v.get(1);
	//	return (x - 500) * (x - 500) + (y - 500) * (y - 500); 
	//	
	//}



}





//// this was the code as of 20 October 2012, before I took out all the "biasValue" stuff
//
//       public class TestFunctions {
//
//
//       	public static final int MIN_OPT = 0;
//       	public static final int MAX_OPT = 1;
//
//       	public static final int RESULTS_VAL_INDEX = 0;
//       	public static final int RESULTS_ERR_INDEX = 1;
//
//
//       	public static final int NUMBER_OPT_FUNCTIONS = 8;
//
//       	// the initialization and search space ranges for all but the two Penalized Functions are taken from:
//       	//  J.J.Liang, A.K.Qin, P.N.Suganthan, and S.Baskar
//       	// "Comprehensive learning particle swarm optimizer for global optimization of multimodal functions"
//       	// IEEE Transactions on Evolutionary Computation, v. 10, no. 3, pp. 281-295, 2006.
//       	//
//       	// the initialization and search space ranges for the two Penalized Functions are taken from:
//       	//  M.E.H.Pederson
//       	// "Good parameters for particle swarm optimization"
//       	// Hvass Labs, Tech. Report No. HL1001, 2010.
//       	
//       	// -------------------------------------------------------------------------------------------------------------------------------
//       	
//       	//  Schwefel Problem 2.26
//       	// 	minimum is 0.0, which occurs at (-420.9687,...,-420.9687)
//       	public static final int SCHWEFEL_FUNCTION_NUM = 0;
//       	private static final int SCHWEFEL_FUNCTION_OPT_TYPE = MIN_OPT;
//       	
//       	private static final double SCHWEFEL_INIT_MIN_VAL = -500.0;
//       	private static final double SCHWEFEL_INIT_MAX_VAL = 500.0;
//       	private static final double SCHWEFEL_INIT_RANGE = SCHWEFEL_INIT_MAX_VAL - SCHWEFEL_INIT_MIN_VAL;
//       	
//       	private static final double SCHWEFEL_SEARCH_SPACE_MIN_VAL = -500.0;
//       	private static final double SCHWEFEL_SEARCH_SPACE_MAX_VAL = 500.0;
//       	private static final double SCHWEFEL_SEARCH_SPACE_RANGE =SCHWEFEL_SEARCH_SPACE_MAX_VAL - SCHWEFEL_SEARCH_SPACE_MIN_VAL;
//       	
//       	private static final double SCHWEFEL_SPEED_MIN_VAL = SCHWEFEL_SEARCH_SPACE_MIN_VAL;
//       	private static final double SCHWEFEL_SPEED_MAX_VAL = SCHWEFEL_SEARCH_SPACE_MAX_VAL;
//       	private static final double SCHWEFEL_SPEED_RANGE = SCHWEFEL_SPEED_MAX_VAL - SCHWEFEL_SPEED_MIN_VAL;
//       	
//       	private static double SCHWEFEL_OPT_VALUE = 0.0;
//       	private static double schwefelBiasValue;
//       	private static double schwefelShiftedOptValue;
//       	private static DoubleVector schwefelShiftVector;
//
//       	private static final double SCHWEFEL_OPT_COORD = -420.9687;
//       	public static final double SCHWEFEL_SHIFT_RANGE = -(SCHWEFEL_SEARCH_SPACE_MIN_VAL - SCHWEFEL_OPT_COORD);
//
//       	// -------------------------------------------------------------------------------------------------------------------------------
//       	
//       	//  Rastrigin Function
//       	// 	minimum is 0.0, which occurs at (0.0,...,0.0)
//       	public static final int RASTRIGIN_FUNCTION_NUM = 1;
//       	private static final int RASTRIGIN_FUNCTION_OPT_TYPE = MIN_OPT;
//       	
//          	private static final double RASTRIGIN_INIT_MIN_VAL = -5.12;
//       	private static final double RASTRIGIN_INIT_MAX_VAL = 2.0;
//       	private static final double RASTRIGIN_INIT_RANGE = RASTRIGIN_INIT_MAX_VAL - RASTRIGIN_INIT_MIN_VAL;
//       	
//       	private static final double RASTRIGIN_SEARCH_SPACE_MIN_VAL = -5.12;
//       	private static final double RASTRIGIN_SEARCH_SPACE_MAX_VAL = 5.12;
//       	private static final double RASTRIGIN_SEARCH_SPACE_RANGE = RASTRIGIN_SEARCH_SPACE_MAX_VAL - RASTRIGIN_SEARCH_SPACE_MIN_VAL;
//       	
//       	private static final double RASTRIGIN_SPEED_MIN_VAL = RASTRIGIN_SEARCH_SPACE_MIN_VAL;
//       	private static final double RASTRIGIN_SPEED_MAX_VAL = RASTRIGIN_SEARCH_SPACE_MAX_VAL;
//       	private static final double RASTRIGIN_SPEED_RANGE = RASTRIGIN_SPEED_MAX_VAL - RASTRIGIN_SPEED_MIN_VAL;
//       	
//       	private static final double RASTRIGIN_OPT_VALUE = 0.0;
//       	private static double rastriginBiasValue;
//       	private static double rastriginShiftedOptValue;
//       	private static DoubleVector rastriginShiftVector;
//
//       	public static final double RASTRIGIN_OPT_COORD = 0.0; 
//       	public static final double RASTRIGIN_SHIFT_RANGE = RASTRIGIN_SEARCH_SPACE_MAX_VAL - RASTRIGIN_OPT_COORD;
//       	
//       	// -------------------------------------------------------------------------------------------------------------------------------
//       	
//       	//  Ackley Function
//       	// 	minimum is 0.0, which occurs at (0.0,...,0.0)
//       	public static final int ACKLEY_FUNCTION_NUM = 2;
//       	private static final int ACKLEY_FUNCTION_OPT_TYPE = MIN_OPT;
//       	
//       	// have also seen (10.0, 20.0) for initialization and (-32.0, 32.0) for search range
//       	private static final double ACKLEY_INIT_MIN_VAL = -32.768; 
//       	private static final double ACKLEY_INIT_MAX_VAL = 16.0;   
////       	private static final double ACKLEY_INIT_MIN_VAL = 15.0; 
////       	private static final double ACKLEY_INIT_MAX_VAL = 30.0;   
//       	private static final double ACKLEY_INIT_RANGE = ACKLEY_INIT_MAX_VAL - ACKLEY_INIT_MIN_VAL;
//       	
//       	private static final double ACKLEY_SEARCH_SPACE_MIN_VAL = -32.768;   
//       	private static final double ACKLEY_SEARCH_SPACE_MAX_VAL = 32.768;   
////       	private static final double ACKLEY_SEARCH_SPACE_MIN_VAL = -30.0;   
////       	private static final double ACKLEY_SEARCH_SPACE_MAX_VAL = 30.0;   
//       	private static final double ACKLEY_SEARCH_SPACE_RANGE = ACKLEY_SEARCH_SPACE_MAX_VAL - ACKLEY_SEARCH_SPACE_MIN_VAL;
//       	
//       	private static final double ACKLEY_SPEED_MIN_VAL = ACKLEY_SEARCH_SPACE_MIN_VAL;
//       	private static final double ACKLEY_SPEED_MAX_VAL = ACKLEY_SEARCH_SPACE_MAX_VAL;
//       	private static final double ACKLEY_SPEED_RANGE = ACKLEY_SPEED_MAX_VAL - ACKLEY_SPEED_MIN_VAL;
//       	
//       	private static final double ACKLEY_OPT_VALUE = 0.0;
//       	private static double ackleyBiasValue;
//       	private static double ackleyShiftedOptValue;
//       	private static DoubleVector ackleyShiftVector;
//       	
//       	public static final double ACKLEY_OPT_COORD = 0.0; 
//       	public static final double ACKLEY_SHIFT_RANGE = ACKLEY_SEARCH_SPACE_MAX_VAL - ACKLEY_OPT_COORD;
//
//       	// -------------------------------------------------------------------------------------------------------------------------------
//       	
//       	// Griewank function
//       	// 	minimum is 0.0, which occurs at (0.0,...,0.0)
//       	public static final int GRIEWANK_FUNCTION_NUM = 3;
//       	private static final int GRIEWANK_FUNCTION_OPT_TYPE = MIN_OPT;
//       	
//       	// have also seen (300.0, 600.0) for initialization and (-600.0, 600.0) for search range
//       	private static final double GRIEWANK_INIT_MIN_VAL = -600.0;
//       	private static final double GRIEWANK_INIT_MAX_VAL = 200.0;
//       	private static final double GRIEWANK_INIT_RANGE = GRIEWANK_INIT_MAX_VAL - GRIEWANK_INIT_MIN_VAL;
//       	
//       	private static final double GRIEWANK_SEARCH_SPACE_MIN_VAL = -600.0;
//       	private static final double GRIEWANK_SEARCH_SPACE_MAX_VAL = 600.0;
//       	private static final double GRIEWANK_SEARCH_SPACE_RANGE = GRIEWANK_SEARCH_SPACE_MAX_VAL -  GRIEWANK_SEARCH_SPACE_MIN_VAL;
//       	
//       	private static final double GRIEWANK_SPEED_MIN_VAL = GRIEWANK_SEARCH_SPACE_MIN_VAL;
//       	private static final double GRIEWANK_SPEED_MAX_VAL = GRIEWANK_SEARCH_SPACE_MAX_VAL;
//       	private static final double GRIEWANK_SPEED_RANGE = GRIEWANK_SPEED_MAX_VAL - GRIEWANK_SPEED_MIN_VAL;
//       	
//       	private static final double GRIEWANK_OPT_VALUE = 0.0;
//       	private static double griewankBiasValue;
//       	private static double griewankShiftedOptValue;
//       	private static DoubleVector griewankShiftVector;
//       	
//       	public static final double GRIEWANK_OPT_COORD = 0.0; 
//       	public static final double GRIEWANK_SHIFT_RANGE = GRIEWANK_SEARCH_SPACE_MAX_VAL - GRIEWANK_OPT_COORD;
//
//       	// -------------------------------------------------------------------------------------------------------------------------------
//       	
//       	// Penalized Function 1
//       	// 	minimum is 0.0, which occurs at (1.0,...,1.0)   // I think this is wrong even though it's what everyone says...
//       	// 	minimum is 0.0, which occurs at (-1.0,...,-1.0)
//       	public static final int PENALIZED_FUNCTION_1_NUM = 4;
//       	private static final int PENALIZED_FUNCTION_1_OPT_TYPE = MIN_OPT;
//       	
//       	private static final double PENALIZED_FUNCTION_1_INIT_MIN_VAL = 5.0;
//       	private static final double PENALIZED_FUNCTION_1_INIT_MAX_VAL = 50.0;
//       	private static final double PENALIZED_FUNCTION_1_INIT_RANGE = PENALIZED_FUNCTION_1_INIT_MAX_VAL - PENALIZED_FUNCTION_1_INIT_MIN_VAL;
//       	
//       	private static final double PENALIZED_FUNCTION_1_SEARCH_SPACE_MIN_VAL = -50.0;
//       	private static final double PENALIZED_FUNCTION_1_SEARCH_SPACE_MAX_VAL = 50.0;
//       	private static final double PENALIZED_FUNCTION_1_SEARCH_SPACE_RANGE = PENALIZED_FUNCTION_1_SEARCH_SPACE_MAX_VAL -  PENALIZED_FUNCTION_1_SEARCH_SPACE_MIN_VAL;
//       	
//       	private static final double PENALIZED_FUNCTION_1_SPEED_MIN_VAL = PENALIZED_FUNCTION_1_SEARCH_SPACE_MIN_VAL;
//       	private static final double PENALIZED_FUNCTION_1_SPEED_MAX_VAL = PENALIZED_FUNCTION_1_SEARCH_SPACE_MAX_VAL;
//       	private static final double PENALIZED_FUNCTION_1_SPEED_RANGE = PENALIZED_FUNCTION_1_SPEED_MAX_VAL - PENALIZED_FUNCTION_1_SPEED_MIN_VAL;
//       	
//       	private static final double PENALIZED_FUNCTION_1_OPT_VALUE = 0.0;
//       	private static double penalizedFunction1BiasValue;
//       	private static double penalizedFunction1ShiftedOptValue;
//       	private static DoubleVector penalizedFunction1ShiftVector;
//       	
//       	public static final double PENALIZED_FUNCTION_1_OPT_COORD = -1.0; 
//       	public static final double PENALIZED_FUNCTION_1_SHIFT_RANGE = -(PENALIZED_FUNCTION_1_SEARCH_SPACE_MIN_VAL - PENALIZED_FUNCTION_1_OPT_COORD);
//
//       	// -------------------------------------------------------------------------------------------------------------------------------
//       	
//       	// Penalized Function 2
//       	// 	minimum is 0.0, which occurs at (1.0,...,1.0)
//       	public static final int PENALIZED_FUNCTION_2_NUM = 5;
//       	private static final int PENALIZED_FUNCTION_2_OPT_TYPE = MIN_OPT;
//       	
//       	private static final double PENALIZED_FUNCTION_2_INIT_MIN_VAL = 5.0;
//       	private static final double PENALIZED_FUNCTION_2_INIT_MAX_VAL = 50.0;
//       	private static final double PENALIZED_FUNCTION_2_INIT_RANGE = PENALIZED_FUNCTION_2_INIT_MAX_VAL - PENALIZED_FUNCTION_2_INIT_MIN_VAL;
//       	
//       	private static final double PENALIZED_FUNCTION_2_SEARCH_SPACE_MIN_VAL = -50.0;
//       	private static final double PENALIZED_FUNCTION_2_SEARCH_SPACE_MAX_VAL = 50.0;
//       	private static final double PENALIZED_FUNCTION_2_SEARCH_SPACE_RANGE = PENALIZED_FUNCTION_2_SEARCH_SPACE_MAX_VAL -  PENALIZED_FUNCTION_2_SEARCH_SPACE_MIN_VAL;
//       	
//       	private static final double PENALIZED_FUNCTION_2_SPEED_MIN_VAL = PENALIZED_FUNCTION_2_SEARCH_SPACE_MIN_VAL;
//       	private static final double PENALIZED_FUNCTION_2_SPEED_MAX_VAL = PENALIZED_FUNCTION_2_SEARCH_SPACE_MAX_VAL;
//       	private static final double PENALIZED_FUNCTION_2_SPEED_RANGE = PENALIZED_FUNCTION_2_SPEED_MAX_VAL - PENALIZED_FUNCTION_2_SPEED_MIN_VAL;
//       	
//       	private static final double PENALIZED_FUNCTION_2_OPT_VALUE = 0.0;
//       	private static double penalizedFunction2BiasValue;
//       	private static double penalizedFunction2ShiftedOptValue;
//       	private static DoubleVector penalizedFunction2ShiftVector;
//
//       	public static final double PENALIZED_FUNCTION_2_OPT_COORD = 1.0; 
//       	public static final double PENALIZED_FUNCTION_2_SHIFT_RANGE = PENALIZED_FUNCTION_2_SEARCH_SPACE_MAX_VAL - PENALIZED_FUNCTION_2_OPT_COORD;
//
//       	// -------------------------------------------------------------------------------------------------------------------------------
//       	
//       	// Sphere function
//       	// 	minimum is 0.0, which occurs at (0.0,...,0.0)
//       	public static final int SPHERE_FUNCTION_NUM = 6;
//       	private static final int SPHERE_FUNCTION_OPT_TYPE = MIN_OPT;
//       	
//       	private static final double SPHERE_INIT_MIN_VAL = -100.0;
//       	private static final double SPHERE_INIT_MAX_VAL = 50.0;
////       	private static final double SPHERE_INIT_MIN_VAL = 50.0;      // from Pedersen paper: used when testing function shift of 25
////       	private static final double SPHERE_INIT_MAX_VAL = 100.0;     // from Pedersen paper: used when testing function shift of 25
//       	private static final double SPHERE_INIT_RANGE = SPHERE_INIT_MAX_VAL - SPHERE_INIT_MIN_VAL;
//       	
//       	private static final double SPHERE_SEARCH_SPACE_MIN_VAL = -100.0;
//       	private static final double SPHERE_SEARCH_SPACE_MAX_VAL = 100.0;
////       	private static final double SPHERE_SEARCH_SPACE_MIN_VAL = -100.0;
////       	private static final double SPHERE_SEARCH_SPACE_MAX_VAL = 100.0;
//       	private static final double SPHERE_SEARCH_SPACE_RANGE = SPHERE_SEARCH_SPACE_MAX_VAL - SPHERE_SEARCH_SPACE_MIN_VAL;
//       	
//       	private static final double SPHERE_SPEED_MIN_VAL = SPHERE_SEARCH_SPACE_MIN_VAL;
//       	private static final double SPHERE_SPEED_MAX_VAL = SPHERE_SEARCH_SPACE_MAX_VAL;
//       	private static final double SPHERE_SPEED_RANGE = SPHERE_SPEED_MAX_VAL - SPHERE_SPEED_MIN_VAL;
//       	
//       	private static double SPHERE_OPT_VALUE = 0.0;
//       	private static double sphereBiasValue;
//       	private static double sphereShiftedOptValue;
//       	private static DoubleVector sphereShiftVector;
//
//       	public static final double SPHERE_OPT_COORD = 0.0; 
//       	public static final double SPHERE_SHIFT_RANGE = SPHERE_SEARCH_SPACE_MAX_VAL - SPHERE_OPT_COORD;
//
//       	// -------------------------------------------------------------------------------------------------------------------------------
//       	
//       	//  Rosenbrock Function
//       	// 	minimum is 0.0, which occurs at (1.0,...,1.0)
//       	public static final int ROSENBROCK_FUNCTION_NUM = 7;
//       	private static final int ROSENBROCK_FUNCTION_OPT_TYPE = MIN_OPT;
//       	
//       	private static final double ROSENBROCK_INIT_MIN_VAL = -2.048;
//       	private static final double ROSENBROCK_INIT_MAX_VAL = 2.048;
////       	private static final double ROSENBROCK_INIT_MIN_VAL = 15.0;        // from Pedersen paper: used when testing function shift of 25
////       	private static final double ROSENBROCK_INIT_MAX_VAL = 30.0;        // from Pedersen paper: used when testing function shift of 25
//       	private static final double ROSENBROCK_INIT_RANGE = ROSENBROCK_INIT_MAX_VAL - ROSENBROCK_INIT_MIN_VAL;
//       	
//       	private static final double ROSENBROCK_SEARCH_SPACE_MIN_VAL = -2.048;
//       	private static final double ROSENBROCK_SEARCH_SPACE_MAX_VAL = 2.048;
////       	private static final double ROSENBROCK_SEARCH_SPACE_MIN_VAL = -100.0;     // from Pedersen paper: used when testing function shift of 25
////       	private static final double ROSENBROCK_SEARCH_SPACE_MAX_VAL = 100.0;      // from Pedersen paper: used when testing function shift of 25
//       	private static final double ROSENBROCK_SEARCH_SPACE_RANGE =ROSENBROCK_SEARCH_SPACE_MAX_VAL - ROSENBROCK_SEARCH_SPACE_MIN_VAL;
//       	
//       	private static final double ROSENBROCK_SPEED_MIN_VAL = ROSENBROCK_SEARCH_SPACE_MIN_VAL;
//       	private static final double ROSENBROCK_SPEED_MAX_VAL = ROSENBROCK_SEARCH_SPACE_MAX_VAL;
//       	private static final double ROSENBROCK_SPEED_RANGE = ROSENBROCK_SPEED_MAX_VAL - ROSENBROCK_SPEED_MIN_VAL;
//       	
//       	private static final double ROSENBROCK_OPT_VALUE = 0.0;
//       	private static double rosenbrockBiasValue;
//       	private static double rosenbrockShiftedOptValue;
//       	private static DoubleVector rosenbrockShiftVector;
//       	
//       	public static final double ROSENBROCK_OPT_COORD = 1.0; 
//       	public static final double ROSENBROCK_SHIFT_RANGE = ROSENBROCK_SEARCH_SPACE_MAX_VAL - ROSENBROCK_OPT_COORD;
//
//
//       	// -------------------------------------------------------------------------------------------------------------------------------
//       	
//
//       	private static final int[] OPT_TYPES = 
//       	{ SCHWEFEL_FUNCTION_OPT_TYPE, RASTRIGIN_FUNCTION_OPT_TYPE, ACKLEY_FUNCTION_OPT_TYPE, GRIEWANK_FUNCTION_OPT_TYPE, 
//       		PENALIZED_FUNCTION_1_OPT_TYPE, PENALIZED_FUNCTION_2_OPT_TYPE, SPHERE_FUNCTION_OPT_TYPE, ROSENBROCK_FUNCTION_OPT_TYPE};
//
//       	
//       	public static final double[] INIT_MIN_VALS  = 
//       	{ SCHWEFEL_INIT_MIN_VAL, RASTRIGIN_INIT_MIN_VAL, ACKLEY_INIT_MIN_VAL, GRIEWANK_INIT_MIN_VAL, 
//       		PENALIZED_FUNCTION_1_INIT_MIN_VAL, PENALIZED_FUNCTION_2_INIT_MIN_VAL, SPHERE_INIT_MIN_VAL, ROSENBROCK_INIT_MIN_VAL };
////       	private static final double[] INIT_MAX_VALS  = 
////       	{ SCHWEFEL_INIT_MAX_VAL, RASTRIGIN_INIT_MAX_VAL, ACKLEY_INIT_MAX_VAL, GRIEWANK_INIT_MAX_VAL, 
////       		PENALIZED_FUNCTION_1_INIT_MAX_VAL, PENALIZED_FUNCTION_2_INIT_MAX_VAL, SPHERE_INIT_MAX_VAL, ROSENBROCK_INIT_MAX_VAL };
//       	public static final double[] INIT_RANGES  = 
//       	{ SCHWEFEL_INIT_RANGE, RASTRIGIN_INIT_RANGE, ACKLEY_INIT_RANGE, GRIEWANK_INIT_RANGE, 
//       		PENALIZED_FUNCTION_1_INIT_RANGE, PENALIZED_FUNCTION_2_INIT_RANGE, SPHERE_INIT_RANGE, ROSENBROCK_INIT_RANGE };
//
//       	
////       	private static final double[] SEARCH_SPACE_MIN_VALS  = 
////       	{ SCHWEFEL_SEARCH_SPACE_MIN_VAL, RASTRIGIN_SEARCH_SPACE_MIN_VAL, ACKLEY_SEARCH_SPACE_MIN_VAL, GRIEWANK_SEARCH_SPACE_MIN_VAL, 
////       		PENALIZED_FUNCTION_1_SEARCH_SPACE_MIN_VAL, PENALIZED_FUNCTION_2_SEARCH_SPACE_MIN_VAL, SPHERE_SEARCH_SPACE_MIN_VAL, ROSENBROCK_SEARCH_SPACE_MIN_VAL };
//       	public static final double[] SEARCH_SPACE_MAX_VALS  = 
//       	{ SCHWEFEL_SEARCH_SPACE_MAX_VAL, RASTRIGIN_SEARCH_SPACE_MAX_VAL, ACKLEY_SEARCH_SPACE_MAX_VAL, GRIEWANK_SEARCH_SPACE_MAX_VAL, 
//       		PENALIZED_FUNCTION_1_SEARCH_SPACE_MAX_VAL, PENALIZED_FUNCTION_2_SEARCH_SPACE_MAX_VAL, SPHERE_SEARCH_SPACE_MAX_VAL, ROSENBROCK_SEARCH_SPACE_MAX_VAL };
//       	
//       	public static final double[] SEARCH_SPACE_RANGES  = 
//       	{ SCHWEFEL_SEARCH_SPACE_RANGE, RASTRIGIN_SEARCH_SPACE_RANGE, ACKLEY_SEARCH_SPACE_RANGE, GRIEWANK_SEARCH_SPACE_RANGE, 
//       		PENALIZED_FUNCTION_1_SEARCH_SPACE_RANGE, PENALIZED_FUNCTION_2_SEARCH_SPACE_RANGE, SPHERE_SEARCH_SPACE_RANGE, ROSENBROCK_SEARCH_SPACE_RANGE };
//
//       	
//       	public static final double[] SPEED_MIN_VALS  = 
//       	{ SCHWEFEL_SPEED_MIN_VAL, RASTRIGIN_SPEED_MIN_VAL, ACKLEY_SPEED_MIN_VAL, GRIEWANK_SPEED_MIN_VAL, 
//       		PENALIZED_FUNCTION_1_SPEED_MIN_VAL, PENALIZED_FUNCTION_2_SPEED_MIN_VAL, SPHERE_SPEED_MIN_VAL, ROSENBROCK_SPEED_MIN_VAL };
//       	
//       	public static final double[] SPEED_MAX_VALS  = 
//       	{ SCHWEFEL_SPEED_MAX_VAL, RASTRIGIN_SPEED_MAX_VAL, ACKLEY_SPEED_MAX_VAL, GRIEWANK_SPEED_MAX_VAL, 
//       		PENALIZED_FUNCTION_1_SPEED_MAX_VAL, PENALIZED_FUNCTION_2_SPEED_MAX_VAL, SPHERE_SPEED_MAX_VAL, ROSENBROCK_SPEED_MAX_VAL };
//       	
//       	public static final double[] SPEED_RANGES  = 
//       	{ SCHWEFEL_SPEED_RANGE, RASTRIGIN_SPEED_RANGE, ACKLEY_SPEED_RANGE, GRIEWANK_SPEED_RANGE, 
//       		PENALIZED_FUNCTION_1_SPEED_RANGE, PENALIZED_FUNCTION_2_SPEED_RANGE, SPHERE_SPEED_RANGE, ROSENBROCK_SPEED_RANGE };
//       	
//       	public static final double[] OPT_COORD  = 
//       	{ SCHWEFEL_OPT_COORD, RASTRIGIN_OPT_COORD, ACKLEY_OPT_COORD, GRIEWANK_OPT_COORD, 
//       		PENALIZED_FUNCTION_1_OPT_COORD, PENALIZED_FUNCTION_2_OPT_COORD, SPHERE_OPT_COORD, ROSENBROCK_OPT_COORD };
//
//       	public static final double[] SHIFT_RANGE  = 
//       	{ SCHWEFEL_SHIFT_RANGE, RASTRIGIN_SHIFT_RANGE, ACKLEY_SHIFT_RANGE, GRIEWANK_SHIFT_RANGE, 
//       		PENALIZED_FUNCTION_1_SHIFT_RANGE, PENALIZED_FUNCTION_2_SHIFT_RANGE, SPHERE_SHIFT_RANGE, ROSENBROCK_SHIFT_RANGE };
//
//       	public static final double[] OPT_VALUES  = 
//       	{ SCHWEFEL_OPT_VALUE, RASTRIGIN_OPT_VALUE, ACKLEY_OPT_VALUE, GRIEWANK_OPT_VALUE, 
//       		PENALIZED_FUNCTION_1_OPT_VALUE, PENALIZED_FUNCTION_2_OPT_VALUE, SPHERE_OPT_VALUE, ROSENBROCK_OPT_VALUE };
//
//
//       	
//       	
//       	// -------------------------------------------------------------------------------------------------------------------------------
//
//
//       	
//       	public static void setTestFunction(int functionNum, DoubleVector shiftVector, double inBiasValue, boolean useIncomingBias) {
//
//       		if (functionNum == TestFunctions.SCHWEFEL_FUNCTION_NUM) {
//       			schwefelShiftVector = shiftVector.getCopy();
//       			schwefelBiasValue = useIncomingBias? inBiasValue: evalSchwefelNoShiftNoBias(shiftVector);
//       			schwefelShiftedOptValue = SCHWEFEL_OPT_VALUE + schwefelBiasValue;				
//       		}
//       		else if (functionNum == TestFunctions.ROSENBROCK_FUNCTION_NUM) {
//       			rosenbrockShiftVector = shiftVector.getCopy();
//       			rosenbrockBiasValue = useIncomingBias? inBiasValue: evalRosenbrockNoShiftNoBias(shiftVector);
//       			rosenbrockShiftedOptValue = ROSENBROCK_OPT_VALUE + rosenbrockBiasValue;				
//       		}	
//       		else if (functionNum == TestFunctions.RASTRIGIN_FUNCTION_NUM) {
//       			rastriginShiftVector = shiftVector.getCopy();
//       			rastriginBiasValue = useIncomingBias? inBiasValue: evalRastriginNoShiftNoBias(shiftVector);
//       			rastriginShiftedOptValue = RASTRIGIN_OPT_VALUE + rastriginBiasValue;					
//       		}	
//       		else if (functionNum == TestFunctions.ACKLEY_FUNCTION_NUM) {
//       			ackleyShiftVector = shiftVector.getCopy();
//       			ackleyBiasValue = useIncomingBias? inBiasValue: evalAckleyNoShiftNoBias(shiftVector);
//       			ackleyShiftedOptValue = ACKLEY_OPT_VALUE + ackleyBiasValue;	
//       		}	
//       		else if (functionNum == TestFunctions.GRIEWANK_FUNCTION_NUM) {
//       			griewankShiftVector = shiftVector.getCopy();
//       			griewankBiasValue = useIncomingBias? inBiasValue: evalGriewankNoShiftNoBias(shiftVector);
//       			griewankShiftedOptValue = GRIEWANK_OPT_VALUE + griewankBiasValue;						
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_1_NUM) {
//       			penalizedFunction1ShiftVector = shiftVector.getCopy();
//       			penalizedFunction1BiasValue = useIncomingBias? inBiasValue: evalPenalizedFunction1NoShiftNoBias(shiftVector);
//       			penalizedFunction1ShiftedOptValue = PENALIZED_FUNCTION_1_OPT_VALUE + penalizedFunction1BiasValue;							
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_2_NUM) {
//       			penalizedFunction2ShiftVector = shiftVector.getCopy();
//       			penalizedFunction2BiasValue = useIncomingBias? inBiasValue: evalPenalizedFunction2NoShiftNoBias(shiftVector);
//       			penalizedFunction2ShiftedOptValue = PENALIZED_FUNCTION_2_OPT_VALUE + penalizedFunction2BiasValue;							
//       		}
//       		else if (functionNum == TestFunctions.SPHERE_FUNCTION_NUM) {
//       			sphereShiftVector = shiftVector.getCopy();
//       			sphereBiasValue = useIncomingBias? inBiasValue: evalSphereNoShiftNoBias(shiftVector);
//       			sphereShiftedOptValue = SPHERE_OPT_VALUE + sphereBiasValue;							
//       		}
//
//       	}
//       		
//       		
//       	public static double getShiftedOptValue(int functionNum) {
//
//       		double retVal = 0.0;
//       		
//       		if (functionNum == TestFunctions.SCHWEFEL_FUNCTION_NUM) {
//       			retVal = schwefelShiftedOptValue;
//       		}
//       		else if (functionNum == TestFunctions.ROSENBROCK_FUNCTION_NUM) {
//       			retVal = rosenbrockShiftedOptValue;
//       		}	
//       		else if (functionNum == TestFunctions.RASTRIGIN_FUNCTION_NUM) {
//       			retVal = rastriginShiftedOptValue;
//       		}	
//       		else if (functionNum == TestFunctions.ACKLEY_FUNCTION_NUM) {
//       			retVal = ackleyShiftedOptValue;
//       		}	
//       		else if (functionNum == TestFunctions.GRIEWANK_FUNCTION_NUM) {
//       			retVal = griewankShiftedOptValue;
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_1_NUM) {
//       			retVal = penalizedFunction1ShiftedOptValue;
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_2_NUM) {
//       			retVal = penalizedFunction2ShiftedOptValue;
//       		}
//       		else if (functionNum == TestFunctions.SPHERE_FUNCTION_NUM) {
//       			retVal = sphereShiftedOptValue;
//       		}
//
//       		return retVal;
//
//       	}
//
//
//       	public static DoubleVector getShiftedOptVector(int functionNum, int numDimensions) {
//
//       		DoubleVector retVal = new DoubleVector(numDimensions, OPT_COORD[functionNum]);
//       		
//       		if (functionNum == TestFunctions.SCHWEFEL_FUNCTION_NUM) {
//       			retVal.add(schwefelShiftVector);
//       		}
//       		else if (functionNum == TestFunctions.ROSENBROCK_FUNCTION_NUM) {
//       			retVal.add(rosenbrockShiftVector);
//       		}	
//       		else if (functionNum == TestFunctions.RASTRIGIN_FUNCTION_NUM) {
//       			retVal.add(rastriginShiftVector);
//       		}	
//       		else if (functionNum == TestFunctions.ACKLEY_FUNCTION_NUM) {
//       			retVal.add(ackleyShiftVector);
//       		}	
//       		else if (functionNum == TestFunctions.GRIEWANK_FUNCTION_NUM) {
//       			retVal.add(griewankShiftVector);
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_1_NUM) {
//       			retVal.add(penalizedFunction1ShiftVector);
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_2_NUM) {
//       			retVal.add(penalizedFunction2ShiftVector);
//       		}
//       		else if (functionNum == TestFunctions.SPHERE_FUNCTION_NUM) {
//       			retVal.add(sphereShiftVector);
//       		}
//
//       		return retVal;
//
//       	}
//
//       	public static DoubleVector getShiftVector(int functionNum) {
//
//       		DoubleVector retVal = null;
//       		
//       		if (functionNum == TestFunctions.SCHWEFEL_FUNCTION_NUM) {
//       			retVal = schwefelShiftVector;
//       		}
//       		else if (functionNum == TestFunctions.ROSENBROCK_FUNCTION_NUM) {
//       			retVal = rosenbrockShiftVector;
//       		}	
//       		else if (functionNum == TestFunctions.RASTRIGIN_FUNCTION_NUM) {
//       			retVal = rastriginShiftVector;
//       		}	
//       		else if (functionNum == TestFunctions.ACKLEY_FUNCTION_NUM) {
//       			retVal = ackleyShiftVector;
//       		}	
//       		else if (functionNum == TestFunctions.GRIEWANK_FUNCTION_NUM) {
//       			retVal = griewankShiftVector;
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_1_NUM) {
//       			retVal = penalizedFunction1ShiftVector;
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_2_NUM) {
//       			retVal = penalizedFunction2ShiftVector;
//       		}
//       		else if (functionNum == TestFunctions.SPHERE_FUNCTION_NUM) {
//       			retVal = sphereShiftVector;
//       		}
//
//       		return retVal.getCopy();
//
//       	}
//
//
//       	public static int getOptType(int functionNum) {
//       		return OPT_TYPES[functionNum];
//       	}
//
//
//       	public static double[] evalWithError(DoubleVector v, int functionNum) {
//
//       		++PSO.numFunctionEvaluations;
////       		System.out.println("                in TF  numFunctionEvaluations = " + PSO.numFunctionEvaluations);
//
//       		double[] retArray = new double[2];
//
//       		if (functionNum == TestFunctions.SCHWEFEL_FUNCTION_NUM) {
//       			retArray = TestFunctions.evalSchwefelWithError(v);
//       		}
//       		else if (functionNum == TestFunctions.ROSENBROCK_FUNCTION_NUM) {
//       			retArray = TestFunctions.evalRosenbrockWithError(v);
//       		}	
//       		else if (functionNum == TestFunctions.RASTRIGIN_FUNCTION_NUM) {
//       			retArray = TestFunctions.evalRastriginWithError(v);
//       		}	
//       		else if (functionNum == TestFunctions.ACKLEY_FUNCTION_NUM) {
//       			retArray = TestFunctions.evalAckleyWithError(v);
//       		}	
//       		else if (functionNum == TestFunctions.GRIEWANK_FUNCTION_NUM) {
//       			retArray = TestFunctions.evalGriewankWithError(v);
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_1_NUM) {
//       			retArray = TestFunctions.evalPenalizedFunction1WithError(v);
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_2_NUM) {
//       			retArray = TestFunctions.evalPenalizedFunction2WithError(v);
//       		}
//       		else if (functionNum == TestFunctions.SPHERE_FUNCTION_NUM) {
//       			retArray = TestFunctions.evalSphereWithError(v);
//       		}	
//
//       		return retArray;
//
//       	}
//
//       	// does not increment PSO.numFunctionEvaluations 
//       	public static double[] evalWithErrorDoNotCountFE (DoubleVector v, int functionNum) {
//
//
//       		double[] retArray = new double[2];
//
//       		if (functionNum == TestFunctions.SCHWEFEL_FUNCTION_NUM) {
//       			retArray = TestFunctions.evalSchwefelWithError(v);
//       		}
//       		else if (functionNum == TestFunctions.ROSENBROCK_FUNCTION_NUM) {
//       			retArray = TestFunctions.evalRosenbrockWithError(v);
//       		}	
//       		else if (functionNum == TestFunctions.RASTRIGIN_FUNCTION_NUM) {
//       			retArray = TestFunctions.evalRastriginWithError(v);
//       		}	
//       		else if (functionNum == TestFunctions.ACKLEY_FUNCTION_NUM) {
//       			retArray = TestFunctions.evalAckleyWithError(v);
//       		}	
//       		else if (functionNum == TestFunctions.GRIEWANK_FUNCTION_NUM) {
//       			retArray = TestFunctions.evalGriewankWithError(v);
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_1_NUM) {
//       			retArray = TestFunctions.evalPenalizedFunction1WithError(v);
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_2_NUM) {
//       			retArray = TestFunctions.evalPenalizedFunction2WithError(v);
//       		}
//       		else if (functionNum == TestFunctions.SPHERE_FUNCTION_NUM) {
//       			retArray = TestFunctions.evalSphereWithError(v);
//       		}	
//
//       		return retArray;
//
//       	}
//
//       	
////       	private static double eval(DoubleVector v, int functionNum) {
//       //
////       		double retVal = 0.0;
////       		
////       		if (functionNum == TestFunctions.SCHWEFEL_FUNCTION_NUM) {
////       			retVal = TestFunctions.evalSchwefel(v);
////       		}
////       		else if (functionNum == TestFunctions.ROSENBROCK_FUNCTION_NUM) {
////       			retVal = TestFunctions.evalRosenbrock(v);
////       		}	
////       		else if (functionNum == TestFunctions.RASTRIGIN_FUNCTION_NUM) {
////       			retVal = TestFunctions.evalRastrigin(v);
////       		}	
////       		else if (functionNum == TestFunctions.ACKLEY_FUNCTION_NUM) {
////       			retVal = TestFunctions.evalAckley(v);
////       		}	
////       		else if (functionNum == TestFunctions.GRIEWANK_FUNCTION_NUM) {
////       			retVal = TestFunctions.evalGriewank(v);
////       		}	
////       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_1_NUM) {
////       			retVal = TestFunctions.evalPenalizedFunction1(v);
////       		}	
////       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_2_NUM) {
////       			retVal = TestFunctions.evalPenalizedFunction2(v);
////       		}
////       		else if (functionNum == TestFunctions.SPHERE_FUNCTION_NUM) {
////       			retVal = TestFunctions.evalSphere(v);
////       		}	
//       //
////       		return retVal;
//       //
////       	}
//
//
//       	//  Schwefel Problem 2.26
//       	// 	minimum is 0.0, which occurs at (-420.9687,...,-420.9687)
//
//       	public static double[] evalSchwefelWithError (DoubleVector v) {
//
//       		double value = evalSchwefel(v);
//       		double error = value - schwefelShiftedOptValue;
//       		double[] retArray = new double[2];
//       		retArray[RESULTS_VAL_INDEX] = value;
//       		retArray[RESULTS_ERR_INDEX] = error;
//       		return retArray;
//
//       	}
//       	
//
//       	public static double evalSchwefel (DoubleVector y) {
//
//       		DoubleVector v = DoubleVector.sub(y, schwefelShiftVector);
//       		
//       		double sum = 0;
//       		for(int i = 0 ; i < v.size() ; ++i) {
//       			double xi = v.get(i);
//       			sum += xi * Math.sin(Math.sqrt(Math.abs(xi)));
//       		}
//
//       		
//       		return 418.9829 * v.size() + sum + schwefelBiasValue;
//
//       	}
//
//       	
//       	public static double evalSchwefelNoShiftNoBias (DoubleVector v) {
//
//       		double sum = 0;
//       		for(int i = 0 ; i < v.size() ; ++i) {
//       			double xi = v.get(i);
//       			sum += xi * Math.sin(Math.sqrt(Math.abs(xi)));
//       		}
//
//       		return 418.9829 * v.size() + sum;
//
//       	}
//
//
//       	
//       	
//       	
//       	//  Rosenbrock Function
//       	// 	minimum is 0.0, which occurs at (1.0,...,1.0)
//
//       	public static double[] evalRosenbrockWithError (DoubleVector v) {
//
//       		double value = evalRosenbrock(v);
//       		double error = value - rosenbrockShiftedOptValue;
//       		double[] retArray = new double[2];
//       		retArray[RESULTS_VAL_INDEX] = value;
//       		retArray[RESULTS_ERR_INDEX] = error;
//       		return retArray;
//
//       	}
//
//       	public static double evalRosenbrock (DoubleVector y) {
//
//       		DoubleVector v = DoubleVector.sub(y, rosenbrockShiftVector);
//
//       		double retVal = 0;
//       		// NOTE:  LAST DIMENSION SHOULD *NOT* BE INCLUDED IN CALCULATION
//       		for(int i = 0 ; i < v.size() - 1 ; ++i) {
//       			double xi = v.get(i);
//       			double xiPlusOne = v.get(i+1);
//       			retVal += 100.0 * Math.pow(xiPlusOne - xi*xi, 2.0) + Math.pow(xi-1.0, 2.0);
//       		}
//
//       		return retVal + rosenbrockBiasValue;
//
//       	}
//
//       	
//       	public static double evalRosenbrockNoShiftNoBias (DoubleVector v) {
//
//       		double retVal = 0;
//       		// NOTE:  LAST DIMENSION SHOULD *NOT* BE INCLUDED IN CALCULATION
//       		for(int i = 0 ; i < v.size() - 1 ; ++i) {
//       			double xi = v.get(i);
//       			double xiPlusOne = v.get(i+1);
//       			retVal += 100.0 * Math.pow(xiPlusOne - xi*xi, 2.0) + Math.pow(xi-1.0, 2.0);
//       		}
//
//       		return retVal;
//
//       	}
//
//
//       	
//       	
//       	//  Rastrigin Function
//       	// 	minimum is 0.0, which occurs at (0.0,...,0.0)
//
//       	public static double[] evalRastriginWithError (DoubleVector v) {
//
//       		double value = evalRastrigin(v);
//       		double error = value - rastriginShiftedOptValue;
//       		double[] retArray = new double[2];
//       		retArray[RESULTS_VAL_INDEX] = value;
//       		retArray[RESULTS_ERR_INDEX] = error;
//       		return retArray;
//
//       	}
//
//       	public static double evalRastrigin (DoubleVector y) {
//
//       		DoubleVector v = DoubleVector.sub(y, rastriginShiftVector);
//
//       		double retVal = 0;
//       		for(int i = 0 ; i < v.size() ; ++i) {
//       			double xi = v.get(i);
//       			retVal += xi*xi - 10.0*Math.cos(2.0*Math.PI*xi) + 10.0;
//       		}
//
//       		return retVal + rastriginBiasValue;
//
//       	}
//
//
//       	public static double evalRastriginNoShiftNoBias (DoubleVector v) {
//
//       		double retVal = 0;
//       		for(int i = 0 ; i < v.size() ; ++i) {
//       			double xi = v.get(i);
//       			retVal += xi*xi - 10.0*Math.cos(2.0*Math.PI*xi) + 10.0;
//       		}
//
//       		return retVal;
//
//       	}
//       	
//       	
//       	
//
//       	//  Ackley Function
//       	// 	minimum is 0.0, which occurs at (0.0,...,0.0)
//
//       	public static double[] evalAckleyWithError (DoubleVector v) {
//
//       		double value = evalAckley(v);
//       		double error = value - ackleyShiftedOptValue;
//       		double[] retArray = new double[2];
//       		retArray[RESULTS_VAL_INDEX] = value;
//       		retArray[RESULTS_ERR_INDEX] = error;
//       		return retArray;
//
//       	}
//
//       	public static double evalAckley (DoubleVector y) {
//
//       		DoubleVector v = DoubleVector.sub(y, ackleyShiftVector);
//
//       		double firstSum = 0.0;
//       		double secondSum = 0.0;
//       		for(int i = 0 ; i < v.size() ; ++i) {
//       			double xi = v.get(i);
//       			firstSum += xi * xi;
//       			secondSum += Math.cos(2.0*Math.PI*xi);
//       		}
//
//       		return -20.0 * Math.exp(-0.2 * Math.sqrt(firstSum/v.size())) - 
//       			Math.exp(secondSum/v.size()) + 20.0 + Math.E + ackleyBiasValue;
//
//       	}	
//
//       	
//       	public static double evalAckleyNoShiftNoBias (DoubleVector v) {
//       		
//       		double firstSum = 0.0;
//       		double secondSum = 0.0;
//       		for(int i = 0 ; i < v.size() ; ++i) {
//       			double xi = v.get(i);
//       			firstSum += xi * xi;
//       			secondSum += Math.cos(2.0*Math.PI*xi);
//       		}
//
//       		return -20.0 * Math.exp(-0.2 * Math.sqrt(firstSum/v.size())) - 
//       			Math.exp(secondSum/v.size()) + 20.0 + Math.E;
//
//       	}	
//
//       	
//       	
//
//       	// Griewank function
//       	// 	minimum is 0.0, which occurs at (0.0,...,0.0)
//
//       	public static double[] evalGriewankWithError (DoubleVector v) {
//
//       		double value = evalGriewank(v);
//       		double error = value - griewankShiftedOptValue;
//       		double[] retArray = new double[2];
//       		retArray[RESULTS_VAL_INDEX] = value;
//       		retArray[RESULTS_ERR_INDEX] = error;
//       		return retArray;
//
//       	}
//
//       	public static double evalGriewank (DoubleVector y) {
//
//       		DoubleVector v = DoubleVector.sub(y, griewankShiftVector);
//
//       		double sumSquares = 0.0;
//       		double productCos = 1.0;
//       		for(int i = 0 ; i < v.size() ; ++i) {
//       			double xi = v.get(i);
//       			sumSquares += xi * xi;
//       			productCos *= Math.cos(xi/Math.sqrt(i+1));
//       		}
//
//       		return sumSquares/4000.0 - productCos + 1.0 + griewankBiasValue;
//
//       	}	
//
//       	public static double evalGriewankNoShiftNoBias (DoubleVector v) {
//
//       		double sumSquares = 0.0;
//       		double productCos = 1.0;
//       		for(int i = 0 ; i < v.size() ; ++i) {
//       			double xi = v.get(i);
//       			sumSquares += xi * xi;
//       			productCos *= Math.cos(xi/Math.sqrt(i+1));
//       		}
//
//       		return sumSquares/4000.0 - productCos + 1.0;
//
//       	}	
//
//
//
//       	// Penalized Function 1
//       	// 	minimum is 0.0, which occurs at (1.0,...,1.0)
//       	
//       	public static double[] evalPenalizedFunction1WithError (DoubleVector v) {
//
//       		double value = evalPenalizedFunction1(v);
//       		double error = value - penalizedFunction1ShiftedOptValue;
//       		double[] retArray = new double[2];
//       		retArray[RESULTS_VAL_INDEX] = value;
//       		retArray[RESULTS_ERR_INDEX] = error;
//       		return retArray;
//
//       	}
//
//
//       	
//       	public static double evalPenalizedFunction1 (DoubleVector y) {
//
//       		DoubleVector x = DoubleVector.sub(y, penalizedFunction1ShiftVector);
//
//       		int xSize = x.size();
//
//       		double firstSum = 0.0;
//       		for(int i = 0 ; i < xSize - 1 ; ++i) {
//       			double yiAti = yi(x.get(i));
//       			double sinTerm = Math.sin(Math.PI*yi(x.get(i+1)));
//       			firstSum += (yiAti-1.0) * (yiAti-1.0) * (1.0 + 10.0*sinTerm*sinTerm);
//       		}
//
//       		double uSum = 0.0;
//       		for(int i = 0 ; i < xSize ; ++i) {
//       			uSum += u(x.get(i), 10.0, 100.0, 4.0);
//       		}
//
//       		double firstSinTerm = Math.sin(Math.PI*yi(x.get(0)));
//       		double yiAtnMinus1 = yi(x.get(xSize-1));
//       		
//       		// some papers have PI/30.0 (first coefficient), but more papers have PI/n  where n = number of dimensions
//       		return Math.PI/xSize * (10.0*firstSinTerm*firstSinTerm + firstSum + (yiAtnMinus1-1.0)*(yiAtnMinus1-1.0)) + uSum + penalizedFunction1BiasValue;
//
//       	}	
//
//       	
//       	public static double evalPenalizedFunction1NoShiftNoBias (DoubleVector x) {
//
//       		int xSize = x.size();
//
//       		double firstSum = 0.0;
//       		for(int i = 0 ; i < xSize - 1 ; ++i) {
//       			double yiAti = yi(x.get(i));
//       			double sinTerm = Math.sin(Math.PI*yi(x.get(i+1)));
//       			firstSum += (yiAti-1.0) * (yiAti-1.0) * (1.0 + 10.0*sinTerm*sinTerm);
//       		}
//
//       		double uSum = 0.0;
//       		for(int i = 0 ; i < xSize ; ++i) {
//       			uSum += u(x.get(i), 10.0, 100.0, 4.0);
//       		}
//
//       		double firstSinTerm = Math.sin(Math.PI*yi(x.get(0)));
//       		double yiAtnMinus1 = yi(x.get(xSize-1));
//       		
//       		// some papers have PI/30.0 (first coefficient), but more papers have PI/n  where n = number of dimensions
//       		return Math.PI/xSize * (10.0*firstSinTerm*firstSinTerm + firstSum + (yiAtnMinus1-1.0)*(yiAtnMinus1-1.0)) + uSum;
//
//       	}	
//
//       	
//       	// needed in Penalized Function 1
//       	private static double yi(double xi) {
//       		return 1.0 + (xi + 1.0) / 4.0;
//       	}
//
//
//
//
//       	// Penalized Function 2
//       	// 	minimum is 0.0, which occurs at (1.0,...,1.0)
//
//       	public static double[] evalPenalizedFunction2WithError (DoubleVector v) {
//
//       		double value = evalPenalizedFunction2(v);
//       		double error = value - penalizedFunction2ShiftedOptValue;
//       		double[] retArray = new double[2];
//       		retArray[RESULTS_VAL_INDEX] = value;
//       		retArray[RESULTS_ERR_INDEX] = error;
//       		return retArray;
//
//       	}
//
//
//       	
//       	public static double evalPenalizedFunction2 (DoubleVector y) {
//
//       		DoubleVector x = DoubleVector.sub(y, penalizedFunction2ShiftVector);
//
//       		int xSize = x.size();
//
//       		double firstSum = 0.0;
//       		for(int i = 0 ; i < xSize - 1 ; ++i) {
//       			double xi = x.get(i);
//       			double sinTerm = Math.sin(3.0*Math.PI*x.get(i+1));
//       			firstSum += (xi-1.0)*(xi-1.0) * (1.0 + sinTerm*sinTerm);
//       		}
//
//       		double uSum = 0.0;
//       		for(int i = 0 ; i < xSize ; ++i) {
//       			uSum += u(x.get(i), 5.0, 100.0, 4.0);
//       		}
//
//       		double firstSinTerm = Math.sin(3.0*Math.PI*x.get(0));
//       		double lastXi = x.get(xSize-1);
//       		double lastSinTerm = Math.sin(2.0*Math.PI*lastXi);
//       		return 0.1 * (firstSinTerm*firstSinTerm + firstSum + (lastXi-1.0)*(lastXi-1.0) * (1.0 + lastSinTerm*lastSinTerm)) + uSum + penalizedFunction2BiasValue;
//
//       	}	
//
//       	public static double evalPenalizedFunction2NoShiftNoBias (DoubleVector x) {
//
//       		int xSize = x.size();
//
//       		double firstSum = 0.0;
//       		for(int i = 0 ; i < xSize - 1 ; ++i) {
//       			double xi = x.get(i);
//       			double sinTerm = Math.sin(3.0*Math.PI*x.get(i+1));
//       			firstSum += (xi-1.0)*(xi-1.0) * (1.0 + sinTerm*sinTerm);
//       		}
//
//       		double uSum = 0.0;
//       		for(int i = 0 ; i < xSize ; ++i) {
//       			uSum += u(x.get(i), 5.0, 100.0, 4.0);
//       		}
//
//       		double firstSinTerm = Math.sin(3.0*Math.PI*x.get(0));
//       		double lastXi = x.get(xSize-1);
//       		double lastSinTerm = Math.sin(2.0*Math.PI*lastXi);
//       		return 0.1 * (firstSinTerm*firstSinTerm + firstSum + (lastXi-1.0)*(lastXi-1.0) * (1.0 + lastSinTerm*lastSinTerm)) + uSum;
//
//       	}	
//
//
//       	// needed in Penalized Function 1 and Penalized Function 2
//       	private static double u(double xi, double a, double k, double m) {
//
//       		double retVal = 0;
//       		if (xi > a) {
//       			retVal = k * Math.pow(xi-a, m);
//       		}
//       		else if (xi <= a && xi >= -a) {
//       			retVal = 0.0;
//       		}
//       		else {  // xi < -a
//       			retVal = k * Math.pow(-xi-a, m);
//       		}
//
//       		return retVal;
//
//       	}
//
//
//       	
//       	// Sphere function
//       	// 	minimum is 0.0, which occurs at (0.0,...,0.0)
//       	
//       	public static double[] evalSphereWithError (DoubleVector v) {
//
//       		double value = evalSphere(v);
//       		double error = value - sphereShiftedOptValue;
//       		double[] retArray = new double[2];
//       		retArray[RESULTS_VAL_INDEX] = value;
//       		retArray[RESULTS_ERR_INDEX] = error;
//       		return retArray;
//
//       	}
//
//       	public static double evalSphere (DoubleVector y) {
//
//       		DoubleVector v = DoubleVector.sub(y, sphereShiftVector);
//
//       		double sumSquares = 0.0;
//       		for(int i = 0 ; i < v.size() ; ++i) {
//       			double xi = v.get(i);
//       			sumSquares += xi * xi;
//       		}
//
//       		return sumSquares + sphereBiasValue;
//
//       	}	
//
//       	public static double evalSphereNoShiftNoBias (DoubleVector v) {
//
//       		double sumSquares = 0.0;
//       		for(int i = 0 ; i < v.size() ; ++i) {
//       			double xi = v.get(i);
//       			sumSquares += xi * xi;
//       		}
//
//       		return sumSquares;
//
//       	}	
//
//
//       	
//       	public static String getFunctionName(int functionNum) {
//
//       		String retName = "";
//       		
//       		if (functionNum == TestFunctions.SCHWEFEL_FUNCTION_NUM) {
//       			retName = "SCHWEFEL_FUNCTION";
//       		}
//       		else if (functionNum == TestFunctions.ROSENBROCK_FUNCTION_NUM) {
//       			retName = "ROSENBROCK_FUNCTION";
//       		}	
//       		else if (functionNum == TestFunctions.RASTRIGIN_FUNCTION_NUM) {
//       			retName = "RASTRIGIN_FUNCTION";
//       		}	
//       		else if (functionNum == TestFunctions.ACKLEY_FUNCTION_NUM) {
//       			retName = "ACKLEY_FUNCTION";
//       		}	
//       		else if (functionNum == TestFunctions.GRIEWANK_FUNCTION_NUM) {
//       			retName = "GRIEWANK_FUNCTION";
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_1_NUM) {
//       			retName = "PENALIZED_FUNCTION_1";
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_2_NUM) {
//       			retName = "PENALIZED_FUNCTION_2";
//       		}
//       		else if (functionNum == TestFunctions.SPHERE_FUNCTION_NUM) {
//       			retName = "SPHERE_FUNCTION";
//       		}	
//       		
//       		return retName;
//
//       	}
//
//
//       	public static void printlnFunctionName(int functionNum) {
//
//       		if (functionNum == TestFunctions.SCHWEFEL_FUNCTION_NUM) {
//       			System.out.println("SCHWEFEL_FUNCTION");
//       		}
//       		else if (functionNum == TestFunctions.ROSENBROCK_FUNCTION_NUM) {
//       			System.out.println("ROSENBROCK_FUNCTION");
//       		}	
//       		else if (functionNum == TestFunctions.RASTRIGIN_FUNCTION_NUM) {
//       			System.out.println("RASTRIGIN_FUNCTION");
//       		}	
//       		else if (functionNum == TestFunctions.ACKLEY_FUNCTION_NUM) {
//       			System.out.println("ACKLEY_FUNCTION");
//       		}	
//       		else if (functionNum == TestFunctions.GRIEWANK_FUNCTION_NUM) {
//       			System.out.println("GRIEWANK_FUNCTION");
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_1_NUM) {
//       			System.out.println("PENALIZED_FUNCTION_1");
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_2_NUM) {
//       			System.out.println("PENALIZED_FUNCTION_2");
//       		}
//       		else if (functionNum == TestFunctions.SPHERE_FUNCTION_NUM) {
//       			System.out.println("SPHERE_FUNCTION");
//       		}	
//
//       	}
//
//       	public static void printFunctionName(int functionNum) {
//
//       		if (functionNum == TestFunctions.SCHWEFEL_FUNCTION_NUM) {
//       			System.out.print("SCHWEFEL_FUNCTION");
//       		}
//       		else if (functionNum == TestFunctions.ROSENBROCK_FUNCTION_NUM) {
//       			System.out.print("ROSENBROCK_FUNCTION");
//       		}	
//       		else if (functionNum == TestFunctions.RASTRIGIN_FUNCTION_NUM) {
//       			System.out.print("RASTRIGIN_FUNCTION");
//       		}	
//       		else if (functionNum == TestFunctions.ACKLEY_FUNCTION_NUM) {
//       			System.out.print("ACKLEY_FUNCTION");
//       		}	
//       		else if (functionNum == TestFunctions.GRIEWANK_FUNCTION_NUM) {
//       			System.out.print("GRIEWANK_FUNCTION");
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_1_NUM) {
//       			System.out.print("PENALIZED_FUNCTION_1");
//       		}	
//       		else if (functionNum == TestFunctions.PENALIZED_FUNCTION_2_NUM) {
//       			System.out.print("PENALIZED_FUNCTION_2");
//       		}
//       		else if (functionNum == TestFunctions.SPHERE_FUNCTION_NUM) {
//       			System.out.print("SPHERE_FUNCTION");
//       		}	
//
//       	}
//
//
//       	//	public static double eval (DoubleVector v) {
//       	//	
//       	//	double x = v.get(0);
//       	//	double y = v.get(1);
//       	//	return x * x + y * y; 
//       	//	
//       	//}
//
//       	//public static double eval (DoubleVector v) {
//       	//	
//       	//	double x = v.get(0);
//       	//	double y = v.get(1);
//       	//	return (x - 500) * (x - 500) + (y - 500) * (y - 500); 
//       	//	
//       	//}
//
//
//
//       }
//
