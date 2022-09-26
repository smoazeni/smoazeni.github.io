package EnergySequential;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;
import java.util.logging.StreamHandler;
import java.util.StringTokenizer;
import java.text.DecimalFormat;

//import LSOC.Exception;
//import LSOC.Main;

import ilog.concert.*;
import ilog.cplex.*;

/**
 * <p>
 * Title: Load & Source Optimization Controller - LSOC
 * </p>
 * <p>
 * Description: Solves an LP or a MIP, using CPLEX.
 * </p>
 * <p>
 * Copyright: Copyright (c) 2011 by The Trustees of Princeton University
 * </p>
 * <p>
 * Company: CASTLE Lab - Princeton University
 * </p>
 * 
 * @author Hugo P. Simao, modified by Hyun Bin (Vince) Jeong
 * @version 1.0, Jan/11
 */

public class MIPSparseSolver {

	// an ID to be associated to this model
	private String modelId = null;
	// the tolerance for integer precision
	private double epsInt = 1e-6;
	// the tolerance for relative gap
	private double epsGap = 1.0 / 10000.;
	// the optimizer time limit, in seconds
	private double timLim = 5000.0;
	// the maximum # of threads to be used
	private int maxThreads = 8;
	// the CPLEX model
	private IloCplex iloModel = null;
	// the LP matrix behind this CPLEX model
	private IloLPMatrix iloLPMatrix = null;
	// the objective function associated to this model, as well as the variables
	// and their coefficients
	private IloObjective iloObjFctn = null;
	private IloNumVar[] iloVarArray = null;
	private double[] coeffArray = null;
	// the array of internal indices of the constraints
	private ArrayList<Integer> rowIdx = null;

	// some useful flags
	private boolean solved = false;
	private int maxSolveTrials = 1;
	// for debugging purposes
	private boolean DEBUG = false;
	public static int exportNumb = 0;
	// for reporting purposes
	private Logger logger = null;
	private FileOutputStream foStream = null;

	private static double plusInf = 1e60;
	private static double minusInf = -1e60;
	private static DecimalFormat decFmt = new DecimalFormat("#,##0.00");

	/**
	 * Returns an instance of this object, initialized with the given arguments.
	 * 
	 * @param c
	 *            double[]
	 * @param A
	 *            ArrayList<ArrayList<Double>>
	 * @param b
	 *            ArrayList<Double>
	 * @param Aeq
	 *            ArrayList<ArrayList<Double>>
	 * @param beq
	 *            ArrayList<Double>
	 * @param lB
	 *            double[]
	 * @param uB
	 *            double[]
	 * @param yIdx
	 *            int[]
	 * @param varIdList
	 *            ArrayList<String>
	 * @param minObj
	 *            boolean
	 */
	public MIPSparseSolver(String modelId_, double[] c,
			ArrayList<ArrayList<Double>> A, ArrayList<Double> b,
			ArrayList<ArrayList<Double>> Aeq, ArrayList<Double> beq,
			double[] lB, double[] uB, int[] yIdx, ArrayList<String> varIdList,
			boolean minObj, double epsInt_, FileOutputStream myFoStream,
			Logger myLogger, boolean DEBUG_EXPORT_LP, boolean debugFlag) {
		DEBUG = debugFlag;
		// create a logger
		if (myLogger == null)
			createLogger();
		else {
			foStream = myFoStream;
			logger = myLogger;
		}

		try {
			long elapTime = System.currentTimeMillis();
			if (DEBUG)
				logger.info("> creating CPLEX model");
			// instantiate the CPLEX model object
			IloCplex tempIloModel = new IloCplex();
			modelId = modelId_;
			epsInt = epsInt_;

			// set the level of verbose of the solver and direct its messaging
			// output
			if (DEBUG) {
				tempIloModel.setOut(foStream);
				tempIloModel.setWarning(foStream);
			} else {
				tempIloModel.setOut(null);
				tempIloModel.setWarning(null);
			}

			// initialize some parameters for the IP solver
			if (yIdx != null && yIdx.length > 0) {
				if (DEBUG)
					logger.info("> setting IP solver parameters");
				tempIloModel.setParam(IloCplex.DoubleParam.EpGap, epsGap);
				tempIloModel.setParam(IloCplex.DoubleParam.EpInt, epsInt);
				tempIloModel.setParam(IloCplex.DoubleParam.TiLim, timLim);
				maxSolveTrials = 2;
			} else
				maxSolveTrials = 1;
			// find out the # of processors available in the running environment
			int numProcessors = Runtime.getRuntime().availableProcessors();
			/*
			 * for temporary debugging purposes only
			 * logger.warning("> max # of 'processors' available: " +
			 * numProcessors); //
			 */
			// limit the # of threads to be used, if applicable
			if (numProcessors > maxThreads) {
				tempIloModel.setParam(IloCplex.IntParam.Threads, maxThreads);
				logger.info("> limited # of threads to be used to "
						+ maxThreads);
				// force "deterministic" parallel mode (when the # of threads is
				// limited, "opportunistic" mode becomes the default)
				tempIloModel.setParam(IloCplex.IntParam.ParallelMode, 1);
			}

			if (DEBUG)
				logger.info("> creating LP matrix");
			// add an LP matrix to the problem
			iloLPMatrix = tempIloModel.addLPMatrix();

			// checking the # of variables and constraints

			int numVars = 0;
			if (c != null)
				numVars = c.length;

			int numConstraints = 0;
			if (A != null) {
				numConstraints += A.size();
				if (b == null) {
					throw new IloException(
							"LP has no vector of RHS values for the inequalities");
				} else if (A.size() != b.size()) {
					throw new IloException(
							"LP has size of vector of RHS values for the inequalities ("
									+ b.size()
									+ ") different from corresponding # of constraints ("
									+ A.size() + ")");
				}
			}
			if (Aeq != null) {
				numConstraints += Aeq.size();
				if (beq == null) {
					throw new IloException(
							"LP has no vector of RHS values for the equalities");
				} else if (Aeq.size() != beq.size()) {
					throw new IloException(
							"LP has size of vector of RHS values for the equalities ("
									+ beq.size()
									+ ") different from corresponding # of constraints ("
									+ Aeq.size() + ")");
				}
			}
			if (numVars == 0) {
				throw new IloException("LP has no variables");
			} else if (numConstraints == 0) {
				throw new IloException("LP has no constraints");
			}

			if (DEBUG)
				logger.info("> creating LP columns");
			// create the set of variables
			IloNumVarType[] iloType = new IloNumVarType[numVars];
			for (int i = 0; i < numVars; i++)
				iloType[i] = IloNumVarType.Float;
			if (yIdx != null && yIdx.length > 0) {
				if (DEBUG)
					logger.info("> designating integer columns");
				// set the binary variables
				for (int ii = 0; ii < yIdx.length; ii++) {
					iloType[yIdx[ii]] = IloNumVarType.Int;
				}
			}
			if (lB == null) {
				// assuming the typical default values for the lower bounds (0)
				lB = new double[numVars];
				for (int i = 0; i < numVars; i++)
					lB[i] = 0;
			}
			if (uB == null) {
				// assuming the typical default values for the upper bounds (max
				// double value)
				uB = new double[numVars];
				for (int i = 0; i < numVars; i++)
					uB[i] = plusInf;
			}
			iloVarArray = tempIloModel.numVarArray(
					tempIloModel.columnArray(iloLPMatrix, numVars), lB, uB,
					iloType);
			if (varIdList != null || DEBUG_EXPORT_LP) {
				// set the variable names for debugging purposes only
				for (int i = 0; i < numVars; i++)
					iloVarArray[i].setName(varIdList.get(i));
			}

			if (DEBUG)
				logger.info("> creating LP objective function");
			// set the objective function
			if (minObj)
				iloObjFctn = tempIloModel.minimize(tempIloModel.scalProd(
						iloVarArray, c));
			else
				iloObjFctn = tempIloModel.maximize(tempIloModel.scalProd(
						iloVarArray, c));
			tempIloModel.add(iloObjFctn);
			coeffArray = c;
			// set a name for the objective function row for debugging purposes
			// only
			if (DEBUG)
				iloObjFctn.setName("TotCost");

			// allocate the array of constraint indices
			rowIdx = new ArrayList<Integer>(numConstraints);

			// create all the equality constraints
			if (Aeq != null && Aeq.size() > 0) {
				if (DEBUG)
					logger.info("> creating LP equality constraints");
				// loop over all equality constraints
				for (int i = 0; i < Aeq.size(); i++) {
					ArrayList<Double> constraint = Aeq.get(i);
					int numNZ = constraint.size() / 2;
					if (numNZ > 0) {
						// allocate the arrays to store the nonzero elements of
						// the row
						int ind[] = new int[numNZ];
						double a[] = new double[numNZ];
						int k = 0;
						for (int j = 0; j < constraint.size(); j = j + 2) {
							ind[k] = (int) (double) constraint.get(j);
							a[k] = constraint.get(j + 1);
							k++;
						}
						// add the equality constraint to the problem
						int conIdx = iloLPMatrix.addRow(beq.get(i), beq.get(i),
								ind, a);
						/*
						 * for temporary debugging purposes only //if (i == 0)
						 * logger
						 * .info("the very first row in LPMatrix has row number "
						 * + conIdx);
						 * //logger.info("C-PLEX is adding a constraint to the "
						 * + conIdx + "-th row."); if (i != conIdx)
						 * logger.info(" i = " + i + " but conIdx = " + conIdx);
						 * //
						 */
						rowIdx.add(conIdx);
					} else {
						rowIdx.add(-1);
						logger.warning("row "
								+ i
								+ " of equality matrix 'Aeq' has all coefficients null; please check!");
					}
				}
			}

			// create all the inequality constraints (less than or equal to)
			if (A != null && A.size() > 0) {
				if (DEBUG)
					logger.info("> creating LP inequality constraints");
				// loop over all inequality constraints
				for (int i = 0; i < A.size(); i++) {
					ArrayList<Double> constraint = A.get(i);

					int numNZ = constraint.size() / 2;
					/*
					 * // find the # of nonzero elements in this row of the
					 * matrix of coefficients int numNZ = 0; for (int j = 0; j <
					 * A[i].length; j++) { if (A[i][j] != 0) numNZ++; }
					 */
					if (numNZ > 0) {
						// allocate the arrays to store the nonzero elements of
						// the row
						int ind[] = new int[numNZ];
						double a[] = new double[numNZ];
						int k = 0;
						for (int j = 0; j < constraint.size(); j = j + 2) {
							ind[k] = (int) (double) constraint.get(j);
							a[k] = constraint.get(j + 1);
							k++;
						}
						// add the inequality constraint to the problem
						int conIdx = iloLPMatrix.addRow(minusInf, b.get(i),
								ind, a);
						rowIdx.add(conIdx);
					} else {
						rowIdx.add(-1);
						logger.warning("row "
								+ i
								+ " of inequality matrix 'A' has all coefficients null; please check!");
					}
				}
			}

			// store the reference to the CPLEX model
			iloModel = tempIloModel;
			solved = false;
			if (DEBUG) {
				logger.info("> total # of vars: "
						+ numVars
						+ " ("
						+ (yIdx != null ? yIdx.length : 0)
						+ " int); # of rows: "
						+ ((A != null ? A.size() : 0) + (Aeq != null ? Aeq
								.size() : 0)));
			}

			// export the model to an external file
			if (varIdList != null && DEBUG_EXPORT_LP) {
				logger.warning("> exporting LP model #" + (++exportNumb)
						+ " to an external file right after instantiating it");
				String fname = modelId + "_LPModel_" + exportNumb + ".lp";
				iloModel.exportModel(fname);
			}
			if (DEBUG) {
				// Main.summarizeMemory();
				elapTime = System.currentTimeMillis() - elapTime;
				logger.info("> elapsed time to generate CPLEX model: "
						+ elapTime / 1000. + " sec");
			}
		} catch (IloException e) {
			logger.warning("Could not create a CPLEX model; caught exception "
					+ e);
			cleanModel();
		}

	}

	/**
	 * Creates a logger for this object.
	 */
	private void createLogger() {
		if (logger == null) {
			logger = Logger.getLogger("MIPSparseSolver");
			logger.setUseParentHandlers(false);
			try {
				foStream = new FileOutputStream("MIPSparseSolver.log");
				LogFormatter formatter = new LogFormatter();
				StreamHandler handler = new StreamHandler(foStream, formatter);
				formatter.addHandler(handler);
				logger.addHandler(handler);
			} catch (FileNotFoundException e) {
				System.out
						.println("Cannot open log file for CPLEX MIPSparseSolver object; please check! Exception: "
								+ e);
				System.out.flush();
			}
			if (DEBUG)
				logger.info("> created logger");
		}
	}

	/** returns whether or not a valid CPLEX model was instantiated */
	public boolean isValid() {
		return (iloModel != null);
	}

	/** return the model ID */
	public String getModelId() {
		return modelId;
	}

	/** return the model integer tolerance */
	public double getEpsInt() {
		return epsInt;
	}

	/** return LP model for the C-PLEX */
	public IloCplex getIloModel() {
		return iloModel;
	}

	/** return LP matrix for the C-PLEX */
	public IloLPMatrix getIloLPMatrix() {
		return iloLPMatrix;
	}

	/** return LP matrix for the C-PLEX */
	public ArrayList<Integer> getRowIdx() {
		return rowIdx;
	}

	/** return LP matrix for the C-PLEX */
	public void setSolved(boolean isSolved) {
		solved = isSolved;
	}

	/** set debug flag */
	public void setDebugFlag(boolean isOn) {
		DEBUG = isOn;
	}

	/**
	 * Solves the LP stored in this object. Returns a true/false flag to tell
	 * whether this problem has been solved.
	 * 
	 * @return boolean
	 */
	public boolean solve() {
		try {
			// solve it now
			int solveTrialsCount = 1;
			while (!solved && solveTrialsCount <= maxSolveTrials) {
				if (DEBUG)
					logger.info("> solving the problem");
				/*
				 * for temporary debugging purposes only // export the LP model
				 * to a file logger.warning("> exporting LP model #" +
				 * (++exportNumb) +
				 * " to an external file right before solving it"); String
				 * afname = modelId + "_LPModel_" + exportNumb + ".lp";
				 * iloModel.exportModel(afname); //
				 */
				solved = iloModel.solve();
				if (!solved) {
					logger.warning("> failed to solve CPLEX model at trial #"
							+ solveTrialsCount + "; CPLEX solution status: "
							+ iloModel.getCplexStatus());
					solveTrialsCount++;
					if (solveTrialsCount <= maxSolveTrials) {
						if (iloModel.getCplexStatus().toString()
								.contains("InfOrUnbd")
								|| iloModel.getCplexStatus().toString()
										.contains("Infeasible")) {
							// export the LP model to a file
							logger.warning("> exporting LP model #"
									+ (++exportNumb)
									+ " to an external file right after failing to solve it");
							String fname = modelId + "_LPModel_" + exportNumb
									+ ".lp";
							iloModel.exportModel(fname);
							// clear the current model
							logger.warning("> removing all modeling objects from the LP model");
							iloModel.clearModel();
							// re-import the LP model
							logger.warning("> re-importing the *same* LP model from the external file; will try to solve it again!");
							iloModel.importModel(fname);
							int numMatrices = 0;
							for (Iterator loop = iloModel.LPMatrixIterator(); loop
									.hasNext(); numMatrices++) {
								iloLPMatrix = (IloLPMatrix) loop.next();
							}
							if (numMatrices > 1) {
								logger.log(
										Level.SEVERE,
										"> after clearing and re-importing LP model, active model has "
												+ numMatrices
												+ " LP matrices; do not know what to do, please check!");
							}
							iloObjFctn = iloModel.getObjective();
							iloVarArray = null;
							// initialize some parameters for the IP solver, if
							// applicable
							if (iloModel.isMIP()) {
								logger.warning("> resetting IP solver parameters");
								iloModel.setParam(IloCplex.DoubleParam.EpGap,
										epsGap);
								iloModel.setParam(IloCplex.DoubleParam.EpInt,
										epsInt);
								iloModel.setParam(IloCplex.DoubleParam.TiLim,
										timLim);
							}
						} else {
							logger.warning("> relaxing MIP solver parameters; will try to solve model again!");
							iloModel.setParam(
									IloCplex.DoubleParam.EpGap,
									5 * iloModel
											.getParam(IloCplex.DoubleParam.EpGap));
							iloModel.setParam(
									IloCplex.DoubleParam.EpInt,
									10 * iloModel
											.getParam(IloCplex.DoubleParam.EpInt));
							// 'reset' the objective function, in order to force
							// a re-optimization
							resetObjectiveFunction();
						}
						// reset the level of verbose of the solver and direct
						// its messaging output
						iloModel.setOut(foStream);
						iloModel.setWarning(foStream);
					} else {
						logger.warning("> exporting LP model #"
								+ (++exportNumb)
								+ " to an external file right before aborting");
						String fname = modelId + "_LPModel_" + exportNumb
								+ ".lp";
						iloModel.exportModel(fname);
						logger.warning("> CPLEX status: "
								+ iloModel.getCplexStatus().toString());
						System.exit(-1);
					}
				}
			}
			// restore the original settings, if applicable
			if (solveTrialsCount > 1) {
				logger.warning("> restoring MIP solver parameters");
				iloModel.setParam(IloCplex.DoubleParam.EpGap, epsGap);
				if (!DEBUG) {
					iloModel.setOut(null);
					iloModel.setWarning(null);
				}
			}
		} catch (IloException e) {
			logger.warning("MIPSparseSolver.solve: caught an IloException: "
					+ e);
		}
		return solved;
	}

	/** Reset the objective function of the active model */
	public void resetObjectiveFunction() {
		try {
			// 'reset' all coefficients of the objective function
			iloModel.setLinearCoefs(iloObjFctn, iloVarArray, coeffArray);
		} catch (IloException e) {
			logger.warning("MIPSparseSolver.resetObjectiveFunction: caught an IloException: "
					+ e);
		}
	}

	/**
	 * Returns a vector with the solution to this problem, if it has been
	 * solved. Otherwise, returns null.
	 * 
	 * @return double[]
	 */
	public double[] retrieveSolution() {
		if (solved) {
			try {
				if (DEBUG) {
					logger.info("> retrieving solution");
					double objFctnVal = iloModel.getObjValue();
					logger.info("  obj fctn value: "
							+ decFmt.format(objFctnVal));
				}
				return iloModel.getValues(iloLPMatrix);
			} catch (IloException e) {
				logger.warning("Caught exception " + e
						+ " while trying to retrieve solution; please check!");
				return null;
			}
		} else
			return null;
	}

	/**
	 * Returns the objective function value associated to the solution to this
	 * problem, if it has been solved. Otherwise, returns 0.
	 * 
	 * @return double
	 */
	public double retrieveObjFctnValue() {
		if (solved) {
			try {
				if (DEBUG)
					logger.info("> retrieving objective function value");
				return iloModel.getObjValue();
			} catch (IloException e) {
				logger.warning("Caught exception "
						+ e
						+ " while trying to retrieve objective function value; please check!");
				return 0;
			}
		} else
			return 0;
	}

	/**
	 * Release all the internal CPLEX objects currently stored in this object.
	 */
	public void cleanModel() {
		try {
			if (isValid()) {
				iloModel.endModel();
				iloModel = null;
			}
			modelId = null;
			epsInt = 1e-6;
			epsGap = 1.0 / 100.;
			timLim = 5000.0;
			iloLPMatrix = null;
			iloObjFctn = null;
			iloVarArray = null;
			coeffArray = null;
			rowIdx = null;
			solved = false;
			maxSolveTrials = 1;
			if (DEBUG)
				logger.info("> ended CPLEX model");
			DEBUG = false;
		} catch (IloException e) {
			logger.warning("Caught exception " + e
					+ " while trying to end CPLEX model; please check!");
		}
	}

	/**
	 * Creates an instance of this object, initialized with the given arguments,
	 * solves it, and returns the vector of solution values for the underlying
	 * variable.
	 * 
	 * Solving min c^t x
	 * 
	 * subject to: A x <= b Aeq x - beq lB <= x <= uB some x may be {0,1}
	 * 
	 * @param c
	 *            double[]
	 * @param A
	 *            ArrayList<ArrayList<Double>> array list of array lists, each
	 *            one corresponding to an inequality constraint each constraint
	 *            array list is entered as a sequence of pairs of doubles,
	 *            corresponding to the sequence of (index of non-zero
	 *            coefficient, actual value of coefficient)
	 * @param b
	 *            ArrayList<Double> array list of actual double values (NOT in
	 *            sparse format)
	 * @param Aeq
	 *            ArrayList<ArrayList<Double>> see explanation for A above and
	 *            apply it to equality constraints
	 * @param beq
	 *            ArrayList<Double>
	 * @param lB
	 *            double[]
	 * @param uB
	 *            double[]
	 * @param yIdx
	 *            int[] // indices of the integer variables; if none, specify a
	 *            null array
	 * @param varIdList
	 *            ArrayList<String>
	 * @param minObj
	 *            boolean
	 * @param epsInt_
	 *            double
	 * @param debugFlag
	 *            boolean
	 * 
	 * @return double[]
	 */
	public static double[] createSolve(double[] c,
			ArrayList<ArrayList<Double>> A, ArrayList<Double> b,
			ArrayList<ArrayList<Double>> Aeq, ArrayList<Double> beq,
			double[] lB, double[] uB, int[] yIdx, ArrayList<String> varIdList,
			boolean minObj, double epsInt_, boolean debugFlag) {

		double[] xSol = null; // the return variable
		boolean DEBUG_EXPORT_LP = false;
		// instantiate the solver
		MIPSparseSolver solver = new MIPSparseSolver("Std", c, A, b, Aeq, beq,
				lB, uB, yIdx, varIdList, minObj, epsInt_, null, null,
				DEBUG_EXPORT_LP, debugFlag);
		// case when a valid CPLEX model was instantiated
		if (solver.isValid()) {
			// try to solve the problem
			if (solver.solve()) {
				// retrieve the solution
				xSol = solver.retrieveSolution();
			}
			// release the solver
			solver.cleanModel();
		}
		/* for temporary debugging purposes only */
		else {
			solver.logger.log(Level.WARNING,
					"Could not instantiate CPLEX object.");
		}
		// */

		return xSol;
	}

	/** Nested class which defines a simple log formatter */
	public class LogFormatter extends SimpleFormatter {
		protected final ArrayList<Handler> handlers = new ArrayList<Handler>(3);;

		/**
		 * @param handler_
		 *            Handler
		 */
		public void addHandler(Handler handler_) {
			handlers.add(handler_);
		}

		/**
		 * @param record_
		 *            LogRecord
		 * @return String
		 */
		@Override
		public synchronized String format(LogRecord record_) {
			record_.getLevel().toString();
			String string = record_.getLevel() + ": " + record_.getMessage()
					+ System.getProperty("line.separator");
			for (int i = 0; i < handlers.size(); i++) {
				Handler handler = (Handler) handlers.get(i);
				handler.flush();
			}
			return string;
		}
	}

	/** Main program: for testing purposes only */
	public static void main(String[] args) {
		// initialize the arrays of a small test MIP
		double[] c = { -9, -5, -6, -4 };
		double[][] A_mat = { { 6, 3, 5, 2 }, { 0, 0, 1, 1 }, { -1, 0, 1, 0 },
				{ 0, -1, 0, 1 } };
		double[] b_vec = { 9, 1, 0, 0 };
		ArrayList<ArrayList<Double>> A = new ArrayList<ArrayList<Double>>();
		ArrayList<Double> b = new ArrayList<Double>();
		for (int i = 0; i < A_mat.length; i++) {
			ArrayList<Double> a_row = new ArrayList<Double>();
			for (int j = 0; j < A_mat[i].length; j++) {
				if (A_mat[i][j] != 0) {
					a_row.add((double) j);
					a_row.add(A_mat[i][j]);
				}
			}
			A.add(a_row);
			b.add(b_vec[i]);
		}
		double[][] Aeq_mat = null;
		ArrayList<ArrayList<Double>> Aeq = null;
		double[] beq_vec = null;
		ArrayList<Double> beq = null;
		double[] lB = null;
		double[] uB = { 1, 1, 1, 1 };
		int[] yIdx = { 0, 1, 3 };
		boolean debug = true;
		// */
		/*
		 * parse the arguments of a larger MIP int numVars = 0; int numIneqCons
		 * = 0; int numEqCons = 0; if (args.length >= 3) { numVars =
		 * Integer.parseInt(args[0]); numIneqCons = Integer.parseInt(args[1]);
		 * numEqCons = Integer.parseInt(args[2]); System.out.println();
		 * System.out.println("MIP dimensions: #vars " + numVars +
		 * ", #IneqCons " + numIneqCons + ", #EqCons " + numEqCons);
		 * System.out.flush(); } else { System.out.println();
		 * System.out.println("Usage:"); System.out.println(
		 * "\tjava [RUNTIME-ARGUMENTS] Solvers.MIPSolver [PARAMETERS]");
		 * System.out.println(); System.out.println("PARAMETERS:");
		 * System.out.println
		 * ("\t<# Vars> <# Ineq Constraints> <# Eq Constraints>");
		 * System.out.flush(); System.exit( -1); } // read-in the data arrays
		 * double[] c = readNumVector("./c.txt", numVars); double[][] A =
		 * readNumMatrix("./A.txt", numIneqCons, numVars); double[] b =
		 * readNumVector("./b.txt", numIneqCons); double[][] Aeq =
		 * readNumMatrix("./Aeq.txt", numEqCons, numVars); double[] beq =
		 * readNumVector("./beq.txt", numEqCons); double[] lB =
		 * readNumVector("./lb.txt", numVars); double[] uB =
		 * readNumVector("./ub.txt", numVars); int[] yIdx =
		 * readBoolVector("./yidx.txt", numVars); boolean debug = true; //
		 */
		/* solve the MIP */
		double[] x = MIPSparseSolver.createSolve(c, A, b, null, null, lB, uB,
				null, null, true, 1e-8, debug);
		// output the solution vector
		System.out.println();
		if (x != null) {
			System.out.println("Writing solution vector onto an output file.");
			writeNumVector("./out_x.txt", x);
			// computing the objective function value of the solution
			double objFctnVal = 0;
			for (int i = 0; i < x.length; i++)
				objFctnVal += c[i] * x[i];
			System.out.println("Obj Fctn Value: " + decFmt.format(objFctnVal));
		} else {
			System.out.println("Could not find a solution. Check log file.");
		}
		System.out.flush();
		System.exit(0);
		// */
	}

	/** Methods to read-in/write-out the several data files */
	private static double[] readNumVector(String fname, int numElem) {
		double[] vector = new double[numElem]; // the return variable
		int i = -1;
		try {
			// setup the buffered reader to read in the file one record at a
			// time
			File filePtr = new File(fname);
			FileReader fileReader = new FileReader(filePtr);
			BufferedReader bufReader = new BufferedReader(fileReader);
			// read the file, record by record
			String nextRec = bufReader.readLine();
			i = 0;
			while (nextRec != null && i < numElem) {
				// parse the double value
				vector[i] = Double.parseDouble(nextRec.trim());
				// get the next record
				nextRec = bufReader.readLine();
				i++;
			}
			// close the input files
			bufReader.close();
			fileReader.close();
		} catch (Exception e) {
			System.out.println("Caught exception " + e + " while reading rec #"
					+ i + " of file " + fname);
			System.out.println("Exiting!");
			System.out.flush();
			System.exit(-1);
		}
		return vector;
	}

	private static int[] readBoolVector(String fname, int numElem) {
		int[] vector = null; // the return variable
		int i = -1;
		try {
			// setup the buffered reader to read in the file one record at a
			// time
			File filePtr = new File(fname);
			FileReader fileReader = new FileReader(filePtr);
			BufferedReader bufReader = new BufferedReader(fileReader);
			boolean[] origVector = new boolean[numElem]; // the return variable
			// read the file, record by record
			String nextRec = bufReader.readLine();
			i = 0;
			int actNumElem = 0;
			while (nextRec != null && i < numElem) {
				// parse the boolean value
				origVector[i] = new Boolean(nextRec.trim()).booleanValue();
				if (origVector[i])
					actNumElem++;
				// get the next record
				nextRec = bufReader.readLine();
				i++;
			}
			// create the actual desired vector of indices
			vector = new int[actNumElem];
			int ii = 0;
			for (i = 0; i < numElem; i++) {
				if (origVector[i])
					vector[ii] = i;
				ii++;
			}
			// close the input files
			bufReader.close();
			fileReader.close();
		} catch (Exception e) {
			System.out.println("Caught exception " + e + " while reading rec #"
					+ i + " of file " + fname);
			System.out.println("Exiting!");
			System.out.flush();
			System.exit(-1);
		}
		return vector;
	}

	private static double[][] readNumMatrix(String fname, int numRows,
			int numCols) {
		double[][] matrix = new double[numRows][numCols]; // the return variable
		int i = -1;
		int j = -1;
		try {
			// setup the buffered reader to read in the file one record at a
			// time
			File filePtr = new File(fname);
			FileReader fileReader = new FileReader(filePtr);
			BufferedReader bufReader = new BufferedReader(fileReader);
			// read the file, record by record
			String nextRec = bufReader.readLine();
			i = 0;
			while (nextRec != null && i < numRows) {
				// break the record into tokens
				StringTokenizer tokens = new StringTokenizer(nextRec);
				j = 0;
				while (tokens.hasMoreTokens() && j < numCols) {
					matrix[i][j] = Double
							.parseDouble(tokens.nextToken().trim());
					j++;
				}
				// get the next record
				nextRec = bufReader.readLine();
				i++;
			}
			// close the input readers
			bufReader.close();
			fileReader.close();
		} catch (Exception e) {
			System.out.println("Caught exception " + e
					+ " while reading cell (" + i + ',' + j + ") of file "
					+ fname);
			System.out.println("Exiting!");
			System.out.flush();
			System.exit(-1);
		}
		return matrix;
	}

	private static void writeNumVector(String fname, double[] vector) {
		try {
			// setup the writers
			File filePtr = new File(fname);
			FileWriter fileWriter = new FileWriter(filePtr);
			BufferedWriter bufWriter = new BufferedWriter(fileWriter);
			PrintWriter prnWriter = new PrintWriter(bufWriter);
			// loop over all elements in the vector
			for (int i = 0; i < vector.length; i++) {
				prnWriter.println(decFmt.format(vector[i]));
			}
			// close the writers
			prnWriter.close();
			bufWriter.close();
			fileWriter.close();
		} catch (Exception e) {
			System.out.println("Caught exception " + e + " while writing file "
					+ fname);
			System.out.flush();
		}
	}

}
