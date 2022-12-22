package EnergySequential;

//import ilog.concert.*;
//import ilog.cplex.*;
//import java.io.*;
//import java.util.*;   //Random is defined in the "java.util" library package, Also Vector is defined here import java.util.Vector;
import java.lang.Math;
//import mpi.*;
//import matlabcontrol.*;
import java.io.PrintWriter;
import java.io.FileNotFoundException;

public class PatternSearchSolver {

	private double[] arg_opt = null; // output

	public PatternSearchSolver(double[] arg0, double[][] PMat, double[][] DMat,
			double[][] EMat, double etaD, double etaC, double Rcap,
			double DeltaRC, double DeltaRD, double RInitial, int K)
			throws Exception {
		
		int n = arg0.length;
		int M = PMat.length;
		int T = PMat[0].length;
		double[] argk = new double[n];
		double[] argk1Mat = new double[n];
		double[] argk2Mat = new double[n];
		double ExpansionParameter = 2; // a default value is 2
		double ContractionParameter = 0.5;
		double[] Gamma0 = new double[n]; // initial step length (Note that this
											// value should be chosen very
		// carefully if the problem is unconstrained. Because with such a
		// fixed step size and fixed direction, the arriving point may be
		// infeasible.)
		for (int elli = 0; elli < n; elli++) {
			Gamma0[elli] = 0.5;
		}
		double argError = 1;
		double[] Gammak = new double[n];
		double SufficientDecreaseFunction;
		int IterationCounter;
		int[] TempVec = new int[n];
		// -----------------------Computing Expected Values from
		// Simulations--------------------
		double[] expectedPrices = new double[T + 1];
		double[] expectedDemands = new double[T + 1];
		double[] expectedWinds = new double[T + 1];
		double[] sumvecPrice = new double[T];
		double[] sumvecDemand = new double[T];
		double[] sumvecEnergy = new double[T];
		for (int ellt = 0; ellt < T; ellt++) {
			sumvecPrice[ellt] = 0;
			sumvecDemand[ellt] = 0;
			sumvecEnergy[ellt] = 0;
			for (int ellM = 0; ellM < M; ellM++) {
				sumvecPrice[ellt] = sumvecPrice[ellt] + PMat[ellM][ellt];
				sumvecDemand[ellt] = sumvecDemand[ellt] + DMat[ellM][ellt];
				sumvecEnergy[ellt] = sumvecEnergy[ellt] + EMat[ellM][ellt];
			}
			expectedPrices[ellt] = sumvecPrice[ellt] / M;
			expectedDemands[ellt] = sumvecDemand[ellt] / M;
			expectedWinds[ellt] = sumvecEnergy[ellt] / M;
		}
		expectedPrices[T] = expectedPrices[T - 1];
		expectedDemands[T] = expectedDemands[T - 1];
		expectedWinds[T] = expectedWinds[T - 1];
		
        //when Phi=0
		/*for (int ellt = 0; ellt < T; ellt++) {
			expectedWinds[ellt] = (0.5 * 0.45 * 1.225 * (Math.PI * 50 * 50) * (Math
					.pow(1.4781, 6)
					+ 15
					* (Math.pow(1.4781, 4))
					* (Math.pow(0.4020, 2))
					+ 45
					* (Math.pow(1.4781, 2))
					* (Math.pow(0.4020, 4)) + 15 * (Math.pow(0.4020, 6))));
		}
		expectedPrices[T] = expectedPrices[T - 1];
		expectedWinds[T] = expectedWinds[T - 1];*/

		try {

			PrintWriter writerPP = new PrintWriter("MySequentialTxtFiles/MySequential_Pattern.txt");
			writerPP.println("---------------------------------------------");
			writerPP.println("Expected Prices over T Periods: ");
			for (int ellt = 0; ellt < T + 1; ellt++) {
				writerPP .println(expectedPrices[ellt] + "");
			}
			writerPP .println("Expected Demands over T Periods: ");
			for (int ellt = 0; ellt < T + 1; ellt++) {
				writerPP .println(expectedDemands[ellt] + "");
			}
			writerPP .println("Expected Energy over T Periods: ");
			for (int ellt = 0; ellt < T + 1; ellt++) {
				writerPP .println(expectedWinds[ellt] + "");
			}
			writerPP .flush();

			for (int elli = 0; elli < n; elli++) {
				Gammak[elli] = Gamma0[elli];
			}
			double[] PMatVec = new double[M * T];
			double[] DMatVec = new double[M * T];
			double[] EMatVec = new double[M * T];
			for (int ellt = 0; ellt < T; ellt++) {
				for (int ellm = 0; ellm < M; ellm++) {
					PMatVec[M * ellt + ellm] = PMat[ellm][ellt];
					DMatVec[M * ellt + ellm] = DMat[ellm][ellt];
					EMatVec[M * ellt + ellm] = EMat[ellm][ellt];
				}
			}
			// -----------------Initialization---------------------------
			IterationCounter = 0;
			for (int elli = 0; elli < n; elli++) {
				argk[elli] = arg0[elli];
			}
			// ----------------------------------------------------------
			// ------------Using Matlab Control Interface----------------

		//	try {
				/*MatlabProxyFactory factory = new MatlabProxyFactory();
				MatlabProxy proxy = factory.getProxy();

				proxy.setVariable("argk", argk);
				proxy.setVariable("PMatVec", PMatVec);
				proxy.setVariable("DMatVec", DMatVec);
				proxy.setVariable("EMatVec", EMatVec);
				proxy.setVariable("ExpectedPrices", expectedPrices);
				proxy.setVariable("ExpectedDemands", expectedDemands);
				proxy.setVariable("ExpectedWind", expectedWinds);
				proxy.setVariable("M", M);
				proxy.setVariable("T", T);
				proxy.setVariable("etaD", etaD);
				proxy.setVariable("etaC", etaC);
				proxy.setVariable("Rcap", Rcap);
				proxy.setVariable("DeltaRC", DeltaRC);
				proxy.setVariable("DeltaRD", DeltaRD);
				proxy.setVariable("RInitial", RInitial);
				proxy.eval("f_argk=f(argk,  PMatVec, DMatVec, EMatVec, M, T, ExpectedPrices,ExpectedDemands,ExpectedWind, etaD,etaC,Rcap,DeltaRC,DeltaRD,RInitial )");
				double fmatlab_valuek = ((double[]) proxy.getVariable("f_argk"))[0];
				writerPP .println("Objective Value at current point from matlab: "
						+ fmatlab_valuek);
                */

				double Norm_Gamma = 0.5;

				while ((argError > 10e-3 || Norm_Gamma > 10e-4) && (15 > IterationCounter)) {
					writerPP .println("Iteration Number: " + IterationCounter);
					writerPP .flush();
					writerPP .println("Error:  " + argError);
					writerPP .flush();
					writerPP .println("Norm Gamma:  " + Norm_Gamma);
					writerPP .flush();

					int SuccessAtPointK = 0;
					int i = 0;
					long startf = System.nanoTime(); 
					fBlackBox FValueK = new fBlackBox(argk, PMat, DMat,
							EMat, M, T, expectedPrices, expectedDemands,
							expectedWinds, etaD, etaC, Rcap, DeltaRC,
							DeltaRD, RInitial, K);
					double f_valuek = FValueK.getObjective_velue();
					double secsfvalue = 1.0 * (System.nanoTime() - startf) / 1000000000.0;
					writerPP.println("Evaluate the function f: "+ String.format("%8.10f secs!", secsfvalue));
					writerPP.println("Objective Value at current point: "
							+ f_valuek);
					writerPP.flush();
					
					while (i < n && SuccessAtPointK == 0) {
						writerPP .println("Dimension:  " + i);
						writerPP .flush();
						writerPP .println("Gamma_i:  " + Gammak[i]);
						writerPP .flush();
						// DMat=[eye(n,n) -eye(n,n) ones(n,1)]; %initial
						// direction set
						for (int ell = 0; ell < n; ell++) {
							TempVec[ell] = 0;
						}
						TempVec[i] = 1;
						writerPP .println("Vector of Directions");
						writerPP .flush();
						for (int ell = 0; ell < n; ell++) {
							writerPP .println(TempVec[ell] + "");
						}
						writerPP .println("Current Point");
						writerPP .flush();
						for (int ell = 0; ell < n; ell++) {
							writerPP .println(argk[ell] + "");
						}
						Norm_Gamma = Math.max(Norm_Gamma, Gammak[i]);
						for (int ell = 0; ell < n; ell++) {
							argk1Mat[ell] = argk[ell] + Gammak[i]
									* TempVec[ell];
							if (argk1Mat[ell] < 0) {
								argk1Mat[ell] = 0;
							}
							if (argk1Mat[ell] > 2) {
								argk1Mat[ell] = 2;
							}
							argk2Mat[ell] = argk[ell] - Gammak[i]
									* TempVec[ell];
							if (argk2Mat[ell] < 0) {
								argk2Mat[ell] = 0;
							}
							if (argk2Mat[ell] > 2) {
								argk2Mat[ell] = 2;
							}
						}
						writerPP .println("Positive Direction Point");
						for (int ell = 0; ell < n; ell++) {
							writerPP .println(argk1Mat[ell] + "");
							writerPP .flush();
						}
						writerPP .println("Negative Direction Point");
						for (int ell = 0; ell < n; ell++) {
							writerPP .println(argk2Mat[ell] + "");
							writerPP .flush();
						}

						/*
						 * proxy.setVariable("argk", argk);
						 * proxy.setVariable("argk1Mat", argk1Mat);
						 * proxy.setVariable("argk2Mat", argk2Mat);
						 * proxy.setVariable("PMatVec", PMatVec);
						 * proxy.setVariable("DMatVec", DMatVec);
						 * proxy.setVariable("EMatVec", EMatVec);
						 * proxy.setVariable("ExpectedPrices", expectedPrices);
						 * proxy.setVariable("ExpectedDemands",
						 * expectedDemands); proxy.setVariable("ExpectedWind",
						 * expectedWinds); proxy.setVariable("M", M);
						 * proxy.setVariable("T", T);
						 * proxy.setVariable("etaD",etaD);
						 * proxy.setVariable("etaC",etaC);
						 * proxy.setVariable("Rcap",Rcap );
						 * proxy.setVariable("DeltaRC",DeltaRC );
						 * proxy.setVariable("DeltaRD",DeltaRD );
						 * proxy.setVariable("RInitial",RInitial); proxy.eval(
						 * "f_argk=f(argk,  PMatVec, DMatVec, EMatVec, M, T, ExpectedPrices,ExpectedDemands,ExpectedWind, etaD,etaC,Rcap,DeltaRC,DeltaRD,RInitial )"
						 * ); proxy.eval(
						 * "f_argk1 = f(argk1Mat, PMatVec, DMatVec, EMatVec, M, T, ExpectedPrices,ExpectedDemands,ExpectedWind, etaD,etaC,Rcap,DeltaRC,DeltaRD,RInitial )"
						 * ); proxy.eval(
						 * "f_argk2 = f(argk2Mat, PMatVec, DMatVec, EMatVec, M, T, ExpectedPrices,ExpectedDemands,ExpectedWind, etaD,etaC,Rcap,DeltaRC,DeltaRD,RInitial )"
						 * ); double fmatlab_valuek = ((double[])
						 * proxy.getVariable("f_argk"))[0]; double
						 * fmatlab_valuek1 = ((double[])
						 * proxy.getVariable("f_argk1"))[0]; double
						 * fmatlab_valuek2 = ((double[])
						 * proxy.getVariable("f_argk2"))[0];
						 * writerPP .println("Objective Value at current point: "
						 * +fmatlab_valuek);
						 * writerPP .println("Objective Value at current point: "
						 * +fmatlab_valuek1);
						 * writerPP .println("Objective Value at current point: "
						 * +fmatlab_valuek2);
						 */

						long startfk1 = System.nanoTime(); 
						fBlackBox FValueK1 = new fBlackBox(argk1Mat, PMat,
								DMat, EMat, M, T, expectedPrices,
								expectedDemands, expectedWinds, etaD, etaC,
								Rcap, DeltaRC, DeltaRD, RInitial, K);
						double f_valuek1 = FValueK1.getObjective_velue();
						double secsfvaluek1 = 1.0 * (System.nanoTime() - startfk1) / 1000000000.0;
						
						writerPP.println("Evaluate the function fk1: "+ String.format("%8.10f secs!", secsfvaluek1));
						
						writerPP.println("Difference Objective Value at point k1: "
								+ (f_valuek1));
						writerPP.flush();
						
						long startfk2 = System.nanoTime(); 
						fBlackBox FValueK2 = new fBlackBox(argk2Mat, PMat,
								DMat, EMat, M, T, expectedPrices,
								expectedDemands, expectedWinds, etaD, etaC,
								Rcap, DeltaRC, DeltaRD, RInitial, K);
						double f_valuek2 = FValueK2.getObjective_velue();
						double secsfvaluek2 = 1.0 * (System.nanoTime() - startfk2) / 1000000000.0;
						writerPP.println("Evaluate the function fk2: "+ String.format("%8.10f secs!", secsfvaluek2));
						
						
						
						writerPP.println("Difference Objective Value at point k2: "
								+ (f_valuek2));
						writerPP.flush();


						SufficientDecreaseFunction = 0.1;// (Math.abs(f_valuek))*(Math.pow(Gammak,1.5));
						writerPP.println("***Sufficient Decrease Value: "+ SufficientDecreaseFunction);
						writerPP.flush();
						if ((f_valuek1 <= (f_valuek - SufficientDecreaseFunction))
								&& f_valuek1 <= f_valuek2) {
							argError = 0;
							for (int j = 0; j < n; j++) {
								argError += (argk1Mat[j] - argk[j])
										* (argk1Mat[j] - argk[j]);
								argk[j] = argk1Mat[j];
							}
							writerPP .println("***Success at positive direction*** "
									+ i);
							SuccessAtPointK = 1;
							Gammak[i] = ExpansionParameter * Gammak[i];
						}
						if ((f_valuek2 <= (f_valuek - SufficientDecreaseFunction))
								&& f_valuek1 >= f_valuek2) {
							argError = 0;
							for (int j = 0; j < n; j++) {
								argError += (argk2Mat[j] - argk[j])
										* (argk2Mat[j] - argk[j]);
								argk[j] = argk2Mat[j];
							}
							writerPP .println("***Success at negative direction*** "
									+ i);
							SuccessAtPointK = 1;
							Gammak[i] = ExpansionParameter * Gammak[i];
						}
						if (f_valuek1 > (f_valuek - SufficientDecreaseFunction)
								&& (f_valuek2 > (f_valuek - SufficientDecreaseFunction))) {
							Gammak[i] = ContractionParameter * Gammak[i];
						}

						i = i + 1;
					}
					// if (SuccessAtPointK==0 && i==n){for (int
					// elli=0;elli<n;elli++){Gammak[i]=ContractionParameter*Gammak[i];}}
					IterationCounter = IterationCounter + 1;
					writerPP .println("Iteration Number: " + IterationCounter);
				}
				arg_opt = argk;

			/*} catch (MatlabConnectionException e) {
				writerPP.println("Caught exception " + e);
			} catch (MatlabInvocationException e) {
				writerPP.println("Caught exception " + e);
			}*/
			writerPP.close();
		} catch (FileNotFoundException e) {
			System.out.println("Caught exception " + e);
		}
		
	}

	public double[] getarg_opt() {
		return arg_opt;
	}
}
