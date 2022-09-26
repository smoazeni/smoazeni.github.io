package EnergySequential; //This is the main class in the package. It includes main().

import java.util.*; //Random is defined in the "java.util" library package
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.lang.Math;
//import ilog.concert.*;
//import ilog.cplex.*;
//import java.io.*;

public class SequentialAlgoCode {
	private static int numberOfSimulations = 100;
	private static int numberOfPeriods = 4;
	private static int numberOfBasisFunctions = 1;
	private static double deltaT = 1;
	private static double currentDemand = 100000;
	private static double currentPrice = 100;
	private static double etaD = 0.85;
	private static double etaC = 0.75;
	private static double Rcap = 40000;
	private static double DeltaRC = 1;
	private static double DeltaRD = 1;
	private static double RInitial = 0.5;
	static double[][] simulatedPrices = new double[numberOfSimulations][numberOfPeriods];
	static double[][] simulatedDemands = new double[numberOfSimulations][numberOfPeriods];
	static double[][] simulatedWinds = new double[numberOfSimulations][numberOfPeriods];

	// Begin Generating Winds
	public static double[] getWindsRand() { // It generates a sample path
		double sigmaEps = 0.4020; // Standard Devaision
		double phi1 = 0.7633;
		double expectedSqrtWt = 1.4781;
		double[] yWind = new double[numberOfPeriods];
		double[] wWind = new double[numberOfPeriods];
		double[] eWind = new double[numberOfPeriods];
		wWind[0] = 4;
		yWind[0] = Math.sqrt(wWind[0]) - expectedSqrtWt;
		eWind[0] = 0.5 * 0.45 * 1.225 * (Math.PI * Math.pow(50, 2))
				* (Math.pow(wWind[0], 3)) * deltaT;
		Random rand = new Random(System.nanoTime()); /*
													 * An instance of this class
													 * is used to generate a
													 * stream of pseudorandom
													 * numbers.
													 */
		for (int i = 1; i < numberOfPeriods; i++) {
			yWind[i] = phi1 * yWind[i - 1] + sigmaEps * rand.nextGaussian();
			wWind[i] = Math.pow(yWind[i] + expectedSqrtWt, 2);
			eWind[i] = 0.5 * 0.45 * 1.225 * (Math.PI * Math.pow(50, 2))
					* (Math.pow(wWind[i], 3)) * deltaT;
		}
		return eWind;
	}

	// End Generating Winds
	// Begin Generating Prices------------------------------------------------
	public static double[] getPricesRand() {
		double phiPrice = 1;// 1;
		double sigmaPrice = 40;// 8;
		double[] dailyPrice = new double[numberOfPeriods];
		double[] hourlyPrice = new double[numberOfPeriods];
		Random rand1 = new Random(System.nanoTime());
		double[] dsPrice = new double[numberOfPeriods];
		dsPrice[0] = currentPrice;
		double[] dPrice = new double[numberOfPeriods];
		dPrice[0] = dsPrice[0] + 30 + 180;
		for (int ellt = 1; ellt < numberOfPeriods; ellt++) {
			if (ellt >= 1 & ellt <= 120) {
				dailyPrice[ellt] = 180;
			} else if (ellt > 120 & ellt <= 168) {
				dailyPrice[ellt] = 50;
			}
			// double hour_of_day = ellt - 24 * Math.floor(ellt / 24);
			/*
			if (0 <= hour_of_day & hour_of_day <= 5) {
				hourlyPrice[ellt] = 30;
			} else if (6 <= hour_of_day & hour_of_day <= 11) {
				hourlyPrice[ellt] = 40;
			} else if (12 <= hour_of_day & hour_of_day <= 17) {
				hourlyPrice[ellt] = 60;
			} else if (18 <= hour_of_day & hour_of_day <= 23) {
				hourlyPrice[ellt] = 100;
			}
			*/
			
				hourlyPrice[0] = 30;
				hourlyPrice[1] = 40;
				hourlyPrice[2] = 60;
				hourlyPrice[3] = 100;
			
			
			
			
			if (hourlyPrice[ellt] == 0) {
				System.out.println("Ops! Zero Demand!");
			}
			dsPrice[ellt] = phiPrice * dsPrice[ellt - 1] + sigmaPrice
					* rand1.nextGaussian();
			dPrice[ellt] = dsPrice[ellt] + dailyPrice[ellt]+hourlyPrice[ellt];
			if (dPrice[ellt] < 0) {
				dPrice[ellt] = 0;
			}
		}

		/*---------Mean Reversion Jump Diffussion Price Model----------------
		double seasonalPrice=0;
		double lambdaP=0.3;
		double muP=4;
		double sigmaP=0.6;
		double sigmaJ=0.3; // Standard Devaision
		double[] epsJump=new double[numberOfPeriods];
		double pJump=0.06;
		double[] JJump=new double[numberOfPeriods];	
		double helpVarP;
		double[] PPrice=new double[numberOfPeriods];
		double[] yPrice=new double[numberOfPeriods];
		PPrice[0]=currentPrice;
		yPrice[0]=Math.log(currentPrice);
		Random rand1 = new Random(System.nanoTime());
		Random rand2 = new Random(System.nanoTime());
		Random rand3 = new Random(System.nanoTime());
		for (int i = 1; i < numberOfPeriods; i++) {
			JJump[i]= 0;
			epsJump[i]=sigmaJ*rand1.nextGaussian();   // This message takes no parameters and returns a random number from a normal distribution with mean 0 and standard deviation 1.	
			helpVarP=rand2.nextDouble();  //To generate a random real number uniformly distributed between 0 and 1, use the "nextDouble" message. This message takes no parameters.
			if (helpVarP<pJump){JJump[i]=  epsJump[i];}
			yPrice[i]=yPrice[i-1]+lambdaP*(muP-yPrice[i-1])*deltaT+sigmaP*Math.sqrt(deltaT)*rand3.nextGaussian()+JJump[i];
			PPrice[i]=Math.exp(yPrice[i]+seasonalPrice);
		                                          }
		-------------------------------------------------------*/
		return dPrice;
	}

	// End Generating Prices---------------------------------------------------
	// Begin Generating Demands
	public static double[] getDemandsRand() {
		double phiDemand = 0.96;// 0.96;//0.9636;
		double sigmaDemand = 4500;// 914870;
		double[] dailyDemand = new double[numberOfPeriods];
		double[] hourlyDemand = new double[numberOfPeriods];
		Random rand1 = new Random(System.nanoTime());
		double[] dsDemand = new double[numberOfPeriods];
		dsDemand[0] = currentDemand;
		double[] dDemand = new double[numberOfPeriods];
		dDemand[0] = dsDemand[0] + 30000 + 90000;
		for (int ellt = 1; ellt < numberOfPeriods; ellt++) {
			if (ellt >= 1 & ellt <= 120) {
				dailyDemand[ellt] = 90000;
			} else if (ellt > 120 & ellt <= 168) {
				dailyDemand[ellt] = 30000;
			}
			//double hour_of_day = ellt - 24 * Math.floor(ellt / 24);
			/*
			if (0 <= hour_of_day & hour_of_day <= 5) {
				hourlyDemand[ellt] = 20000;
			} else if (6 <= hour_of_day & hour_of_day <= 11) {
				hourlyDemand[ellt] = 30000;
			} else if (12 <= hour_of_day & hour_of_day <= 17) {
				hourlyDemand[ellt] = 50000;
			} else if (18 <= hour_of_day & hour_of_day <= 23) {
				hourlyDemand[ellt] = 90000;
			}
			*/
			
		
				hourlyDemand[0] = 20000;
		
				hourlyDemand[1] = 30000;
		
				hourlyDemand[2] = 50000;
		
				hourlyDemand[3] = 90000;
		
			
			
			if (hourlyDemand[ellt] == 0) {
				System.out.println("Ops! Zero Demand!");
			}
			dsDemand[ellt] = phiDemand * dsDemand[ellt - 1] + sigmaDemand
					* rand1.nextGaussian();
			dDemand[ellt] = dsDemand[ellt] + dailyDemand[ellt]
					+ hourlyDemand[ellt];
			if (dDemand[ellt] < 0) {
				dDemand[ellt] = 0;
			}
		}

		return dDemand;
	}

	// End Generating Demands
	// Begin main method
	public static void main(String args[]) throws Exception{
		long start = System.nanoTime(); 
		for (int ellj = 0; ellj < numberOfSimulations; ellj++) {
			simulatedPrices[ellj] = getPricesRand();
			simulatedDemands[ellj] = getDemandsRand();
			simulatedWinds[ellj] = getWindsRand();
		}
		double secs1 = 1.0 * (System.nanoTime() - start) / 1000000000.0;
		
		try {
			PrintWriter writer = new PrintWriter("MySequentialTxtFiles/MySimulated_Input_Data.txt");
			
			writer.println("Data Generation Time: "+ String.format("%8.10f secs!", secs1));
			writer.flush();
			// --------------------------Print Prices---------------------
			writer.println("---------------------------------------------");
			writer.println("Prices over numberOfPeriods Periods:");
			for (int ellM = 0; ellM < numberOfSimulations; ellM++) {
				writer.println("Simulatin Path:  " + ellM);
				for (int ellt = 0; ellt < numberOfPeriods; ellt++) {
					writer.println(simulatedPrices[ellM][ellt] + "");
				}
			}
			// --------------------------Print Demand-----------------------
			writer.println("---------------------------------------------");
			writer.println("Demands over numberOfPeriods Periods:");
			for (int ellM = 0; ellM < numberOfSimulations; ellM++) {
				writer.println("Simulatin Path:  " + ellM);
				for (int ellt = 0; ellt < numberOfPeriods; ellt++) {
					writer.println(simulatedDemands[ellM][ellt] + "");
				}
			}
			// --------------------------Print Wind Energy-------------------
			writer.println("---------------------------------------------");
			writer.println("Winds over numberOfPeriods Periods:");
			for (int ellM = 0; ellM < numberOfSimulations; ellM++) {
				writer.println("Simulatin Path:   " + ellM);
				for (int ellt = 0; ellt < numberOfPeriods; ellt++) {
					writer.println(simulatedWinds[ellM][ellt] + "");
				}
			}
			// --------------------------------------------------------------
			double[] InitialSol = new double[numberOfBasisFunctions
					* (numberOfPeriods - 1)];
			writer.println("---------------------------------------------");
			writer.println("Initial Theta Point:");
			writer.flush();
			for (int ellT = 0; ellT < numberOfPeriods - 1; ellT++) {
				for (int ellK = 0; ellK < numberOfBasisFunctions; ellK++) {
					InitialSol[ellT * numberOfBasisFunctions + ellK] = 1.5;
					writer.println(""
							+ InitialSol[ellT * numberOfBasisFunctions + ellK]);
					writer.flush();
				}
			}

			long startPattern = System.nanoTime(); 
			PatternSearchSolver PSS = new PatternSearchSolver(InitialSol,
					simulatedPrices, simulatedDemands, simulatedWinds, etaD,
					etaC, Rcap, DeltaRC, DeltaRD, RInitial,
					numberOfBasisFunctions);
			double[] Arg_opt = PSS.getarg_opt();
			if (Arg_opt != null) {
				writer.println("---------------------------------------------");
				writer.println("Final Optimal Solution:");
				writer.flush();
				for (int ell = 0; ell < InitialSol.length; ell++) {
					writer.println(Arg_opt[ell] + "");
					writer.flush();
				}
			} else {
				writer.println("Output optimal solution is null!");
				writer.flush();
			}
			double secs = 1.0 * (System.nanoTime() - startPattern) / 1000000000.0;
			writer.println("Done in: "+ String.format("%8.10f secs!", secs));
			writer.flush();
			writer.close();
		} catch (FileNotFoundException e) {
			System.out.println("Caught exception " + e);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
}