package EnergySequential; //This is the main class in the package. It includes main().

import java.util.*; //Random is defined in the "java.util" library package
import mpi.*;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.lang.Math;
import java.lang.reflect.Field;
import java.net.InetAddress;
import java.net.UnknownHostException;
//import ilog.concert.*;
//import ilog.cplex.*;
//import java.io.*;

public class ParallelAlgoCode {
	private static int numberOfSimulations = 1;
	private static int numberOfPeriods = 6;// 168;
	private static int numberOfBasisFunctions = 1;
	private static int n = numberOfBasisFunctions * (numberOfPeriods - 1);
	private static double deltaT = 1;
	private static double currentDemand = 100000;
	private static double currentPrice = 100;
	private static double etaD = 1;// 0.85;
	private static double etaC = 1;// 0.75;
	private static double Rcap = 40000;
	private static double DeltaRC = 1;
	private static double DeltaRD = 1;
	private static double RInitial = 0.5;
	static double[][] simulatedPrices = new double[numberOfSimulations][numberOfPeriods];
	static double[][] simulatedDemands = new double[numberOfSimulations][numberOfPeriods];
	static double[][] simulatedWinds = new double[numberOfSimulations][numberOfPeriods];
	public static final int COMPUTE = 0;
	public static final int END = -1;
	public static final int N = 10; // number of nodes
	public static double[] end = new double[1];
	private static double[] argk = new double[n];
	private static double[] InitialSol = new double[n];

	// Begin Generating Winds
	public static double[] getWindsRand() { // It generates a sample path
		double sigmaEps = 0.4020; // Standard Devaision
		double phi1 = 0;// 0.7633;
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
		double phiPrice = 0;// 1;
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
			double hour_of_day = ellt - 24 * Math.floor(ellt / 24);
			if (0 <= hour_of_day & hour_of_day <= 5) {
				hourlyPrice[ellt] = 30;
			} else if (6 <= hour_of_day & hour_of_day <= 11) {
				hourlyPrice[ellt] = 40;
			} else if (12 <= hour_of_day & hour_of_day <= 17) {
				hourlyPrice[ellt] = 60;
			} else if (18 <= hour_of_day & hour_of_day <= 23) {
				hourlyPrice[ellt] = 100;
			}
			if (hourlyPrice[ellt] == 0) {
				System.out.println("Ops! Zero Demand!");
			}
			dsPrice[ellt] = phiPrice * dsPrice[ellt - 1] + sigmaPrice
					* rand1.nextGaussian();
			dPrice[ellt] = dsPrice[ellt] + 250;// dailyPrice[ellt]+hourlyPrice[ellt];
			if (dPrice[ellt] < 0) {
				dPrice[ellt] = 0;
			}
		}
		return dPrice;
	}

	// End Generating Prices---------------------------------------------------
	// Begin Generating Demands
	public static double[] getDemandsRand() {
		double phiDemand = 0;// 0.96;//0.9636;
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
			double hour_of_day = ellt - 24 * Math.floor(ellt / 24);
			if (0 <= hour_of_day & hour_of_day <= 5) {
				hourlyDemand[ellt] = 20000;
			} else if (6 <= hour_of_day & hour_of_day <= 11) {
				hourlyDemand[ellt] = 30000;
			} else if (12 <= hour_of_day & hour_of_day <= 17) {
				hourlyDemand[ellt] = 50000;
			} else if (18 <= hour_of_day & hour_of_day <= 23) {
				hourlyDemand[ellt] = 90000;
			}
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
	public static void main(String args[]) throws Exception {
		MPI.Init(args);
		/*
		 * identify the host where this object is running try { String hostName
		 * = null; InetAddress addr = InetAddress.getLocalHost(); hostName =
		 * addr.getHostName(); if (hostName.contains("moazeni"))
		 * System.setProperty( "java.library.path",
		 * "C:/Program Files/ibm/ILOG/CPLEX_Studio124/cplex/bin/x64_win64" );
		 * else System.setProperty( "java.library.path",
		 * "/opt/ilog/cplex/bin/x86-64_sles10_4.1" ); } catch
		 * (UnknownHostException e) {
		 * System.out.println("Could not identify host machine.  Exiting!");
		 * System.exit(-1); }
		 * 
		 * Field fieldSysPath = ClassLoader.class.getDeclaredField( "sys_paths"
		 * ); fieldSysPath.setAccessible( true ); fieldSysPath.set( null, null
		 * ); //
		 */

		long start = System.nanoTime();

		int me = MPI.COMM_WORLD.Rank();
		int size = MPI.COMM_WORLD.Size();
		// ------------------------------------------------------
		// creating a stop checker
		File file = new File("DFO/MyParalleTxtFiles/control.txt");
		if (!file.exists()) {
			file.createNewFile();
		}
		StopChecker checker = new StopChecker(file);
		checker.start();
		// ------------------------------------------------------

		try {
			if (me == 0) {

				PrintWriter writerS = new PrintWriter(
						"DFO/MyParallelTxtFiles/MyParallel_Simulated_Input_Data.txt");
				PrintWriter writerP = new PrintWriter(
						"DFO/MyParallelTxtFiles/MyParallel_PatternSearchPerformance.txt");

				// identify the host where this object is running
				String hostName = null;
				try {
					InetAddress addr = InetAddress.getLocalHost();
					hostName = addr.getHostName();
					writerP.println("> running on " + hostName);
				} catch (UnknownHostException e) {
					writerP.println("Could not identify host machine.  Exiting!");
					System.exit(-1);
				}

				for (int ellj = 0; ellj < numberOfSimulations; ellj++) {
					simulatedPrices[ellj] = getPricesRand();
					simulatedDemands[ellj] = getDemandsRand();
					simulatedWinds[ellj] = getWindsRand();
				}

				double secs1 = 1.0 * (System.nanoTime() - start) / 1000000000.0;
				writerS.println("Data Generation Time: "
						+ String.format("%8.10f secs!", secs1));
				writerS.flush();
				// --------------------------Print Prices---------------------
				writerS.println("---------------------------------------------");
				writerS.println("Prices over numberOfPeriods Periods:");
				for (int ellM = 0; ellM < numberOfSimulations; ellM++) {
					writerS.println("Simulatin Path:  " + ellM);
					for (int ellt = 0; ellt < numberOfPeriods; ellt++) {
						writerS.println(simulatedPrices[ellM][ellt] + "");
					}
				}
				// --------------------------Print Demand-----------------------
				writerS.println("---------------------------------------------");
				writerS.println("Demands over numberOfPeriods Periods:");
				for (int ellM = 0; ellM < numberOfSimulations; ellM++) {
					writerS.println("Simulatin Path:  " + ellM);
					for (int ellt = 0; ellt < numberOfPeriods; ellt++) {
						writerS.println(simulatedDemands[ellM][ellt] + "");
					}
				}
				// --------------------------Print Wind
				// Energy-------------------
				writerS.println("---------------------------------------------");
				writerS.println("Winds over numberOfPeriods Periods:");
				for (int ellM = 0; ellM < numberOfSimulations; ellM++) {
					writerS.println("Simulatin Path:   " + ellM);
					for (int ellt = 0; ellt < numberOfPeriods; ellt++) {
						writerS.println(simulatedWinds[ellM][ellt] + "");
					}
				}

				double[] expectedPrices = new double[numberOfPeriods + 1];
				double[] expectedDemands = new double[numberOfPeriods + 1];
				double[] expectedWinds = new double[numberOfPeriods + 1];
				double[] sumvecPrice = new double[numberOfPeriods];
				double[] sumvecDemand = new double[numberOfPeriods];
				double[] sumvecEnergy = new double[numberOfPeriods];
				for (int ellt = 0; ellt < numberOfPeriods; ellt++) {
					sumvecPrice[ellt] = 0;
					sumvecDemand[ellt] = 0;
					sumvecEnergy[ellt] = 0;
					for (int ellM = 0; ellM < numberOfSimulations; ellM++) {
						sumvecPrice[ellt] = sumvecPrice[ellt]
								+ simulatedPrices[ellM][ellt];
						sumvecDemand[ellt] = sumvecDemand[ellt]
								+ simulatedDemands[ellM][ellt];
						sumvecEnergy[ellt] = sumvecEnergy[ellt]
								+ simulatedWinds[ellM][ellt];
					}
					expectedPrices[ellt] = sumvecPrice[ellt]
							/ numberOfSimulations;
					expectedDemands[ellt] = sumvecDemand[ellt]
							/ numberOfSimulations;
					expectedWinds[ellt] = sumvecEnergy[ellt]
							/ numberOfSimulations;
				}
				expectedPrices[numberOfPeriods] = expectedPrices[numberOfPeriods - 1];
				expectedDemands[numberOfPeriods] = expectedDemands[numberOfPeriods - 1];
				expectedWinds[numberOfPeriods] = expectedWinds[numberOfPeriods - 1];

				for (int ellt = 0; ellt < numberOfPeriods; ellt++) {
					expectedPrices[ellt] = 250;
					expectedWinds[ellt] = (0.5 * 0.45 * 1.225 * (Math.PI * 50 * 50) * (Math
							.pow(1.4781, 6)
							+ 15
							* (Math.pow(1.4781, 4))
							* (Math.pow(0.4020, 2))
							+ 45
							* (Math.pow(1.4781, 2)) * (Math.pow(0.4020, 4)) + 15 * (Math
							.pow(0.4020, 6))));
				}
				expectedPrices[numberOfPeriods] = expectedPrices[numberOfPeriods - 1];
				expectedWinds[numberOfPeriods] = expectedWinds[numberOfPeriods - 1];

				writerS.println("---------------------------------------------");
				writerS.println("Expected Prices over numberOfPeriods Periods: ");
				for (int ellt = 0; ellt < numberOfPeriods + 1; ellt++) {
					writerS.println(expectedPrices[ellt] + "");
				}
				writerS.println("Expected Demands over numberOfPeriods Periods: ");
				for (int ellt = 0; ellt < numberOfPeriods + 1; ellt++) {
					writerS.println(expectedDemands[ellt] + "");
				}
				writerS.println("Expected Energy over numberOfPeriods Periods: ");
				for (int ellt = 0; ellt < numberOfPeriods + 1; ellt++) {
					writerS.println(expectedWinds[ellt] + "");
				}
				writerS.flush();

				for (int inode = 1; inode < size; inode++) {
					Request[] sendRequests = new Request[3 * numberOfSimulations + 3];
					int counter = 0;
					for (int ellj = 0; ellj < numberOfSimulations; ellj++) {
						sendRequests[counter] = MPI.COMM_WORLD.Isend(
								simulatedPrices[ellj], 0,
								simulatedPrices[ellj].length, MPI.DOUBLE,
								inode, 3 * ellj);
						counter++;
						sendRequests[counter] = MPI.COMM_WORLD.Isend(
								simulatedDemands[ellj], 0,
								simulatedDemands[ellj].length, MPI.DOUBLE,
								inode, 3 * ellj + 1);
						counter++;
						sendRequests[counter] = MPI.COMM_WORLD.Isend(
								simulatedWinds[ellj], 0,
								simulatedWinds[ellj].length, MPI.DOUBLE, inode,
								3 * ellj + 2);
						counter++;
					}
					sendRequests[counter] = MPI.COMM_WORLD.Isend(
							expectedPrices, 0, expectedPrices.length,
							MPI.DOUBLE, inode, 30);
					counter++;
					sendRequests[counter] = MPI.COMM_WORLD.Isend(
							expectedDemands, 0, expectedDemands.length,
							MPI.DOUBLE, inode, 31);
					counter++;
					sendRequests[counter] = MPI.COMM_WORLD.Isend(expectedWinds,
							0, expectedWinds.length, MPI.DOUBLE, inode, 32);
					counter++;
					Request.Waitall(sendRequests);
					writerP.println("hi " + counter + " node: " + 1
							+ " sendRequests.size: " + sendRequests.length);
					writerP.flush();
				}
				// --------------------------------------------------------------
				writerS.println("---------------------------------------------");
				writerS.println("Initial Theta Point:");
				writerS.flush();
				for (int ellT = 0; ellT < numberOfPeriods - 1; ellT++) {
					for (int ellK = 0; ellK < numberOfBasisFunctions; ellK++) {
						InitialSol[ellT * numberOfBasisFunctions + ellK] = 1.5;
						writerS.println(""
								+ InitialSol[ellT * numberOfBasisFunctions
										+ ellK]);
						writerS.flush();
					}
				}

				writerP.println("1-preparing to sent x. . .");
				writerP.flush();

				// **************************************************************************
				double[] f_valuek1 = new double[size - 1];
				double[] f_valuek2 = new double[size - 1];
				double[] arg_opt = new double[n];
				double[] argk1Mat = new double[n];
				double[] argk2Mat = new double[n];
				double ExpansionParameter = 2; // a default value is 2
				double ContractionParameter = 0.5;
				double[] Gammak = new double[n];
				double SufficientDecreaseFunction;
				int IterationCounter;
				int[] TempVec = new int[n];
				double[] Gamma0 = new double[n];
				double f_valuek = 0;
				writerP.println("2-preparing to sent x. . .");
				writerP.flush();

				for (int elli = 0; elli < n; elli++) {
					Gamma0[elli] = 0.5;
					Gammak[elli] = Gamma0[elli];
				}

				writerP.println("3-preparing to sent x. . .");
				writerP.flush();

				IterationCounter = 0;
				for (int elli = 0; elli < n; elli++) {
					argk[elli] = InitialSol[elli];
				}

				writerP.println("4-preparing to sent x. . .");
				writerP.flush();

				double argError = 1;
				double Norm_Gamma = 0.5;

				while ((argError > 10e-3 || Norm_Gamma > 10e-10)
						&& (100 > IterationCounter)) {
				//	InputStream inputfileHi = new FileInputStream("Hi.txt");
					//System.out.println("file: " + inputfileHi.read());
					//if ((inputfileHi.read()) != -1) {

						writerP.println("Iteration Number: " + IterationCounter);
						writerP.flush();
						writerP.println("Error:  " + argError);
						writerP.flush();
						writerP.println("Norm Gamma:  " + Norm_Gamma);
						writerP.flush();

						writerP.println("Current Point");
						writerP.flush();
						for (int ell = 0; ell < n; ell++) {
							writerP.println(argk[ell] + "");
						}
						writerP.flush();

						try {

							fBlackBox FValueK = new fBlackBox(argk,
									simulatedPrices, simulatedDemands,
									simulatedWinds, numberOfSimulations,
									numberOfPeriods, expectedPrices,
									expectedDemands, expectedWinds, etaD, etaC,
									Rcap, DeltaRC, DeltaRD, RInitial,
									numberOfBasisFunctions);

							writerP.println("getting obj. . .");
							writerP.flush();

							f_valuek = FValueK.getObjective_velue();

							writerP.println("f_valuek in me node: " + f_valuek);
							writerP.flush();

						} catch (java.lang.UnsatisfiedLinkError e) {
							writerP.println("I am in main exception: " + e);
						}

						Request sent[] = new Request[3 * (size - 1)];

						writerP.println("Value of size" + size);
						writerP.flush();
						writerP.println("" + Gammak.length);
						writerP.flush();

						for (int inode = 1; inode < size; inode++) {

							writerP.println("Dimension:  " + inode
									+ "Gamma_i:  " + Gammak[inode - 1]);
							writerP.flush();
							for (int ell = 0; ell < n; ell++) {
								TempVec[ell] = 0;
							}
							TempVec[inode - 1] = 1;
							writerP.println("Direction");
							writerP.flush();
							for (int ell = 0; ell < n; ell++) {
								writerP.println(TempVec[ell] + "");
							}
							Norm_Gamma = Math
									.max(Norm_Gamma, Gammak[inode - 1]);

							for (int ell = 0; ell < n; ell++) {
								argk1Mat[ell] = argk[ell] + Gammak[inode - 1]
										* TempVec[ell];
								if (argk1Mat[ell] < 0) {
									argk1Mat[ell] = 0;
								}
								if (argk1Mat[ell] > 2) {
									argk1Mat[ell] = 2;
								}
								argk2Mat[ell] = argk[ell] - Gammak[inode - 1]
										* TempVec[ell];
								if (argk2Mat[ell] < 0) {
									argk2Mat[ell] = 0;
								}
								if (argk2Mat[ell] > 2) {
									argk2Mat[ell] = 2;
								}
							}

							for (int ell = 0; ell < n; ell++) {
								writerP.println(argk1Mat[ell] + "");
								writerP.flush();
							}

							for (int ell = 0; ell < n; ell++) {
								writerP.println(argk2Mat[ell] + "");
								writerP.flush();
							}

							writerP.println("sent x to node: " + inode
									+ ". . .");
							writerP.flush();

							sent[3 * (inode - 1)] = MPI.COMM_WORLD.Isend(argk,
									0, argk.length, MPI.DOUBLE, inode, -50);
							sent[3 * (inode - 1) + 1] = MPI.COMM_WORLD.Isend(
									argk1Mat, 0, argk1Mat.length, MPI.DOUBLE,
									inode, -100);
							sent[3 * (inode - 1) + 2] = MPI.COMM_WORLD.Isend(
									argk2Mat, 0, argk2Mat.length, MPI.DOUBLE,
									inode, -200);

						}

						writerP.println("sent x");
						writerP.flush();
						// Request.Waitall(sent);

						writerP.println("all sent!");
						writerP.flush();

						SufficientDecreaseFunction = 0.1;// (Math.abs(f_valuek))*(Math.pow(Gammak,1.5));
						Request[] resp = new Request[size - 1];

						int Tag_f = -1;
						double Benchmark_f = f_valuek
								- SufficientDecreaseFunction;
						writerP.println("Current Benchmark is: " + Benchmark_f);
						writerP.flush();

						double[] result_vector_node1 = new double[2];
						resp[0] = MPI.COMM_WORLD
								.Irecv(result_vector_node1, 0,
										result_vector_node1.length, MPI.DOUBLE,
										1, -300);
						double[] result_vector_node2 = new double[2];
						resp[1] = MPI.COMM_WORLD
								.Irecv(result_vector_node2, 0,
										result_vector_node2.length, MPI.DOUBLE,
										2, -300);
						double[] result_vector_node3 = new double[2];
						resp[2] = MPI.COMM_WORLD
								.Irecv(result_vector_node3, 0,
										result_vector_node3.length, MPI.DOUBLE,
										3, -300);
						// Request.Waitany(resp);
						Request.Waitall(resp);
						writerP.println("done by receiving!");
						writerP.flush();
						f_valuek1[0] = result_vector_node1[0];
						f_valuek2[0] = result_vector_node1[1];
						f_valuek1[1] = result_vector_node2[0];
						f_valuek2[1] = result_vector_node2[1];
						f_valuek1[2] = result_vector_node3[0];
						f_valuek2[2] = result_vector_node3[1];

						for (int inode = 1; inode < size; inode++) {
							writerP.println("f_valuek1 received from node "
									+ inode + " is " + f_valuek1[inode - 1]);
							writerP.flush();
							writerP.println("f_valuek2 received from node "
									+ inode + " is " + f_valuek2[inode - 1]);
							writerP.flush();
							if (f_valuek1[inode - 1] < Benchmark_f) {
								Tag_f = inode;
								Benchmark_f = f_valuek1[inode - 1];
							} else if (f_valuek2[inode - 1] < Benchmark_f) {
								Tag_f = inode;
								Benchmark_f = f_valuek2[inode - 1];
							}
						}
						writerP.println("Tag_f: " + Tag_f);
						writerP.flush();
						writerP.println("Benchmark_f: " + Benchmark_f);
						writerP.flush();

						if (Tag_f != -1) {

							for (int ell = 0; ell < n; ell++) {
								TempVec[ell] = 0;
							}
							TempVec[Tag_f - 1] = 1;
							for (int ell = 0; ell < n; ell++) {
								argk1Mat[ell] = argk[ell] + Gammak[Tag_f - 1]
										* TempVec[ell];
								if (argk1Mat[ell] < 0) {
									argk1Mat[ell] = 0;
								}
								if (argk1Mat[ell] > 2) {
									argk1Mat[ell] = 2;
								}
								argk2Mat[ell] = argk[ell] - Gammak[Tag_f - 1]
										* TempVec[ell];
								if (argk2Mat[ell] < 0) {
									argk2Mat[ell] = 0;
								}
								if (argk2Mat[ell] > 2) {
									argk2Mat[ell] = 2;
								}
							}
							if ((f_valuek1[Tag_f - 1] <= (f_valuek - SufficientDecreaseFunction))
									&& f_valuek1[Tag_f - 1] <= f_valuek2[Tag_f - 1]) {
								argError = 0;
								for (int j = 0; j < n; j++) {
									argError += (argk1Mat[j] - argk[j])
											* (argk1Mat[j] - argk[j]);
									argk[j] = argk1Mat[j];
								}
								writerP.println("***Success at positive direction*** "
										+ Tag_f);
								Gammak[Tag_f - 1] = ExpansionParameter
										* Gammak[Tag_f - 1];
							}
							if ((f_valuek2[Tag_f - 1] <= (f_valuek - SufficientDecreaseFunction))
									&& f_valuek1[Tag_f - 1] >= f_valuek2[Tag_f - 1]) {
								argError = 0;
								for (int j = 0; j < n; j++) {
									argError += (argk2Mat[j] - argk[j])
											* (argk2Mat[j] - argk[j]);
									argk[j] = argk2Mat[j];
								}
								writerP.println("***Success at negative direction*** "
										+ Tag_f);
								Gammak[Tag_f - 1] = ExpansionParameter
										* Gammak[Tag_f - 1];
							}
						} else if (Tag_f == -1) {
							for (int elli = 0; elli < n; elli++) {
								Gammak[elli] = ContractionParameter
										* Gammak[elli];
							}
						}

						// try {
						// } catch (java.lang.ArrayIndexOutOfBoundsException e)
						// {
						// System.out.println(" caught an exc. at main node: " +
						// me);
						// System.out.flush();
						// System.exit(-1);
						// }

						IterationCounter = IterationCounter + 1;
						writerP.println("Iteration Counter:  "
								+ IterationCounter);
						writerP.flush();
				/*	} else if (inputfileHi.read() == -1) {
						System.out.println("Aborted!");
						writerS.flush();
						Request freq[] = new Request[size - 1];
						for (int inode = 1; inode < size; inode++) {
							freq[inode - 1] = MPI.COMM_WORLD.Isend(InitialSol,
									0, InitialSol.length, MPI.DOUBLE, inode,
									END);
						}
						Request.Waitall(freq);

						writerS.close();
						writerP.close();
						System.exit(-1);
					}*/

				}// End While

				for (int elln = 0; elln < n; elln++) {
					arg_opt[elln] = argk[elln];
				}

				if (arg_opt != null) {
					writerP.println("---------------------------------------------");
					writerP.println("Final Optimal Solution:");
					writerP.flush();
					for (int ell = 0; ell < n; ell++) {
						writerP.println(arg_opt[ell] + "");
						writerP.flush();
					}
				} else {
					writerP.println("Output optimal solution is null!");
					writerP.flush();
				}

				// ***************************************************************************
				double secs = 1.0 * (System.nanoTime() - start) / 1000000000.0;
				writerS.println("Done Everything Within: "
						+ String.format("%8.10f secs!", secs));
				writerS.flush();
				Request freq[] = new Request[size - 1];
				for (int inode = 1; inode < size; inode++) {
					freq[inode - 1] = MPI.COMM_WORLD.Isend(InitialSol, 0,
							InitialSol.length, MPI.DOUBLE, inode, END);
				}
				Request.Waitall(freq);

				writerS.close();
				writerP.close();

			} else {// ///////////////////////////////////////////////////////////////////////////
				PrintWriter writer_Client = new PrintWriter("DFO/MyParallelTxtFiles/Client_output_"
						+ me + ".txt");
				// identify the host where this object is running
				String hostName = null;
				try {
					InetAddress addr = InetAddress.getLocalHost();
					hostName = addr.getHostName();
					writer_Client.println("> running on " + hostName);
				} catch (UnknownHostException e) {
					writer_Client
							.println("Could not identify host machine.  Exiting!");
					System.exit(-1);
				}

				double[][] simulatedPrices = new double[numberOfSimulations][numberOfPeriods];
				double[][] simulatedDemands = new double[numberOfSimulations][numberOfPeriods];
				double[][] simulatedWinds = new double[numberOfSimulations][numberOfPeriods];
				double[] expectedPrices = new double[numberOfPeriods + 1];
				double[] expectedDemands = new double[numberOfPeriods + 1];
				double[] expectedWinds = new double[numberOfPeriods + 1];
				int counter = 0;
				Request[] respNode = new Request[3 * numberOfSimulations + 3];
				for (int ellj = 0; ellj < numberOfSimulations; ellj++) {
					respNode[counter] = MPI.COMM_WORLD.Irecv(
							simulatedPrices[ellj], 0,
							simulatedPrices[ellj].length, MPI.DOUBLE, 0,
							3 * ellj);
					counter++;
					respNode[counter] = MPI.COMM_WORLD.Irecv(
							simulatedDemands[ellj], 0,
							simulatedDemands[ellj].length, MPI.DOUBLE, 0,
							3 * ellj + 1);
					counter++;
					respNode[counter] = MPI.COMM_WORLD.Irecv(
							simulatedWinds[ellj], 0,
							simulatedWinds[ellj].length, MPI.DOUBLE, 0,
							3 * ellj + 2);
					counter++;
				}
				respNode[counter] = MPI.COMM_WORLD.Irecv(expectedPrices, 0,
						expectedPrices.length, MPI.DOUBLE, 0, 30);
				counter++;
				respNode[counter] = MPI.COMM_WORLD.Irecv(expectedDemands, 0,
						expectedDemands.length, MPI.DOUBLE, 0, 31);
				counter++;
				respNode[counter] = MPI.COMM_WORLD.Irecv(expectedWinds, 0,
						expectedWinds.length, MPI.DOUBLE, 0, 32);

				Request.Waitall(respNode);
				// --------------------------Print Prices---------------------
				writer_Client
						.println("---------------------------------------------");
				writer_Client.println("Prices over numberOfPeriods Periods:");
				for (int ellM = 0; ellM < numberOfSimulations; ellM++) {
					writer_Client.println("Simulatin Path:  " + ellM);
					for (int ellt = 0; ellt < numberOfPeriods; ellt++) {
						writer_Client.println(simulatedPrices[ellM][ellt] + "");
					}
				}
				// --------------------------Print Demand-----------------------
				writer_Client
						.println("---------------------------------------------");
				writer_Client.println("Demands over numberOfPeriods Periods:");
				for (int ellM = 0; ellM < numberOfSimulations; ellM++) {
					writer_Client.println("Simulatin Path:  " + ellM);
					for (int ellt = 0; ellt < numberOfPeriods; ellt++) {
						writer_Client
								.println(simulatedDemands[ellM][ellt] + "");
					}
				}
				// --------------------------Print Wind
				// Energy-------------------
				writer_Client
						.println("---------------------------------------------");
				writer_Client.println("Winds over numberOfPeriods Periods:");
				for (int ellM = 0; ellM < numberOfSimulations; ellM++) {
					writer_Client.println("Simulatin Path:   " + ellM);
					for (int ellt = 0; ellt < numberOfPeriods; ellt++) {
						writer_Client.println(simulatedWinds[ellM][ellt] + "");
					}
				}

				writer_Client
						.println("---------------------------------------------");
				writer_Client
						.println("Expected Prices over numberOfPeriods Periods: ");
				for (int ellt = 0; ellt < numberOfPeriods + 1; ellt++) {
					writer_Client.println(expectedPrices[ellt] + "");
				}
				writer_Client
						.println("Expected Demands over numberOfPeriods Periods: ");
				for (int ellt = 0; ellt < numberOfPeriods + 1; ellt++) {
					writer_Client.println(expectedDemands[ellt] + "");
				}
				writer_Client
						.println("Expected Energy over numberOfPeriods Periods: ");
				for (int ellt = 0; ellt < numberOfPeriods + 1; ellt++) {
					writer_Client.println(expectedWinds[ellt] + "");
				}
				writer_Client
						.println("---------------------------------------------");
				writer_Client.flush();

				while (true) {// //////////////////////////////////////////////////////////////////////////////////
					Status stat = null;
					try {
						double recArr[] = new double[argk.length];
						stat = MPI.COMM_WORLD.Recv(recArr, 0, recArr.length,
								MPI.DOUBLE, 0, MPI.ANY_TAG);

						if (stat.tag != END && stat.tag < 0) {
							writer_Client.println("We are in!");
							writer_Client.flush();

							double[] argk1Mat_client = new double[n];
							double[] argk2Mat_client = new double[n];
							Request[] resv_x = new Request[2];
							resv_x[0] = MPI.COMM_WORLD.Irecv(argk1Mat_client,
									0, argk1Mat_client.length, MPI.DOUBLE, 0,
									-100);
							resv_x[1] = MPI.COMM_WORLD.Irecv(argk2Mat_client,
									0, argk2Mat_client.length, MPI.DOUBLE, 0,
									-200);
							Request.Waitall(resv_x);
							writer_Client.println("argk1Mat_client");
							writer_Client.flush();
							for (int elln = 0; elln < n; elln++) {
								writer_Client.println(argk1Mat_client[elln]
										+ "");
								writer_Client.flush();
							}
							writer_Client.println("argk2Mat_client");
							writer_Client.flush();
							for (int elln = 0; elln < n; elln++) {
								writer_Client.println(argk2Mat_client[elln]
										+ "");
								writer_Client.flush();
							}

							try {

								fBlackBox FValueK1_client = new fBlackBox(
										argk2Mat_client, simulatedPrices,
										simulatedDemands, simulatedWinds,
										numberOfSimulations, numberOfPeriods,
										expectedPrices, expectedDemands,
										expectedWinds, etaD, etaC, Rcap,
										DeltaRC, DeltaRD, RInitial,
										numberOfBasisFunctions);

								double f_valuek1_client = FValueK1_client
										.getObjective_velue();
								writer_Client.println("f_valuek1:  "
										+ f_valuek1_client);
								writer_Client.flush();

								fBlackBox FValueK2_client = new fBlackBox(
										argk2Mat_client, simulatedPrices,
										simulatedDemands, simulatedWinds,
										numberOfSimulations, numberOfPeriods,
										expectedPrices, expectedDemands,
										expectedWinds, etaD, etaC, Rcap,
										DeltaRC, DeltaRD, RInitial,
										numberOfBasisFunctions);

								double f_valuek2_client = FValueK2_client
										.getObjective_velue();

								writer_Client.println("f_valuek2:  "
										+ f_valuek2_client);
								writer_Client.flush();

								Request[] send_optVal = new Request[1];
								double res[] = new double[2];
								res[0] = f_valuek1_client;
								res[1] = f_valuek2_client;

								send_optVal[0] = MPI.COMM_WORLD.Isend(res, 0,
										res.length, MPI.DOUBLE, 0, -300);
								Request.Waitall(send_optVal);

							} catch (java.lang.UnsatisfiedLinkError e) {
								writer_Client.println("I am in client node: "
										+ me + " exception: " + e);
							}

							writer_Client.println("Done by sending!");
							writer_Client.flush();
						}
					} catch (java.lang.ArrayIndexOutOfBoundsException e) {
						System.out.println(" catch in client");
						System.out.flush();
						System.exit(-1);
					}

					if (stat.tag == END) {
						java.net.InetAddress localMachine = java.net.InetAddress
								.getLocalHost();
						String hname = localMachine.getHostName();
						String message = "Message from: "
								+ String.format("%3d done. hostname: %s!", me,
										hname);
						writer_Client.println(message);
						writer_Client.flush();
						//-----------------------------------------------------
						if (checker.isAlive()) {
							checker.end();
						}
						//-----------------------------------------------------
						break;
					}
				}
				writer_Client.close();
			}
		} catch (FileNotFoundException e) {
			System.out.println("Caught exception " + e);
		} catch (Exception e) {
			System.out.println("Caught exception " + e);
		}
		MPI.Finalize();
		System.exit(-1);
	}
}