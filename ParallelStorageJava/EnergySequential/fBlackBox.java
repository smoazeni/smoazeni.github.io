package EnergySequential;

import java.util.ArrayList;
import mpi.*;
//import java.io.*;
//import java.util.*;   //Random is defined in the "java.util" library package, Also Vector is defined here import java.util.Vector;
//import java.lang.Math;
//import matlabcontrol.*;
//import ilog.concert.*;
//import ilog.cplex.*;

public class fBlackBox {
	
	private double ObjectiveValue=10000000;

	public fBlackBox(double[] thetaVec, double[][] PMat, double[][] DMat,double[][] EMat, int M, int T, double[] ExpectedP, double[] ExpectedD, double[] ExpectedE, double etaD, double etaC, double Rcap, double DeltaRC, double DeltaRD, double RInitial, int K){
		double[] x_sol_opt= new double[4];
		double[] xcplex=new double[4];
		double[] lB={0,0,0,0};
		double R_Final=0;
		double[] linearTerm_Objective=new double[4];
		double risk_aversion_parameter=0;
		double risk;
		boolean debug = true;
		double[] V_vec = {1, 1,(-etaD),(-1)};
	    double[] U_vec={etaC, 1/etaD, -1, -etaC};
        double[] bvec_t=new double[6];
	    double[][] theta=new double[T-1][K];
	    for (int ellt=0; ellt<T-1; ellt++){
		    for (int ellk=0; ellk<K; ellk++){
	             theta[ellt][ellk]=thetaVec[ellt*K+ellk];
	             System.out.println("Theta Value:    "+theta[ellt][ellk]);
	             
		     }
	    }

	
	    
	    double[][] A_mat = {{-etaC,-(1/etaD),1,etaC}, {0,-(1/etaD),1,0}, {etaC,0,0,-etaC}, {etaC,1/etaD,-1,-etaC},{0, 1,0,0},{0, 0,0,1}};//6*4
	    ArrayList<ArrayList<Double>> A = new ArrayList<ArrayList<Double>>();
	    for (int i = 0; i < A_mat.length; i++) {
			ArrayList<Double> a_row = new ArrayList<Double>();
			for (int j = 0; j < A_mat[i].length; j++) {
				if (A_mat[i][j] != 0) {
					a_row.add((double)j);
					a_row.add(A_mat[i][j]);
				}
			}
			A.add(a_row);
		}
  
	    double[][] R_TransitionFunction=new double[M][T+1];
	    for (int ellm=0; ellm<M; ellm++){
	    	R_TransitionFunction[ellm][0]=RInitial;
	    	for (int ellt=1; ellt<T+1; ellt++){R_TransitionFunction[ellm][ellt]=0;}
	                                    }
	    double sum=0;

	   
	    for (int elli=0; elli<M; elli++){
	        for (int ellt=0; ellt<T-1; ellt++){
	        	double alpha_t=(etaC*(EMat[elli][ellt]-Math.min(EMat[elli][ellt],DMat[elli][ellt]))-(1/etaD)*(DMat[elli][ellt]-Math.min(EMat[elli][ellt],DMat[elli][ellt])))/Rcap;
	        	bvec_t[0]=Rcap*R_TransitionFunction[elli][ellt]+alpha_t*Rcap;
	        	bvec_t[1]=DeltaRD*Rcap-(1/etaD)*(DMat[elli][ellt]-Math.min(EMat[elli][ellt],DMat[elli][ellt]));
	        	bvec_t[2]=DeltaRC*Rcap-etaC*(EMat[elli][ellt]-Math.min(EMat[elli][ellt],DMat[elli][ellt]));
	        	bvec_t[3]=(1-R_TransitionFunction[elli][ellt])*Rcap-alpha_t*Rcap;
	        	bvec_t[4]=DMat[elli][ellt]-Math.min(EMat[elli][ellt],DMat[elli][ellt]);	
	        	bvec_t[5]=EMat[elli][ellt]-Math.min(EMat[elli][ellt],DMat[elli][ellt]);		

                ArrayList<Double> b = new ArrayList<Double>();
                for (int i = 0; i < A_mat.length; i++) {b.add(bvec_t[i]);}
	            		
	            //HMatrix=2*(theta(t,2)/(Rcap^2))*(v*v');    %Hessian Matrix 4*4
                /*double coefficient=theta[ellt][1]+theta[ellt][3]*ExpectedE[ellt+1]+
	                               theta[ellt][4]*ExpectedD[ellt+1]+theta[ellt][5]*ExpectedP[ellt+1]+
	                               theta[ellt][6]*ExpectedP[ellt+1]*ExpectedD[ellt+1]+
	                               theta[ellt][7]*ExpectedP[ellt+1]*ExpectedE[ellt+1]+
	                               theta[ellt][8]*ExpectedD[ellt+1]*ExpectedE[ellt+1]+
	                               theta[ellt][9]*ExpectedD[ellt+1]*ExpectedE[ellt+1]*ExpectedP[ellt+1];*/
 	    
                
                
                //System.out.println("---------------------------------------------");
                //System.out.println("Number of Period: "+ellt);
                //System.out.println("Number of Simulation: "+elli);
                //System.out.println("Theta Value:    "+theta[ellt][0]);
                //System.out.println("Expected Price at next period:    "+ExpectedP[ellt+1]);
                //System.out.println("Spot Price of this path at this period:    "+PMat[elli][ellt]);
                //System.out.println("Difference of Expected Price and Price:    "+(PMat[elli][ellt]-ExpectedP[ellt+1]));
                //System.out.println("Difference of Expected Price and Theta Price:    "+(PMat[elli][ellt]-theta[ellt][0]*ExpectedP[ellt+1]));
                
                
        	    for (int ellj=0;ellj<4;ellj++){
        	    linearTerm_Objective[ellj]=PMat[elli][ellt]*V_vec[ellj]-theta[ellt][0]*etaD*ExpectedP[ellt+1]*U_vec[ellj];
        	    //System.out.println("Linear Term Objective:  "+linearTerm_Objective[ellj]);
        	    }
        	    
        	    
        	    //long cplex_time = System.nanoTime(); 
        	    //MIPSparseSolver MS_LP;
        	    xcplex=MIPSparseSolver.createSolve(linearTerm_Objective, A, b, null, null, lB, null, null, null, true, 1e-8, debug);
        	   // double secs_cplex = 1.0 * (System.nanoTime() - cplex_time) / 1000000000.0;
    			//System.out.println("Done in: "+ String.format("%8.10f secs!",secs_cplex ));
  
        	    for (int i = 0; i < 4; i++) {x_sol_opt[i]=xcplex[i];}
	            R_TransitionFunction[elli][ellt+1]=R_TransitionFunction[elli][ellt]+alpha_t+(1/Rcap)* U_vec[0]*x_sol_opt[0]+(1/Rcap)* U_vec[1]*x_sol_opt[1]+(1/Rcap)* U_vec[2]*x_sol_opt[2]+(1/Rcap)* U_vec[3]*x_sol_opt[3];
	 	        sum=sum+PMat[elli][ellt]*(x_sol_opt[0]+x_sol_opt[1]-etaD*x_sol_opt[2]-x_sol_opt[3]);
                //for (int ellj=0;ellj<4;ellj++){System.out.println(x_sol_opt[ellj]+"");}
                //System.out.println("R_TransitFunction:     "+ R_TransitionFunction[elli][ellt+1]);
            }
	        x_sol_opt[0]=0;
	        x_sol_opt[1]=Math.max(0,DMat[elli][T-1]-Math.min(DMat[elli][T-1],EMat[elli][T-1])-etaD*Rcap*Math.min(R_TransitionFunction[elli][T-1],DeltaRD));
	        x_sol_opt[2]=(1/etaD)*Math.max(0,etaD*Math.min(R_TransitionFunction[elli][T-1],DeltaRD)*Rcap-(DMat[elli][T-1]-Math.min(DMat[elli][T-1],EMat[elli][T-1])));
	        x_sol_opt[3]=EMat[elli][T-1]-Math.min(EMat[elli][T-1],DMat[elli][T-1]);
	        R_Final=R_TransitionFunction[elli][T-1]+((EMat[elli][T-1]-Math.min(EMat[elli][T-1],DMat[elli][T-1]))-(DMat[elli][T-1]-Math.min(EMat[elli][T-1],DMat[elli][T-1])))/Rcap+(1/Rcap)*(etaC*x_sol_opt[0]+(1/etaD)*x_sol_opt[1]-x_sol_opt[2]-etaC*x_sol_opt[3]);
	        //System.out.println("Final Amount in Battery:   "+R_Final);
	        sum=sum+PMat[elli][T-1]*(x_sol_opt[0]+x_sol_opt[1]-etaD*x_sol_opt[2]-x_sol_opt[3]);
	    }

	    risk=0;
	    ObjectiveValue=(sum/M)+risk_aversion_parameter*risk;

}
	public double getObjective_velue(){return ObjectiveValue;}
}
