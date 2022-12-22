package EnergySequential;

import java.util.ArrayList;
import mpi.*;
//import java.io.*;
//import java.util.*;   //Random is defined in the "java.util" library package, Also Vector is defined here import java.util.Vector;
//import java.lang.Math;
//import matlabcontrol.*;
//import ilog.concert.*;
//import ilog.cplex.*;

public class fBlackBoxLinear {
	
	private double ObjectiveValue=10000000;

	public fBlackBoxLinear(double[] thetaVec, double[][] PMat, double[][] DMat,double[][] EMat, int M, int T, double[] ExpectedP, double[] ExpectedD, double[] ExpectedE, double etaD, double etaC, double Rcap, double DeltaRC, double DeltaRD, double RInitial, int K, int Tag_GR_RG){
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

	
	   
  
	    double[][] R_TransitionFunction=new double[M][T+1];
	    for (int ellm=0; ellm<M; ellm++){
	    	R_TransitionFunction[ellm][0]=RInitial;
	    	for (int ellt=1; ellt<T+1; ellt++){R_TransitionFunction[ellm][ellt]=0;}
	                                    }
	    double sum=0;

	   
	    for (int elli=0; elli<M; elli++){
	        for (int ellt=0; ellt<T-1; ellt++){
	        	
	        	x_sol_opt[0]=thetaVec[0]*DMat[elli][ellt]+thetaVec[1]*EMat[elli][ellt]+
	        			     thetaVec[2]*R_TransitionFunction[elli][ellt]+thetaVec[3];
	        	x_sol_opt[1]=thetaVec[4]*DMat[elli][ellt]+thetaVec[5]*EMat[elli][ellt]+
       			             thetaVec[6]*R_TransitionFunction[elli][ellt]+thetaVec[7];
	        	if (Tag_GR_RG==0){
	        	x_sol_opt[2]=thetaVec[8]*DMat[elli][ellt]+thetaVec[9]*EMat[elli][ellt]+
	       			         thetaVec[10]*R_TransitionFunction[elli][ellt]+thetaVec[11];
	        	x_sol_opt[3]=0;
	        	}else if (Tag_GR_RG==1){
		        	x_sol_opt[2]=0;
		        	x_sol_opt[3]=thetaVec[8]*DMat[elli][ellt]+thetaVec[9]*EMat[elli][ellt]+
		       			         thetaVec[10]*R_TransitionFunction[elli][ellt]+thetaVec[11];
		        	}
	        	double alpha_t=(etaC*(EMat[elli][ellt]-Math.min(EMat[elli][ellt],DMat[elli][ellt]))-(1/etaD)*(DMat[elli][ellt]-Math.min(EMat[elli][ellt],DMat[elli][ellt])))/Rcap;
	            R_TransitionFunction[elli][ellt+1]=R_TransitionFunction[elli][ellt]+alpha_t+(1/Rcap)* U_vec[0]*x_sol_opt[0]+(1/Rcap)* U_vec[1]*x_sol_opt[1]+(1/Rcap)* U_vec[2]*x_sol_opt[2]+(1/Rcap)* U_vec[3]*x_sol_opt[3];
	 	        sum=sum+PMat[elli][ellt]*(x_sol_opt[0]+x_sol_opt[1]-etaD*x_sol_opt[2]-x_sol_opt[3]);
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
