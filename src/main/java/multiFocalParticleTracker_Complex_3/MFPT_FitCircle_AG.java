package multiFocalParticleTracker_Complex_3;

import org.apache.commons.math3.linear.*;

public class MFPT_FitCircle_AG{
	/*
	for a better understanding, the following link is recommended:
	https://meshlogic.github.io/posts/jupyter/curve-fitting/fitting-a-circle-to-cluster-of-3d-points/
	*/
//	public void run(String arg){// for testing
//		double r=20;
//		double [] x=new double [360];
//		double [] y=new double [360];
//		double [] intensity=new double [360];
//		for(int i=0;i<360; i++){
//			x[i]=3+r*Math.cos(Math.toRadians(i))+2*(1-2*Math.random());
//			y[i]=15+r*Math.sin(Math.toRadians(i))+2*(1-2*Math.random());
//			intensity[i]=100*Math.random();
//		}
//	        
//		double [] result1=getRadiusAndR2(x,y,intensity, 3.0,15.0);
//		IJ.log(String.format("radius is %.2f, R2 is %.2f", result1[0], result1[1]));
//		
//		double [] result2=getCenterRadiusAndR2(x,y,intensity);
//		IJ.log(String.format("center is [%.2f, %.2f], radius is %.2f, R2 is %.2f", result2[0], result2[1],result2[2],result2[3]));
//		
//	}

	/**
	 * Fit a circle with defined centre to an intensity matrix and output the radius of the circle and the goodness of the fit (R^2)
	 * Code developed by An Gong and Luis Alvarez, 01.07.2019
	 * */
	public static double [] getRadiusAndR2(double [] x, double [] y, double [] intensity, double centerX, double centerY){
		//[weight]*[(x-centerX)^2+(y-centerY)^2]=[weight]C, where C=radius^2 and we take intensity^2 as the weight. 
		int n=x.length;
		double [][] arrayA=new double [n][1];
		double [] arrayB=new double [n];
		double meanB=0;

		double weight;
		for(int i=0;i<n;i++){
			weight=intensity[i]*intensity[i];
			arrayA[i][0]=weight;
			arrayB[i]=weight*((x[i]-centerX)*(x[i]-centerX)+(y[i]-centerY)*(y[i]-centerY));
			meanB+=arrayB[i];
		}
		meanB/=n;
		RealMatrix matA=new BlockRealMatrix(arrayA);
		RealVector vecB = new ArrayRealVector(arrayB);
 		DecompositionSolver solver = new SingularValueDecomposition(matA).getSolver();
 		RealVector vecC = solver.solve(vecB);
 		double [] arrayC=vecC.toArray();
		double radius=Math.sqrt(arrayC[0]);

		//Here is calcualting R2, using the definition from wiki entry of "Coefficient of determination" 
		double SSTot=0;
		double SSRes=0;
		for(int i=0; i<n; i++){
			weight=intensity[i]*intensity[i];
			SSTot+=(arrayB[i]-meanB)*(arrayB[i]-meanB);
			SSRes+=(weight*arrayC[0]-arrayB[i])*(weight*arrayC[0]-arrayB[i]);
		}
		//IJ.log("SSRes is "+SSRes+", SSTot is "+SSTot);
		double R2=1;
		if(SSTot/meanB>1e-8){  //when SSTot is zero, there might be some calculation error from the digital precision issue 
			R2=1-SSRes/SSTot;
		}
		double [] result=new double [2];
		result[0]=radius;
		result[1]=R2;
		return result;
		
	}

	/**
	 * Fit a circle to an intensity matrix and output the x centre, the y centre, the radius of the circle, and the goodness of the fit (R^2)
	 * Code developed by An Gong and Luis Alvarez, 01.07.2019
	 * */
	public static double [] getCenterRadiusAndR2(double [] x, double [] y, double [] intensity){
		/* The circle equation of  (x-centerX)^2+(y-centerY)^2=radius^2 can be rewrite as:
		C0*(2*x) + C1*(2*y) + C2 = x^2 + y^2
		where C0 = centerX;
		C1 = centerY;
		C2 = radius^2 - centerX^2 - centerY^2;
		We look then for parameters C0, C1, and C2 such that our experimental
		the linear fit is contsructing into the maxtix form of [weight][2x, 2y, 1]*[C0, C1, C2]'=[weight][x^2+y^2] with intensity^2 as the weight
		*/
		int n=x.length;
		double [][] arrayA=new double [n][3];
		//double [][] arrayWeight=new double[n][n];
		double [] arrayB=new double [n];
		double meanB=0;
		double weight;
		for(int i=0;i<n;i++){
			weight=intensity[i]*intensity[i];
			arrayA[i][0]=2*x[i]*weight;
			arrayA[i][1]=2*y[i]*weight;
			arrayA[i][2]=1*weight;
			arrayB[i]=(x[i]*x[i]+y[i]*y[i])*weight;
			meanB+=arrayB[i];		
		}
		meanB/=n;
 		RealVector vecB = new ArrayRealVector(arrayB);
 		RealMatrix matA=new BlockRealMatrix(arrayA);
 		
 		DecompositionSolver solver = new SingularValueDecomposition(matA).getSolver();
 		RealVector vecC = solver.solve(vecB);
 		double [] arrayC=vecC.toArray();
 		
 		//Here is calcualting R2, using the definition from wiki entry of "Coefficient of determination" 
		double SSTot=0;
		double SSRes=0;
		double res;
		for(int i=0; i<n; i++){
			weight=intensity[i]*intensity[i];
			SSTot+=(arrayB[i]-meanB)*(arrayB[i]-meanB);
			res=2*x[i]*weight*arrayC[0]+2*y[i]*weight*arrayC[1]+weight*arrayC[2]-arrayB[i];
			SSRes+=res*res;
		}
		//IJ.log("SSRes is "+SSRes+", SSTot is "+SSTot);
		double R2=1;
		if(SSTot/meanB>1e-8){  //when SSTot is too small, there might be some calculation error from the digital precision issue 
			R2=1-SSRes/SSTot;
		}
		
		double [] result=new double [4];
		result[0]=arrayC[0];
		result[1]=arrayC[1];
		result[2]=Math.sqrt(arrayC[0]*arrayC[0]+arrayC[1]*arrayC[1]+arrayC[2]);
		result[3]=R2;
		return result;
		
	}

}
