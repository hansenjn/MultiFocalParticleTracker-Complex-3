/***===============================================================================
 
 https://github.com/hansenjn/MultiFocalParticleTracker-Complex-3, Version v0.0.8

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation (http://www.gnu.org/licenses/gpl.txt )

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 See the GNU General Public License for more details.
 
 @authors Jan N Hansen, Luis Alvarez, An Gong
 @copyright (C) 2018-2020: Jan N Hansen, Luis Alvarez, An Gong
   
 For any questions please feel free to contact me (jan.hansen@uni-bonn.de).

==============================================================================**/

package multiFocalParticleTracker_Complex_3;

public class MFPTCPoint3{
	private double x,y,correctedZFact,correctedZDiff,z,zByRadiusMethod2;
	private int t;
	String selectedZCombo = "";
	double [] selZComboAsArray = new double [0];
	double [] correctionFactors = new double [0];
	double [] correctionDifference = new double [0];
	String involvedPlanes = "";
	int nrOfInvolvedPlanes;
	
	private double [] radiusCentred; //plane
	private double [] fitX; //plane
	private double [] fitY; //plane
	private double [] r2; //plane
	private double [] r2Centred; //plane
	
	public MFPTCPoint3(double px, double py, int frame, int nSlices){
		x = px;
		y = py;		
		t = frame;
		
		correctedZFact = Double.NaN;
		correctedZDiff = Double.NaN;
		z = Double.NaN;
		zByRadiusMethod2 = Double.NaN;
		nrOfInvolvedPlanes = 0;
		
		fitX = new double [nSlices];
		fitY = new double [nSlices];
		r2 = new double [nSlices];
		r2Centred = new double [nSlices];
		radiusCentred = new double [nSlices];
		
		for(int s = 0; s < nSlices; s++){
			fitX [s] = Double.NEGATIVE_INFINITY;
			fitY [s] = Double.NEGATIVE_INFINITY;
			r2 [s] = Double.NEGATIVE_INFINITY;
			r2Centred [s] = Double.NEGATIVE_INFINITY;
			radiusCentred [s] = Double.NEGATIVE_INFINITY;
		}	
	}	
	
	public double X(){
		return x;
		
	}
	
	public double Y(){
		return y;
		
	}
	
	public int T(){
		return t;
		
	}
	
	public double Z(){
		return z;
		
	}
		
	public double zAvgCorrectedFact(){
		return correctedZFact;
	}
	
	public double zAvgCorrectedDiff(){
		return correctedZDiff;
	}
	
	public double zByRadiusCentMethod2(){
		return zByRadiusMethod2;
	}
	
	public int getNrOfInvolvedPlanes(){
		return nrOfInvolvedPlanes;
		
	}
			
	void setZ(double newZ){
		z = newZ;
	}
	
	void setCorrectedZAvgFact (double newZ){
		correctedZFact = newZ;
	}
	
	void setCorrectedZAvgDiff (double newZ){
		correctedZDiff = newZ;
	}
	
	void setZByRadiusMethod2(double newZ){
		zByRadiusMethod2 = newZ;
	}
		
	void setSelectedZCombo(String combo){
		selectedZCombo = "" + combo;
	}
	
	void setSelectedZComboAsArray(double [] selZComboArray){
		selZComboAsArray = new double [selZComboArray.length];
		for(int i = 0; i < selZComboArray.length; i++){
			selZComboAsArray [i] = selZComboArray [i]; 
		}
	}
	
	void setInvolvedPlanes(String combo){
		involvedPlanes = "" + combo;
	}
	
	void setNrOfInvolvedPlanes(int nr){
		nrOfInvolvedPlanes = nr;
		
	}
	
	public boolean zCorrDefFact(){
		if(Double.isNaN(correctedZFact)){
			return false;
		}
		return true;
	}
	
	public boolean zCorrDefDiff(){
		if(Double.isNaN(correctedZDiff)){
			return false;
		}
		return true;
	}
	
	public boolean zDef(){
		if(Double.isNaN(z)){
			return false;
		}
		return true;
	}
	
	public double radiusCentred (int slice){
		if(radiusCentred [slice]!=Double.NEGATIVE_INFINITY){
			return radiusCentred [slice];
		}else{
			return 0.0;
		}	
	}
	
	void setRadiusCentred (double radius, int slice){
		radiusCentred [slice] = radius;
	}
	
	public double fitX(int slice){
		if(fitX[slice]!=Double.NEGATIVE_INFINITY){
			return fitX[slice];
		}else{
			return 0.0;
		}	
	}
	
	void setFitX(double fitValue, int slice){
		fitX[slice] = fitValue;
	}
	
	public double fitY(int slice){
		if(fitY[slice]!=Double.NEGATIVE_INFINITY){
			return fitY[slice];
		}else{
			return 0.0;
		}	
	}
	
	void setFitY(double fitValue, int slice){
		fitY[slice] = fitValue;
	}
	
	public double r2(int slice){
		if(r2[slice]!=Double.NEGATIVE_INFINITY){
			return r2[slice];
		}else{
			return 0.0;
		}	
	}
	
	void setR2(double goodness, int slice){
		r2[slice] = goodness;
	}
	
	public double r2Centred(int slice){
		if(r2Centred[slice]!=Double.NEGATIVE_INFINITY){
			return r2Centred[slice];
		}else{
			return 0.0;
		}	
	}
	
	void setR2Centred(double goodness, int slice){
		r2Centred[slice] = goodness;
	}
	
	void setCorrectionFactors(double [] corrFactors){
		correctionFactors = corrFactors;
	}
	
	void setCorrectionDifference(double [] corrDiff){
		correctionDifference = corrDiff;
	}
}
