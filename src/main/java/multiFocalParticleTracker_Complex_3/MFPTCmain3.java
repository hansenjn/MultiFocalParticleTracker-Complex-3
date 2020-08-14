/***===============================================================================
 
 https://github.com/hansenjn/MultiFocalParticleTracker-Complex-3, Version v0.1.0

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

import java.awt.Color;
import java.awt.event.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;

import javax.swing.UIManager;

import ij.*;
import ij.gui.*;
import ij.io.*;
import ij.measure.*;
import ij.measure.SplineFitter;
import ij.plugin.*;
import ij.text.TextPanel;
import multiFocalParticleTracker_Complex_3.MFPT_FitCircle_AG;
import multiFocalParticleTracker_Complex_3.jnhsupport.*;

public class MFPTCmain3 implements PlugIn, Measurements{
	//Name
		public static final String PLUGINNAME = "MultiFocalParticleTracker-Complex-3";
		public static final String PLUGINVERSION = "v0.1.0";
		
		double xyCal = 0.34375;	//0.34375 for 32x, 0.55 for 20x		
		int maxRadius = 10;
		boolean savePlot = false;
		double widthThresholdMin = 1.5, widthThresholdMax = 4.5;
		int upscaling = 10;
		
		ProgressDialog progress;
		boolean done = false;
		
		double [][] correctionFactors, correctionDifferences;
		int refPlane;
		
		double binningFactorForCorrFactors = 10.0;
		
	@Override
	public void run(String arg) {
		/**&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		GenericDialog
		&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
		
		GenericDialog gd = new GenericDialog("mulit-focal analysis");		
//		setInsets(top, left, bottom)
		gd.setInsets(0,0,0);	gd.addMessage(PLUGINNAME + ", version " + PLUGINVERSION + " (\u00a9 2017-" + constants.dateY.format(new Date()) + ")", constants.Head1);
					
		gd.setInsets(10,0,0);	gd.addNumericField("xy calibration [um]", xyCal, 5);
		gd.setInsets(10,0,0);	gd.addNumericField("Radius threshold [um]: min | max", widthThresholdMin, 5);
		gd.setInsets(-23,55,0);	gd.addNumericField("", widthThresholdMax, 5);		
		gd.setInsets(10,0,0);	gd.addNumericField("upscaling of LUT", upscaling, 0);
		gd.setInsets(10,0,0);	gd.addNumericField("maximum radius for xy fitting (px)", maxRadius, 0);
		gd.setInsets(10, 0, 0);	gd.addCheckbox("Save Likelihood Plot", savePlot);
		gd.setInsets(10, 0, 0);	gd.addCheckbox("Verify calibration mode", false);
						
		gd.showDialog();
			 	
	 	xyCal = gd.getNextNumber();	 	
		widthThresholdMin = (double) gd.getNextNumber();
		widthThresholdMax = (double) gd.getNextNumber();
		upscaling = (int)gd.getNextNumber();
		maxRadius = (int)gd.getNextNumber();		
		savePlot = gd.getNextBoolean();
		boolean verify = gd.getNextBoolean();
		
		if (gd.wasCanceled())return;
	  	
  		/**&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		Initiate multi task management
		&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
		try{
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		}catch(Exception e){}
		
		//select tasks
		OpenFilesDialog od = new OpenFilesDialog ();
		od.setLocation(0,0);
		od.setVisible(true);		
		od.addWindowListener(new java.awt.event.WindowAdapter() {
	        public void windowClosing(WindowEvent winEvt) {
	        	//Analysis canceled
	        	return;
	        }
	    });	
		while(od.done==false){
			 try{
				 Thread.currentThread().sleep(50);
		     }catch(Exception e){
		     }
		}
		
		int tasks = od.filesToOpen.size();
		String [] name = new String [tasks];
		String [] dir = new String [tasks];
		boolean tasksSuccessfull [] = new boolean [tasks];
		for(int task = 0; task < tasks; task++){
			name [task] = od.filesToOpen.get(task).getName();
			dir [task] = od.filesToOpen.get(task).getParent() + System.getProperty("file.separator");
			tasksSuccessfull [task] = false;
		}	
				
		//start progress dialog
		progress = new ProgressDialog(name, tasks);
		progress.setVisible(true);
		progress.addWindowListener(new java.awt.event.WindowAdapter() {
	        public void windowClosing(WindowEvent winEvt) {
	        	progress.stopProcessing();
	        	if(done==false){
	        		IJ.error("Script stopped...");
	        	}       	
	        	System.gc();
	        	return;
	        }
		});	
		
		/**&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		Processing
		&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
//		String homePath = FileSystemView.getFileSystemView().getHomeDirectory().getAbsolutePath();
		
		// initialize
		Date startDate;
		TextPanel tp, tp2;	
		String nameImage = "";
    	String dirImage = "";
    	String nameLUT = "";
    	String dirLUT = "";   
    	String nameLUTSD = "";
    	String dirLUTSD = "";    	
		double rawLUT [][], LUT [][], LUTSD [][], LUTPrecision [][];
		
				
		{
			// get image location
	    	OpenDialog oImage;	    	
    		oImage = new OpenDialog("Open corresponding image", null);
    		progress.replaceBarText("Open corresponding image");
    		nameImage = oImage.getFileName();
	    	dirImage = oImage.getDirectory();
	    	
	    	// get look up table location
	    	OpenDialog oLUT;	    	
    		oLUT = new OpenDialog("Open corresponding LUT file (.txt)", null);
    		progress.replaceBarText("Open corresponding LUT file (.txt)");
    		nameLUT = oLUT.getFileName();
	    	dirLUT = oLUT.getDirectory();
	    	
	    	// get sd of look up table location
	    	OpenDialog oLUTSD;	    	
    		oLUTSD = new OpenDialog("Open corresponding LUT standard deviation file (.txt)", null);
    		progress.replaceBarText("Open corresponding LUT standard deviation file (.txt)");
    		nameLUTSD = oLUTSD.getFileName();
	    	dirLUTSD = oLUTSD.getDirectory();
		}
		System.gc();
		
		//Initialize and open basic files
		ImagePlus imp;		
    	imp = IJ.openImage(dirImage+nameImage);	
    	imp.getCalibration().pixelHeight = xyCal;
    	imp.getCalibration().pixelWidth = xyCal;	
    	rawLUT = this.readLUT(dirLUT + nameLUT);
		LUT = this.readUpscaledLUT(dirLUT + nameLUT, upscaling);
		LUTSD = this.readUpscaledLUT(dirLUTSD + nameLUTSD, upscaling);
		LUTPrecision = getLUTPrecision(LUT, LUTSD);
		
		int lutCenterPos [] = new int [LUT.length-1];
		double LutMin, LutTemp; int LutCt;

		//Initialize for initial params
			int binnedDistanceMatrix [][] = this.getBinnedDistanceMatrix();
			double uncalibDistanceMatrix [][] = this.getUncalibDistanceMatrix();
					
		//Determine center line
		for(int s = 1; s < LUT.length; s++){
			LutMin = Double.POSITIVE_INFINITY;
			for(int i = 2; i < LUT[0].length-2; i++){
				LutTemp = 0.0;
				LutCt = 0;
				//smooth (window size 5)
				for(int ii = i-2; ii < i+3; ii++){
					if(LUT[s][ii] != Double.NEGATIVE_INFINITY){
						LutTemp += LUT[s][ii];
						LutCt++;
					}
				}	
				if(LutCt == 5){
					LutTemp /= LutCt;
					if(LutTemp < LutMin){
						LutMin = LutTemp;
						lutCenterPos [s-1] = i;
					}
				}
			}
			if(lutCenterPos [s-1] == -1){
				IJ.error("LUT is not in a correct format - center line is missing!");
			}else{
				progress.notifyMessage("LUT center pos for " + s + ": " + lutCenterPos [s-1], ProgressDialog.LOG);
			}
		}
						
		//processing
		tasking: for(int task = 0; task < tasks; task++){
			progress.updateBarText("in progress...");
			startDate = new Date();
			progress.clearLog();
			
			running: while(true){
				//Check if image is processable
				if(imp.getBitDepth()==24){
		    		progress.notifyMessage("ERROR: Multi focal analysis cannot be used for RGB images!", ProgressDialog.ERROR);
			  		break running;
		    	}				
				
				//determine save path
		    	String savePath;
		    	if(name[task].contains(".")){
					savePath = name[task].substring(0,name[task].lastIndexOf("."));
				}else{
					savePath = name[task];
				}	
		    	savePath = dir [task] + savePath + "_CMFPT" + "_" + constants.df0.format(upscaling)+ "_" + constants.df0.format(maxRadius);
		    	
		    	//Save raw LUT and interpolated LUT
		    	plotXY2DArray(LUT, "Raw LUT", "Position", "Fit Width", savePath + "_rawLUT", false);
		    	output2DArray(LUT, "SplineFitted LUT", name[task], savePath + "_iLUT");
		    	plotXY2DArray(LUT, "Spline Interpolated LUT", "Position", "Fit Width", savePath + "_iLUT", false);
		    	plotXY2DArray(LUTSD, "Spline Interpolated SD LUT", "Position", "Standard Deviation of Fit Width", savePath + "_iLUTSD", false);
		    	plotXY2DArray(LUTPrecision, "Precision of Spline Interpolated LUT", "Position", "z-precision (µm)", savePath + "_iLUTP", false);

		    	/**
		    	 * ANALYSIS
		    	 * */
				{			
					//read x,y, points
					ArrayList<MFPTCPoint3> track;
					if(verify){
						track = readTrackListInVerificationMode(dir[task]+name[task],imp.getNSlices(),imp.getNFrames());
					}else{
						track = readTrackList(dir[task]+name[task],imp.getNSlices());
					}					
					if(track.equals(null))	break running;
					  	
					//Determine z
					{
						//Constants for processing
						String out;
						double [][] radiiNew; //[nr][0 = value, 1 = slicePos]
						double [][] radiusZ;
						double [][] radiusZPrecision;
						boolean [] keepPlane;
						boolean keep1, keep2;
						int selCombo, minDifferencePos;
						double finalZ, minSD, SD, minDifference, differenceTemp;
						String selectedZCombo;
						String selectedZComboPrecision;
						double selZCombo [] = new double [imp.getNSlices()];
						double selZComboPrecision [] = new double [imp.getNSlices()];
						double positions [] = new double [maxRadius*2+1];
						for(int i = 0; i < maxRadius*2+1; i++){
							positions [i] = i - maxRadius;
						}
						
						double pxi;
						double pyi;
						double sigma = 0.0;
						double intensities [][] = new double [1][1];
						double [] result;
						int frame;
						
						//z method one
						double [][] zValues, zValuesPrecision;
						int bestPrecZ;
						double minPrec;
						int nrOfCombinations;
						int comboCounter;
						int wZPsCounter;
						int arraySize, wzArraySize;
						int nrOfValues;
						
						//method two;
						double stepSizeInLut;
						double [][] comparPlot = new double [2][LUT[0].length];
						double [][] diffPlot = new double [0][0]; 
						double LUTValue;
						double minValue;
						int index; 
						int nrOfIncludedPlanes;
						String includedPlanes;
											
						//Defining most central plane as reference;
						refPlane = (int) (((double)LUT.length-2.0)/2.0);
						int refIndex;
						double minDist;
						double [] transientCorrFactors, transientCorrDifferences;
						progress.notifyMessage("Reference plane for correction factor estimation: " + refPlane, ProgressDialog.LOG);
												
						for(int i = 0; i < track.size(); i++){
							try{
								/**
								 * Radial fitting and testing radial fit
								 * */
								{	
									pxi = track.get(i).X();
									pyi = track.get(i).Y();
									frame = track.get(i).T();
//									progress.updateBarText("analyzing particle " + constants.df6US.format(track.get(i).X()) 
//											+ " Y " + constants.df6US.format(track.get(i).Y())
//											+ " ID " + constants.df0.format(i)
//											+ " T " + constants.df0.format(frame));
//									skippedPlanes = "";
									for(int s = 0; s < imp.getNSlices(); s++){
										intensities = new double [1][1];
										try{
											//GET INTENSITIES	
											intensities = getIntensitiesInRadius(imp, pxi, pyi, frame, s, binnedDistanceMatrix, uncalibDistanceMatrix, maxRadius);
										}catch(Exception e){
											out = "";
											for(int err = 0; err < e.getStackTrace().length; err++){
												out += " \n " + e.getStackTrace()[err].toString();
											}
											progress.notifyMessage("Exception in forming intensity profile " + constants.df6US.format(track.get(i).X()) 
												+ " Y " + constants.df6US.format(track.get(i).Y())
												+ " ID " + constants.df0.format(i)
												+ " T " + constants.df0.format(frame)
												+ " S " + constants.df0.format(s) 
												+ "Cause: " + e.getCause()+ "\n"	
												+ "Throwable: " + e.getMessage() + "\n"	
												+ out, 
												ProgressDialog.LOG);
//											progress.notifyMessage("sigma " + sigma, ProgressDialog.LOG);
										}
										if(intensities.length == 1){
//											progress.notifyMessage("Intensities could not be grabbed for " + constants.df6US.format(track.get(i).X()) 
//											+ " Y " + constants.df6US.format(track.get(i).Y())
//											+ " ID " + constants.df0.format(i)
//											+ " T " + constants.df0.format(frame)
//											+ " S " + constants.df0.format(s) 
//											+ "Sigma: " + sigma, ProgressDialog.LOG);
											continue;
										}
										try{
											result = MFPT_FitCircle_AG.getCenterRadiusAndR2(intensities [0], intensities [1], intensities [2]);
											track.get(i).setFitX(result[0], s);
											track.get(i).setFitY(result[1], s);
											track.get(i).setRadiusCentred(result[2], s);
											track.get(i).setR2Centred(result[3], s);							
										}catch(Exception e){
											out = "";
											for(int err = 0; err < e.getStackTrace().length; err++){
												out += " \n " + e.getStackTrace()[err].toString();
											}
											progress.notifyMessage("Exception in fitting (with free centre) for particle X " 
												+ constants.df6US.format(track.get(i).X()) 
												+ " Y " + constants.df6US.format(track.get(i).Y())
												+ " ID " + constants.df0.format(i)
												+ " T " + constants.df0.format(frame)
												+ " S " + constants.df0.format(s) 
												+ "Cause: " + e.getCause()+ "\n"	
												+ "Throwable: " + e.getMessage() + "\n"	
												+ out, 
												ProgressDialog.LOG);
											progress.notifyMessage("sigma " + sigma, ProgressDialog.LOG);
											for(int abc = 0; abc < intensities[0].length; abc++){
												progress.notifyMessage(abc + ": " + intensities[0][abc] + " / "
														+ intensities[1][abc] + " / " + intensities[2][abc], ProgressDialog.LOG);	
											}
										}
									}															
									progress.addToBar(0.8/(double)track.size());
//									if(skippedPlanes.length()>0){
//										progress.notifyMessage("particle " + constants.df0.format(i+1)
//											+ " - planes " + skippedPlanes + "skipped - shift undefined", 
//											multi_focal_pt_complex.jnhsupport.ProgressDialog.LOG);
//									}								
								}
								
								/**
								 * Determining Z
								 * */
								nrOfValues = 0;
								for(int s = 0; s < imp.getNSlices(); s++){
									if(track.get(i).r2Centred(s) > 0.8
											&& track.get(i).radiusCentred(s) >= widthThresholdMin
											&& track.get(i).radiusCentred(s) <= widthThresholdMax){
										nrOfValues++;
									}
								}
								
								if(nrOfValues>1){
									radiiNew = new double [nrOfValues][2];
									nrOfValues = 0;
									for(int s = 0; s < imp.getNSlices(); s++){
										if(track.get(i).r2Centred(s) > 0.80
												&& track.get(i).radiusCentred(s) >= widthThresholdMin
												&& track.get(i).radiusCentred(s) <= widthThresholdMax){
											radiiNew [nrOfValues][0] = track.get(i).radiusCentred(s);
											radiiNew [nrOfValues][1] = s+1;
											nrOfValues++;
										}
									}
									/**
									 * determine width Z values - method 1 
									 * */
									{
										/**
										 * Searching LUT for best matches in planes
										 * */
										radiusZ = new double [nrOfValues][2];
										radiusZPrecision = new double [nrOfValues][2];
										keepPlane = new boolean [nrOfValues];
										for(int wz = 0; wz < nrOfValues; wz++){
											minDifference = Double.POSITIVE_INFINITY;
											minDifferencePos = -1;
											for(int pos = lutCenterPos[(int)radiiNew[wz][1]-1]; pos >= 0; pos--){
												if(LUT[(int)radiiNew[wz][1]][pos] == Double.NEGATIVE_INFINITY)	continue;
												differenceTemp = Math.abs(radiiNew[wz][0] - LUT[(int)radiiNew[wz][1]][pos]);
												if(differenceTemp < minDifference){
													minDifference = differenceTemp;
													minDifferencePos = pos;
												}
											}
//											radiusZ[wz][0] = slicePosition [(int)radiiNew[wz][1]-1] + LUT[0][minDifferencePos];
											radiusZ[wz][0] = LUT[0][minDifferencePos];
											radiusZPrecision[wz][0] = LUTPrecision[(int)radiiNew[wz][1]][minDifferencePos];
											
											keep1 = checkMinPositionInLUT((int)radiiNew[wz][1],minDifferencePos, LUT, radiiNew[wz][0]);
											
											minDifference = Double.POSITIVE_INFINITY;
											minDifferencePos = -1;
											for(int pos = lutCenterPos[(int)radiiNew[wz][1]-1]; pos < LUT[0].length; pos++){
												if(LUT[(int)radiiNew[wz][1]][pos] == Double.NEGATIVE_INFINITY)	continue;
												differenceTemp = Math.abs(radiiNew[wz][0] - LUT[(int)radiiNew[wz][1]][pos]);
												if(differenceTemp < minDifference){
													minDifference = differenceTemp;
													minDifferencePos = pos;
												}
											}
//											widthsZ[wz][1] = slicePosition [(int)widthsNew[wz][1]-1] + LUT[0][minDifferencePos];
											radiusZ[wz][1] = LUT[0][minDifferencePos];
											radiusZPrecision[wz][1] = LUTPrecision[(int)radiiNew[wz][1]][minDifferencePos];
											keep2 = checkMinPositionInLUT((int)radiiNew[wz][1],minDifferencePos, LUT, radiiNew[wz][0]);
											if(keep1 || keep2){
												keepPlane [wz] = true;
//												IJ.log("Kept plane " + wz + " " + keep1 + " - " + keep2);
											}else{
//												IJ.log("Skipped plane " + wz + " " + keep1 + " - " + keep2);
												keepPlane [wz] = false;
											}
										}
										
										/**
										 * ESTIMATING CORRECTION FACTORS
										 * */
										{
											refIndex = -1;
											transientCorrFactors = new double [LUT.length-1];
											Arrays.fill(transientCorrFactors,Double.NaN);
											transientCorrDifferences = new double [LUT.length-1];
											Arrays.fill(transientCorrDifferences,Double.NaN);
											for(int wz = 0; wz < nrOfValues; wz++){
												if(radiiNew [wz][1]-1 == refPlane){
													refIndex = wz;
												}
											}
											if(refIndex != -1){
												for(int wz = 0; wz < nrOfValues; wz++){
													if(wz == refIndex){
														transientCorrFactors [(int)radiiNew[wz][1]-1] = 1.0;
														transientCorrDifferences [(int)radiiNew[wz][1]-1] = 0.0;
														continue;
													}
													minDist = Double.POSITIVE_INFINITY;
													for(int i1 = 0; i1 < 2; i1++){
														for(int i2 = 0; i2 < 2; i2++){
															if(Math.abs(radiusZ[refIndex][i1]-radiusZ[wz][i2]) < minDist){
																minDist = Math.abs(radiusZ[refIndex][i1]-radiusZ[wz][i2]);
																transientCorrFactors [(int)radiiNew[wz][1]-1] = radiusZ[wz][i2]/radiusZ[refIndex][i1];
																transientCorrDifferences [(int)radiiNew[wz][1]-1] = radiusZ[wz][i2]-radiusZ[refIndex][i1];
															}
														}
													}													
												}												
											}
											track.get(i).setCorrectionFactors(transientCorrFactors);
											track.get(i).setCorrectionDifference(transientCorrDifferences);
										}
									}
									
									//determine array size
									arraySize = 0;
									wzArraySize = 0;
									
									for(int wzi = 0; wzi < keepPlane.length; wzi++){
										if(keepPlane [wzi]){
											arraySize++;
											wzArraySize++;
										}
									}
//									finalZIncl = widthsZ.length;
									
									/**
									 * Width Method 1 -> Combining
									 * Also including method for best precision selection from version v0.1.0 on
									 * */
									if(wzArraySize > 0 && arraySize > 1){
										//calculate nr of combinations
										nrOfCombinations = (int) Math.pow(2, wzArraySize);
										zValues = new double [nrOfCombinations][arraySize];
										zValuesPrecision = new double [nrOfCombinations][arraySize];
										
										//create all combinations
										comboCounter = 0;
										wZPsCounter = 0;
										for(int wZPs = 0; wZPs < radiusZ.length; wZPs++){
											if(!keepPlane[wZPs])	continue;
											comboCounter = 0;
											for(int nrC = 0; nrC < nrOfCombinations/Math.pow(2,(wZPsCounter+1)); nrC++){					
												for(int w = 0; w < Math.pow(2,wZPsCounter); w++){
//														combinations [comboCounter] += "0";
													zValues [comboCounter][wZPsCounter] = radiusZ [wZPs][0];	//minus
													zValuesPrecision [comboCounter][wZPsCounter] = radiusZPrecision [wZPs][0];	//minus
													comboCounter++;
												}
												for(int w = 0; w < Math.pow(2,wZPsCounter); w++){
//														combinations [comboCounter] += "1";
													zValues [comboCounter][wZPsCounter] = radiusZ [wZPs][1];	//plus
													zValuesPrecision [comboCounter][wZPsCounter] = radiusZPrecision [wZPs][1];	//plus
													comboCounter++;
												}
											}
											wZPsCounter ++;
										}									
										{
//												IJ.log("ID " + i + " comboCounter " + comboCounter + " arraylength " + nrOfCombinations);
//												String output = "";
//												for(int iii = 0; iii < nrOfCombinations; iii++){
//													output = "ID " + i + "combo " + iii + "";
//													for(int jjj = 0; jjj < arraySize; jjj++){
//														output += " " + zValues [iii][jjj];
//													}		
//													IJ.log(output);
//												}
										}									
										
										minSD = Double.MAX_VALUE;
										selCombo = -1;
										minPrec = Double.POSITIVE_INFINITY;
										bestPrecZ = -1;
//											IJ.log("combinations");
										for(int nrC = 0; nrC < nrOfCombinations; nrC++){
//												IJ.log("combination " + combinations [nrC]);
//												IJ.log("zvPre " + zValues [nrC][arraySize-1]);
											
											//calculate SD
											SD = tools.getSD(zValues[nrC]);
											if(SD < minSD){
												minSD = SD;
												selCombo = nrC;
											}
										}
										
										try{											
//											finalZ = tools.getMedian(zValues[selCombo]);
											finalZ = tools.getAverage(zValues[selCombo]);
											if(finalZ == Double.NEGATIVE_INFINITY)	IJ.log("finalZ is negative infinity");
											if(finalZ == Double.POSITIVE_INFINITY)	IJ.log("finalZ is positive infinity");
											track.get(i).setZ(finalZ);
											
											//select z with best precision
											for(int pr = 0; pr < zValues[selCombo].length; pr++){
												if(zValuesPrecision[selCombo][pr] < minPrec) {
													minPrec = zValuesPrecision[selCombo][pr];
													bestPrecZ = pr;
												}
											}
											track.get(i).setZByPrecise(zValues[selCombo][bestPrecZ]);											
											
										}catch(Exception e){
											out = "";
											for(int err = 0; err < e.getStackTrace().length; err++){
												out += " \n " + e.getStackTrace()[err].toString();
											}
											progress.notifyMessage("Error combo: " + out +" \n"
													+ ("ID " + i + "nrOfC " + nrOfCombinations + " - selCombo " + selCombo
															+ " - minSD " + minSD + " - widthZ[0][0] " + (radiusZ[0][0]) + " - widthZ[0][1]" + (radiusZ[0][1]) 
															+ " - comboCounter " + comboCounter + " - widthZLength " + radiusZ.length
															+ " - arraysize " + arraySize + " - min precision " + minPrec + " - best precision index " + bestPrecZ), 
											multiFocalParticleTracker_Complex_3.jnhsupport.ProgressDialog.LOG);
										}
										
//											IJ.log("ID " + i + "selected Combo: " + selCombo);
										selectedZCombo = "";
										selectedZComboPrecision = "";
										Arrays.fill(selZCombo,Double.NaN);
										Arrays.fill(selZComboPrecision,Double.NaN);
										for(int r = 0; r < zValues[selCombo].length; r++){
											selectedZCombo += (", " + zValues[selCombo][r]);
											selZCombo [(int)(radiiNew[r][1]) - 1] = zValues[selCombo][r];
											selectedZComboPrecision += (", " + zValuesPrecision[selCombo][r]);
											selZComboPrecision [(int)(radiiNew[r][1]) - 1] = zValuesPrecision[selCombo][r];
										}
										if(selectedZCombo.length()>0)	selectedZCombo = selectedZCombo.substring(2);
										if(selectedZComboPrecision.length()>0)	selectedZComboPrecision = selectedZComboPrecision.substring(2);
										track.get(i).setSelectedZCombo(selectedZCombo);
										track.get(i).setSelectedZComboAsArray(selZCombo);
										track.get(i).setSelectedZComboPrecisions(selectedZComboPrecision);
										track.get(i).setSelectedZComboAsArrayPrecisions(selZComboPrecision);
									}else{
//										track.get(i).setZByWidthsMethod2(Double.NaN);
//										progress.notifyMessage("X " + constants.df0.format(track.get(i).X()) 
//										+ " Y " + constants.df0.format(track.get(i).Y())
//										+ " ID " + i + " - to few appropriate minima for width measurement (" + wzArraySize
//										+ ")",ProgressDialog.LOG);
									}							
									
									/**
									 * determine Z - method 2 = LIKELIHOOD
									 * */
									nrOfIncludedPlanes = wzArraySize;
									if(nrOfIncludedPlanes > 1){										
										track.get(i).setNrOfInvolvedPlanes(nrOfIncludedPlanes);
										stepSizeInLut = Math.abs(LUT[0][0] - LUT[0][1]);
										Arrays.fill(comparPlot[0], 1.0);
										Arrays.fill(comparPlot[1], 0.0);
										if(savePlot)	diffPlot = new double [radiiNew.length+1][LUT[0].length];
										LUTValue = Double.NaN;;
										includedPlanes = "";									
										for(int wz = 0; wz < radiiNew.length; wz++){
											includedPlanes += constants.df0.format(radiiNew[wz][1]) + "	";
											if(savePlot)	Arrays.fill(diffPlot[wz], Double.NaN);
											for(int pos = 0; pos < LUT[0].length; pos++){
												if(LUT[(int)radiiNew[wz][1]][pos] != Double.NEGATIVE_INFINITY){
													LUTValue = LUT[(int)radiiNew[wz][1]][pos];
												}else{											
													searchInLut: for(int lutsearch = 0; lutsearch < LUT[0].length; lutsearch++){
														if(pos+lutsearch<LUT[0].length){
															if(LUT[(int)radiiNew[wz][1]][pos+lutsearch] != Double.NEGATIVE_INFINITY){
																LUTValue = LUT[(int)radiiNew[wz][1]][pos+lutsearch];
//																IJ.log("lutvalue " + LUTValue);
																break searchInLut;
															}
														}
														if(pos-lutsearch>=0){
															if(LUT[(int)radiiNew[wz][1]][pos-lutsearch] != Double.NEGATIVE_INFINITY){
																LUTValue = LUT[(int)radiiNew[wz][1]][pos-lutsearch];
//																IJ.log("lutvalue " + LUTValue);
																break searchInLut;
															}
														}
													}												
												}
												if(!Double.isNaN(LUTValue)){
													differenceTemp = Math.abs(radiiNew[wz][0] - LUTValue);
													if(savePlot)	diffPlot [wz][pos] = differenceTemp;													
													comparPlot [0][pos] *= (1.0+differenceTemp);
													comparPlot [1][pos] += 1.0;
													LUTValue = Double.NaN;
												}else{
//													IJ.log("DOUBLE NAN " + LUT[(int)widthsNew[wz][1]][pos]);
												}
											}
										}									
										includedPlanes.substring(0, includedPlanes.length()-1);
										while(includedPlanes.split("	",-1).length < 3){
											includedPlanes += "	";											
										}
										track.get(i).setInvolvedPlanes(includedPlanes);									
										
										minValue = Double.POSITIVE_INFINITY;
										index = -1; 
										for(int j = 0; j < comparPlot[0].length; j++){
//											if(comparPlot[1][j] > 2.0){
												if(comparPlot[0][j] < minValue){
													if(j-1>0 && comparPlot[0][j-1] > comparPlot[0][j]){
														if(j+1<comparPlot[0].length && comparPlot[0][j+1] > comparPlot[0][j]){
															minValue = comparPlot [0][j];
															index = j;
														}
													}
												}
//											}
										}
										if(index > 0 && index < comparPlot[0].length){
											/**
											 * determine Z - method 2 - save
											 * */
																						
											//Tracking in method 0.0.2/0.0.3
											track.get(i).setZByRadiusMethod2(index*stepSizeInLut + LUT[0][0]);
											
											/**
											 * Output Plots
											 * */
											if(savePlot){
												for(int dp = 0; dp < diffPlot[0].length; dp++){
													diffPlot[radiiNew.length][dp] = comparPlot[0][dp];
												}
												output1DArray(comparPlot[0], name [task], nameImage, savePath 
														+ "_X_" + constants.df0.format(track.get(i).X()) 
														+ "_Y_" + constants.df0.format(track.get(i).Y())
														+ "_ID_" + i 
														+ ".txt");
												plot2DArray(comparPlot, 
														"X" + constants.df0.format(track.get(i).X()) 
															+ " Y" + constants.df0.format(track.get(i).Y())
															+ " ID" + i, 
														"x",
														"y", 
														savePath
															+ "_X_" + constants.df0.format(track.get(i).X()) 
															+ "_Y_" + constants.df0.format(track.get(i).Y())
															+ "_ID_" + i,
														true);
												
												output2DArray(diffPlot, name [task], nameImage, savePath
														+ "_X_" + constants.df0.format(track.get(i).X()) 
														+ "_Y_" + constants.df0.format(track.get(i).Y())
														+ "_ID_" + i + "_differences.txt");
												
												plot2DArray(diffPlot, 
														"X" + constants.df0.format(track.get(i).X()) 
															+ " Y" + constants.df0.format(track.get(i).Y())
															+ " ID" + i, 
														"x",
														"y", 
														savePath
															+ "_X_" + constants.df0.format(track.get(i).X()) 
															+ "_Y_" + constants.df0.format(track.get(i).Y())
															+ "_ID_" + i + "_differences",
														true);
											}											
										}else{
//											track.get(i).setZByWidthsMethod2(Double.NaN);
//											IJ.log("ID " + i + " not determined: min pos " + index);
										}
									}else{
//										track.get(i).setZByWidthsMethod2(Double.NaN);
//										out = "";
//										for(int kp = 0; kp < keepPlane.length; kp++){
//											if(keepPlane[kp])
//											out += "" + (kp+1) + " ";
//										}
//										progress.notifyMessage("X " + constants.df0.format(track.get(i).X()) 
//										+ " Y " + constants.df0.format(track.get(i).Y())
//										+ " ID " + i + " - to few planes available (" + track.get(i).nrOfInvolvedPlanes + " planes kept: " + out + ")",ProgressDialog.LOG);
										
									}							
								}
							}catch(Exception e){
								out = "";
								for(int err = 0; err < e.getStackTrace().length; err++){
									out += " \n " + e.getStackTrace()[err].toString();
								}
								progress.notifyMessage("Exception for particle X " + constants.df6US.format(track.get(i).X()) 
									+ " Y " + constants.df6US.format(track.get(i).Y())
									+ " ID " + constants.df0.format(i) 
									+ "Cause: " + e.getCause()+ "\n"	
									+ "Throwable: " + e.getMessage() + "\n"	
									+ out, 
									multiFocalParticleTracker_Complex_3.jnhsupport.ProgressDialog.LOG);
							}
							
						}			
						progress.addToBar(0.1);
						
						/**
						 * Estimate Correction Factors
						 * */
						int LUTCorrBins = (int)(LUT[0].length/(double)upscaling/binningFactorForCorrFactors)+1;
						correctionFactors = new double [LUT.length-1][LUTCorrBins];
						correctionDifferences = new double [LUT.length-1][LUTCorrBins];
						double LUTOffset = LUT[0][0], LUTStepsize = LUT [0][1] - LUT [0][0], LUTMax = LUT[0][LUT[0].length-1], LUTSpan = LUTMax-LUTOffset;
						progress.notifyMessage("LUT offset " + constants.df6US.format(LUTOffset) 
							+ " LUT max " + constants.df6US.format(LUTMax)  
							+ " LUT span " + constants.df6US.format(LUTSpan)  
							+ "	LUT stepsize " + constants.df6US.format(LUTStepsize)
							+ " LUTCorrBins " + constants.df0.format(LUTCorrBins),
							ProgressDialog.LOG);
						{
							double transientCorrectionFactors [][] = new double [LUTCorrBins][track.size()];
							double transientCorrectionDifferences [][] = new double [LUTCorrBins][track.size()];
							int transientCorrectionFactorsCounter [] = new int [LUTCorrBins];
							int index2;
							for(int s = 0; s < LUT.length-1; s++){
								for(int cb = 0; cb < transientCorrectionFactors.length; cb++){
									Arrays.fill(transientCorrectionFactors [cb],Double.POSITIVE_INFINITY);
									Arrays.fill(transientCorrectionDifferences [cb],Double.POSITIVE_INFINITY);
								}
								Arrays.fill(transientCorrectionFactorsCounter, 0);
								
								for(int i = 0; i < track.size(); i++){
									if(track.get(i).correctionFactors.length > 0
											&& !Double.isNaN(track.get(i).correctionFactors[s])
											&& track.get(i).selZComboAsArray.length > 0
											&& !Double.isNaN(track.get(i).selZComboAsArray[s])
											&& track.get(i).correctionDifference.length > 0){
										index2 = (int)((track.get(i).selZComboAsArray[s]-LUTOffset)/LUTStepsize/LUTSpan/(double)upscaling/binningFactorForCorrFactors*(LUTCorrBins-1.0));
										transientCorrectionFactors [index2][i] = track.get(i).correctionFactors[s];
										transientCorrectionDifferences [index2][i] = track.get(i).correctionDifference[s];
										transientCorrectionFactorsCounter [index2]++;
									}
								}
								for(int cb = 0; cb < transientCorrectionFactors.length; cb++){
									Arrays.sort(transientCorrectionFactors[cb]);
									Arrays.sort(transientCorrectionDifferences[cb]);
									if(transientCorrectionFactorsCounter[cb]>0){
										correctionFactors [s][cb] = tools.getAverageOfRange(transientCorrectionFactors[cb],0,transientCorrectionFactorsCounter[cb]-1);
										correctionDifferences [s][cb] = tools.getAverageOfRange(transientCorrectionDifferences[cb],0,transientCorrectionFactorsCounter[cb]-1);
									}else{
										correctionFactors [s][cb] = Double.NaN;
										correctionDifferences [s][cb] = Double.NaN;
									}									
								}								
//								output1DArray(correctionFactors[s], name [task], nameImage, savePath + "_CFact_p" + s);
							}
							
							//Plot correction maps
							plot2DArray(correctionFactors, "Correction Factors", "Position in LUT","Correction Factor", savePath + "_CFact", true);
							output2DArray(correctionFactors, name [task], nameImage, savePath + "_CFact_l");
							plot2DArray(correctionDifferences, "Correction Differences", "Position in LUT","Correction Differences", savePath + "_CDiff", true);
							output2DArray(correctionDifferences, name [task], nameImage, savePath + "_CDiff_l");
						}
						
						/**
						 * Apply Correction Factors
						 * */
						{
							double estimFactor, estimDiff;
							int estimCounter;
							int index2;
							for(int i = 0; i < track.size(); i++){
								estimFactor = 0.0;
								estimDiff = 0.0;
								estimCounter = 0;
								for(int s = 0; s < LUT.length-1; s++){	
									if(track.get(i).selZComboAsArray.length > 0
											&& !Double.isNaN(track.get(i).selZComboAsArray[s])){
										index2 = (int)((track.get(i).selZComboAsArray[s]-LUTOffset)/LUTStepsize/LUTSpan/(double)upscaling/binningFactorForCorrFactors*(LUTCorrBins-1.0));
										if(!Double.isNaN(correctionFactors[s][index2])){
											estimFactor += track.get(i).selZComboAsArray[s] / correctionFactors[s][index2];
											estimDiff += track.get(i).selZComboAsArray[s] - correctionDifferences[s][index2];
											estimCounter ++;
										}
									}
								}
								if(estimCounter!=0){
									track.get(i).setCorrectedZAvgFact(estimFactor/estimCounter);
									track.get(i).setCorrectedZAvgDiff(estimDiff/estimCounter);
								}
								
							}
						}
						
						
						
					}
					System.gc();
					
					//save into file
					  	tp = new TextPanel("Results");
					  	tp2 = new TextPanel("Point list");
					  	
					  	tp.append("Saving date:	" + constants.dateTab.format(new Date()) + "	Analysis started:	" + constants.dateTab.format(startDate));
						tp.append("Processed file:	" + name [task]);
						tp.append("Processed image:	" + nameImage);
						tp.append("LUT file:	" + nameLUT + "	LUT SD file:	" + nameLUTSD);
						tp.append("");
						tp.append("Settings: ");
						tp.append("	" + "xy calibration [um]:	" + constants.df6US.format(xyCal));
						tp.append("	" + "Reference plane for correction-factor estimation (>= 1 and  <= nr of planes:	" + constants.df0.format(refPlane));
						tp.append("	" + "scaling of lut:	" + constants.df0.format(upscaling));
						tp.append("	" + "maximum radius size:	" + constants.df0.format(maxRadius));
						tp.append("	" + "width threshold min:	" + constants.df3US.format(widthThresholdMin));
						tp.append("	" + "width threshold max:	" + constants.df3US.format(widthThresholdMax));
						tp.append("");
						
						String appendTxt = "";
						appendTxt = "index	T	X (µm)	Y (µm)	Z average (µm)	factor-corr. Z avg (µm)	diff-corr. Z avg (µm)	Z from precisest plane (µm)"
								+ "	" + "radius plane 1 (µm)	radius plane 2 (µm)	radius plane 3 (µm)	radius plane 4 (µm)	chosen combination"
								+ "	" + "radius-derived Z (µm) (method 2)" 
								+ "	" + "number of used planes (method 2)"
								+ "	" + "IDs of used planes (method 2)			"
								+ "	" + "z plane 1 (µm)	z plane 2 (µm)	z plane 3 (µm)	z plane 4 (µm)"
								+ "	" + "z precision plane 1 (µm)	z precision plane 2 (µm)	z precision plane 3 (µm)	z precision plane 4 (µm)";
						tp.append(appendTxt);
						tp2.append(appendTxt);
						
						for(int i = 0; i < track.size(); i++){
							appendTxt = "" + constants.df0.format(i);
							appendTxt += "	" + constants.df6US.format(track.get(i).T());
							appendTxt += "	" + constants.df6US.format(track.get(i).X());
							appendTxt += "	" + constants.df6US.format(track.get(i).Y());
							
							appendTxt += "	";
							if(track.get(i).zDef()){
								appendTxt += constants.df6US.format(track.get(i).Z());
							}
														
							appendTxt += "	";
							if(track.get(i).zCorrDefFact()){
								appendTxt += constants.df6US.format(track.get(i).zAvgCorrectedFact());
							}
							
							appendTxt += "	";
							if(track.get(i).zCorrDefDiff()){
								appendTxt += constants.df6US.format(track.get(i).zAvgCorrectedDiff());
							}
							
							appendTxt += "	";
							if(track.get(i).zPreciseDef()){
								appendTxt += constants.df6US.format(track.get(i).zByPreciseMethod());
							}
							
							for(int s = 0; s < imp.getNSlices(); s++){
								appendTxt += "	";
								if(track.get(i).radiusCentred(s) != 0.0){
									appendTxt += constants.df6US.format(track.get(i).radiusCentred(s));
								}
							}
							
							appendTxt += "	";
							if(!track.get(i).selectedZCombo.equals("")){
								appendTxt += track.get(i).selectedZCombo;
							}
							
							appendTxt += "	";
							if(!Double.isNaN(track.get(i).zByRadiusCentMethod2())){
								appendTxt += constants.df6US.format(track.get(i).zByRadiusCentMethod2());
							}
							
							appendTxt += "	";
							if(!Double.isNaN(track.get(i).zByRadiusCentMethod2())){
								appendTxt += constants.df0.format(track.get(i).getNrOfInvolvedPlanes());
							}
							
							appendTxt += "	";
							if(!Double.isNaN(track.get(i).zByRadiusCentMethod2())){
								appendTxt += track.get(i).involvedPlanes;
							}else{
								appendTxt += "			";
							}
							
							
							if(track.get(i).selZComboAsArray.length!=0){
								for(int j = 0; j < track.get(i).selZComboAsArray.length; j++){
									appendTxt += "	";
									if(!Double.isNaN(track.get(i).selZComboAsArray[j]))	appendTxt += constants.df6US.format(track.get(i).selZComboAsArray[j]);	
								}
								for(int j = 0; j < imp.getNSlices() - track.get(i).selZComboAsArray.length; j++) {
									appendTxt += "	";
								}
							}else{
								for(int j = 0; j < imp.getNSlices(); j++){
									appendTxt += "	";										
								}
							}
							
							if(track.get(i).selZComboAsArrayPrecision.length!=0){
								for(int j = 0; j < track.get(i).selZComboAsArrayPrecision.length; j++){
									appendTxt += "	";
									if(!Double.isNaN(track.get(i).selZComboAsArrayPrecision[j]))	appendTxt += constants.df6US.format(track.get(i).selZComboAsArrayPrecision[j]);	
								}
								for(int j = 0; j < imp.getNSlices() - track.get(i).selZComboAsArray.length; j++) {
									appendTxt += "	";
								}
							}else{
								for(int j = 0; j < imp.getNSlices(); j++){
									appendTxt += "	";										
								}
							}
							
							tp.append(appendTxt);
							tp2.append(appendTxt);
							progress.addToBar(0.1/(double)track.size());
						}					
						tp.append("");
						
						addFooter(tp);	
					  	tp.saveAs(savePath + ".txt");
//					  	addFooter(tp2);
					  	tp2.saveAs(savePath + "_s.txt");
					  	
					//plot
//						plotTrace(imp, track, slicesSorted, savePath + ".tif");			
				}
				
				//save progress dialog log file
				  	progress.saveLog(savePath + "_log.txt");
				  
			  	//finish progress dialog				  	
				  	progress.setBar(1.0);
				  	done = true;
				  	tasksSuccessfull [task] = true;
				  	break running;		  	
			}//(end runnning)
			
			if(progress.isStopped()) break tasking;
			progress.moveTask(task);			
		}
		imp.changes = false;
	  	imp.close();
	}
	
	static void addFooter (TextPanel tp){
		tp.append("");
		tp.append("Datafile was generated by '"+PLUGINNAME+"', (\u00a9 2017-" + constants.dateY.format(new Date()) + ": Jan N Hansen (jan.hansen@uni-bonn.de, https://github.com/hansenjn/MultiFocalParticleTracker-Complex-3))");
		tp.append("The plug-in '"+PLUGINNAME+"' is distributed in the hope that it will be useful,"
				+ " but WITHOUT ANY WARRANTY; without even the implied warranty of"
				+ " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.");		
		tp.append("Plug-in version:	"+PLUGINVERSION);	  	
	}
	
	private ArrayList<MFPTCPoint3> readTrackList(String filePath, int nSlices){
		ArrayList<MFPTCPoint3> track = new ArrayList<MFPTCPoint3>();
		track.ensureCapacity(400);
		{
			//read information about width max and min
			String px, py, frame;
			try {
				FileReader fr = new FileReader(new File(filePath));
				BufferedReader br = new BufferedReader(fr);
				String line = "";							
				reading: while(true){
					try{
						line = br.readLine();	
						if(line.equals(null)){
							progress.updateBarText("track of " + track.size() + " points detected - now analyzing each point...");
							break reading;
						}
					}catch(Exception e){
						progress.updateBarText("track of " + track.size() + " points detected - now analyzing each point...");
						break reading;
					}					
					
					//check if comma instead of points and if so replace "," by "."
					frame = (line.substring(0,line.indexOf("	")));
					px = (line.substring(line.indexOf("	")+1,line.lastIndexOf("	")));
					py = (line.substring(line.lastIndexOf("	")+1));
					if(frame.contains(",") && !frame.contains("."))	frame = frame.replace(",", ".");
					if(px.contains(",") && !px.contains("."))	px = px.replace(",", ".");
					if(py.contains(",") && !py.contains("."))	py = py.replace(",", ".");
					
					track.add(new MFPTCPoint3(Double.parseDouble(px),Double.parseDouble(py),Integer.parseInt(frame), nSlices));							
				}					
				br.close();
				fr.close();
			}catch (IOException e) {
				IJ.error("Problem with text loading");
				e.printStackTrace();
				return null;
			}
		}
		track.trimToSize();
		return track;
	}
		
	private ArrayList<MFPTCPoint3> readTrackListInVerificationMode(String filePath, int nSlices, int nFrames){
		ArrayList<MFPTCPoint3> track = new ArrayList<MFPTCPoint3>();
		track.ensureCapacity(400);
		{
			//read information about width max and min
			String px, py;
			try {
				FileReader fr = new FileReader(new File(filePath));
				BufferedReader br = new BufferedReader(fr);
				String line = "";							
				reading: while(true){
					try{
						line = br.readLine();	
						if(line.equals(null)){
							progress.updateBarText("track of " + track.size() + " points detected");
							break reading;
						}
					}catch(Exception e){
						progress.updateBarText("track of " + track.size() + " points detected");
						break reading;
					}
					
					//check if comma instead of points and if so replace "," by "."
					px = (line.substring(0,line.lastIndexOf("	")));
					py = (line.substring(line.lastIndexOf("	")+1));
					if(px.contains(",") && !px.contains("."))	px = px.replace(",", ".");
					if(py.contains(",") && !py.contains("."))	py = py.replace(",", ".");
					for(int t = 0; t < nFrames; t++){
						track.add(new MFPTCPoint3(Double.parseDouble(px),Double.parseDouble(py),t,nSlices));	
					}											
				}					
				br.close();
				fr.close();
			}catch (IOException e) {
				IJ.error("Problem with text loading");
				e.printStackTrace();
				return null;
			}
		}
		track.trimToSize();
		return track;
	}
	
	/**
	 * Used as template LUT until 17.12.2019
	 * */
	private double [][] readLUT (String filePath){
		ArrayList<String> lines = new ArrayList<String>();
		lines.ensureCapacity(400);
		String line = "";
		{
			//read information from file
			try {
				FileReader fr = new FileReader(new File(filePath));
				BufferedReader br = new BufferedReader(fr);
				reading: while(true){
					try{
						line = br.readLine();	
						if(line.equals(null)){
							progress.updateBarText("Detected " + lines.size() + " lines in LUT!");
							progress.notifyMessage("Detected " + lines.size() + " lines in LUT!", ProgressDialog.LOG);
							break reading;
						}else{
							lines.add(line);
						}
					}catch(Exception e){
//						String out = "";
//						for(int err = 0; err < e.getStackTrace().length; err++){
//							out += " \n " + e.getStackTrace()[err].toString();
//						}
//						progress.notifyMessage("LUT error: " + out, ProgressDialog.ERROR);
						progress.updateBarText("Detected " + lines.size() + " lines in LUT!");
						progress.notifyMessage("Detected " + lines.size() + " lines in LUT!", ProgressDialog.LOG);
						break reading;
					}					
				}					
				br.close();
				fr.close();
			}catch (IOException e) {
				IJ.error("Problem with text loading");
				e.printStackTrace();
				return null;
			}
		}
		lines.trimToSize();
		
		double LUT [][] = new double [5][lines.size()];
		{
			//create LUT array
			String temp = "";
			for(int i = 0; i < lines.size(); i++){
				line = lines.get(i);
				for(int j = 0; j < 4; j ++){
					temp = (line.substring(0,line.indexOf("	")));
					line = line.substring(line.indexOf("	")+1);
					if(temp.contains(",") && !temp.contains("."))	temp = temp.replace(",", ".");
					if(temp=="" || temp.isEmpty()){
						LUT [j][i] = Double.NEGATIVE_INFINITY;
					}else{
						LUT [j][i] = Double.parseDouble(temp);
					}
//					IJ.log("LUT " + i + " - " + j + ": " + LUT[j][i]);
				}
				if(line.contains(",") && !line.contains("."))	line = line.replace(",", ".");
				if(line=="" || line.isEmpty()){
					LUT [4][i] = Double.NEGATIVE_INFINITY;
				}else{
					LUT [4][i] = Double.parseDouble(line);
				}				
//				IJ.log("LUT " + i + " - 4: " + LUT[4][i]);
			}
		}
		lines.clear();
		lines = null;
		System.gc();
		return LUT;
	}
	
	/**
	 * Used from 17.12.2019
	 * */
	private double [][] readUpscaledLUT (String filePath, int scaleFactor){
		ArrayList<String> lines = new ArrayList<String>();
		lines.ensureCapacity(400);
		String line = "";
		{
			//read information from file
			try {
				FileReader fr = new FileReader(new File(filePath));
				BufferedReader br = new BufferedReader(fr);
				reading: while(true){
					try{
						line = br.readLine();	
						if(line.equals(null)){
							progress.updateBarText("Detected " + lines.size() + " lines in LUT!");
							progress.notifyMessage("Detected " + lines.size() + " lines in LUT!", ProgressDialog.LOG);
							break reading;
						}else{
							lines.add(line);
						}
					}catch(Exception e){
//						String out = "";
//						for(int err = 0; err < e.getStackTrace().length; err++){
//							out += " \n " + e.getStackTrace()[err].toString();
//						}
//						progress.notifyMessage("LUT error: " + out, ProgressDialog.ERROR);
						progress.updateBarText("Detected " + lines.size() + " lines in LUT!");
						progress.notifyMessage("Detected " + lines.size() + " lines in LUT!", ProgressDialog.LOG);
						break reading;
					}					
				}					
				br.close();
				fr.close();
			}catch (IOException e) {
				IJ.error("Problem with text loading");
				e.printStackTrace();
				return null;
			}
		}
		lines.trimToSize();
		
		double LUT [][] = new double [5][lines.size()];
		{
			//create LUT array
			String temp = "";
			for(int i = 0; i < lines.size(); i++){
				line = lines.get(i);
				for(int j = 0; j < 4; j ++){
					temp = (line.substring(0,line.indexOf("	")));
					line = line.substring(line.indexOf("	")+1);
					if(temp.contains(",") && !temp.contains("."))	temp = temp.replace(",", ".");
					if(temp=="" || temp.isEmpty()){
						LUT [j][i] = Double.NaN;
					}else{
						LUT [j][i] = Double.parseDouble(temp);
					}
//					IJ.log("LUT " + i + " - " + j + ": " + LUT[j][i]);
				}
				if(line.contains(",") && !line.contains("."))	line = line.replace(",", ".");
				if(line=="" || line.isEmpty()){
					LUT [4][i] = Double.NaN;
				}else{
					LUT [4][i] = Double.parseDouble(line);
				}				
//				IJ.log("LUT " + i + " - 4: " + LUT[4][i]);
			}
		}
		lines.clear();
		lines = null;
		System.gc();
		
		LUT = upscaleLutBySplineFitter(LUT, scaleFactor);
				
		return LUT;
	}
	
	/**
	 * Upscale LUT using a Spline Fit
	 * */
	private static double [][] upscaleLutBySplineFitter(double inLUT [][], int scaleFactor){		
		double LUT [][] = new double [inLUT.length][inLUT[0].length*scaleFactor-scaleFactor];
		double pos;
		SplineFitter splF;
		ArrayList<float []> points = new ArrayList<float []>(inLUT[0].length);
		for(int p = 1; p < inLUT.length; p++){
			points.clear();
			points.ensureCapacity(inLUT[0].length);
			for(int i = 0; i < inLUT [p].length; i++){
				if(!Double.isNaN(inLUT[p][i])){
					points.add(new float [] {(float) inLUT[0][i], (float) inLUT[p][i]});
				}
			}
			{
				float flLUT [][] = new float [2][points.size()];
				for(int i = 0; i < flLUT[0].length; i++){
					flLUT [0][i] = points.get(i) [0];
					flLUT [1][i] = points.get(i) [1];
				}			
				splF = new SplineFitter(flLUT[0], flLUT[1], flLUT[0].length);
			}
			
			
			Arrays.fill(LUT[p], Double.NEGATIVE_INFINITY);
			for(int i = 0; i < inLUT[p].length-1; i++){
				if(Double.isNaN(inLUT[p][i]))	continue;
				if(Double.isNaN(inLUT[p][i+1]))	continue;
				
				for(int j = 0; j < scaleFactor; j++){
					pos = inLUT[0][i] + (inLUT[0][i+1]-inLUT[0][i]) * ((double)j/(double)scaleFactor);
					LUT [0][i*scaleFactor+j] = pos;
					LUT [p][i*scaleFactor+j] = splF.evalSpline(pos);
				}
			}
			if(!Double.isNaN(inLUT [p][inLUT[p].length-1]))	LUT [p][LUT[p].length-1] = inLUT [p][inLUT[p].length-1];
		}		
		
		return LUT;		
	}
		
	private boolean checkMinPositionInLUT(int column, int row, double [][] LUT, double widthValue){
		double value = Math.abs(LUT[column][row]-widthValue);
		//check up
		int counter = 0;
		for(int add = 0; add < 6; add++){
			if(row+add<LUT[column].length){
				if(LUT[column][row+add]!=Double.NEGATIVE_INFINITY){
					if(Math.abs(LUT[column][row+add]-widthValue)>value){
						counter++;
					}					
				}
			}			
		}
		if(counter<5){
//			IJ.log("column " +  column + ", row " + row + ", value " + value);
			return false;
		}
		//check down
		counter = 0;
		for(int add = 0; add < 6; add++){
			if(row-add>=0){
				if(LUT[column][row-add]!=Double.NEGATIVE_INFINITY){
					if(Math.abs(LUT[column][row-add]-widthValue)>value){
						counter++;
					}					
				}
			}
		}
		if(counter<5){
//			IJ.log("column " +  column + ", row " + row + ", value " + value);
			return false;
		}
		return true;
	}
	
	/**
	 * Calculates the precision of the LUT
	 * (New from version 0.1.0 on)
	 * @author Luis Alvarez, translated from Matlab to Java by Jan N. Hansen
	 * @param LUT = the LUT array
	 * @param sdLUT = the array describing the standard deviation of the LUT, must be same dimension as the LUT array.
	 * */
	private double [][] getLUTPrecision (double LUT [][], double sdLUT [][]){
		double [][] precisionLUT = new double [LUT.length][LUT[0].length];
		if(LUT[0].length != sdLUT[0].length) {
			progress.notifyMessage("ERROR: LUT and LUT sd files do not contain same number of lines. Cannot calculate precisions.", ProgressDialog.ERROR);
			return precisionLUT;
		}
		if(LUT.length != sdLUT.length) {
			progress.notifyMessage("ERROR: LUT and LUT sd files do not contain same number of columns. Cannot calculate precisions.", ProgressDialog.ERROR);
			return precisionLUT;
		}
		for(int i = 0; i < precisionLUT[0].length; i++) {
			if(LUT[0][i] != sdLUT[0][i]) {
				progress.notifyMessage("ERROR: LUT positions do not match sdLUT positions. Cannot calculate precisions.", ProgressDialog.ERROR);
				return precisionLUT;
			}
			precisionLUT[0][i] = LUT [0][i];
		}
		
		try {
			//Derive derivative of the LUT
			for(int i = 1; i < precisionLUT.length; i++) {
				precisionLUT [i][0] = (LUT [i][1] - LUT [i][0]) / (LUT [0][1] - LUT [0][0]);
				for(int j = 1; j < precisionLUT[i].length-1; j++) {
					precisionLUT [i][j] = (LUT [i][j+1] - LUT [i][j]) / (LUT [0][j+1] - LUT [0][j]);
					precisionLUT [i][j] += (LUT [i][j] - LUT [i][j-1]) / (LUT [0][j] - LUT [0][j-1]);
					precisionLUT [i][j] /= 2.0;
				}
				precisionLUT [i][precisionLUT[i].length-1] = (LUT [i][LUT[i].length-1] - LUT [i][LUT[i].length-2]) 
						/ (LUT [0][LUT[i].length-1] - LUT [0][LUT[i].length-2]);
			}
		}catch(Exception e) {
			String out = "";
			for(int err = 0; err < e.getStackTrace().length; err++){
				out += " \n " + e.getStackTrace()[err].toString();
			}
			progress.notifyMessage("ERROR in LUT interpretation: " + out, ProgressDialog.ERROR);
			return precisionLUT;
		}	
		
		try {
			//Convert derivative to precision by multiplication of the absolute derivative with the standard deviation of the LUT
			for(int i = 1; i < precisionLUT.length; i++) {
				for(int j = 0; j < precisionLUT[i].length; j++) {
					precisionLUT [i][j] = Math.abs(precisionLUT[i][j])*sdLUT[i][j];
				}
			}
		}catch(Exception e) {
			String out = "";
			for(int err = 0; err < e.getStackTrace().length; err++){
				out += " \n " + e.getStackTrace()[err].toString();
			}
			progress.notifyMessage("ERROR in LUT-SD interpretation: " + out, ProgressDialog.ERROR);
			return precisionLUT;
		}
			
		return precisionLUT;
	}
	
	private static void output2DArray(double [][] array, String name, String nameImage, String savePath){
		TextPanel tp = new TextPanel("Results");
	  	tp.append("Processed file:	" + name);
		tp.append("Processed image:	" + nameImage);
		tp.append("");
		tp.append("ARRAY:");
		String appendTxt;
		for(int i = 0; i < array.length; i++){
			appendTxt = "";
			for(int j = 0; j < array[0].length; j++){
				if(j!=0)	appendTxt += "	";
				if(!Double.isNaN(array[i][j]) && !(array[i][j]==Double.NEGATIVE_INFINITY)){
					appendTxt += constants.df6US.format(array [i][j]);
				}
			}			
			tp.append(appendTxt);
		}					
		tp.append("");		
		addFooter(tp);	
	  	tp.saveAs(savePath + ".txt");
	  	tp = null;
	  	System.gc();
	}
	
	/**
	 * 1st dimension > different graphs
	 * 2nd dimension > y points
	 * */
	private static void plot2DArray(double [][] array, String label, String xLabel, String yLabel, String savePath, boolean logarithmic){
		double [] xValues = new double [array[0].length];
		for(int i = 0; i < xValues.length; i++)	xValues [i] = i;
		double yMax = Double.NEGATIVE_INFINITY, yMin = 0.0;
		double max, min;
		for(int i = 0; i < array.length; i++){
			max = tools.getMaximum(array[i]);
			if(yMax < max) yMax = max;
			min = tools.getMinimum(array[i]);
			if(yMin > min) yMin = min;
		}
		Color c;
		Plot p;
		ImagePlus pImp;
		String legend = "";
		PlotWindow.noGridLines = true;
		
		p = new Plot(label, xLabel, yLabel);
		p.setAxisYLog(logarithmic);
		p.setSize(600, 400);
		p.setLimits(0, xValues.length-1, yMin, yMax);		
		for(int i = 0; i < array.length; i++){
			c = new Color(54+(int)(i/(double)array.length*200.0), 54+(int)(i/(double)array.length*200.0), 54+(int)(i/(double)array.length*200.0));
			p.setColor(c);
			p.addPoints(xValues,array[i],PlotWindow.LINE);
			legend += "" + i;
			legend += "\n";
		}
		p.addLegend(legend);
		p.setLimitsToFit(true);
		pImp = p.makeHighResolution("plot",1,true,false);
		IJ.saveAs(pImp,"PNG",savePath + ".png");
		pImp.changes = false;
		pImp.close();
		p.dispose();			
	  	System.gc();
	}
	
	/**
	 * 1st dimension > different graphs and in index 0 y axis
	 * 2nd dimension > y points
	 * */
	private static void plotXY2DArray(double [][] array, String label, String xLabel, String yLabel, String savePath, boolean logarithmic){
		double [] xValues = new double [array[0].length];
		for(int i = 0; i < xValues.length; i++)	xValues [i] = array[0][i];
		double xMin = tools.getMinimum(array[0]), xMax = tools.getMaximum(array[0]), yMax = Double.NEGATIVE_INFINITY;
		double max;
		for(int i = 1; i < array.length; i++){
			max = tools.getMaximum(array[i]);
			if(yMax < max) yMax = max;			
		}
		Color c;
		Plot p;
		ImagePlus pImp;
		String legend = "";
		PlotWindow.noGridLines = true;
		
		p = new Plot(label, xLabel, yLabel);
		p.setAxisYLog(logarithmic);
		p.setSize(600, 400);
		p.setLimits(xMin, xMax, 0.0, yMax);		
		for(int i = 1; i < array.length; i++){
			c = new Color(54+(int)(i/(double)array.length*200.0), 54+(int)(i/(double)array.length*200.0), 54+(int)(i/(double)array.length*200.0));
			p.setColor(c);
			p.addPoints(xValues,array[i],PlotWindow.LINE);
			legend += "" + i;
			legend += "\n";
		}
		p.addLegend(legend);
		p.setLimitsToFit(true);
		pImp = p.makeHighResolution("plot",1,true,false);
		IJ.saveAs(pImp,"PNG",savePath + ".png");
		pImp.changes = false;
		pImp.close();
		p.dispose();			
	  	System.gc();
	}
	
	private static void output1DArray(double [] array, String name, String nameImage, String savePath){
		TextPanel tp = new TextPanel("Results");
	  	tp.append("Processed file:	" + name);
		tp.append("Processed image:	" + nameImage);
		tp.append("");
		tp.append("ARRAY:");
		String appendTxt;
		for(int i = 0; i < array.length; i++){ 
			appendTxt = "";
			if(!Double.isNaN(array[i])) appendTxt = constants.df6US.format(array [i]);
			tp.append(appendTxt);
		}					
		tp.append("");		
		addFooter(tp);	
	  	tp.saveAs(savePath + ".txt");
	  	tp = null;
	  	System.gc();
	}
	
	/**
	 * New from 2019-04-22
	 * */
	private int [][] getBinnedDistanceMatrix(){
		int matrix [][] = new int [(maxRadius*2+1)][(maxRadius*2+1)];
		
		double distance;
		for(int x = 0; x < matrix.length; x++){
			for(int y = 0; y < matrix[x].length; y++){
				distance = Math.sqrt((x-maxRadius)*(x-maxRadius)+(y-maxRadius)*(y-maxRadius));
//				if(distance<=maxRadius){
					matrix [x][y] = (int)Math.round(distance);
//				}				
			}
		}
		
		return matrix;	
//		return tools.getMedianOfRange(values, firstIndex, lastIndex);
	}
	
	/**
	 * New from 2019-04-22
	 * */
	private double [][] getUncalibDistanceMatrix(){
		double matrix [][] = new double [(maxRadius*2+1)][(maxRadius*2+1)];		
		double distance;
		for(int x = 0; x < matrix.length; x++){
			for(int y = 0; y < matrix[x].length; y++){
				distance = Math.sqrt(((double)x-(double)maxRadius)*((double)x-(double)maxRadius)+((double)y-(double)maxRadius)*((double)y-(double)maxRadius));
//				if(distance<=maxRadius){
					matrix [x][y] = distance;
//				}				
			}
		}
		
		return matrix;	
//		return tools.getMedianOfRange(values, firstIndex, lastIndex);
	}

	/**
	 * Developed from 02.07. on
	 * @param centreX calibrated x centre position
	 * @param centreY calibrated y centre position
	 * @param slice >= 0 and < nrOfSlices
	 * @param frame >= 0 and < or nrOfFrames
	 * @returns 2D array: first dimension: 0 = x coord (calibrated), 1 = y coord (calibrated), 2 = intensity; 
	 * second dimension: different pixels
	 * */
	private double [][] getIntensitiesInRadius(ImagePlus imp, double centreX, double centreY, int frame, int slice, 
			int [][] binnedDistanceMatrix, double [][] uncalibDistanceMatrix,
			double uncalibRadius) {
		int px = (int)Math.round(centreX/imp.getCalibration().pixelWidth) - maxRadius;
		int py = (int)Math.round(centreY/imp.getCalibration().pixelHeight) - maxRadius;
		if(px >= (double)imp.getWidth()-1 || px < 0
				||py >= (double)imp.getHeight()-1 || py < 0){
			return new double [1][1];
		}	
		
		int counter = 0;
		for(int x = 0; x < binnedDistanceMatrix.length; x++){
			for(int y = 0; y < binnedDistanceMatrix[x].length; y++){
				if(binnedDistanceMatrix[x][y] <= maxRadius){
//					if(uncalibDistanceMatrix[x][y] <= uncalibRadius){	//From v 0.0.8 on commented
						if(px+x < imp.getWidth() && py+y < imp.getHeight()){
							counter++;							
						}	
//					}										
				}				
			}
		}	
		
		double intensities [][] = new double [3][counter]; //x, y, z
		Arrays.fill(intensities[0], 0.0);
		Arrays.fill(intensities[1], 0.0);
		Arrays.fill(intensities[2], 0.0);
			
		counter = 0;
		for(int x = 0; x < binnedDistanceMatrix.length; x++){
			for(int y = 0; y < binnedDistanceMatrix[x].length; y++){
				if(binnedDistanceMatrix[x][y] <= maxRadius){
//					if(uncalibDistanceMatrix[x][y] <= uncalibRadius){		//From v 0.0.8 on commented		
						if(px+x < imp.getWidth() && py+y < imp.getHeight()){
							intensities[0][counter] = (double)(px+x)*imp.getCalibration().pixelWidth;
							intensities[1][counter] = (double)(py+y)*imp.getCalibration().pixelHeight;
							intensities[2][counter] = imp.getStack().getVoxel(px+x, py+y,imp.getStackIndex(1, slice+1, frame+1)-1);
							counter++;							
						}	
//					}							
				}				
			}
		}
		return intensities;		
	}
	
	
}