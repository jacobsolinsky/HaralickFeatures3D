/*
 * To the extent possible under law, the ImageJ developers have waived
 * all copyright and related or neighboring rights to this tutorial code.
 *
 * See the CC0 1.0 Universal license for details:
 *     https://creativecommons.org/publicdomain/zero/1.0/
 */

import java.lang.Math;
import static java.lang.Math.sin;

import java.io.IOException;
import java.io.PrintWriter;

import static java.lang.Math.cos;
import static java.lang.Math.PI;


import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.plugin.*;
import ij.process.*;
import ij.gui.*;
import ij.*;
import ij.measure.Calibration;
import java.util.Arrays;
import java.util.HashMap;


import net.imglib2.img.ImagePlusAdapter;
import net.imagej.ImageJ;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.UnsignedShortType;
import net.imglib2.type.Type;
import net.imglib2.Cursor;
import net.imglib2.view.Views;
import net.imglib2.view.IterableRandomAccessibleInterval;
import net.imglib2.RandomAccess;

import org.scijava.Context;


public class Haralick_Features_3D implements PlugInFilter {
	ImagePlus img;
	int dataSlice = 0;
	int dataDim = 2;
	int maskSlice = 1;
	double baseLength = 0.5;
	int distanceSampleDepth = 10;
	int angleSampleDepth = 8;
	int binno = 20;
	Calibration imgCalibration;
	
	public int setup(String arg, ImagePlus img){
		this.img = img;
		this.imgCalibration = img.getCalibration();

		IJ.register(Haralick_Features_3D.class);
		return DOES_16+NO_CHANGES;
	}

	@Override
	public void run(ImageProcessor ip) {
		String units = imgCalibration.getUnits();
		GenericDialog gd1 = new GenericDialog("Select Algorithm Version and parameters");
	    gd1.addNumericField("Feature size in " + units + ": ", baseLength, 0);
	    gd1.addNumericField("Distance Search Depth", distanceSampleDepth, 0);
	    gd1.addNumericField("Bin number", binno, 0);
	    gd1.addNumericField("Angle Search Depth", angleSampleDepth, 0);
		gd1.showDialog();
	    if (gd1.wasCanceled()) return;
	    this.baseLength = gd1.getNextNumber();
	    this.distanceSampleDepth = (int) gd1.getNextNumber();
	    this.binno = (int) gd1.getNextNumber();
	    this.angleSampleDepth = (int) gd1.getNextNumber();
		circumventJava();
	}
	
	public<T extends RealType< T > & Type< T > > void circumventJava() {
        final Context context = (Context) IJ.runPlugIn("org.scijava.Context", "");
        final ImageJ ij = new ImageJ(context);
        
		Img<UnsignedShortType> img3d = ImagePlusAdapter.wrapNumeric(this.img);
		ImageJFunctions.show(img3d);
		Img<UnsignedShortType> binnedImg = binning(img3d, binno);
		ImageJFunctions.show(binnedImg);
		double[][][] angleSpace = sampleAngles(angleSampleDepth);
		double[][][][] distanceArr = sampleAnglesAndDistances(angleSpace, distanceSampleDepth, baseLength);
		int[][][][] pixelOffsetArr = pixelOffsets(distanceArr);
		int[][][][][] glcmArr = greyLevelCorrelationMatrices(pixelOffsetArr, binnedImg);
		HashMap<double[], int[][]> glcmDict = linkGlcmToParameters(distanceArr, glcmArr);
		HashMap<double[], Double> corrDict = linkCorrelations(glcmDict);
		PrintWriter writer = null;
		try {
			writer = new PrintWriter("C:/Users/van Berlo Lab/Desktop/Jacob/glcm.json", "UTF-8");
			writer.println("[");
			for(double[] ptd: glcmDict.keySet()) {
				int[][] glcm = glcmDict.get(ptd);
				writer.println(
						"{" + 
						"\"phi\": " + Double.toString(ptd[0]) + "," +
						"\"theta\": " + Double.toString(ptd[1]) + "," +
						"\"distance\": " + Double.toString(ptd[2]) + "," +
						"\"glcm\": " + Arrays.deepToString(glcm) +
						"},"
						);
			}
			writer.println("]");
		}
		catch (IOException ex) {
			IJ.log("printing failure?");
		}
		finally {
			writer.close();
		}
		writer = null;
		try {
			writer = new PrintWriter("C:/Users/van Berlo Lab/Desktop/Jacob/corr.json", "UTF-8");
			writer.println("[");
			for(double[] ptd: corrDict.keySet()) {
				double corr = corrDict.get(ptd);
				writer.println(
						"{" + 
						"\"phi\": " + Double.toString(ptd[0]) + "," +
						"\"theta\": " + Double.toString(ptd[1]) + "," +
						"\"distance\": " + Double.toString(ptd[2]) + "," +
						"\"glcm\": " + Double.toString(corr) +
						"},"
						);
			}
			writer.println("]");
		}
		catch (IOException ex) {
			IJ.log("printing failure?");
		}
		finally {
			writer.close();
		}
		writer = null;
		try {
			writer = new PrintWriter("C:/Users/van Berlo Lab/Desktop/Jacob/spheres.csv", "UTF-8");
			writer.println("phi, theta, distance, x, y, z");
			for (int i = 0; i < distanceArr.length; i++) {
				double[][][] iSubArr = distanceArr[i];
				for (int j = 0; j < iSubArr.length; j++) {
					double[][] ijSubArr = iSubArr[j];
					for (int k = 0; k < ijSubArr.length; k ++) {
						double[] ijkSubArr = ijSubArr[k];
						writer.println(
								Double.toString(ijkSubArr[0])  + "," +
								Double.toString(ijkSubArr[1])  + "," +
								Double.toString(ijkSubArr[2])  + "," +
								Double.toString(ijkSubArr[3])  + "," +
								Double.toString(ijkSubArr[4])  + "," + 
								Double.toString(ijkSubArr[5])
								);
					}
				}
			}
		}
		catch (IOException ex) {
			
		}
		finally {
			writer.close();
		}
	}
	
    public Img<UnsignedShortType>
    binning(Img<UnsignedShortType> inputfull, final int binno){
    		
            // create a cursor for the image (the order does not matter)
    		IterableRandomAccessibleInterval<UnsignedShortType> input = 
    				new IterableRandomAccessibleInterval<UnsignedShortType>(
    						Views.hyperSlice(inputfull, dataDim, dataSlice));
    		RandomAccess<UnsignedShortType> maskPointer = 
    				Views.hyperSlice(inputfull, dataDim, maskSlice).randomAccess();
            Cursor<UnsignedShortType> inputCursor = input.cursor();
    		int min = 0;
    		int max = 0;
    		UnsignedShortType type;
     
            // initialize min and max with the first image value
            boolean loopInit = true;
            int totalPixels = 0;
            // loop over the rest of the data and determine min and max value in the masked area
            long[] pos = new long[input.numDimensions()];
            while ( inputCursor.hasNext() ){
            	inputCursor.fwd();
            	inputCursor.localize(pos);
            	maskPointer.setPosition(pos);
            	if (maskPointer.get().getInteger() > 0) {
            		if (loopInit) {
	            		min = inputCursor.get().getInteger();
	            		max = inputCursor.get().getInteger();
	            		loopInit = false;
            		}
                    type = inputCursor.get();
                    totalPixels += 1;
         
                    if ( type.getInteger() < min ) {
                        min = type.getInteger();
                    }
         
                    if ( type.getInteger() > max ) {
                        max = type.getInteger();
                    }
            	}
            	
                // we need this type more than once
            	

            }
            
            //Create histogram of image pixel intensities in the masked area
            int[] histogram = new int[max - min + 1];
            inputCursor = input.cursor();
            pos = new long[input.numDimensions()];
            while ( inputCursor.hasNext() ){
            	inputCursor.fwd();
            	inputCursor.localize(pos);
            	maskPointer.setPosition(pos);
            	if (maskPointer.get().getInteger() > 0) {
	                type = inputCursor.get();
	                histogram[type.getInteger() - min] += 1;
            	}
            }
            
            
            //Calculates number of histogram bins
            int[] binThreshes = new int[binno];
            int binThresh = totalPixels / binno;
            int runningTotal = 0;
            int runningBinTotal = 1;
            for (int i = 0; i < histogram.length; i++) {
            	runningTotal += histogram[i];
            	if (runningTotal > (runningBinTotal * binThresh)) {
            		binThreshes[runningBinTotal] = i;
            		runningBinTotal += 1;
            		if (runningBinTotal == binno) {break;}
            	}
            }
            
            //Assigns pixels in copied image to bins

            Img<UnsignedShortType> retval = copyImage(inputfull);
    		IterableRandomAccessibleInterval<UnsignedShortType> retvalInput = 
    				new IterableRandomAccessibleInterval<UnsignedShortType>(
    						Views.hyperSlice(retval, dataDim, dataSlice));
            Cursor<UnsignedShortType> retvalCursor = retvalInput.cursor();
            IJ.log(Arrays.toString(binThreshes));
            while (retvalCursor.hasNext()) {
            	retvalCursor.fwd();
            	int compInt = retvalCursor.get().getInteger();
            	for (int i = 0; i < binno; i++) {
            		if (i == binno - 1) {
            			retvalCursor.get().set(i);
            			break;
            		}
            		else if (compInt > binThreshes[i]) {continue;}
            		else {
            			retvalCursor.get().set(i);
            			break;
            		}
            	}
            }
            
            return retval;
        }
    public < T extends Type< T > > Img< T > copyImage( final Img< T > input ){
        // create a new Image with the same properties
        // note that the input provides the size for the new image as it implements
        // the Interval interface
        Img< T > output = input.factory().create( input );
 
        // create a cursor for both images
        Cursor< T > cursorInput = input.cursor();
        Cursor< T > cursorOutput = output.cursor();
 
        // iterate over the input
        while ( cursorInput.hasNext()){
            // move both cursors forward by one pixel
            cursorInput.fwd();
            cursorOutput.fwd();
 
            // set the value of this pixel of the output image to the same as the input,
            // every Type supports T.set( T type )
            cursorOutput.get().set( cursorInput.get() );
        }
 
        // return the copy
        return output;
    }
    
    public static double[][][] sampleAngles(int binno){
    	double[][][] retval = new double[binno][][];
    	double phiInc = PI / (double) binno;
    	for (int i = 0; i < binno; i++) {
    		double phi = i * phiInc;
    		int thetaSamples = (int) (sin(phi) * (double) binno);
    		if (thetaSamples < 1) {thetaSamples = 1;}
    		retval[i] = new double[thetaSamples][2];
    		double thetaInc = PI / (double) thetaSamples;
    		for (int j = 0; j < thetaSamples; j++) {
    			retval[i][j] = new double[] {phi, (double) j * thetaInc};
    		}
    		
    	}
    	return retval;
    }
    public static double[][][][] sampleAnglesAndDistances(double[][][] angleSpace, int distanceSampleDepth, double distanceInc){
 	   double[][][][] retval = new double[angleSpace.length][][][];
 	   double[] distances = new double[distanceSampleDepth];
 	   for (int i = 0; i < distanceSampleDepth; i++) {
 		   distances[i] = (double) i * distanceInc;
 	   }
 	   for (int phiI = 0; phiI < angleSpace.length; phiI ++) {
 		   double[][] thetasAtPhi = angleSpace[phiI];
 		   retval[phiI] = new double[thetasAtPhi.length][][];
 		   double[][][] phiSubArr = retval[phiI];
 		   for (int thetaI = 0; thetaI < thetasAtPhi.length; thetaI ++) {
 			   double[] phiTheta = thetasAtPhi[thetaI];
 			   double phi = phiTheta[0];
 			   double theta = phiTheta[1];
 			   phiSubArr[thetaI] = new double[distanceSampleDepth][];
 			   double[][] phiThetaSubArr = phiSubArr[thetaI];
 			   for (int distanceI = 0; distanceI < distanceSampleDepth; distanceI ++) {
 				   double distance = distances[distanceI];
 				   double x = distance * cos(theta) * sin(phi);
 				   double y = distance * sin(theta) * sin(phi);
 				   double z = distance * cos(phi);
 				   phiThetaSubArr[distanceI] = new double[] {phi, theta, distance, x, y, z};
 			   }
 		   }
 	   }
 	   return retval;
    }
    public int[][][][] pixelOffsets(double[][][][] distanceArr){
    	int[][][][] retval = new int[distanceArr.length][][][];
    	for (int phiI = 0; phiI < distanceArr.length; phiI++ ) {
    		double[][][] phiSubArr = distanceArr[phiI];
    		retval[phiI] = new int[phiSubArr.length][][];
    		int[][][] phiSubIntArr = retval[phiI];
    		for (int thetaI = 0; thetaI < phiSubArr.length; thetaI ++){
    			double[][] phiThetaSubArr = phiSubArr[thetaI];
    			phiSubIntArr[thetaI] = new int[phiThetaSubArr.length][];
    			int[][] phiThetaSubIntArr = phiSubIntArr[thetaI];
    			for (int distanceI = 0; distanceI < phiThetaSubArr.length; distanceI++) {
    				double[] doubleTriplet = phiThetaSubArr[distanceI];
    				phiThetaSubIntArr[distanceI] = new int[3];
    				double x = doubleTriplet[3];
    				double y = doubleTriplet[4];
    				double z = doubleTriplet[5];
    				phiThetaSubIntArr[distanceI] = new int[] {
							(int) Math.round(imgCalibration.getRawX(x)),
							(int) Math.round(imgCalibration.getRawY(y)),
							(int) Math.round(imgCalibration.getRawZ(z))
					};
    			}
    		}
    	}
    	return retval;
    }
    public int[][][][][] greyLevelCorrelationMatrices(
    		int[][][][] pixelOffsets, Img<UnsignedShortType> binImg){
    	IJ.log("gothere");
    	
    	int[][][][][] retval = new int[pixelOffsets.length][][][][];
		IterableRandomAccessibleInterval<UnsignedShortType> input = 
				new IterableRandomAccessibleInterval<UnsignedShortType>(
						Views.hyperSlice(binImg, dataDim, dataSlice));
		
		RandomAccess<UnsignedShortType> maskPointer = 
				Views.extendValue(
						Views.hyperSlice(binImg, dataDim, maskSlice),
						new UnsignedShortType( 0 )).randomAccess();
		
		RandomAccess<UnsignedShortType> neighborPointer = 
				Views.hyperSlice(binImg, dataDim, dataSlice).randomAccess();
		
        Cursor<UnsignedShortType> inputCursor = input.cursor();
        
    	//Initializes the array of glcms
    	for (int phiI = 0; phiI < pixelOffsets.length; phiI++ ) {
    		int[][][] phiSubArr = pixelOffsets[phiI];
    		retval[phiI] = new int[phiSubArr.length][][][];
    		int[][][][] phiSubIntArr = retval[phiI];
    		for (int thetaI = 0; thetaI < phiSubArr.length; thetaI ++){
    			int[][] phiThetaSubArr = phiSubArr[thetaI];
    			phiSubIntArr[thetaI] = new int[phiThetaSubArr.length][][];
    			int[][][] phiThetaSubIntArr = phiSubIntArr[thetaI];
    			for (int distanceI = 0; distanceI < phiThetaSubArr.length; distanceI++) {
    				phiThetaSubIntArr[distanceI] = new int[binno][binno];
    			}
    		}
    	}
        
        long iterCount = 0;
        long[] pos = new long[input.numDimensions()];
        while ( inputCursor.hasNext() ){
        	iterCount ++;
        	if (iterCount % 1000 == 0) {
        		IJ.log(Long.toString(iterCount));
        	}
        	inputCursor.fwd();
        	inputCursor.localize(pos);
        	maskPointer.setPosition(pos);
        	//Fills the array of grey level correlation matrices
        	if (maskPointer.get().getInteger() > 0) {
	        	for (int phiI = 0; phiI < pixelOffsets.length; phiI++ ) {
	        		int[][][] phiSubArr = pixelOffsets[phiI];
	        		int[][][][] phiSubIntArr = retval[phiI];
	        		for (int thetaI = 0; thetaI < phiSubArr.length; thetaI ++){
	        			int[][] phiThetaSubArr = phiSubArr[thetaI];
	        			int[][][] phiThetaSubIntArr = phiSubIntArr[thetaI];
	        			for (int distanceI = 0; distanceI < phiThetaSubArr.length; distanceI++) {
	        				int[] offset = phiThetaSubArr[distanceI];
	        				long[] absolute = new long[3];
	        				int[][] glcm = phiThetaSubIntArr[distanceI];
	        				absolute[0] = pos[0] + (long) offset[0];
	        				absolute[1] = pos[1] + (long) offset[1];
	        				absolute[2] = pos[2] + (long) offset[2];
	                		for (int i = 0; i < pos.length; i++) {
	                			absolute[i] = pos[i] + (long) offset[i];
	                		}
	                		maskPointer.setPosition(absolute);
	                		if (maskPointer.get().getInteger() > 0) {
	                			int comparor = inputCursor.get().getInteger();
	                			neighborPointer.setPosition(absolute);
	                			int comparee = neighborPointer.get().getInteger();
	                			glcm[comparor][comparee] += 1;
	                		}
	        			}
	        		}
	        	}
	        }
        }
    	return retval;
    }
    public HashMap<double[], int[][]> linkGlcmToParameters(double[][][][] distanceArr, int[][][][][] glcmArr){
    	HashMap<double[], int[][]> retval = new HashMap<double[], int[][]>();
    	for (int phiI = 0; phiI < distanceArr.length; phiI++ ) {
    		double[][][] phiSubArr = distanceArr[phiI];
    		int[][][][] phiGlcmSubArr = glcmArr[phiI];
    		for (int thetaI = 0; thetaI < phiSubArr.length; thetaI ++){
    			double[][] phiThetaSubArr = phiSubArr[thetaI];
    			int[][][] phiThetaGlcmSubArr = phiGlcmSubArr[thetaI];
    			for (int distanceI = 0; distanceI < phiThetaSubArr.length; distanceI++) {
    				int[][] glcm = phiThetaGlcmSubArr[distanceI];
    				double[] doubleTriplet = phiThetaSubArr[distanceI];
    				double phi = doubleTriplet[0];
    				double theta = doubleTriplet[1];
    				double distance = doubleTriplet[2];
    				retval.put( new double[] {phi, theta, distance}, glcm);
    			}
    		}
    	}
    	return retval;
    }
    
    public HashMap<double[], Double> linkCorrelations(HashMap<double[], int[][]> glcmDict){
    	HashMap<double[], Double> retval = new HashMap<double[], Double>();
		for(double[] ptd: glcmDict.keySet()) {
			int[][] glcm = glcmDict.get(ptd);
			double mux = 0, muy = 0;
			int cardinality = 0;
			for (int x = 0; x < glcm.length; x++) {
				int[] glcmSub = glcm[x];
				for (int y = 0; y < glcmSub.length; y++) {
					int entry = glcm[x][y];
					cardinality += entry;
					mux += x * entry;
					muy += y * entry;
				}
			}
			mux = mux / (double) cardinality;
			muy = muy / (double) cardinality;
			double sdx = 0, sdy = 0, cov = 0;
			for (int x = 0; x < glcm.length; x++) {
				int[] glcmSub = glcm[x];
				for (int y = 0; y < glcmSub.length; y++) {
					int entry = glcm[x][y];
					double xdif = (double) x - mux;
					double ydif = (double) y - muy;
					sdx += (double) entry * Math.pow(xdif, 2);
					sdy += (double) entry * Math.pow(ydif, 2);
					cov += (double) entry * (xdif * ydif);
				}
			}
			sdx = Math.sqrt(sdx);
			sdy = Math.sqrt(sdy);
			retval.put(ptd, cov / (sdx * sdy));
		}		
    	return retval;
    }

	
}
