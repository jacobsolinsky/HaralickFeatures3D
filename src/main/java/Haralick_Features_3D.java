/*
 * To the extent possible under law, the ImageJ developers have waived
 * all copyright and related or neighboring rights to this tutorial code.
 *
 * See the CC0 1.0 Universal license for details:
 *     https://creativecommons.org/publicdomain/zero/1.0/
 */

import java.util.stream.Stream;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.plugin.*;
import ij.process.*;
import ij.gui.*;
import ij.*;
import java.util.Arrays;

import net.imagej.ImageJ;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.integer.UnsignedShortType;
import net.imglib2.util.Fraction;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.Cursor;

import org.scijava.Context;
import org.scijava.table.DefaultGenericTable;
import org.scijava.table.GenericTable;
public class Haralick_Features_3D implements PlugInFilter {
	ImagePlus img;
	
	public int setup(String arg, ImagePlus img){
		this.img = img;
		IJ.register(Haralick_Features_3D.class);
		return DOES_16+NO_CHANGES;
	}

	@Override
	public void run(ImageProcessor ip) {
		// ask user for number of rows & columns

		ImageStack stack = this.img.getStack();
		int width = this.img.getWidth();
		int height = this.img.getHeight();
		int depth = stack.size();
		short[] pixels = new short[width * height * depth];

		Img<UnsignedShortType> img3d = new ArrayImgFactory<UnsignedShortType>()
				.create(new long[] {width,  height,  depth}, new UnsignedShortType()
				); 
		Cursor<UnsignedShortType> cursor = img3d.cursor();
		cursor.fwd();
		while (cursor.hasNext()) {
			int[] posar = new int[] {0,0,0};
			
			cursor.localize(posar);
			UnsignedShortType thispixel = cursor.get();
			thispixel.set((int) stack.getVoxel(posar[0], posar[1], posar[2]));
			cursor.fwd();
		}
		
		final GenericDialog gd = new GenericDialog("Display a Table2");
		gd.addNumericField("Rows",  (int) pixels[179], 0);
		gd.addNumericField("Columns", 10, 0);
		gd.showDialog();
		if (gd.wasCanceled()) return;
		final int rowCount = (int) gd.getNextNumber();
		final int colCount = (int) gd.getNextNumber();
		ImageJFunctions.show(img3d);
		displayTable(10, colCount);
	}

	private void displayTable(final int rowCount, final int colCount) {
		// retrieve the ImageJ application context
		final Context context = (Context) IJ.runPlugIn("org.scijava.Context", "");
		final ImageJ ij = new ImageJ(context);

		// create a spreadsheet
		final GenericTable spreadsheet =
			new DefaultGenericTable(colCount, rowCount);

		// display the spreadsheet
		ij.ui().show("Spreadsheet", spreadsheet);
	}

	/** Tests the plugin. */
	public static void main(final String... args) {
		final ImageJ ij = new ImageJ();
		ij.launch(args);

		IJ.runPlugIn(Haralick_Features_3D.class.getName(), "");
	}
}
