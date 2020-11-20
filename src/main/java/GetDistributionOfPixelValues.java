
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.text.DecimalFormat;
import java.util.UUID;

import gb.esac.binner.Binner;
import gb.esac.io.AsciiDataFileWriter;
import gb.esac.tools.BasicStats;
import hep.aida.IAnalysisFactory;
import hep.aida.IFitFactory;
import hep.aida.IFitResult;
import hep.aida.IFitter;
import hep.aida.IFunction;
import hep.aida.IPlotter;
import hep.aida.IPlotterFactory;
import hep.aida.ref.histogram.FixedAxis;
import hep.aida.ref.histogram.Histogram1D;
import nom.tam.fits.Fits;
import nom.tam.fits.FitsException;
import nom.tam.fits.ImageHDU;
import nom.tam.util.ArrayFuncs;
import nom.tam.util.BufferedDataInputStream;
import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;


public class GetDistributionOfPixelValues {

    private static DecimalFormat decimal = new DecimalFormat("0.0000");   
    private static String filename;
    private static int extension;
    public static void main(String[] args) throws Exception {
	configureLogger();
	handleArgs(args);
	double stdDev = getStdDevOfPixelValuesDistribution(filename, extension);
	System.out.println(stdDev);
    }

    private static double getStdDevOfPixelValuesDistribution(String filename, int ext) throws Exception {
	// Get the data
	Fits fitsImage  = openFits(filename);
	ImageHDU hdu = (ImageHDU) fitsImage.getHDU(ext);
	float[][] pixelValues_2d = (float[][]) hdu.getKernel();
	float[] pixelValues = (float[]) ArrayFuncs.flatten(pixelValues_2d);
	// Make preliminary histo
	logger.info("Building preliminary histo:");
	logger.info("(excluding NaN and 0.0 values)");
	logger.info("  min = -3");
	logger.info("  max = 3");
	logger.info("  nbins = 1000");
	FixedAxis axis1 = new FixedAxis(1000, -3, 3);
	Histogram1D histo1 = new Histogram1D("", "Pixel Values", axis1);
	double weight = 0;
	for ( int i=0; i < pixelValues.length; i++ )  {
	    if ( Double.isNaN(pixelValues[i]) || pixelValues[i] == 0.0 ) {
		weight = 0;
	    }
	    else {
		weight = 1;
	    }
	    histo1.fill(pixelValues[i], weight);
	}
	// Get mean and rms and do the histo again
	double avg = histo1.mean();
	double rms = Math.abs(histo1.rms());
	int entries = histo1.entries();
	logger.info("  meam = "+avg);
	logger.info("  rms = "+rms);
	logger.info("  entries ="+entries);
	logger.info("Building final histo:");
	logger.info("(excluding NaN and 0.0 values)");
	double min = avg - 3.*rms;
	double max = avg + 3.*rms;
	int nbins = 1000;
	logger.info("  min = "+min+" (mean - 3*rms)");
	logger.info("  max = "+max+" (mean + 3*rms)");
	logger.info("  nbins = 1000");
	FixedAxis axis2 = new FixedAxis(nbins, min, max);
	Histogram1D histo2 = new Histogram1D("", "Pixel Values", axis2);
	weight = 0;
	for ( int i=0; i < pixelValues.length; i++ )  {
	    if ( Double.isNaN(pixelValues[i]) || pixelValues[i] == 0.0 ) {
		weight = 0;
	    }
	    else {
		weight = 1;
	    }
	    histo2.fill(pixelValues[i], weight);
	}
	avg = histo2.mean();
	rms = Math.abs(histo2.rms());
	entries = histo2.entries();
	logger.info("  meam = "+avg);
	logger.info("  rms = "+rms);
	logger.info("  entries ="+entries);
	logger.info("Fitting:");
	// Fit
	double amplitude = entries*axis2.binWidth(0)/Math.sqrt(2*Math.PI*rms*rms);
	double[] initial = new double[] {amplitude, avg, rms};
	IAnalysisFactory af = IAnalysisFactory.create();
	IFitFactory fitf = af.createFitFactory();
	IFitter fitter = fitf.createFitter("Chi2");
	IFitResult fitResult = fitter.fit(histo2, "g", initial);
	IFunction fittedFunction = fitResult.fittedFunction();
	double mean = fitResult.fittedParameter("mean");
	double sigma = fitResult.fittedParameter("sigma");
	String[] paramNames = fitResult.fittedParameterNames();
	double[] paramValues = fitResult.fittedParameters();
	logger.info("Fit status = "+fitResult.fitStatus());
	logger.info("Fit quality = "+fitResult.quality());
	logger.info("Fit results:");
	for ( int i=0; i < paramNames.length; i++ ) { logger.info("  "+paramNames[i]+" = "+Math.abs(paramValues[i])); }
	double[] y = new double[nbins];
	double sigmaSqrd = Math.pow(sigma, 2);
	for ( int i=0; i < nbins; i++ ) {
	    double x = axis2.binCenter(i);
	    y[i] = fittedFunction.value(new double[] {x}); 
	}
	// Write
	String histoName = "histoOfPixelValues.qdp";
	AsciiDataFileWriter out = new AsciiDataFileWriter(histoName);
	out.writeHisto(Binner.makeHisto(pixelValues, min, max, nbins), y, "Pixel Value");
	// Show
	// IPlotterFactory plotf = af.createPlotterFactory();
	// IPlotter plotter = plotf.create("Plot");
	// plotter.region(0).plot(histo);
	// plotter.region(0).plot(fitResult.fittedFunction());
	// plotter.show();
	return Math.abs(sigma);
    }

    //  Logger
    private static Logger logger  = Logger.getLogger(GetDistributionOfPixelValues.class);
    private static File loggerFile;
    private static void configureLogger() throws IOException {
	String loggerFilename= "logger.config";
	InputStream log = ClassLoader.getSystemResourceAsStream(loggerFilename);
	UUID uuid = UUID.randomUUID();
	String homeDir = System.getProperty("user.home");
	loggerFilename = new String(homeDir+File.pathSeparator+"logger.config_"+uuid.toString());
	loggerFile = new File(loggerFilename);
	loggerFile.deleteOnExit();
	inputStreamToFile(log, loggerFilename);
        PropertyConfigurator.configure(loggerFilename);
    }
    private static void inputStreamToFile(InputStream io, String fileName) throws IOException {       
	FileOutputStream fos = new FileOutputStream(fileName);
	byte[] buf = new byte[256];
	int read = 0;
	while ((read = io.read(buf)) > 0) {
	    fos.write(buf, 0, read);
	}
	fos.flush();
	fos.close();
    }

    // Arguments
    private static void handleArgs(String[] args) {
	if ( args.length != 2 ) {
	    logger.error("Usage: java -jar GetDistributionOfPixelValues.jar filename extension");
	    System.exit(-1);
	}
	filename = args[0];
	extension = (int) Integer.valueOf(args[1]);
    }

    // isGzipped
    public static boolean isGzipped(String fileName) throws IOException {
	return isGzipped(new File(fileName));
    }    
    public static boolean isGzipped(File file) throws IOException {
	InputStream in = new FileInputStream(file);
	int magic1 = in.read();
	int magic2 = in.read();
	in.close();
	return (magic1 == 0037 && magic2 == 0213);
    }

    // openFits
    public static Fits openFits(String filename) throws IOException, FitsException {
	return openFits(new File(filename));
    }
    public static Fits openFits(File file) throws IOException, FitsException {
	boolean isGzipped = isGzipped(file);
	BufferedDataInputStream dis = new BufferedDataInputStream(new FileInputStream(file));
	Fits fitsFile = new Fits(dis, isGzipped);
	return fitsFile;
    }

}
