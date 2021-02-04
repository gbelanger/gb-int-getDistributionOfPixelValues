import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.text.DecimalFormat;
import java.util.UUID;

import hep.aida.IAnalysisFactory;
import hep.aida.IFitFactory;
import hep.aida.IFitResult;
import hep.aida.IFitter;
import hep.aida.IAxis;
import hep.aida.ITree;
import hep.aida.IFunction;
import hep.aida.IPlotter;
import hep.aida.IPlotterFactory;
import hep.aida.IHistogramFactory;
import hep.aida.IHistogram1D;
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
    for ( int i=0; i < pixelValues.length; i++ ) {
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
    double[] function = new double[nbins];
    double sigmaSqrd = Math.pow(sigma, 2);
    for ( int i=0; i < nbins; i++ ) {
      double x = axis2.binCenter(i);
      function[i] = fittedFunction.value(new double[] {x});
    }
    // Write
    writeHisto("histoOfPixelValues.qdp", pixelValues, min, max, nbins, function);
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

  // writeHisto
  private static void writeHisto(final String filename, final float[] values, final double xmin, final double xmax, final int nbins, final double[] function) throws IOException {
    IHistogram1D histo = makeHisto(values, xmin, xmax, nbins);
    String xLabel = "Pixel Value";
    String yLabel = "Entries per bin";
    boolean showStats = true;
    String[] header = makeHistoHeader((Histogram1D)histo, xLabel, yLabel, showStats);
    double[][] data = getData(histo);
    printToFile(filename, header, data[0], data[1], function);
  }

  private static IHistogram1D makeHisto(final float[] data, final double xmin, final double xmax, final int nBins) {
    IAnalysisFactory af = IAnalysisFactory.create();
    ITree tree = af.createTreeFactory().create();
    IHistogramFactory hf = af.createHistogramFactory(tree);
    IHistogram1D histo = hf.createHistogram1D("histo", nBins, xmin, xmax);
    double weight = 0;
    for ( int i=0; i < data.length; i++ ) {
      if ( Double.isNaN(new Double(data[i])) || data[i] == 0.0 ) {
        weight = 0;
      }
      else {
        weight = 1;
      }
      histo.fill(data[i], weight);
    }
    return histo;
  }

  private static String[] makeHistoHeader(final Histogram1D histo, final String xLabel, final String yLabel, final boolean showStats) {
    double[] yMinMax = calculateYMinYMax(histo);
    double yMin = yMinMax[0];
    double yMax = yMinMax[1];
    return makeHistoHeader(histo, yMin, yMax, xLabel, yLabel, showStats);
  }

  private static String[] makeHistoHeader(final Histogram1D histo, final double yMin, final double yMax, final String xLabel, final String yLabel, final boolean showStats) {
    IAxis axis = histo.axis();
    int nBins = axis.bins();
    double binWidth = axis.binWidth(0);
    double lowerEdge = axis.binLowerEdge(0);
    double upperEdge = axis.binUpperEdge(nBins-1);
    double xRange = upperEdge - lowerEdge;
    double nMajorDivs = xRange/(2*binWidth);
    String xMinStr = decimal.format(lowerEdge);
    String xMaxStr = decimal.format(upperEdge);
    String yMinStr = decimal.format(yMin);
    String yMaxStr = decimal.format(yMax);
    String[] header = null;
    if ( showStats == true ) {
      int entries = histo.entries();
      double mean = histo.mean();
      double rms = histo.rms();
      String meanStr = decimal.format(mean);
      String rmsStr = decimal.format(rms);
      header = new String[] {
        "DEV /XS",
        "READ 1",
        "LAB T", "LAB F",
        "TIME OFF",
        "LINE STEP",
        "LINE ON 3",
        "LW 4", "CS 1.5",
        "LAB 1 VPOS 0.76 0.8 \"Entries = "+entries+"\" JUST RIGHT CS 1",
        "LAB 2 VPOS 0.76 0.77 \"Mean = "+meanStr+"\" JUST RIGHT CS 1",
        "LAB 3 VPOS 0.76 0.74 \"RMS = "+rmsStr+"\" JUST RIGHT CS 1",
        "LAB X "+xLabel,
        "LAB Y "+yLabel,
        "VIEW 0.2 0.1 0.8 0.9",
        "R X "+xMinStr+" "+xMaxStr,
        "R Y "+yMinStr+" "+yMaxStr,
        //"GRID X "+nMajorDivs+",2",
        "!"
      };
    }
    else {
      header = new String[] {
        "DEV /XS",
        "READ 1",
        "LAB T", "LAB F",
        "TIME OFF",
        "LINE STEP",
        "LINE ON 3",
        "LW 4", "CS 1.5",
        "LAB X "+xLabel,
        "LAB Y "+yLabel,
        "VIEW 0.2 0.1 0.8 0.9",
        "R X "+xMinStr+" "+xMaxStr,
        "R Y "+yMinStr+" "+yMaxStr,
        //"GRID X "+nMajorDivs+",2",
        "!"
      };
    }
    return header;
  }

  public static double[][] getData(final Histogram1D histo) {
    IAxis axis = histo.axis();
    int nBins = axis.bins();
    double[] binHeights = new double[nBins];
    double[] binCentres = new double[nBins];
    for ( int i=0; i < nBins; i++ ) {
      binHeights[i] = histo.binHeight(i);
      binCentres[i] = axis.binCenter(i);
    }
    return new double[][] {binCentres, binHeights};
  }

  public static double[][] getData(final IHistogram1D iHisto) {
    return getData((Histogram1D) iHisto);
  }

  private static double[] calculateYMinYMax(final Histogram1D histo) {
    double maxBinHeight = histo.maxBinHeight();
    double minBinHeight = histo.minBinHeight();
    double margin = 0.05*maxBinHeight;
    double max = maxBinHeight + margin;
    double yMin = minBinHeight - margin;
    yMin = Math.max(0, yMin-margin);
    double yMax = max;
    return new double[] {yMin, yMax};
  }

  private static void printToFile(final String filename, final String[] header, final double[] binCentres, final double[] binHeights, final double[] function) throws IOException {
    int bufferSize = 256000;
    PrintWriter printWriter = new PrintWriter(new BufferedWriter(new FileWriter(filename), bufferSize));
    for ( int i=0; i < header.length; i++ ) {
      printWriter.println(header[i]);
    }
    for ( int i=0; i < binCentres.length; i++ ) {
      printWriter.println((binCentres[i]) +"\t"+ (binHeights[i]) +"\t"+ function[i]);
    }
    printWriter.close();
  }


}
