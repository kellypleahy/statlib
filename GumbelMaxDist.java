/*
 * GumbelDist.java
 *
 * Created on October 15, 2002, 10:34 AM
 */

package statlib;

/**
 *
 * @author  KLeahy
 */
public class GumbelMaxDist implements Distribution {
  
  private final static Double[][] paramRanges 
    = {{null, null}, {new Double(0), null}};
  private final static String distrName = "Gumbel";
  private final static int paramCount = 2;
  private final static String[] paramNames = {"Mu", "Sigma"};

  private double[] paramValues = new double[2];
  
  private final static double euler = 0.577215664901532860606512;
  private final static double sqrt6dpi = Math.sqrt(6) / Math.PI;

  /** Creates a new instance of GumbelDist */
  public GumbelMaxDist(double mu, double sigma) {
    paramValues[0] = mu;
    paramValues[1] = sigma;
  }
  
  public GumbelMaxDist(FrequencyDist d) {
    paramValues = estimate(d);
  }
  
  public static double[] estimate(Distribution d) {
    double[] rv = new double[2];
    if(d instanceof FrequencyDist) {
      FrequencyDist d2 = (FrequencyDist)d;
      rv[1] = Math.sqrt(d2.getSampleVar()) * sqrt6dpi;
    } else {
      rv[1] = d.getStdDev().doubleValue() * sqrt6dpi;
    }
    rv[0] = d.getMean().doubleValue() - rv[1] * euler;
    return rv;
  }
  
  /** Convolve the distribution N times with itself (if implemented)
   */
  public Distribution convolve(int N) {
    return null;
  }
  
  /** Get the specified central moment from the distribution.
   */
  public Double getCentralMoment(int i) {
    return null;
  }
  
  /** Get the probability that a variable with this distribution lies below or
   * equal to the specified value.
   */
  public Double getCumulativeProb(double v) {
    double F = Math.exp(-Math.exp((paramValues[0] - v)/paramValues[1]));
    return new Double(F);
  }

  /**
   * Get the value probabilility function (p.d.f.) evaluated at the given point.
   */
  public Double getProbability(double v) {
    // TODO: don't remember how to do this!
    return null;
  }

  /** Get an instance string for the distribution (this would be for example
   * "Normal(0, 1)" for the standard normal distribution).
   */
  public String getDistributionInstance() {
    return getDistributionName() + "(" + paramValues[0] + ", " + paramValues[1]
      + ")";
  }
  
  /** Get the standard name of the distribution.
   */
  public String getDistributionName() {
    return distrName;
  }
  
  /** Get the mth central moment of the distribution limited to the range
   *  between a and b
   */
  public Double getLimitedCentralMoment(int m, Double a, Double b) {
    return null;
  }
  
  /** Get the first m central moments of the distribution limited to the range
   *  between a and b
   */
  public Double[] getLimitedCentralMoments(int m, Double a, Double b) {
    return null;
  }
  
  /** Get the mth raw moment of the distribution limited to the range
   *  between a and b.
   */
  public Double getLimitedRawMoment(int m, Double a, Double b) {
    return null;
  }
  
  /** Get the first m raw moments of the distribution limited to the range
   *  between a and b.
   */
  public Double[] getLimitedRawMoments(int m, Double a, Double b) {
    return null;
  }
  
  /** Get the mean of the distribution.
   */
  public Double getMean() {
    return new Double(paramValues[0] + paramValues[1] * euler);
  }
  
  /** Get the number of parameters required by the distribution.
   */
  public int getParameterCount() {
    return paramCount;
  }
  
  /** Get the name of the ith parameter required by the distribution.
   */
  public String getParameterName(int i) {
    if(i < 0 || i >= paramCount) return null;
    else return paramNames[i];
  }
  
  /** Get the allowed parameter value ranges for the ith parameter.
   */
  public Double[] getParameterRange(int i) {
    if(i < 0 || i >= paramCount) return null;
    else return paramRanges[i];
  }
  
  /** Get the value of the ith parameter of this instance of the distribution.
   */
  public double getParameterValue(int i) {
    if(i < 0 || i >= paramCount) return 0;
    else return paramValues[i];
  }
  
  /** Get the parameter values of this instance of the distribution.
   */
  public double[] getParameterValues() {
    return paramValues;
  }
  
  /** Get the specified raw moment from the distribution.
   */
  public Double getRawMoment(int i) {
    if(i < 0) return null;
    else if(i == 0) return new Double(1.0);
    else if(i == 1) return getMean();
    else return null;
  }
  
  /** Get the standard deviation of the distribution.
   */
  public Double getStdDev() {
    return new Double(paramValues[1] / sqrt6dpi);
  }
  
  /** Sample the distribution at the points specified.
   */
  public Distribution sample(double[] points) {
    return null;
  }
  
  /** Sample the distribution at evenly spaced intervals, using the min value
   * and max value specified to make the correct number of buckets.
   */
  public Distribution sampleBuckets(double minValue, double maxValue, int nBuckets) {
    return null;
  }
  
  /** Sample the distribution at evenly spaced intervals.
   */
  public Distribution sampleStepped(double minValue, double stepSize, int nSteps) {
    return null;
  }
  
  /** Generate random sample from distribution (using Uniform random number
   * generator supplied, or java built-in (Math.random()) if null supplied).
   * @param n The number of values to be simulated (optimized for > 1)
   * @param rand The random number generator to use (if non-null) or use
   *             Math.random() if not supplied (null).
   * @return An array of length n with the simulated values  */
  public double[] simulateValues(int n, IUniformRandom rand) {
    return null;
  }
  
  /** Truncate the distribution using the min and max values specified (pass null
   *   for unbounded).
   * @param minValue the minimum value to use for truncation (or null for none)
   * @param maxValue the maximum value to use for truncation (or null for none)
   * @param redistribute if true, redistribute the probability on the truncated
   *                     tails to the distribution evenly (as if these
   *                     observations aren't possible), otherwise allocate the
   *                     probability to the tails as point masses.
   * @return Null if the distribution does not have a truncated form
   *     (or it's not implemented).
   */
  public Distribution truncate(Double minValue, Double maxValue, boolean redistribute) {
    return null;
  }
  
  /** Calculate the quantile of the distribution at probability p
   */
  public Double getQuantile(double p) {
    if(p <= 0 || p >= 1) return null;
    else return new Double(-Math.log(-Math.log(p)));
  }
  
  public String toString() {
    return getDistributionInstance();
  }
}
