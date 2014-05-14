package statlib;

import java.util.Arrays;
import java.util.Vector;

public interface Distribution {
  
  public interface WeightedLazyFit {
    public void Initialize();
    public Distribution Finalize();
    public void AddPoint(double weight, double point);
  }
  
  public interface WeightedLazyFitFactory {
    public WeightedLazyFit CreateInstance();
  }
  
  public static class TruncatedDistribution implements Distribution {
    private Double min, max;
    private boolean redistribute;
    private double minFreq, maxFreq;
    private Vector moments;
    private Distribution srcDist;
    
    public TruncatedDistribution(Distribution srcDist, Double min, Double max, 
      boolean redistribute)
    {
      this.srcDist = srcDist;
      this.min = min;
      this.max = max;
      this.redistribute = redistribute;
      minFreq = srcDist.getCumulativeProb(min.doubleValue()).doubleValue();
      maxFreq = 1 - srcDist.getCumulativeProb(max.doubleValue()).doubleValue();
      moments = new Vector();
    }
    
    public Double getCumulativeProb(double v) {
      double m, M;
      m = (min == null) ? Double.NEGATIVE_INFINITY : min.doubleValue();
      M = (max == null) ? Double.POSITIVE_INFINITY : max.doubleValue();
      if(v < m)
        return new Double(0.0);
      else if(v > M)
        return new Double(1.0);
      else if(redistribute) 
        return new Double((srcDist.getCumulativeProb(v).doubleValue() - minFreq)
          / (1 - maxFreq - minFreq));
      else
        return srcDist.getCumulativeProb(v);
    }
    
    /**
     * Get the value probabilility function (p.d.f.) evaluated at the given point.
     */
    public Double getProbability(double v) {
      double m, M;
      m = (min == null) ? Double.NEGATIVE_INFINITY : min.doubleValue();
      M = (max == null) ? Double.POSITIVE_INFINITY : max.doubleValue();
      if(v < m || v > M)
        return new Double(0);
      else if(!redistribute)
        return srcDist.getProbability(v);
      else {
        Double d = srcDist.getProbability(v);
        if(d != null)
          return new Double(d.doubleValue() / (1 - maxFreq - minFreq));
        else
          return null;
      }
    }

    /** Get the specified central moment from the distribution.
     */
    public Double getCentralMoment(int i) {
      return null;
    }
    
    /** Get an instance string for the distribution (this would be for example
     * "Normal(0, 1)" for the standard normal distribution).
     */
    public String getDistributionInstance() {
      return "Truncated(" + min + ", " + max + ") " 
        + srcDist.getDistributionInstance();
    }
    
    /** Get the standard name of the distribution.
     */
    public String getDistributionName() {
      return "Truncated " + srcDist.getDistributionName();
    }
    
    /** Get the mean of the distribution.
     */
    public Double getMean() {
      return getRawMoment(1);
    }
    
    /** Get the number of parameters required by the distribution.
     */
    public int getParameterCount() {
      return srcDist.getParameterCount();
    }
    
    /** Get the name of the ith parameter required by the distribution.
     */
    public String getParameterName(int i) {
      return srcDist.getParameterName(i);
    }
    
    /** Get the allowed parameter value ranges for the ith parameter.
     */
    public Double[] getParameterRange(int i) {
      return srcDist.getParameterRange(i);
    }
    
    /** Get the value of the ith parameter of this instance of the distribution.
     */
    public double getParameterValue(int i) {
      return srcDist.getParameterValue(i);
    }
    
    /** Get the parameter values of this instance of the distribution.
     */
    public double[] getParameterValues() {
      return srcDist.getParameterValues();
    }
    
    /** Get the specified raw moment from the distribution.
     */
    public Double getRawMoment(int i) {
      int sz = moments.size();
      if(sz >= i)
        return (Double)moments.get(i-1);
      else {
        for(int r=sz; r < i; r++) {
          double m = srcDist.getLimitedRawMoment(i, min, max).doubleValue();
          if(redistribute) m /= (1 - maxFreq - minFreq);
          moments.add(new Double(m));
        }
        return (Double)moments.get(i-1);
      }
    }
    
    /** Get the standard deviation of the distribution.
     */
    public Double getStdDev() {
      Double m = getRawMoment(1);
      Double m2 = getRawMoment(2);
      if(m == null || m2 == null) return null;
      double mv = m.doubleValue();
      return new Double(m2.doubleValue() - mv * mv);
    }
    
    /** Sample the distribution at the points specified.
     */
    public Distribution sample(double[] points) {
      Distribution rv = srcDist.sample(points);
      if(rv != null)
        return rv.truncate(min, max, redistribute);
      else
        return null;
    }
    
    /** Sample the distribution at evenly spaced intervals, using the min value
     * and max value specified to make the correct number of buckets.
     */
    public Distribution sampleBuckets(double minValue, double maxValue, 
      int nBuckets) 
    {
      Distribution rv = srcDist.sampleBuckets(minValue, maxValue, nBuckets);
      if(rv != null)
        return rv.truncate(min, max, redistribute);
      else
        return null;
    }
    
    /** Sample the distribution at evenly spaced intervals.
     */
    public Distribution sampleStepped(double minValue, double stepSize, 
      int nSteps) 
    {
      Distribution rv = srcDist.sampleStepped(minValue, stepSize, nSteps);
      if(rv != null)
        return rv.truncate(min, max, redistribute);
      else
        return null;
    }
    
    /** Generate random sample from distribution (using Uniform random number
     * generator supplied, or java built-in (Math.random()) if null supplied).
     * @param n The number of values to be simulated (optimized for > 1)
     * @param rand The random number generator to use (if non-null) or use
     *             Math.random() if not supplied (null).
     * @return An array of length n with the simulated values  */
    public double[] simulateValues(int n, IUniformRandom rand) {
      double m, M;
      m = (min == null) ? Double.NEGATIVE_INFINITY : min.doubleValue();
      M = (max == null) ? Double.POSITIVE_INFINITY : max.doubleValue();
      double[] rv = srcDist.simulateValues(n, rand);
      if(rv != null)
        if(redistribute)
          for(int i=0; i<rv.length; i++)
            while(rv[i] < m || rv[i] > M)
              rv[i] = srcDist.simulateValues(1, rand)[0];
        else
          for(int i=0; i<rv.length; i++)
            if(rv[i] < m)
              rv[i] = m;
            else if(rv[i] > M)
              rv[i] = M;
      return rv;
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
    public Distribution truncate(Double minValue, Double maxValue, 
      boolean redistribute) 
    {
      return new TruncatedDistribution(this, minValue, maxValue, redistribute);
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
    
    /** Get the mth central moment of the distribution limited to the range
     *  between a and b
     */
    public Double[] getLimitedCentralMoments(int m, Double a, Double b) {
      return null;
    }
    
    /** Get the mth central moment of the distribution limited to the range
     *  between a and b
     */
    public Double getLimitedCentralMoment(int m, Double a, Double b) {
      return null;
    }
    
    /** Convolve the distribution N times with itself (if implemented)
     */
    public Distribution convolve(int N) {
      return null;
    }
    
    /** Calculate the quantile of the distribution at probability p
     */
    public Double getQuantile(double p) {
      return null;
    }
    
  }

  
  public interface IUniformRandom {
    public void setSeed(double seed);
    public double getNext();
    public void reset();
  }
  
  /**
   * This is a class that implements some functionality usable for most
   *   distributions in a generic way, so as to limit code duplication.
   **/
  public static class StdImpl {
    /**
     * A standard implementation of the IUniformRandom interface that uses Java's
     *   Math.random() to generate the random numbers
     **/
    private static class JavaUniformRandom implements IUniformRandom {
      public void setSeed(double seed) {
        // this can't be implemented for java.Math
      }
      public double getNext() {
        return Math.random();
      }
      public void reset() {
        // this can't be implemented for java.Math
      }
    }
  
    /**
     * The Java standard Math.random() based random number generator for use
     *   in calls to simulateValues()
     **/
    public static IUniformRandom rand = new JavaUniformRandom();
        
    /**
     * A standard implementation of the sample() method for all distributions
     *   to use in lieu of a special method for doing sampling.
     **/
    public static Distribution sample(Distribution src, double[] points) {
      int nSteps = points.length;
      double[] prob = new double[nSteps];
      double cumprob = 0.0;
      double[] mypoints = (double[])points.clone();

      Arrays.sort(mypoints);

      for(int i=0; i<nSteps-1; i++) {
        prob[i] = src.getCumulativeProb(points[i]).doubleValue() - cumprob;
        cumprob += prob[i];
      }

      prob[nSteps-1] = 1 - cumprob;
      return new FrequencyDist(mypoints, prob);
    }
    
    /**
     * A standard implementation of the sampleStepped() method for all 
     *   distributions to use in lieu of a special method for doing sampling.
     **/
    public static Distribution sampleStepped(Distribution src, double minValue, 
      double stepSize, int nSteps)
    {
      double pts[] = new double[nSteps];
      double left = minValue;
      for(int i=0; i<nSteps; i++, left += stepSize)
        pts[i] = left;
      return sample(src, pts);
    }
    
    /**
     * A standard implementation of the sampleBuckets() method for all 
     *   distributions to use in lieu of a special method for doing sampling.
     **/
    public static Distribution sampleBuckets(Distribution src, double minValue, 
      double maxValue, int nBuckets)
    {
      double sz = (maxValue - minValue) / (nBuckets - 1);
      return sampleStepped(src, minValue, sz, nBuckets);
    }
    
    /**
     * A standard implementation of the algorithm to convert from raw moments
     *   to central moments given a list of raw moments is stored up through
     *   the moment requested
     **/
    public static Double centralMoment(Vector moments, int n) {
      // this is based on the formula relating raw moments to central moments
      //   see for example Abramowitz, Stegun: Handbook of Math. functions
      //    ISBN 0-486-61272-4 formula #26.1.14 (page 928).
      if(moments.size() < n) return null;
      double mean = ((Double)moments.get(0)).doubleValue();
      double mul = 1.0, s = 0.0;
      for(int k = 0; k < n; k++) {
        s += mul * ((Double)moments.get(n - k - 1)).doubleValue();
        mul *= -1.0 * mean * (n - k) / (k + 1.0);
      }
      s += mul;
      return new Double(s);
    }
    
    /**
     * Convolve the distribution specified by v with itself N times and return
     *   the result.  The points implied in v must be equally spaced.  The 
     *   parameter v is the frequencies of the points and should be normalized.
     **/
    public static double[] convolve(double[] v, int N) {
      return null;
    }
  }    
  
  /**
   * Get the number of parameters required by the distribution.
   */
  public int getParameterCount();

  /**
   * Get the name of the ith parameter required by the distribution.
   */
  public String getParameterName(int i);
  
  /**
   * Get the allowed parameter value ranges for the ith parameter.
   */
  public Double[] getParameterRange(int i);
  
  /**
   * Get the standard name of the distribution.
   */
  public String getDistributionName();
  
  /**
   * Get an instance string for the distribution (this would be for example 
   * "Normal(0, 1)" for the standard normal distribution).
   */
  public String getDistributionInstance();
  
  /**
   * Get the value of the ith parameter of this instance of the distribution.
   */
  public double getParameterValue(int i);
  
  /**
   * Get the parameter values of this instance of the distribution.
   */
  public double[] getParameterValues();
  
  /**
   * Truncate the distribution using the min and max values specified (pass null
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
  public Distribution truncate(Double minValue, Double maxValue, 
    boolean redistribute);
  
  /**
   * Sample the distribution at evenly spaced intervals.
   */
  public Distribution sampleStepped(double minValue, double stepSize, 
    int nSteps);
  
  /**
   * Sample the distribution at evenly spaced intervals, using the min value
   * and max value specified to make the correct number of buckets.
   */
  public Distribution sampleBuckets(double minValue, double maxValue,
    int nBuckets);
  
  /**
   * Sample the distribution at the points specified.
   */
  public Distribution sample(double[] points);
  
  /**
   * Get the probability that a variable with this distribution lies below or 
   * equal to the specified value.
   */
  public Double getCumulativeProb(double v);
  
  /**
   * Get the value probabilility function (p.d.f.) evaluated at the given point.
   */
  public Double getProbability(double v);
  
  /**
   * Get the mean of the distribution.
   */
  public Double getMean();
  
  /**
   * Get the standard deviation of the distribution.
   */
  public Double getStdDev();
  
  /**
   * Get the specified raw moment from the distribution.
   */
  public Double getRawMoment(int i);
  
  /**
   * Get the specified central moment from the distribution.
   */
  public Double getCentralMoment(int i);
  
  /** Generate random sample from distribution (using Uniform random number
   * generator supplied, or java built-in (Math.random()) if null supplied).
   * @param n The number of values to be simulated (optimized for > 1)
   * @param rand The random number generator to use (if non-null) or use 
   *             Math.random() if not supplied (null).
   * @return An array of length n with the simulated values */
  public double[] simulateValues(int n, IUniformRandom rand);
  
  /**
   * Get the first m raw moments of the distribution limited to the range
   *  between a and b.
   */
  public Double[] getLimitedRawMoments(int m, Double a, Double b);
  
  /**
   * Get the mth raw moment of the distribution limited to the range
   *  between a and b.
   */
  public Double getLimitedRawMoment(int m, Double a, Double b);
  
  /**
   * Get the first m central moments of the distribution limited to the range
   *  between a and b
   */
  public Double[] getLimitedCentralMoments(int m, Double a, Double b);

  /**
   * Get the mth central moment of the distribution limited to the range
   *  between a and b
   */
  public Double getLimitedCentralMoment(int m, Double a, Double b);
  
  /**
   * Convolve the distribution N times with itself (if implemented)
   */
  public Distribution convolve(int N);
  
  /**
   * Calculate the quantile of the distribution at probability p
   */
  public Double getQuantile(double p);

}
