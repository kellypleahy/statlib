package statlib;

import java.util.*;

public class ExponentialDist implements Distribution {
  
  private final static Double[][] paramRanges = {{new Double(0.0), null}};
  private final static String distrName = "Exponential";
  private final static int paramCount = 1;
  private final static String[] paramNames = {"Mean"};
  
  private double[] paramValues = new double[1];
  private Vector moments = new Vector();
  
  public ExponentialDist(double mean) {
    paramValues[0] = mean;
  }
  
  /**
   * Sample the distribution at evenly spaced intervals, using the min value
   * and max value specified to make the correct number of buckets.
   */
  public Distribution sampleBuckets(double minValue, double maxValue, 
    int nBuckets) 
  {
    return StdImpl.sampleBuckets(this, minValue, maxValue, nBuckets);
  }
  
  /**
   * Sample the distribution at evenly spaced intervals.
   */
  public Distribution sampleStepped(double minValue, double stepSize, 
    int nSteps) 
  {
    return StdImpl.sampleStepped(this, minValue, stepSize, nSteps);
  }
  
  /**
   * Sample the distribution at the points specified.
   */
  public Distribution sample(double[] points) {
    return StdImpl.sample(this, points);
  }
  
  /**
   * Get the allowed parameter value ranges for the ith parameter.
   */
  public Double[] getParameterRange(int i) {
    if(i >= paramCount || i < 0)
      return null;
    else
      return paramRanges[i];
  }
  
  /**
   * Get the standard name of the distribution.
   */
  public String getDistributionName() {
    return distrName;
  }
  
  /**
   * Truncate the distribution using the min and max values specified (pass null
   *  for unbounded).
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
  
  /**
   * Get the standard deviation of the distribution.
   */
  public Double getStdDev() {
    return new Double(paramValues[0]);
  }
  
  /**
   * Get an instance string for the distribution (this would be for example
   * "Normal(0, 1)" for the standard normal distribution).
   */
  public String getDistributionInstance() {
    return distrName + "(" + paramValues[0] + ")";
  }
  
  /**
   * Get the probability that a variable with this distribution lies below or
   * equal to the specified value.
   */
  public Double getCumulativeProb(double v) {
    if(v <= 0)
      return new Double(0);
    else
      return new Double(1 - Math.exp(-v / paramValues[0]));
  }
  
  /**
   * Get the value probabilility function (p.d.f.) evaluated at the given point.
   */
  public Double getProbability(double v) {
    if(v <= 0)
      return new Double(0);
    else
      return new Double(Math.exp(-v / paramValues[0]) / v);
  }

  /**
   * Get the specified central moment from the distribution.
   */
  public Double getCentralMoment(int i) {
    if(i == 0) // E[(x - mx)^0] = E[1] = 1
      return new Double(1.0);
    else if(i == 1) // E[(x - mx)^1] = E[x-mx] = E[x] - E[mx] = E[x] - E[x] = 0
      return new Double(0.0);
    else {
      // calculate the central moment -- this is hard!
      // first, force calculation of raw moments
      getRawMoment(i);
      return StdImpl.centralMoment(moments, i);
    }
  }
  
  /**
   * Get the name of the ith parameter required by the distribution.
   */
  public String getParameterName(int i) {
    if(i >= paramCount || i < 0)
      return null;
    else
      return paramNames[i];
  }
  
  /**
   * Get the value of the ith parameter of this instance of the distribution.
   */
  public double getParameterValue(int i) {
    if(i >= paramCount || i < 0)
      return 0.0;
    else
      return paramValues[i];
  }
  
  /**
   * Get the number of parameters required by the distribution.
   */
  public int getParameterCount() {
    return paramCount;
  }
  
  /**
   * Get the parameter values of this instance of the distribution.
   */
  public double[] getParameterValues() {
    return paramValues;
  }
  
  /** Generate random sample from distribution (using Uniform random number
   * generator supplied, or java built-in (Math.random()) if null supplied).
   * @param n The number of values to be simulated (optimized for > 1)
   * @param rand The random number generator to use (if non-null) or use
   *             Math.random() if not supplied (null).
   * @return An array of length n with the simulated values */
  public double[] simulateValues(int n, IUniformRandom rand) {
    if(n <= 0)
      return null;
    else {
      IUniformRandom randL = (rand == null) ? StdImpl.rand : rand;
      double[] v = new double[n];
      if(n < 30) {
        // generate using "slow" way
        for(int i=0; i<n; i++)
          v[i] = - Math.log(randL.getNext()) * paramValues[0];
      } else {
        // generate using the marginally faster way
        double[] U = new double[n+1];
        double S = 1.0;
        for(int i=1; i<n; i++)
          S *= (U[i] = randL.getNext());
        double t = -Math.log(S);
        U[0] = 0.0; U[n] = t;
        Arrays.sort(U, 1, n-1);
        for(int i=0; i<n; i++)
          v[i] = t * (U[i+1] - U[i]);
      }
      return v;
    }
  }
  
  /**
   * Get the mean of the distribution.
   */
  public Double getMean() {
    return new Double(paramValues[0]);
  }
  
  /**
   * Get the specified raw moment from the distribution.
   */
  public Double getRawMoment(int i) {
    int sz = moments.size();
    if(i == 0)
      return new Double(1.0);
    else if(i == 1)
      return new Double(paramValues[0]);
    else if(i <= sz)
      return (Double)moments.get(i-1);
    else {
      double m = paramValues[0];
      if(sz == 0) {
        moments.add(new Double(m));
        sz++;
      }
      double a = ((Double)moments.get(sz-1)).doubleValue();
      for(int j = sz + 1; j <= i; j++) {
        a *= m * j;
        moments.add(new Double(a));
      }
      return getRawMoment(i);
    }
  }
  
  /*
  private double limRawMoment(int m, double lim) {
    if(m == 0) return 1.0;
    if(lim == Double.POSITIVE_INFINITY)
      return getRawMoment(m).doubleValue();
    else if(lim <= 0)
      return 0.0;
    double lm = Math.pow(lim, m);
    double mean = paramValues[0];
    double s = 0;
    double t = lm * lim / m;
    for(int r=0,j=m; r<m; r++,j--)
      s += (t *= j * mean / lim);
    s += t * mean / lim;
    s *= Math.exp(-lim / mean);
    s += lm * (1 - getCumulativeProb(lim).doubleValue());
    return s;
    //return limRawMoments(m, lim)[m];
  }
  
  private double[] limRawMoments(int m, double lim) {
    if(m < 0) return null;
    double rv[] = new double[m+1];
    rv[0] = getCumulativeProb(lim).doubleValue();
    if(m > 0) {
      if(lim == Double.POSITIVE_INFINITY) {
        getRawMoment(m);
        for(int i=1; i<=m; i++) 
          rv[i] = ((Double)moments.get(i-1)).doubleValue();
      } else if(lim <= 0) {
        for(int i=1; i<=m; i++) rv[i] = 0;
      } else {
        double mean = paramValues[0];
        double ex = Math.exp(-lim / mean);
        rv[0] = 1 - ex;
        for(int i=1; i<=m; i++) rv[i] = i * mean * rv[i-1] - (ex *= lim);
        double am = 1.0, q = 1 - rv[0];
        for(int i=1; i<=m; i++) rv[i] += (am *= lim) * q;
      }
    }
    return rv;
  }
  */

  private double limMoment(int m, double lim) {
    if(m==0) return 1.0;
    double apk = 1 + m;
    double ldb = lim / paramValues[0];
    double bk = Math.pow(paramValues[0], m);
    double lk = Math.pow(lim, m);
    return bk * SpecialFunc.gamma(apk) * SpecialFunc.incGamma(apk, ldb) 
      + lk * Math.exp(-lim / paramValues[0]);
  }
  
  /** Get the mth raw moment of the distribution limited to the range
   *  between a and b.
   */
  public Double getLimitedRawMoment(int m, Double a, Double b) {
    double lower, upper;
    lower = (a == null) ? 0 
      : limMoment(m, a.doubleValue()) - Math.pow(a.doubleValue(), m);
    upper = (b == null) ? getRawMoment(m).doubleValue()
      : limMoment(m, b.doubleValue());
    return new Double(upper - lower);
  }
  
  /** Get the first m raw moments of the distribution limited to the range
   *  between a and b.
   */
  public Double[] getLimitedRawMoments(int m, Double a, Double b) {
    Double[] rv = new Double[m+1];
    for(int i=0; i<=m; i++)
      rv[i] = getLimitedRawMoment(i, a, b);
    return rv;
  }
  
  /** Get the mth central moment of the distribution limited to the range
   *  between a and b
   */
  public Double getLimitedCentralMoment(int m, Double a, Double b) {
    Vector v = new Vector(m+1);
    v.addAll(Arrays.asList(getLimitedRawMoments(m, a, b)));
    v.remove(0);
    return StdImpl.centralMoment(v, m);
  }
  
  /** Get the first m central moments of the distribution limited to the range
   *  between a and b
   */
  public Double[] getLimitedCentralMoments(int m, Double a, Double b) {
    Vector v = new Vector(m+1);
    v.addAll(Arrays.asList(getLimitedRawMoments(m, a, b)));
    v.remove(0);
    Double[] rv = new Double[m+1];
    for(int i=0; i<=m; i++)
      rv[i] = StdImpl.centralMoment(v, i);
    return rv;
  }
  
  /** Convolve the distribution N times with itself (if implemented)
   */
  public Distribution convolve(int N) {
    return new GammaDist(N, paramValues[0]);
  }
  
  /** Calculate the quantile of the distribution at probability p
   */
  public Double getQuantile(double p) {
    return new Double(-paramValues[0] * Math.log(1-p));
  }
  
  public String toString() {
    return getDistributionInstance();
  }
}

  
