/*
 * GroupedDataDist.java
 *
 * Created on September 17, 2002, 10:25 PM
 */

package statlib;

import java.util.*;

/**
 *
 * @author  KLeahy
 */
public class GroupedDataDist implements Distribution {
  private double freq[];
  private double lo, hi, totalfreq;
  private int N;
  
  private Vector moments;
  
  /** Creates a new instance of GroupedDataDist */
  public GroupedDataDist(double lo, double hi, int n, double freq[]) {
    assert (n == freq.length);
    this.lo = lo;
    this.hi = hi;
    this.freq = (double[])freq.clone();
    moments = new Vector();
    N = n;
    totalfreq = 0.0;
    for(int i=0; i<N; i++) totalfreq += freq[i];
  }
  
  /** get the bucket to which this value belongs (-2 for below min, -1 for
   *    above max, -3 for error (buckets aren't mutually exhaustive))
   **/
  private int getFreqBucket(double v) {
    double l = lo;
    final double s = (hi - lo) / (N + 1);
    if(v <= lo) return -2;
    else if(v > hi) return -1;
    else for(int i=0; i<N; i++, l+=s) if(l < v && v <= l+s) return i;
    return -3;
  }
  
  /** Convolve the distribution N times with itself (if implemented)
   */
  public Distribution convolve(int N) {
    int l = freq.length * N;
    {
      int j = (int)(Math.log(l) / Math.log(2.0) + 1);
      l = 1 << j;
    }
    
    double newhi = (l - (freq.length * N)) * (hi - lo) / this.N + N * hi;
    
    double[] fqs = new double[l];
    for(int i=0; i<freq.length; i++) fqs[i] = freq[i];
    
    double[] imag = new double[fqs.length];
    double[] fqout = new double[fqs.length], imout = new double[fqs.length];
    try {
      Fourier.FFT(fqs, imag, fqout, imout);
      Fourier.impower(fqout, imout, N);
      Fourier.iFFT(fqout, imout, fqs, imag);
    } catch(Exception e) {
      System.err.println("FFT threw exception: " + e.getMessage());
      e.printStackTrace();
    }
    
    return new GroupedDataDist(lo * N, newhi, l, fqs);
  }
  
  /** Get the specified central moment from the distribution.
   */
  public Double getCentralMoment(int i) {
    if(i == 0) return new Double(1);
    else if(i == 1) return new Double(0);
    int s = moments.size() + 1;
    for(int j=s; j<=i; j++)
      moments.add(getRawMoment(j));
    return StdImpl.centralMoment(moments, i);
  }
  
  /** Get the probability that a variable with this distribution lies below or
   * equal to the specified value.
   */
  public Double getCumulativeProb(double v) {
    if(v <= lo)
      return new Double(0);
    else if(v > hi)
      return new Double(1);
    else {
      double f1 = 0.0, f2 = 0.0;
      double m = lo;
      final double s = (hi - lo) / (N + 1);
      for(int i=0; i<N; i++, m+=s) {
        if(m <= v)
          f2 += freq[i];
        if(m + s <= v)
          f1 += freq[i];
        if(m < v && v <= m + s)
          return new Double(((m + s - v) * f1 + (v - m) * f2) / totalfreq / s);
      }
    }
    return null;
  }
  
  /**
   * Get the value probabilility function (p.d.f.) evaluated at the given point.
   */
  public Double getProbability(double v) {
    // TODO: don't know how to do this yet!
    return null;
  }

  /** Get an instance string for the distribution (this would be for example
   * "Normal(0, 1)" for the standard normal distribution).
   */
  public String getDistributionInstance() {
    return getDistributionName() + " (" + lo + ", " + hi + ", " + N + ")";
  }
  
  /** Get the standard name of the distribution.
   */
  public String getDistributionName() {
    return "Grouped Data";
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
  
  /** Get the mth raw moment of the distribution limited to the range
   *  between a and b.
   */
  public Double getLimitedRawMoment(int m, Double a, Double b) {
    final double s = (hi - lo) / (N + 1);
    double la, lak, lb, lbk;
    double v = 0.0, l = lo, h = lo + s;
    
    // lower limit for calculations
    if(a == null) {
      la = Double.NEGATIVE_INFINITY;
      lak = 0;
    } else {
      la = a.doubleValue();
      lak = Math.pow(la, m);
    }
    // upper limit for calculations
    if(b == null) {
      lb = Double.POSITIVE_INFINITY;
      lbk = 0;
    } else {
      lb = b.doubleValue();
      lbk = Math.pow(lb, m);
    }

    for(int i=0; i<N; i++, l+=s) {
      double mv;
      if(l+s <= la) { // lower limit above interval
        mv = s * lak;
      } else if(l <= la) { // lower limit in interval
        if(l+s <= lb)
          // upper limit above interval
          mv = ((Math.pow(l+s, m+1) - la * lak) / (m+1) + lak * (la - l));
        else
          // upper limit in interval
          mv = (lak * (la - l) 
            + (lbk * lb - lak * la) / (m+1) + lbk * (l+s - lb));
      } else if(l < lb) { // interval starts between limits
        if(l+s <= lb)
          // upper limit above interval
          mv = (Math.pow(l+s, m+1) - Math.pow(l, m+1)) / (m+1);
        else
          // upper limit in interval
          mv = ((lbk * lb - Math.pow(l, m+1)) / (m+1) + lbk * (l+s - lb));
      } else if(l >= lb) { // interval above upper limit
        mv = s * lbk;
      } else {
        assert (false);
        mv = 0;
      }
      m += mv * freq[i];
    }
    return new Double(m / (s * totalfreq));
  }
  
  /** Get the first m raw moments of the distribution limited to the range
   *  between a and b.
   */
  public Double[] getLimitedRawMoments(int m, Double a, Double b) {
    // KPL: change this later to be more optimized for multiple moments!!!
    Double rv[] = new Double[m+1];
    for(int i=0; i<=m; i++)
      rv[i] = getLimitedRawMoment(i, a, b);
    return rv;
  }
  
  /** Get the mean of the distribution.
   */
  public Double getMean() {
    return getRawMoment(1);
  }
  
  /** Get the number of parameters required by the distribution.
   */
  public int getParameterCount() {
    return 0;
  }
  
  /** Get the name of the ith parameter required by the distribution.
   */
  public String getParameterName(int i) {
    return null;
  }
  
  /** Get the allowed parameter value ranges for the ith parameter.
   */
  public Double[] getParameterRange(int i) {
    return null;
  }
  
  /** Get the value of the ith parameter of this instance of the distribution.
   */
  public double getParameterValue(int i) {
    return 0.0;
  }
  
  /** Get the parameter values of this instance of the distribution.
   */
  public double[] getParameterValues() {
    return null;
  }
  
  /** Get the specified raw moment from the distribution.
   */
  public Double getRawMoment(int i) {
    final double s = (hi - lo) / (N + 1);
    double m = 0.0, l = lo, h = lo + s;
    for(int j=0; j<N; i++, l+=s)
      m += freq[j] * (Math.pow(l+s, i+1) - Math.pow(l, i+1)) / s;
    return new Double(m / (totalfreq * (i+1)));
  }
  
  /** Get the standard deviation of the distribution.
   */
  public Double getStdDev() {
    Double e2 = getRawMoment(2), e1 = getRawMoment(1);
    if(e2 != null && e1 != null) {
      double e = e1.doubleValue();
      return new Double(Math.sqrt(e2.doubleValue() - e * e));
    } else {
      return null;
    }
  }
  
  /** Sample the distribution at the points specified.
   */
  public Distribution sample(double[] points) {
    return StdImpl.sample(this, points);
  }
  
  /** Sample the distribution at evenly spaced intervals, using the min value
   * and max value specified to make the correct number of buckets.
   */
  public Distribution sampleBuckets(double minValue, double maxValue, 
    int nBuckets) 
  {
    return StdImpl.sampleBuckets(this, minValue, maxValue, nBuckets);
  }
  
  /** Sample the distribution at evenly spaced intervals.
   */
  public Distribution sampleStepped(double minValue, double stepSize, 
    int nSteps) 
  {
    return StdImpl.sampleStepped(this, minValue, stepSize, nSteps);
  }
  
  private double simulateValue(double rand, double step) {
    double s = 0.0;
    if(rand == 0) return lo;
    else if(rand == 1) return hi;
    for(int i=0; i<freq.length; i++, s+=freq[i]) {
      if(s + freq[i] >= rand && s < rand)
        return lo + step * (i + (rand - s) / freq[i]);
    }
    return hi;
  }
  
  /** Generate random sample from distribution (using Uniform random number
   * generator supplied, or java built-in (Math.random()) if null supplied).
   * @param n The number of values to be simulated (optimized for > 1)
   * @param rand The random number generator to use (if non-null) or use
   *             Math.random() if not supplied (null).
   * @return An array of length n with the simulated values  */
  public double[] simulateValues(int n, IUniformRandom rand) {
    double step = (hi - lo) / (N+1);
    IUniformRandom r = rand == null ? StdImpl.rand : rand;
    double[] rv = new double[n];
    for(int i=0; i<n; i++)
      rv[i] = simulateValue(rand.getNext(), step);
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
  
  /** Calculate the quantile of the distribution at probability p
   */
  public Double getQuantile(double p) {
    // Gotta figure out how to do this...  Should be the same as the code in
    //   sample (p is the uniform random variate).
    return null;
  }
  
  public String toString() {
    return getDistributionInstance();
  }
}
