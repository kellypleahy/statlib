// $Id: LogNormalDist.java,v 1.1 2004/11/09 06:23:15 kpl1 Exp $

package statlib;

import java.util.*;

public class LogNormalDist implements Distribution {
  
  private final static class WeightedLazyFit implements Distribution.WeightedLazyFit {
    private double sw, swz, swzsq;
    
    public void Initialize() { sw = swz = swzsq = 0; }
    
    public Distribution Finalize() {
      if(sw == 0) return null;
      else {
        double zbar, ssq;
        zbar = swz / sw;
        ssq = swzsq / sw - zbar * zbar;
        return new LogNormalDist(zbar, Math.sqrt(ssq));
      }
    }
    
    public void AddPoint(double weight, double point) {
      double z = Math.log(point);
      sw += weight;
      swz += weight * z;
      swzsq += weight * z * z;
    }
  }
  
  private final static class __WeightedLazyFitFactory 
    implements Distribution.WeightedLazyFitFactory 
  {
    public Distribution.WeightedLazyFit CreateInstance() {
      return new WeightedLazyFit();
    }
  }
  
  public final static Distribution.WeightedLazyFitFactory 
    WeightedLazyFitFactory; 
  
  private final static Double[][] paramRanges 
    = {{null, null}, {new Double(0), null}};
  private final static String distrName = "LogNormal";
  private final static int paramCount = 2;
  private final static String[] paramNames = {"Mu", "Sigma"};

  private double[] paramValues = new double[2];
  private Vector moments = new Vector();
  
  private static double sqrt2pi = Math.sqrt(2 * Math.PI);
  
  static {
    WeightedLazyFitFactory = new __WeightedLazyFitFactory();
  }

  public LogNormalDist(double mu, double sigma) {
    paramValues[0] = mu;
    paramValues[1] = sigma;
  }
  
  /**
    * Sample the distribution at evenly spaced intervals, using the min value
    * and max value specified to make the correct number of buckets.
    */
  public Distribution sampleBuckets(double minValue, double maxValue, 
				    int nBuckets) {
    return StdImpl.sampleBuckets(this, minValue, maxValue, 
				 nBuckets);
  }
  
  /**
    * Sample the distribution at evenly spaced intervals.
    */
  public Distribution sampleStepped(double minValue, double stepSize, 
				    int nSteps) {
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
      return new TruncatedDistribution(this, minValue, maxValue,
				       redistribute);
    }
  
  /**
    * Get the standard deviation of the distribution.
    */
  public Double getStdDev() {
    getRawMoment(2);
    double m1 = ((Double)moments.get(0)).doubleValue();
    double m2 = ((Double)moments.get(1)).doubleValue();
    return new Double(m2 - m1 * m1);
  }
  
  /**
    * Get an instance string for the distribution (this would be for example
    * "Normal(0, 1)" for the standard normal distribution).
    */
  public String getDistributionInstance() {
    return distrName + "(" + paramValues[0] + ", " + paramValues[1] + ")";
  }
  
  /**
    * Get the probability that a variable with this distribution lies below or
    * equal to the specified value.
    */
  public Double getCumulativeProb(double v) {
    double z = (Math.log(v) - paramValues[0]) / paramValues[1];
    if(z >= +4.0)
      return new Double(1.0);
    else if(z <= -4.0)
      return new Double(0.0);
    else {
      // calculate the cumulative normal probability -- this is hard!
      // see Handbook of Mathematical Functions 
      //  ed. Abramovitz and Stegun, Dover 1972 (ISBN 0-486-61272-4)
      // power series formula 26.2.10
      double x;
      if(z < 0) 
        x = -z;
      else
        x = z;
      
      int k = 1, n = 1;
      double com = x, term = com/k, d = term;
      final double r = sqrt2pi;
      double a = x * x / -2.0, tol = r * 0.000005;
      
      while(Math.abs(term) > tol) {
        com *= a / n++;
        k += 2;
        term = com / k;
        d += term;
      }
      
      d = d / r + 0.5;
      
      if(d > 1.0) d = 1.0;
      else if(d < 0.0) d = 0.0;
      
      if(z < 0)
        return new Double(1 - d);
      else
        return new Double(d);
    }
  }
  
  /**
   * Get the value probabilility function (p.d.f.) evaluated at the given point.
   */
  public Double getProbability(double v) {
    double d = (Math.log(v) - paramValues[0]) / paramValues[1];
    return new Double(Math.exp(-0.5 * d * d) 
      / (sqrt2pi * paramValues[1] * v));
  }

  /**
    * Get the specified central moment from the distribution.
    */
  public Double getCentralMoment(int i) {
    if(i == 0) // E[(x - mx)^0] = E[1] = 1
      return new Double(1.0);
    else if(i == 1) // E[(x - mx)^1] = E[x-mx] = E[x] - E[mx] = E[x] - E[x] = 0
      return new Double(0.0);
    else if(i == 2) // E[(x - mx)^2] = defn of VAR(x)
      return new Double(paramValues[1]);
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

  public double[] simulateValues(int n) {
    return simulateValues(n, null);
  }

  /**
    * Generate random sample from distribution (using Uniform random number
    * generator supplied, or java built-in (Math.random()) if null supplied).
    * @param n The number of values to be simulated (optimized for > 1)
    * @param rand The random number generator to use (if non-null) or use 
    *             Math.random() if not supplied (null).
    * @return An array of length n with the simulated values
    */
  public double[] simulateValues(int n, IUniformRandom rand) {
    if(n <= 0)
      return null;
    else {
      IUniformRandom randL = (rand == null) ? StdImpl.rand : rand;
      double mu = paramValues[0], sigma = paramValues[1];
      double[] d = new double[n];
      if(n % 2 == 1) {
        // generate a single normal variable using rejection method
        //  (see Ross - Simulation 3rd Ed. pg 71 - example 5f)
        double u, y;
        u = randL.getNext();
        y = randL.getNext();
        y = - Math.log(y);
        while (u > Math.exp(-(y - 1)*(y - 1) / 2)) {
          y = randL.getNext();
          y = - Math.log(y);
        }
        u = randL.getNext();
        if(u < 0.5) y = -y;
        d[0] = Math.exp(y * sigma + mu);
      }
      // generate pairs of normals with Box-Muller 
      //  (see Ross - Simulation 3rd Ed. pg 76)
      int i = 1 + (n%2);
      while(i < n) {
        double v1, v2, s;
        do {
          v1 = 2 * randL.getNext() - 1;
          v2 = 2 * randL.getNext() - 1;
          s = v1 * v1 + v2 * v2;
        } while(s > 1);
        double ss = Math.sqrt(-2 * Math.log(s) / s);
        d[i-1] = Math.exp(ss * v1 * sigma + mu);
        d[i] = Math.exp(ss * v2 * sigma + mu);
        i += 2;
      }
      return d;
    }    
  }

  /**
    * Get the mean of the distribution.
    */
  public Double getMean() {
    return getRawMoment(1);
  }

  /**
    * Get the specified raw moment from the distribution.
    */
  public Double getRawMoment(int i) {
    if(i == 0) 
      return new Double(1.0);
    else {
      int sz = moments.size();
      if(sz >= i)
        return (Double)moments.get(i-1);
      else {
        final double mu = paramValues[1], 
              h_sigmasq = 0.5 * paramValues[1] * paramValues[1];
        
        double k = sz, ksq = sz * sz;
        double m = 0;
        for(int r = sz; r<i; r++) {
          m = Math.exp(k * mu + ksq * h_sigmasq);
          moments.add(new Double(m));
          ksq += 2*k+1; k++;
        }
        return new Double(m);
      }    
    }
  }
  
  /**
    * Get the mth raw moment of the distribution limited to the range
    *  between a and b.
    */
  public Double getLimitedRawMoment(int m, Double a, Double b) {
    /*
    double[] v = binCoef(m);
    double mean = paramValues[0], stddev = Math.sqrt(paramValues[1]);
    double av = (a == null) ? Double.NEGATIVE_INFINITY 
      : (a.doubleValue() - mean) / stddev;
    double bv = (b == null) ? Double.POSITIVE_INFINITY
      : (b.doubleValue() - mean) / stddev;
    double[] nm1 = normLimMoment(m, av), nm2 = normLimMoment(m, bv);
    double facUp = 1.0, facDown = 1.0;
    for(int i=0; i<=m; i++, facUp *= stddev, facDown *= mean) {
      v[i] *= (nm2[i] - nm1[i]) * facUp;
      v[m-i] *= facDown;
    }
    double s = 0;
    for(int i=0; i<=m; i++)
      s += v[i];
    double fa = getCumulativeProb(av).doubleValue(),
      fb = getCumulativeProb(bv).doubleValue();
    if(fa > 0) s += Math.pow(av, m) * fa;
    if(fb < 1) s += Math.pow(bv, m) * (1 - fb);
    return new Double(s);
    */
    return null;
  }

  /**
    * Get the first m raw moments of the distribution limited to the range
    *  between a and b.
    */
  public Double[] getLimitedRawMoments(int m, Double a, Double b) {
    /*
    double mean = paramValues[0], stddev = Math.sqrt(paramValues[1]);
    double av = (a == null) ? Double.NEGATIVE_INFINITY 
      : (a.doubleValue() - mean) / stddev;
    double bv = (b == null) ? Double.POSITIVE_INFINITY
      : (b.doubleValue() - mean) / stddev;
    double fa = getCumulativeProb(av).doubleValue(),
      fb = getCumulativeProb(bv).doubleValue();
    double[] nm1 = normLimMoment(m, av), nm2 = normLimMoment(m, bv);
    Double[] r = new Double[m+1];
    for(int k=0; k<=m; k++) {
      double facUp = 1.0, facDown = 1.0;
      double[] v = binCoef(k);
      for(int i=0; i<=k; i++, facUp *= stddev, facDown *= mean) {
        v[i] *= (nm2[i] - nm1[i]) * facUp;
        v[k-i] *= facDown;
      }
      double s = 0;
      for(int i=0; i<=k; i++)
        s += v[i];
      if(fa > 0) s += Math.pow(av, k) * fa;
      if(fb < 1) s += Math.pow(bv, k) * (1 - fb);
      r[k] = new Double(s);
    }
    return r;
    */
    return null;
  }

  /**
    * Get the mth central moment of the distribution limited to the range
    *  between a and b
    */
  public Double[] getLimitedCentralMoments(int m, Double a, Double b) {
    /*
    Vector v = new Vector(m+1);
    v.addAll(Arrays.asList(getLimitedRawMoments(m, a, b)));
    v.remove(0);
    Double[] rv = new Double[m+1];
    for(int i=0; i<=m; i++)
      rv[i] = StdImpl.centralMoment(v, i);
    return rv;
    */
    return null;
  }

  /**
    * Get the mth central moment of the distribution limited to the range
    *  between a and b
    */
  public Double getLimitedCentralMoment(int m, Double a, Double b) {
    /*
    Vector v = new Vector(m+1);
    v.addAll(Arrays.asList(getLimitedRawMoments(m, a, b)));
    v.remove(0);
    return StdImpl.centralMoment(v, m);
    */
    return null;
  }

  /**
    * Convolve the distribution N times with itself (if implemented)
    */
  public Distribution convolve(int N) {
    /*
    return new LogNormalDist(paramValues[0] * N, paramValues[1] * N);
    */
    return null;
  }

  /**
    * Calculate the quantile of the distribution at probability p
    */
  public Double getQuantile(double p) {
    // need NormInv for this.  Will write later.
    return null;
  }

  public String toString() {
    return getDistributionInstance();
  }
}
