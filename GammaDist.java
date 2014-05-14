/*
 * GammaDist.java
 *
 * Created on September 24, 2002, 10:21 AM
 */

package statlib;

import java.util.*;

/**
 *
 * @author  KLeahy
 */
public class GammaDist implements Distribution {
  private final static Double[][] paramRanges 
    = {{new Double(0), null}, {new Double(0), null}};
  private final static String distrName = "Gamma";
  private final static int paramCount = 2;
  private final static String[] paramNames = {"alpha", "beta"};
  
  private double[] paramValues = new double[2];
  private Vector moments = new Vector();
  private double lnGamAlpha;
  private double betaAlpha;
  
  private double limMoment(int m, double lim) {
    if(m==0) return 1.0;
    double apk = paramValues[0] + m;
    double ldb = lim / paramValues[1];
    double bk = Math.pow(paramValues[1], m);
    double lk = Math.pow(lim, m);
    return bk * Math.exp(SpecialFunc.lnGamma(apk) - lnGamAlpha) 
      * SpecialFunc.incGamma(apk, ldb) 
      + lk * (1 - SpecialFunc.incGamma(paramValues[0], ldb));
  }
  
  public final static class MLEFit {
    public double alpha, beta;
    private double xbar, lxbar, sx;
    
    private static final double tol = 5e-8;
    
    public MLEFit(double xbar, double lxbar, double sx) {
      this(xbar, lxbar, sx, false);
    }
    
    public MLEFit(double xbar, double lxbar, double sx, boolean debug) {
      this.xbar = xbar;
      this.lxbar = lxbar;
      this.sx = sx;
      
      perform(debug);
    }
    
    private void perform(boolean debug) {
      double la, lb;
      
      if(debug) {
        System.out.println("perform called with xbar = " + xbar 
          + ", sx = " + sx);
      }
      
      // initial estimates - the method of moments
      la = xbar / (lb = sx * sx / xbar);
      
      String s;
      
      if(debug) 
        System.out.println("la[0] = " + la + ", lb[0] = " + lb);
      
      int iter = 0;
      
      for(;;) {
        iter++;
        if(debug)
          System.out.println("Entering iteration " + iter);
        
        // the derivatives of lnGamma(alpha)
        double psi   = SpecialFunc.Psi(la),
               dpsi  = SpecialFunc.PolyGamma(1, la);
               
        if(debug) {
          System.out.println("psi[" + la + "] = " + psi);
          System.out.println("dpsi[" + la + "] = " + dpsi);
        }                 
        
        // the "helper" factors
        double t1    = Math.log(lb) + psi - lxbar,
               t2    = la * lb - xbar,
               t3inv = 1.0 / (la * dpsi - 1.0);
               
        if(debug) {
          System.out.println("t1[" + iter + "] = " + t1);
          System.out.println("t2[" + iter + "] = " + t2);
          System.out.println("t3[" + iter + "] = " + (1.0 / t3inv));
        }
        
        // the "deltas"
        double da = (t2 / lb - t1 * la) * t3inv;
        double db = (t1 * lb - dpsi * t2) * t3inv;
        
        if(debug) {
          System.out.println("da[" + iter + "] = " + da);
          System.out.println("db[" + iter + "] = " + db);
        }
        
        // the update equations
        la += da;
        lb += db;
        
        if(debug) {
          System.out.println("la[" + iter + "] = " + la);
          System.out.println("lb[" + iter + "] = " + lb);
        }
        
        // the distance moved (in both parameters)
        double norm = Math.sqrt(da * da + db * db);
        
        if(debug)
          System.out.println("norm[" + iter + "] = " + norm);
        
        if(norm < tol || Double.isNaN(norm)) {
          if(debug) System.out.println("*** breaking on iter " + iter);
          break;
        }
        
        if(Double.isInfinite(norm) || la < 0 || lb < 0 || iter > 100000) {
          System.out.println("GammaMLE diverged, reverting to MoM estimate");
          la = xbar / (lb = sx * sx / xbar);
          break;
        }
      }
      
      // save the results to the class members.
      alpha = la;
      beta = lb;
      
      if(debug) {
        System.out.println("final alpha = " + alpha);
        System.out.println("final beta = " + beta);
      }
    }
  }
  
  private final static class WeightedLazyFit 
    implements Distribution.WeightedLazyFit 
  {
    private double sw, swx, swlx, swxsq;
    
    public void Initialize() { sw = swx = swlx = swxsq = 0; }
    
    public Distribution Finalize() {
      if(sw == 0) return null;
      else {
        double xbar, s, lxbar;
        xbar = swx / sw;
        lxbar = swlx / sw;
        s = Math.sqrt(swxsq / sw - xbar * xbar);
        
        return new GammaDist(new MLEFit(xbar, lxbar, s));
      }
    }
    
    public void AddPoint(double weight, double point) {
      sw += weight;
      swx += weight * point;
      swxsq += weight * point * point;
      swlx += weight * Math.log(point);
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
    
  static {
    WeightedLazyFitFactory = new __WeightedLazyFitFactory();
  }
  
  /**
   * Creates a new instance of GammaDist 
   */
  public GammaDist(double alpha, double beta) {
      //    assert (alpha > 0 && beta > 0);
    paramValues[0] = alpha;
    paramValues[1] = beta;
    lnGamAlpha = SpecialFunc.lnGamma(paramValues[0]);
    betaAlpha = Math.pow(paramValues[1], paramValues[0]);
  }
  
  /**
   * Creates a new instance of GammaDist using the MLE fit object.
   */
  public GammaDist(MLEFit mle) {
    this(mle.alpha, mle.beta);
  }  

  /** Convolve the distribution N times with itself (if implemented)
   */
  public Distribution convolve(int N) {
    return new GammaDist(paramValues[0] * N, paramValues[1]);
  }
  
  /** Get the specified central moment from the distribution.
   */
  public Double getCentralMoment(int i) {
    getRawMoment(i);
    return StdImpl.centralMoment(moments, i);
  }
  
  /** Get the probability that a variable with this distribution lies below or
   * equal to the specified value.
   */
  public Double getCumulativeProb(double v) {
    if(v <= 0) return new Double(0);
    return new Double(SpecialFunc.incGamma(paramValues[0], v / paramValues[1]));
  }
  
  /**
   * Get the value probabilility function (p.d.f.) evaluated at the given point.
   */
  public Double getProbability(double v) {
    if(v <= 0 || Double.isInfinite(lnGamAlpha))
      return new Double(0);
    else {
      double r = -lnGamAlpha - paramValues[0] * Math.log(paramValues[1])
        + (paramValues[0] - 1) * Math.log(v) - v / paramValues[1];
      if(Double.isNaN(r)) {
        System.out.println(getDistributionInstance() 
          + "::getProbability(" + v + ") returning NaN!!!!");
        System.exit(1);
      }
      return new Double(Math.exp(r));
    }
  }

  /** Get an instance string for the distribution (this would be for example
   * "Normal(0, 1)" for the standard normal distribution).
   */
  public String getDistributionInstance() {
    return distrName + "(" + paramValues[0] + ", " + paramValues[1] + ")";
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
  
  /** Get the mean of the distribution.
   */
  public Double getMean() {
    return new Double(paramValues[0] * paramValues[1]);
  }
  
  /** Get the number of parameters required by the distribution.
   */
  public int getParameterCount() {
    return paramCount;
  }
  
  /** Get the name of the ith parameter required by the distribution.
   */
  public String getParameterName(int i) {
    return paramNames[i];
  }
  
  /** Get the allowed parameter value ranges for the ith parameter.
   */
  public Double[] getParameterRange(int i) {
    return paramRanges[i];
  }
  
  /** Get the value of the ith parameter of this instance of the distribution.
   */
  public double getParameterValue(int i) {
    if(i < 0 || i >= paramCount) return 0;
    return paramValues[i];
  }
  
  /** Get the parameter values of this instance of the distribution.
   */
  public double[] getParameterValues() {
    return paramValues;
  }
  
  /** Get the specified raw moment from the distribution.
   */
  public Double getRawMoment(int i) {
    if(i == 0) return new Double(1.0);
    int N = moments.size();
    for(int j=N+1; j<=i; j++) {
      moments.add(new Double(Math.exp(Math.log(paramValues[1]) * j 
        + SpecialFunc.lnGamma(paramValues[0] + j) - lnGamAlpha)));
    }
    return (Double)moments.get(i-1);
  }
  
  /** Get the standard deviation of the distribution.
   */
  public Double getStdDev() {
    return new Double(Math.sqrt(paramValues[0]) * paramValues[1]);
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
    return StdImpl.sampleBuckets(this, minValue, maxValue, 
      nBuckets);
  }
  
  /** Sample the distribution at evenly spaced intervals.
   */
  public Distribution sampleStepped(double minValue, double stepSize, 
    int nSteps) 
  {
    return StdImpl.sampleStepped(this, minValue, stepSize, nSteps);
  }
  
  /** Generate random sample from distribution (using Uniform random number
   * generator supplied, or java built-in (Math.random()) if null supplied).
   * @param n The number of values to be simulated (optimized for > 1)
   * @param rand The random number generator to use (if non-null) or use
   *             Math.random() if not supplied (null).
   * @return An array of length n with the simulated values  */
  public double[] simulateValues(int n, IUniformRandom rand) {
    double T, Y, AB = paramValues[0] * paramValues[1];
    IUniformRandom r = (rand == null) ? StdImpl.rand : rand;
    double[] rv = new double[n];
    double U;
    for(int i=0; i<n; i++) {
      do {
        U = r.getNext();
        while(U == 0) U = r.getNext();
        Y = - Math.log(U) * AB;
        double D = Y / AB;
        double am1 = paramValues[0] - 1;
        T = Math.pow(D, am1) * Math.exp((1 - D) * am1);
        U = r.getNext();
      } while(U > T);
      rv[i] = Y;
    }
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
  public Distribution truncate(Double minValue, Double maxValue, boolean redistribute) {
    return new TruncatedDistribution(this, minValue, maxValue, redistribute);
  }
  
  /** Calculate the quantile of the distribution at probability p
   */
  public Double getQuantile(double p) {
    // need inverse IncGamma for this.  Will write later.
    return null;
  }
  
  public String toString() {
    return getDistributionInstance();
  }
  
}
