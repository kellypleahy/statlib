package statlib;

import java.util.*;
import java.io.*;

/**
 * Represents a frequency distribution with a set of user defined points and
 *   frequencies.  Supports non-normalized frequencies and non-unique point
 *   values.  Also supports truncation at endpoints.
 *
 * @author Kelly Leahy 
 *     (<A HREF="mailto:kellyleahy@swbell.net">kellyleahy@swbell.net</A>)
 * @version 1.0
 */

public class FrequencyDist implements Distribution {
  
  private final static String distrName = "Frequency";
  
  private Vector moments = new Vector();
  private double[] points, freq;
  private double totalfreq;
  
  private void sortPoints() {
    class pointData {
      public double p, f;
      public pointData(double p, double f) { this.p = p; this.f = f; }
    }

    int N = freq.length;
    pointData[] oa = new pointData[N];
    for(int i=0; i<N; i++)
      oa[i] = new pointData(points[i], freq[i]);
    Arrays.sort(oa, 
      new Comparator() {
        public int compare(Object object1, Object object2) {
          return Double.compare(((pointData)object1).p, ((pointData)object2).p);
        }
      });
    for(int i=0; i<N; i++) {
      points[i] = oa[i].p;
      freq[i] = oa[i].f;
    }
  }

  /** 
   * Create an instance of a frequency distribution using the specified
   *   parameters for the data points and frequencies.
   * @param points The list of data points associated with the frequencies, need
   *                 not be unique or sorted.
   * @param freq The frequencies associated with the data points, need not be
   *               normalized (sum to 1.0).
   */
  public FrequencyDist(double[] points, double[] freq) {
    this.points = (double[])points.clone();
    this.freq = (double[])freq.clone();
    //    assert (points.length == freq.length && freq.length > 0);
    totalfreq = 0.0;
    for(int i=0; i<freq.length; i++)
      totalfreq += freq[i];
    //    assert (totalfreq != 0.0);
    sortPoints();
  }
  


  /**
   * Sample the distribution at evenly spaced intervals, using the min value
   * and max value specified to make the correct number of buckets.
   */
  public Distribution sampleBuckets(double minValue, double maxValue, int nBuckets) 
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
   * Get the collection of points on which the distribution is based (this is
   *   a copy of the array and won't affect the behavior of the distribution)
   */
  public double[] getPoints() {
    return (double[])points.clone();
  }
  
  /**
   * Get the collection of frequencies associated with the points in the
   *   distribution (this is a copy of the array and won't affect the behavior
   *   of the distribution)
   */
  public double[] getFrequencies() {
    return (double[])freq.clone();
  }
  
  /**
   * Get the allowed parameter value ranges for the ith parameter.
   */
  public Double[] getParameterRange(int i) {
    return null;
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
    int N = freq.length;
    double m = minValue == null ? Double.MIN_VALUE : minValue.doubleValue(),
      M = maxValue == null ? Double.MAX_VALUE : maxValue.doubleValue();
    if(!redistribute) {
      // no redistribution, need to put extra weight on endpoints.
      // this is really easy, since the points need not be unique.
      double[] p = (double[])points.clone(), f = (double[])freq.clone();
      for(int i=0; i<N; i++) {
        double v = p[i];
        if(v < m) v = m;
        else if(v > M) v = M;
        else continue;
        p[i] = v;
      }      
      return new FrequencyDist(p, f);
    } else {
      // need to redistribute the frequency on the tails
      // this is really easy, since the frequencies need not be normalized.
      int c = 0;
      for(int i=0; i<N; i++) if(points[i] >= m && points[i] <= M) c++;
      double[] p = new double[c], f = new double[c];
      for(int i=0, j=0; i<N && j<c; i++)
        if(points[i] >= m && points[i] <= M) {
          p[j] = points[i];
          f[j] = freq[i];
          j++;
        }
      return new FrequencyDist(p, f);
    }
  }
  
  /**
   * Get the standard deviation of the distribution.
   */
  public Double getStdDev() {
    return new Double(Math.sqrt(getCentralMoment(2).doubleValue()));
  }
  
  /**
   * Get an instance string for the distribution (this would be for example
   * "Normal(0, 1)" for the standard normal distribution).
   */
  public String getDistributionInstance() {
    return distrName + "(special)";
  }
  
  /**
   * Get the probability that a variable with this distribution lies below or
   * equal to the specified value.
   */
  public Double getCumulativeProb(double v) {
    double r = 0.0;
    for(int i=0; i<freq.length && points[i] <= v; i++)
      r += freq[i];
    return new Double(r / totalfreq);
  }
  
  /**
   * Get the value probabilility function (p.d.f.) evaluated at the given point.
   */
  public Double getProbability(double v) {
    double r = 0.0;
    for(int i=0; i<freq.length && points[i] <= v; i++)
      if(points[i] == v) r += freq[i];
    return new Double(r / totalfreq);
  }

  /**
   * Get the specified central moment from the distribution.
   */
  public Double getCentralMoment(int i) {
    if(i < 0)
      return null;
    else if(i == 0)
      return new Double(1.0);
    else if(i == 1)
      return new Double(0.0);
    else {
      // first, force calculation of raw moments
      getRawMoment(i);
      return StdImpl.centralMoment(moments, i);
    }
  }
  
  /**
   * Get the name of the ith parameter required by the distribution.
   */
  public String getParameterName(int i) {
    return null;
  }
  
  /**
   * Get the value of the ith parameter of this instance of the distribution.
   */
  public double getParameterValue(int i) {
    return 0.0;
  }
  
  /**
   * Get the number of parameters required by the distribution.
   */
  public int getParameterCount() {
    return 0;
  }
  
  /**
   * Get the parameter values of this instance of the distribution.
   */
  public double[] getParameterValues() {
    return null;
  }
  
  /** Generate random sample from distribution (using Uniform random number
   * generator supplied, or java built-in (Math.random()) if null supplied).
   * @param n The number of values to be simulated (optimized for > 1)
   * @param rand The random number generator to use (if non-null) or use 
   *             Math.random() if not supplied (null).
   * @return An array of length n with the simulated values */
  public double[] simulateValues(int n, IUniformRandom rand) {
    IUniformRandom rgen = 
      rand == null ? StdImpl.rand : rand;
    int c = 0;
    double d[] = new double[n];
    int N = freq.length;
    for(int i=0; i<n; i++) {
      double F = 0.0;
      double r = rgen.getNext();
      int j;
      for(j=0; j<N; j++) {
        F += freq[j];
        if(F / totalfreq >= r) {
          d[i] = points[j];
          break;
        }
      }
      if(j == N) d[i] = points[j-1];
    }
    return d;
  }
  
  /**
   * Get the mean of the distribution.
   */
  public Double getMean() {
    return getRawMoment(1);
  }
  
  /**
   * Get the sample mean of the distribution
   */
  public double getSampleVar() {
    int N = freq.length;
    return getCentralMoment(2).doubleValue() * N / (N - 1);
  }
  
  /**
   * Get the specified raw moment from the distribution.
   */
  public Double getRawMoment(int i) {
    if(i == 0) 
      return new Double(1.0);
    else if(i < 0)
      return null;      
    else if(moments.size() >= i)
      return (Double)moments.get(i-1);
    else {
      int m = moments.size();
      double s[] = new double[i-m];
      for(int j=0; j<freq.length; j++) {
        double v = 1.0;
        for(int k=0; k<i-m; k++) {
          v *= points[j];
          s[k] += v * freq[j];
        }
      }
      for(int k=0; k<i-m; k++)
        moments.add(new Double(s[k] / totalfreq));
      return new Double(s[i-m-1] / totalfreq);
    }
  }
  
  /**
   * Save the values in the distribution to a file (space delimited) with the 
   *   point value followed by the frequency for the point.  Note: frequencies 
   *   are not normalized (do not sum to 1.0).
   */
  public void saveToFile(String fileName) throws IOException {
    FileWriter fw = new FileWriter(fileName);
    BufferedWriter bw = new BufferedWriter(fw);
    for(int i=0; i<freq.length; i++)
      bw.write(points[i] + "," + freq[i] + "\n");
    bw.close();
    fw.close();
  }

  /**
   * Load the values for a new distribution from a file (space delimited with
   *   the first value on each line as the point value, and the second value
   *   on each line as the frequency.  Note: The frequencies need not be
   *   normalized (sum to 1.0) and the point values need not be sorted.
   */
  public static Distribution loadFromFile(String fileName) throws IOException {
    FileReader fr = new FileReader(fileName);
    BufferedReader br = new BufferedReader(fr);
    Vector pts = new Vector(), freqs = new Vector();
    int i=0;
    while(br.ready()) {
      i++;
      String line = br.readLine();
      line = line.trim();
      if(line.length() == 0) continue;
      String pt = line.substring(1, line.indexOf((int)' ')-1);
      String rest = line.substring(pt.length()+1).trim();
      try {
        Double d1 = Double.valueOf(pt);
        Double d2 = Double.valueOf(rest);
        pts.add(d1);
        freqs.add(d2);
      } catch (NumberFormatException e) {
        System.err.println("Skipping line " + i + " due to exception: "
          + e);
        continue;
      }
    }
    int N = pts.size();
    double[] p = new double[N], f = new double[N];
    for(int j=0; j<N; j++) {
      p[j] = ((Double)pts.get(j)).doubleValue();
      f[j] = ((Double)freqs.get(j)).doubleValue();
    }
    return new FrequencyDist(p, f);
  }
  
  /** Get the mth raw moment of the distribution limited to the range
   *  between a and b.
   */
  public Double getLimitedRawMoment(int m, Double a, Double b) {
    double lm, lM;
    if(a == null) lm = Double.NEGATIVE_INFINITY;
    else lm = a.doubleValue();
    if(b == null) lM = Double.POSITIVE_INFINITY;
    else lM = b.doubleValue();
    
    double s = 0;
    for(int i=0; i<freq.length; i++)
      s += freq[i] * Math.pow(Math.min(lM, Math.max(lm, points[i])), m);
    return new Double(s);
  }
  
  /** Get the first m raw moments of the distribution limited to the range
   *  between a and b.
   */
  public Double[] getLimitedRawMoments(int m, Double a, Double b) {
    double lm, lM;
    if(a == null) lm = Double.NEGATIVE_INFINITY;
    else lm = a.doubleValue();
    if(b == null) lM = Double.POSITIVE_INFINITY;
    else lM = b.doubleValue();

    double s[] = new double[m+1];
    for(int i=0; i<freq.length; i++)
      for(int j=0; j<=m; j++)
        s[j] += freq[i] * Math.pow(Math.min(lM, Math.max(lm, points[i])), j);
    Double rv[] = new Double[m+1];
    for(int i=0; i<=m; i++)
      rv[i] = new Double(s[i]);
    return rv;
  }
  
  /** Get the mth central moment of the distribution limited to the range
   *  between a and b
   */
  public Double[] getLimitedCentralMoments(int m, Double a, Double b) {
    Double r[] = getLimitedRawMoments(m, a, b);
    if(r == null) return null;
    Vector v = new Vector(m+1);
    v.addAll(Arrays.asList(r));
    Double rv[] = new Double[m+1];
    for(int i=0; i<=m; i++)
      rv[i] = StdImpl.centralMoment(v, i);
    return rv;
  }
  
  /** Get the mth central moment of the distribution limited to the range
   *  between a and b
   */
  public Double getLimitedCentralMoment(int m, Double a, Double b) {
    Double r[] = getLimitedRawMoments(m, a, b);
    if(r == null) return null;
    Vector v = new Vector(m+1);
    v.addAll(Arrays.asList(r));
    return StdImpl.centralMoment(v, m);
  }
  
    /** Convolve the distribution N times with itself, and put it into 
     * integer-valued buckets
     *
     * Assumes the points in this distributions are already integers
     */
    
    // kpl1: changed to return FrequencyDist rather than Distribution.
    //   it's stupid to return Distribution when this isn't part of the 
    //   interface, since the caller will always know they're dealing with
    //   the specific "FrequencyDist" type.
    public FrequencyDist iConvolve(int N) {
	// here we need to do the fourier transform on the variables.
	// first, figure out how many intervals we need for the FFT
	// get the maximum and minimum values of the distribution
	int M = freq.length;
	double mv = Double.POSITIVE_INFINITY, MV = Double.NEGATIVE_INFINITY;
	for(int i=0; i<M; i++) {
	    double p = points[i];
	    if(p < mv) mv = p;
	    if(p > MV) MV = p;
	}
	// now mv is min value, and MV is max value of distribution.
	
	//find the nearest power of 2 that the new distribution will fit into
	double range = N*MV - mv, min = mv ;
	range = Math.pow(2, (int)(Math.log(range)/Math.log(2)));
	int nBuckets = (int)range;
	double bucketsize = range / nBuckets;

	// now that we have the range, the number of buckets, and the bucketsize
	//   put the distribution into the approximated form
	FrequencyDist d = (FrequencyDist)sampleStepped(min, bucketsize, nBuckets);
	double[] fqs = d.getFrequencies();
	double [] pts = d.getPoints();
	
	// make sure we got what we expected
	//    assert (fqs.length == nBuckets);
	
	double[] imag = new double[fqs.length];
	double[] fqout = new double[fqs.length], imout = new double[fqs.length];
	try {
	    Fourier.FFT(fqs, imag, fqout, imout);
	    Fourier.impower(fqout, imout, N);;
	    Fourier.iFFT(fqout, imout, fqs, imag);;
	} catch(Exception e) {
	    System.err.println("FFT threw exception: " + e.getMessage());
	    e.printStackTrace();
	}
	
	double shift = mv * (N - 1);
	for(int i=0; i<nBuckets; i++){
	    pts[i] += shift;
	}

	return new FrequencyDist(pts, fqs);
    }
  
  

    /** Convolve this distribution with another FrequencyDist
     *
     * Assumes the points in this distributions are already integers
     */
    
    // kpl1: changed to return FrequencyDist rather than Distribution.
    //   it's stupid to return Distribution when this isn't part of the 
    //   interface, since the caller will always know they're dealing with
    //   the specific "FrequencyDist" type.
    public FrequencyDist iConvolve(FrequencyDist other) {
	// here we need to do the fourier transform on the variables.
	// first, figure out how many intervals we need for the FFT
	// get the maximum and minimum values of the distribution
	int thisM = freq.length;
	double thismv = Double.POSITIVE_INFINITY, thisMV = Double.NEGATIVE_INFINITY;
	for(int i=0; i<thisM; i++) {
	    double p = points[i];
	    if(p < thismv) thismv = p;
	    if(p > thisMV) thisMV = p;
	}

	int otherM = other.freq.length;
	double[] otherPts = other.getPoints();
	double othermv = Double.POSITIVE_INFINITY, otherMV = Double.NEGATIVE_INFINITY;
	for(int i=0; i<otherM; i++) {
	    double p = otherPts[i];
	    if(p < othermv) othermv = p;
	    if(p > otherMV) otherMV = p;
	}

	double min = Math.min(thismv, othermv);
	double max = thisMV + otherMV;
	
	//find the nearest power of 2 that the new distribution will fit into
	double range = max - min + 1;
	range = Math.pow(2, (int)Math.ceil((Math.log(range)/Math.log(2))));
	int nBuckets = (int)range;
	double bucketsize = range / nBuckets;
	
	//assert(bucketsize == 1.0);
	
	// now that we have the range, the number of buckets, and the bucketsize
	//   put the distribution into the approximated form
	FrequencyDist d1 = (FrequencyDist)sampleStepped(min, bucketsize, nBuckets);
	double[] fqs1 = d1.getFrequencies();
	double [] pts1 = d1.getPoints();
	FrequencyDist d2 = (FrequencyDist)other.sampleStepped(min, bucketsize, nBuckets);
	double[] fqs2 = d2.getFrequencies();
	double [] pts2 = d2.getPoints();
	
	// make sure we got what we expected
	//    assert (fqs.length == nBuckets);
	
	double[] imag1 = new double[fqs1.length];
	double[] imag2 = new double[fqs2.length];
	double[] fqout1 = new double[fqs1.length], imout1 = new double[fqs1.length];
	double[] fqout2 = new double[fqs2.length], imout2 = new double[fqs2.length];
	double[] fqresult = new double[fqs2.length], imresult = new double[fqs2.length];
	double[] out = new double[fqs2.length];
	try { 
	    Fourier.FFT(fqs1, imag1, fqout1, imout1);
	    Fourier.FFT(fqs2, imag2, fqout2, imout2);
	    Fourier.immult(fqout1, imout1, fqout2, imout2, fqresult, imresult);;
	    Fourier.iFFT(fqresult, imresult, fqs1, imag1);;
	} catch(Exception e) {
	    System.err.println("FFT threw exception: " + e.getMessage());
	    e.printStackTrace();
	}
	
	for(int i=0; i<nBuckets; i++){
	    pts1[i] += min;
	}

	return new FrequencyDist(pts1, fqs1);
    }
  


  /** Convolve the distribution N times with itself (if implemented)
   */
  public Distribution convolve(int N) {
    // here we need to do the fourier transform on the variables.
    // first, figure out how many intervals we need for the FFT
    // get the maximum and minimum values of the distribution
    int M = freq.length;
    double mv = Double.POSITIVE_INFINITY, MV = Double.NEGATIVE_INFINITY;
    for(int i=0; i<M; i++) {
      double p = points[i];
      if(p < mv) mv = p;
      if(p > MV) MV = p;
    }
    // now mv is min value, and MV is max value of distribution.
    
    // for now, let's aim for 65536 buckets in final distribution
    int nBuckets = 65536;
    double range = N*MV - mv, min = mv ;
    double bucketsize = range / nBuckets;
    
    // now that we have the range, the number of buckets, and the bucketsize
    //   put the distribution into the approximated form
    FrequencyDist d = (FrequencyDist)sampleStepped(min, bucketsize, nBuckets);
    double[] fqs = d.getFrequencies();
    double [] pts = d.getPoints();

    // make sure we got what we expected
    //    assert (fqs.length == nBuckets);
    
    double[] imag = new double[fqs.length];
    double[] fqout = new double[fqs.length], imout = new double[fqs.length];
    try {
	Fourier.FFT(fqs, imag, fqout, imout);
	Fourier.impower(fqout, imout, N);;
	Fourier.iFFT(fqout, imout, fqs, imag);;
    } catch(Exception e) {
	System.err.println("FFT threw exception: " + e.getMessage());
	e.printStackTrace();
    }

    double shift = mv * (N - 1);
    for(int i=0; i<nBuckets; i++){
	pts[i] += shift;
    }

    return new FrequencyDist(pts, fqs);
  }
  
  /** Calculate the quantile of the distribution at probability p
   */
  public Double getQuantile(double p) {
    //KPL: Not tested!!!
    double tp = p * totalfreq;
    double F0 = freq[0], F1 = freq[0];
    int N = freq.length;
    for(int i=1; i<N; i++) {
      F1 = F0 + freq[i];
      if(F1 >= tp && F0 < tp) {
        double rv = ((tp-F0) * points[i-1] + (F1-tp) * points[i]) / (F1-F0);
        return new Double(rv);
      }
      F0 = F1;
    }
    // handle the cases that aren't caught by the above code (endpoints)
    //if(tp > totalfreq - freq[N-1])
    //  return new Double(points[N-1]);
    //else
    if(tp < freq[0])
      return new Double(Double.NEGATIVE_INFINITY);
    else
      return null;
  }
  

    //added by dmo1@cs.wustl.edu - a hastily written method for debugging
    public String toString() {
	String ans ="";
	//\n---------------------------------------------------\n";
	//ans += "Values\t\tCumulative\n";
	for(int i=0; i<freq.length; i++)
	    if(freq[i] > .0001) ans += (points[i] + " -> " + freq[i] + "\n");
	//System.out.println(points[i] + " -> " + freq[i] + "\n");
	//print cumulative prob at arbitrary points over somewhat arbitrary range
	//for(int i=0; i<16; i++) {
	//ans+=  i + "\t\t" + getCumulativeProb(i)  + "\n";
	 //}
	
	 //	for(int i=0; i<points.length && points[i] <=2.1; i++) 
	 //  System.out.println(points[i] + " -> " + freq[i]);
	
	 //     for(int i=0; i<points.length; i++) {
	 // 	    if(freq[i] > .001) {
	 // 	    cumProb += freq[i];
	 // 	    ans += points[i] + "\t\t" + freq[i] + "\t\t" + cumProb + "\n";
	// 	    }
	// 	}
	return ans;
    }
}
