/*
 * MixtureFit.java
 *
 * Created on April 2, 2003, 5:05 PM
 */

package statlib;

/**
 *
 * @author  Kelly Leahy
 */
public class MixtureFit {
  
  private Distribution.WeightedLazyFitFactory m_FitFactory; 
  private Distribution.WeightedLazyFit m_Fitters[];
  private Distribution m_Distrs[];
  private double[] m_Props;
  private double[] m_Data;
  private int m_nComp;
  private int m_nDataPoints;
  
  public interface CustomInitializer {
    public boolean init(double data[], int nComp, Distribution distrs[], 
      double props[]);
  }
  
  private CustomInitializer m_CustomInit = null;
  
  private static double MAX_LOG_INT = Math.log(Integer.MAX_VALUE);

  /** Creates a new instance of MixtureFit */
  public MixtureFit(double[] data, int nComp, 
    Distribution.WeightedLazyFitFactory FitFactory) 
  {
    m_FitFactory = FitFactory;
    m_nComp = nComp;
    m_Fitters = new Distribution.WeightedLazyFit[m_nComp];
    m_Distrs = new Distribution[m_nComp];
    m_Props = new double[m_nComp];
    m_Data = (double[])data.clone();
    m_nDataPoints = data.length;
  }
  
  public MixtureFit(double[] data, int nComp, 
    Distribution.WeightedLazyFitFactory FitFactory,
    CustomInitializer CustomInit)
  {
    this(data, nComp, FitFactory);
    // set up the custom initializer
    this.m_CustomInit = CustomInit;
  }
  
  public void InitializeFit() {
    // if we have a custom initializer, try to use it...
    if(m_CustomInit != null) {
      m_Props = new double[m_nComp];
      m_Distrs = new Distribution[m_nComp];
      
      // if it initializes the data, then exit, otherwise, let the other code 
      //   handle it
      if(m_CustomInit.init(m_Data, m_nComp, m_Distrs, m_Props)) {
        // we need to initialize Fitters, since the custom initializer won't.
        for(int i=0; i<m_nComp; i++)
          (m_Fitters[i] = m_FitFactory.CreateInstance()).Initialize();
        return;
      }
    }
    
    // later we'll put exceptions in here. (4 is arbitrary choice)
    if(m_nDataPoints < m_nComp * 4) {
      System.out.println("Houston, we got problems: too few data points for"
        + " this number of components");
      System.exit(1);
    }
    
    // partition the data into different buckets
    int nInClass[] = new int[m_nComp];
    int nRem = m_nDataPoints;
    
    for(int i=0; i<m_nComp-1; i++)
      nRem -= (nInClass[i] = 
        (int)(Math.random() * (nRem - 4 * (m_nComp - i))) + 4);
    nInClass[m_nComp-1] = nRem;
    
    // sanity check, remove later
    int sum = 0;
    for(int i=0; i<m_nComp; i++)
      sum += nInClass[i];
    assert (sum == m_nDataPoints) : "Sum == m_nDataPoints";
    
    // sorted list of points
    double pts[] = (double[])m_Data.clone();
    java.util.Arrays.sort(pts);
    
    // fit the distributions to the sample partition.
    int idx = 0;
    for(int i=0; i<m_nComp; i++) {
      (m_Fitters[i] = m_FitFactory.CreateInstance()).Initialize();
      for(int j=idx; j < idx + nInClass[i]; j++)
        m_Fitters[i].AddPoint(1.0, pts[j]);
      idx += nInClass[i];
      m_Distrs[i] = m_Fitters[i].Finalize();
      m_Props[i] = nInClass[i] / (double)m_nDataPoints;
    }
  }
  
  double Likelihood(int i, double v) {
    Double d = m_Distrs[i].getProbability(v);
    
    //System.out.println("Likelihood = " + d);
    
    if(d == null)
      return 0;
    else
      return d.doubleValue() * m_Props[i];
  }
  
  double Iterate() {
    double weights[] = new double[m_nComp];
    double totweights[] = new double[m_nComp];
    double NLL = 0.0;
    
    //System.out.println("initializing in iterate");
    
    for(int i=0; i<m_nComp; i++)
      m_Fitters[i].Initialize();
    
    //System.out.println("getting likelihoods");
    
    for(int j=0; j<m_nDataPoints; j++) {
      double sw = 0;
      //System.out.println("getting weights (point " + j + ")");
      for(int i=0; i<m_nComp; i++) {
        //System.out.println("getting likelihoods (component " + i + ")");
        sw += (weights[i] = Likelihood(i, m_Data[j]));
      }
      //System.out.println("adding point " + j + "to fitters");
      for(int i=0; i<m_nComp; i++) {
        m_Fitters[i].AddPoint(weights[i] / sw, m_Data[j]);
        totweights[i] += weights[i] / sw;
      }
      
      NLL -= Math.log(sw);
    }
    
    //System.out.println("Getting distributions");
    
    double sw = 0.0;
    for(int i=0; i<m_nComp; i++) {
      //System.out.println("Fitting component " + i);
      m_Distrs[i] = m_Fitters[i].Finalize();
      sw += totweights[i];
    }
    
    for(int i=0; i<m_nComp; i++) {
      m_Props[i] = totweights[i] / sw;
    }
    
    return NLL;
  }
  
  public double[] GetProportions() {
    return (double[])m_Props.clone();
  }
  
  public Distribution[] GetDistributions() {
    return (Distribution[])m_Distrs.clone();
  }
  
  public String DebugInfo(int level) {
    String r = new String();
    
    if(level >= 1)
      r += "nComp = " + m_nComp + "\n";
    if(level >= 2)
      for(int i=0; i<m_nComp; i++)
        r += "prop[" + (i+1) + "] = " + m_Props[i] + "\n";
    if(level >= 3)
      for(int i=0; i<m_nComp; i++)
        r += "distr[" + (i+1) + "] = '" + m_Distrs[i] + "'\n";
    return r;
  }
  
  private int nComb(int n, int r) {
    double s = 0.0;
    for(int i=r+1; i<=n; i++)
      s += Math.log(i);
    for(int i=1; i<=n-r; i++)
      s -= Math.log(i);
    if(s > MAX_LOG_INT)
      return Integer.MAX_VALUE;
    else
      return (int)Math.round(Math.exp(s));
  }
  
  public double Solve(int MaxStarts, int MaxIters, double Tolerance, 
    boolean bDebug) 
  {
    double MinNLL = Double.POSITIVE_INFINITY;
    double MinProps[] = null;
    Distribution MinDistr[] = null;
    int nStarts = 0;
    
    //{
      /* Don't let it start more times than is theoretically necessary
       *   Theoretically, there are a total of nCr(n+r-1, r-1) different ways
       *   we can partition the first n items from the number line into r 
       *   buckets.  However, since we're dealing with free items only, and we
       *   are putting at least Max(M, 0.25 * (1/r)) items in each bucket, we 
       *   need (n-M*r)+r-1 = (n-(M-1)*r-1) for the first nComb parameter.
       *
       * Once we've determined the theoretically maximum number of possible
       *   starting points, we need to know the theoretically maximum number
       *   of random starts necessary to "cover" this space with 95% confidence.
       */
    /* 
      // not currently implemented
      int MaxThStarts = nComb(m_nDataPoints + m_nComp - 1, m_nComp - 1);
      if(MaxStarts > MaxThStarts)
        MaxStarts = MaxThStarts;
    */
    //}
    
    while(nStarts < MaxStarts) {
      int nIters = 0;
      double lastNLL, NLL;
      
      NLL = Double.POSITIVE_INFINITY;
      
      InitializeFit();
      
      if(bDebug) {
        System.out.println("Fit initialized.");
        System.out.println(DebugInfo(4));
      }

      do {
        nIters++;
        //System.out.println("iteration " + nIters);
        lastNLL = NLL;
        NLL = Iterate();
        if(bDebug) {
          System.out.println(DebugInfo(4));
          System.out.println("NLL = " + NLL);
        }
      } while(Math.abs(lastNLL - NLL) > Tolerance && nIters < MaxIters);
      
      if(NLL < MinNLL) {
        // save the information
        MinNLL = NLL;
        MinProps = (double[])m_Props.clone();
        MinDistr = (Distribution[])m_Distrs.clone();
        //if(bDebug) System.out.println(DebugInfo(4));
      }
      nStarts++;
    }
   
    m_Props = MinProps;
    m_Distrs = MinDistr;
    
    if(bDebug) {
      System.out.println("Final results...");
      System.out.println(DebugInfo(4));
    }
    
    return MinNLL;
  }
}
