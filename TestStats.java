/*
 * TestStats.java
 *
 * Created on September 24, 2002, 2:31 PM
 */

package statlib;

import java.io.*;

/**
 *
 * @author  KLeahy
 */
public class TestStats {

    public static void testDist(Distribution d) {
        System.out.println("Testing distribution: " + d.getDistributionInstance());

        System.out.println("Calculating F(-inf): ");
        try {
            System.out.println(d.getCumulativeProb(Double.NEGATIVE_INFINITY));
        }
        catch (Exception e) {
            System.out.println(e);
        }

        System.out.println("Calculating F(+inf): ");
        try {
            System.out.println(d.getCumulativeProb(Double.POSITIVE_INFINITY));
        }
        catch (Exception e) {
            System.out.println(e);
        }

        System.out.println("Calculating mean: ");
        Double m = null;
        try {
            System.out.println(m = d.getMean());
        }
        catch (Exception e) {
            System.out.println(e);
        }

        Double dv = null;
        System.out.println("Calculating stddev: ");
        try {
            System.out.println(dv = d.getStdDev());
        }
        catch (Exception e) {
            System.out.println(e);
        }

        System.out.println("Calculating stddev (2 c.m.): ");
        try {
            System.out.println(Math.sqrt(d.getCentralMoment(2).doubleValue()));
        }
        catch (Exception e) {
            System.out.println(e);
        }

        System.out.println("Calculating 2nd moment (hard way): ");
        try {
            System.out.println(dv.doubleValue() * dv.doubleValue()
                               + m.doubleValue() * m.doubleValue());
        }
        catch (Exception e) {
            System.out.println(e);
        }

        System.out.println("Calculating 2nd moment (func): ");
        try {
            System.out.println(d.getRawMoment(2));
        }
        catch (Exception e) {
            System.out.println(e);
        }

        System.out.println("Calculating 3 limited raw moments (-inf, +inf): ");
        try {
            Double da[] = d.getLimitedRawMoments(3, null, null);
            for (int i = 1; i <= 3; i++) {
                System.out.println("  [" + i + "] = " + da[i]);
            }
        }
        catch (Exception e) {
            System.out.println(e);
        }

        System.out.println("Calculating 3 raw moments: ");
        try {
            for (int i = 1; i <= 3; i++) {
                System.out.println("  [" + i + "] = " + d.getRawMoment(i));
            }
        }
        catch (Exception e) {
            System.out.println(e);
        }

        System.out.println("Calculating quantile and testing: ");
        try {
            Double q = d.getQuantile(0.8);
            if (q != null) {
                System.out.println("quantile(0.8) = " + q +
                                   ", F(quantile(0.8)) = " +
                                   d.getCumulativeProb(q.doubleValue()));
            }
        }
        catch (Exception e) {
            System.out.println(e);
        }
    }

  
  public static double[][] sqrMatrix(double[][] m) {
    int N = m[0].length;
    double[][] r = new double[N][N];
    for(int i=0; i<N; i++)
      for(int j=0; j<N; j++)
        for(int k=0; k<N; k++)
          r[i][j] += m[i][k] * m[k][j];
    return r;
  }
  
  public static double diffMatrix(double[][] m1, double[][] m2) {
    int N = m1[0].length;
    double r = 0.0;
    for(int i=0; i<N; i++)
      for(int j=0; j<N; j++)
        r += m1[i][j] - m2[i][j];
    return r;
  }

  public static void testEM(int nDist, int nPoints, double prop[], 
    double mean[], double var[]) 
  {
    Distribution d[] = new Distribution[nDist];
    double pts[] = new double[nDist], freqs[] = new double[nDist];
    for(int i=0; i<nDist; i++) { 
      d[i] = new NormalDist(mean[i], var[i]);
      pts[i] = i;
      freqs[i] = prop[i];
    }
    
    Distribution fd = new FrequencyDist(pts, freqs);
    
    double data[] = new double[nPoints];
    double sims[] = fd.simulateValues(nPoints, null);
    for(int i=0; i<nPoints; i++)
      data[i] = d[(int)(sims[i])].simulateValues(1, null)[0];
    MixtureFit m = new MixtureFit(data, nDist, 
      NormalDist.WeightedLazyFitFactory);
    m.Solve(nPoints / 2, 1000, 0.01, true);
  }
  
  public static class EM_System_Examples {
    public double[] example_times;
    public int[] operation;
    public int[] state;
    public int[] op_count;
    public double[] logl;
    public EM_System_Examples(double[] t, int[] op, int[] opc, int[] s, double ll[]) {
      example_times = t;
      operation = op;
      op_count = opc;
      state = s;
      logl = ll;
    }
    public void SaveToCSV(String fileName) throws IOException {
      FileWriter fw = new FileWriter(fileName);
      BufferedWriter bw = new BufferedWriter(fw);
      bw.write("Time,Op,State,LL\n");
      for(int i=0; i<example_times.length; i++)
        bw.write(example_times[i] + "," + operation[i] + "," + state[i] 
          + "," + logl[i] + "\n");
      bw.close();
      fw.close();
      //System.exit(0);
    }
  }
  
  public static EM_System_Examples genExamples(int nStates, int nOper, 
    double[][] del, double[][] m, double[][] dev, double[] poper, 
    int start_state, int nExamples) 
  {
    int curr_state = start_state;
    double r[] = new double[nExamples];
    int op[] = new int[nExamples];
    int opc[] = new int[nStates];
    int s[] = new int[nExamples];
    double logl[] = new double[nExamples];
    Distribution[] tx = new Distribution[nStates];
    Distribution[][] dx = new Distribution[nOper][nStates];
    Distribution opx;
    
    System.out.println("In genExamples():1");

    double[][] mr, probM;
    mr = sqrMatrix(del);
    probM = sqrMatrix(mr);
    while(diffMatrix(mr, probM) > 0.01) {
      mr = probM;
      probM = sqrMatrix(mr);
    }
    
    System.out.println("In genExamples():2");

    for(int i=0; i<nStates; i++) {
      double d[] = new double[nStates];
      double f[] = new double[nStates];
      for(int j=0; j<nStates; j++) {
        d[j] = j;
        f[j] = del[i][j];
      }
      tx[i] = new FrequencyDist(d, f);
      for(int j=0; j<nOper; j++) {
        dx[j][i] = new NormalDist(m[j][i], dev[j][i]*dev[j][i]);
      }
    }

    System.out.println("In genExamples():3");

    {
      double d[] = new double[nOper];
      double f[] = new double[nOper];
      for(int j=0; j<nOper; j++) {
        d[j] = j;
        f[j] = poper[j];
      }
      opx = new FrequencyDist(d, f);
    }

    System.out.println("In genExamples():4");

    for(int j=0; j<nExamples; j++) {
      int this_oper = (int)opx.simulateValues(1, null)[0];
      int next_state = (int)tx[curr_state].simulateValues(1, null)[0];
      r[j] = dx[this_oper][curr_state].simulateValues(1, null)[0];
      op[j] = this_oper;
      opc[this_oper]++;
      s[j] = curr_state;
      logl[j] = Math.log(
        dx[this_oper][curr_state].getProbability(r[j]).doubleValue());
      curr_state = next_state;
    }

    System.out.println("In genExamples():5");

    return new EM_System_Examples(r, op, opc, s, logl);
  }
  
  public static void testEM2() {
    EM_System_Examples ex = genExamples(3, 3, 
      new double[][]{{0.4,0.3,0.3},{0.1,0.6,0.3},{0.5,0.3,0.2}},
      new double[][]{{100,300,1000},{50,25,175},{40,100,100}},
      new double[][]{{2.5,7.5,100},{7.5,2.5,25},{7.5,15,15}},
      new double[]{0.2,0.4,0.4},
      0, 10000);

    System.out.println("After genExamples()");
      
    try {
      ex.SaveToCSV("/home/hokey/kpl1/Examples.csv");
    } catch (IOException e) {
      System.out.println(e);
    }

    System.out.println("After dump");
    
    for(int i=0; i<3; i++) {
      double[] data = new double[ex.op_count[i]];
      for(int j=0, k=0; j<10000; j++)
        if(ex.operation[j] == i)
          data[k++] = ex.example_times[j];
      System.out.println("Fitting state " + i);
      MixtureFit m = new MixtureFit(data, 3, 
        NormalDist.WeightedLazyFitFactory);
      m.Solve(1000, 100000, 0.0001, true);
    }
    
  }
  
  
  private static class EMInitThesis implements MixtureFit.CustomInitializer {
    private double p1;
    private Distribution g1, g2;
    
    public EMInitThesis(double p1, Distribution g1, Distribution g2) {
      this.p1 = p1;
      this.g1 = g1;
      this.g2 = g2;
    }
    
    public boolean init(double data[], int nComp, Distribution distrs[], 
      double props[])
    {
      // don't initialize something we don't understand!
      if(nComp != 2) return false;
      
      distrs[0] = g1;
      distrs[1] = g2;
      props[0] = p1;
      props[1] = 1 - p1;
      
      return true;
    }
  }
  
  /**
   * @param args the command line arguments
   */
  public static void main(String[] args) {
      //testDist(new NormalDist(0, 1));
      //testDist(new NormalDist(13, 4));
      //testDist(new GammaDist(1, 12));
      //testDist(new GammaDist(3, 10));
      //testDist(new ExponentialDist(12));
      //testDist(new ChiSquaredDist(8));
      //double[] values = {4,5,6};
      //double[] freqs = {.3,.4,.3};
      //testDist(new FrequencyDist(values, freqs));

      /*
      CharacteristicFunction n[] = new CharacteristicFunction[2], cv;

      n[0] = new CharacteristicFunction.Gaussian(0.0, 1.0);
      n[1] = new CharacteristicFunction.Gaussian(0.0, 3.0);
      cv = new CharacteristicFunction.Convolution(n);

      double[] real = new double[256];
      double[] imag = new double[256];

      Fourier.CharFnToDFT(cv, -16, 16, 256, real, imag);

      for (int i = 0; i < 256; i++) {
          System.out.println("real[" + i + "] = " + real[i] + ", imag[" + i
                             + "] = " + imag[i]);
      }

      double[] realout = new double[256];
      double[] imagout = new double[256];

      Fourier.iFFT(real, imag, realout, imagout);

      System.out.println("iFFT output--------------------");
      for(int i=0; i<256; i++) {
          System.out.println(realout[i]);
      }
      */

     /*
     double[] p = { 1, 2, 6.5 };
     int[] c = { 0, 1, 4 };
     double a = 2, b = -2, g = 6.5, d = -5, h0 = 1.125, h1 = 1.375;
     BucketOptimize bo = new BucketOptimize(p, c, a, b, g, d, h0, h1);
     bo.Solve();

     */

/*
    CharacteristicFunction cmp[] = new CharacteristicFunction[2];
    double prop[] = new double[2];

    cmp[0] = new CharacteristicFunction.Gaussian(200.0, 40.0);
    prop[0] = 0.1;
    cmp[1] = new CharacteristicFunction.Gaussian(20.0, 2.0);
    prop[1] = 0.9;

    CharacteristicFunction m = new CharacteristicFunction.Mixture(
          prop, cmp);

    CharacteristicFunction cvcmp[] = new CharacteristicFunction[3];
    cvcmp[0] = cvcmp[1] = cvcmp[2] = m;
    CharacteristicFunction cv = new CharacteristicFunction.Convolution(cvcmp);

    double real[] = new double[1024], imag[] = new double[1024];
    double realout[] = new double[1024], imagout[] = new double[1024];
    Fourier.CharFnToDFT(cv, 0, 900, 1024, real, imag);
    Fourier.iFFT(real, imag, realout, imagout);
*/

    /*

    System.out.println("iFFT output--------------------");
    for(int i=0; i<realout.length; i++) {
        System.out.println(realout[i]);
    }

    */

/*
    Fourier.CharFnToDFT(cmp[0], 0, 400, 1024, real, imag);
    Fourier.iFFT(real, imag, realout, imagout);

    System.out.println("iFFT output (2)----------------");
    for(int i=0; i<realout.length; i++) {
        System.out.println(realout[i]);
    }
*/

    /*
    
    System.out.println("Starting testEM2()");
    testEM2();
    System.out.println("Finished testEM2()");
    
    */
    
    /*
    GammaDist.MLEFit mle = new GammaDist.MLEFit(103.1843, 4.61960194, 
      18.77356302, true);
    */
    
    /*
    // this is the data from Part C of KellyThesis.
    double points[] = {
      123.48, 81.29, 96.25, 109.84, 79.14, 77.55, 102.44, 84.23, 98.95, 106.65,
      105.73, 106.82, 74.72, 124.39, 123.74, 83.27, 111.3, 93.11, 68.95, 95.89,
      117.12, 81.78, 118.03, 120.58, 113.64, 111.28, 131.75, 117.14, 119.37,
      94.96, 97.22, 90.15, 120.35, 98.6, 138.83, 126.49, 89.25, 94.22, 100.23,
      99.82, 95.07, 101.87, 118.75, 101.4, 102.13, 177.42, 94.46, 63.79, 82.96,
      67.53, 115.22, 96.97, 108.4, 91.42, 95.53, 83.56, 101.04, 99.64, 94.14,
      113.38, 97.88, 108.14, 116.86, 99.88, 130.45, 109.53, 63.34, 99.9, 88.41,
      108.88, 93.36, 124.94, 73.88, 106.74, 101.41, 93.31, 100.39, 128.75,
      105.94, 111.23, 102.51, 73.78, 106.09, 130.43, 91.42, 132.32, 114.3,
      74.48, 59.25, 121.23, 128.39, 104.18, 94.23, 111.33, 95.52, 114.84,
      143.85, 104.43, 112.2, 107.2
    };
    
    Distribution.WeightedLazyFit wlf 
      = GammaDist.WeightedLazyFitFactory.CreateInstance();
    
    wlf.Initialize();
    for(int j=0; j<points.length; j++) {
      wlf.AddPoint(1, points[j]);
    }
    
    Distribution g = wlf.Finalize();
    System.out.println(g);
    */
    
    if(false)
    {
      // this is the data from Part A of KellyThesis.
      double points[] = {
        981.87, 42.22, 40.58, 37.63, 45.34, 38.38, 1040.32, 999.95, 38.88, 36.72, 
        41.14, 39.51, 38.46, 41.86, 39.6, 43.56, 34.62, 37.58, 40.52, 971.61, 
        40.41, 40.95, 994.83, 40.83, 38.87, 38.32, 37.67, 40.11, 41.81, 43.54, 
        36.92, 39.01, 1007.41, 40.11, 981.36, 42.65, 36.81, 987.98, 39.45, 39.37, 
        39.92, 1010.35, 42.58, 39.05, 36.18, 34.21, 39.29, 39.07, 42.24, 40.41, 
        39.39, 989.21, 42.21, 37.92, 37.88, 39.04, 36.56, 38.56, 40.84, 42.09, 
        993.7, 38.83, 37.63, 1002.91, 39.87, 42.35, 42.45, 38.84, 38.67, 43.27, 
        1044.11, 36.06, 39.4, 39.45, 44.11, 40.72, 39.6, 42.22, 39.52, 41.36, 
        1014.24, 37.6, 46.49, 39.71, 36.85, 38.83, 39.45, 40.96, 41.43, 38.67, 
        41.11, 39.88, 43.12, 41.47, 42.04, 37.92, 41.45, 38.88, 39.29, 35.44
      };
      
      // these are the estimates for initialization of Part A (gamma).
      GammaDist.MLEFit mle1, mle2;
      
      mle1 = new GammaDist.MLEFit(38.3052, 3.64493752, 1.35859522);
      mle2 = new GammaDist.MLEFit(310.6080, 4.62455229, 430.93016683);
      
      Distribution g1, g2;
      g1 = new GammaDist(mle1);//308, 0.13); //mle1);
      g2 = new GammaDist(mle2);//2500, 0.4); //mle2);
      
      double p1 = 0.5;
      
      // this is the intialization object for the EM algorithm
      EMInitThesis init = new EMInitThesis(p1, g1, g2);
      
      MixtureFit mf = new MixtureFit(points, 2, GammaDist.WeightedLazyFitFactory,
        init);
        
      mf.Solve(1, 1000, 1e-7, true);
    }

    if(false)
    {
      // this is the data from Part B of KellyThesis
      double points[] = {
        44.57, 55.27, 49.82, 1031.44, 44.69, 48.78, 49.86, 38.42, 42.87, 50.18, 
        1078.53, 48.66, 49.6, 46.43, 64.17, 38.98, 59.57, 56.12, 43.22, 61.83, 
        52.37, 65.61, 34.28, 46.81, 66.6, 38.85, 51.44, 46.52, 49.66, 48.26, 
        1049.41, 39.15, 41.17, 47.89, 1074.21, 39.19, 41.92, 58.37, 39.37, 49.42, 
        45.98, 55.95, 57.11, 58.32, 41.16, 49.13, 55.52, 51.26, 42.94, 69.36, 
        46, 62.72, 45.33, 38.69, 54.11, 1102.68, 46.51, 51.96, 46.06, 50.7, 
        50.12, 42.6, 56.96, 1002.17, 43.36, 42.42, 47.68, 1041.85, 1054.31, 68.4, 
        41.98, 32.64, 49.32, 1111.46, 42.7, 48.82, 62.89, 25.16, 46.79, 52.37, 
        30.97, 42.52, 44.55, 35.96, 57.39, 54.51, 1034.04, 38.17, 56.9, 49.13, 
        56.67, 49.23, 55.61, 57.94, 53.53, 52.05, 1048.8, 52.19, 45.19, 50.82
      };
      
      // these are the estimates for initialization of Part B (gamma).
      GammaDist.MLEFit mle1, mle2;
      
      mle1 = new GammaDist.MLEFit(43.2008, 3.75753905, 5.22430583);
      mle2 = new GammaDist.MLEFit(276.5026, 4.67256012, 414.87919156);
      
      Distribution g1, g2;
      g1 = new GammaDist(mle1);
      g2 = new GammaDist(mle2);
      
      double p1 = 0.5;
      
      // this is the intialization object for the EM algorithm
      EMInitThesis init = new EMInitThesis(p1, g1, g2);
      
      MixtureFit mf = new MixtureFit(points, 2, GammaDist.WeightedLazyFitFactory,
        init);
        
      mf.Solve(1, 1000, 1e-7, true);
    }
    
    if(true)
    {
      double points[] = {  19,  21,  22,  23,  26,  28,  29,  32 };
      double freqs[] =  { 0.1, 0.1, 0.1, 0.2, 0.1, 0.1, 0.2, 0.1 };
      
      FrequencyDist f = new FrequencyDist(points, freqs);
      FrequencyDist f2 = f.iConvolve(2);
      
      try {f2.saveToFile("convolve.txt"); } catch(IOException e) {};      
    }
  }
}
