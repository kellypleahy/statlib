/*
 * TestStats.java
 *
 * Created on September 24, 2002, 2:31 PM
 */

package statlib;

import java.io.*;
import java.util.*;

/**
 *
 * @author  KLeahy
 */
public class Thesis {

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
  
  public static void main(String[] args) {
    // this is the data from Part A of KellyThesis.
    double pointsA[] = {
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
    
    Arrays.sort(pointsA);
    
    Distribution.WeightedLazyFit wlfA1
      = LogNormalDist.WeightedLazyFitFactory.CreateInstance();
    Distribution.WeightedLazyFit wlfA2
      = LogNormalDist.WeightedLazyFitFactory.CreateInstance();
      
    wlfA1.Initialize();
    wlfA2.Initialize();
    
    for(int i=0; i<50; i++) {
      wlfA1.AddPoint(1, pointsA[i]);
      wlfA2.AddPoint(1, pointsA[i+50]);
    }
    
    Distribution d1A, d2A;
    d1A = wlfA1.Finalize();
    d2A = wlfA2.Finalize();
    
    double p1A = 0.5;
    
    // this is the intialization object for the EM algorithm
    EMInitThesis initA = new EMInitThesis(p1A, d1A, d2A);
    
    MixtureFit mfA = new MixtureFit(pointsA, 2, 
      LogNormalDist.WeightedLazyFitFactory, initA);
      
    double nllA = mfA.Solve(1, 1000, 1e-7, false);
    
    System.out.println("Part A fit with NLL = " + nllA);
    double[] pA = mfA.GetProportions();
    Distribution[] dlA = mfA.GetDistributions();
    
    for(int i=0; i<pA.length; i++) {
      System.out.println("  *** p: " + pA[i] + ", d: " + dlA[i]);
    }

    // this is the data from Part B of KellyThesis
    double pointsB[] = {
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
    
    Arrays.sort(pointsB);
    
    Distribution.WeightedLazyFit wlfB1
      = LogNormalDist.WeightedLazyFitFactory.CreateInstance();
    Distribution.WeightedLazyFit wlfB2
      = LogNormalDist.WeightedLazyFitFactory.CreateInstance();
      
    wlfB1.Initialize();
    wlfB2.Initialize();
    
    for(int i=0; i<50; i++) {
      wlfB1.AddPoint(1, pointsB[i]);
      wlfB2.AddPoint(1, pointsB[i+50]);
    }
    
    Distribution d1B, d2B;
    d1B = wlfB1.Finalize();
    d2B = wlfB2.Finalize();
    
    double p1B = 0.5;
    
    // this is the intialization object for the EM algorithm
    EMInitThesis initB = new EMInitThesis(p1B, d1B, d2B);
    
    MixtureFit mfB = new MixtureFit(pointsB, 2, 
      NormalDist.WeightedLazyFitFactory, initB);
      
    double nllB = mfB.Solve(1, 1000, 1e-7, false);
    
    System.out.println("Part B fit with NLL = " + nllB);
    double[] pB = mfB.GetProportions();
    Distribution[] dlB = mfB.GetDistributions();
    
    for(int i=0; i<pB.length; i++) {
      System.out.println("  *** p: " + pB[i] + ", d: " + dlB[i]);
    }
    
    // this is the data from Part C of KellyThesis.
    double pointsC[] = {
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
    
    Distribution.WeightedLazyFit wlfC 
      = GammaDist.WeightedLazyFitFactory.CreateInstance();
    
    wlfC.Initialize();
    for(int j=0; j<pointsC.length; j++) {
      wlfC.AddPoint(1, pointsC[j]);
    }
    
    Distribution dC = wlfC.Finalize();
    System.out.println("Part C distribution: " + dC.toString());
    
    // this is the data from Part D of KellyThesis.
    double pointsD[] = {
     144.99, 154.91, 155.65, 160.73, 147.14, 155.62, 160.86, 136.52, 144.89, 
     141.68, 145.86, 135.00, 143.86, 162.19, 153.74, 128.37, 145.36, 155.16, 
     146.22, 138.51, 141.33, 134.89, 136.30, 156.48, 144.59, 163.42, 171.11,  
     151.18, 143.77, 169.67, 131.35, 132.28, 165.20, 165.85, 150.78, 123.22, 
     171.79, 118.63, 191.85, 157.02, 127.77, 183.10, 148.33, 147.24, 158.75, 
     174.79, 143.45, 132.53, 156.57, 150.92, 150.76, 153.47, 162.40, 141.98, 
     170.43, 148.11, 136.64, 152.55, 153.98, 151.12, 165.80, 185.25, 138.88, 
     144.39, 174.64, 153.53, 131.53, 150.28, 152.60, 144.14, 130.72, 138.72, 
     132.73, 142.23, 145.57, 141.51, 136.06, 130.99, 149.20, 146.78, 151.26, 
     176.97, 118.60, 122.39, 139.39, 149.36, 139.44, 155.89, 125.63, 144.23, 
     145.78, 134.14, 140.03, 152.39, 137.64, 186.46, 149.79, 144.52, 159.40, 
     144.02
    };
    
    Distribution.WeightedLazyFit wlfD 
      = LogNormalDist.WeightedLazyFitFactory.CreateInstance();
    
    wlfD.Initialize();
    for(int j=0; j<pointsD.length; j++) {
      wlfD.AddPoint(1, pointsD[j]);
    }
    
    Distribution dD = wlfD.Finalize();
    System.out.println("Part D distribution: " + dD.toString());

    double alphaC = dC.getParameterValue(0), 
            betaC = dC.getParameterValue(1);
    CharacteristicFunction C = new CharacteristicFunction.Gamma(alphaC, betaC);

    double muD = dD.getParameterValue(0),
        sigmaD = Math.sqrt(dD.getParameterValue(1));
    CharacteristicFunction D = new CharacteristicFunction.Gaussian(muD, sigmaD);
  }
}
