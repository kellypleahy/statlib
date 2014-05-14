/*
 * SpecialFunc.java
 *
 * Created on March 28, 2003, 3:03 PM
 */

package statlib;

/**
 *
 * @author  Kelly Leahy
 * @version 1.0
 */
public class SpecialFunc {
  private final static double _Epsilon       =  2.2204460492503131e-016;
  private final static double _LogMin        = -7.0839641853226408e+002;
  private final static double _LogMax        =  7.0978271289338397e+002;
  private final static double _SqrtMin       =  1.4916681462400413e-154;
  private final static double _SqrtMax       =  1.3407807929942596e+154;
  private final static double _PosZero       = +0.0000000000000000e+000;
  
  private static class ChebyshevPoly {
    private final double Series[];   // Chebyshev series polynomials
    private final int Order;         // order of expansion
    private final double a, b;       // interval
    
    private final double _h2;        // 2 / (b-a)
    private final double _bpa;       // (b+a)
  
    public ChebyshevPoly(double[] series, int order, double a, double b) 
    {
      this.Series = (double[])series.clone();
      this.Order = order;
      this.a = a;
      this.b = b;
      
      this._h2 = (b - a) / 2;
      this._bpa = (b + a);
    }
    
    private double evaluate(double x) {
      double d = 0.0, dd = 0.0;
      
      final double y = (2.0 * x - _bpa) * _h2;
      
      for(int i=Order; i>0; i--) {
        double t = d;
        d = (y * d) - dd + Series[i];
        dd = t;
      }
      
      return -dd + 0.5 * (y * d + Series[0]);
    }
  }
  
  private final static ChebyshevPoly ChebPsi1
    = new ChebyshevPoly(
       new double[]{
        -0.038057080835217922,
         0.491415393029387130,
        -0.056815747821244730,
         0.008357821225914313,
        -0.001333232857994342,
         0.000220313287069308,
        -0.000037040238178456,
         0.000006283793654854,
        -0.000001071263908506,
         0.000000183128394654,
        -0.000000031353509361,
         0.000000005372808776,
        -0.000000000921168141,
         0.000000000157981265,
        -0.000000000027098646,
         0.000000000004648722,
        -0.000000000000797527,
         0.000000000000136827,
        -0.000000000000023475,
         0.000000000000004027,
        -0.000000000000000691,
         0.000000000000000118,
        -0.000000000000000020
      },
      22,
      -1.0, 1.0
     );

  private final static ChebyshevPoly ChebPsi2
    = new ChebyshevPoly(
       new double[]{
        -0.0204749044678185,
        -0.0101801271534859,
         0.0000559718725387,
        -0.0000012917176570,
         0.0000000572858606,
        -0.0000000038213539,
         0.0000000003397434,
        -0.0000000000374838,
         0.0000000000048990,
        -0.0000000000007344,
         0.0000000000001233,
        -0.0000000000000228,
         0.0000000000000045,
        -0.0000000000000009,
         0.0000000000000002,
        -0.0000000000000000
      },
      15,
      -1.0, 1.0
     );
  
  /**
   * Calculate <code>Y * exp(X)</code> while gracefully handling overflow and
   *   underflow.
   */
  
  public static double YExpX(double x, double y) {
    double ay, r;
    
    ay = Math.abs(y);
    if(y == 0.0) 
      r = 0.0;
    else if((x < 0.5 * _LogMax) && (x > 0.5 * _LogMin)
      && ((ay < 0.8 * _SqrtMax) && (ay > 1.2 * _SqrtMin)))
      r = y * Math.exp(x);
    else {
      double ly, lnr;
      ly = Math.log(y);
      lnr = x + ly;
      
      // handle overflow / underflow...
      if(lnr > _LogMax - 0.01)
        r = Double.POSITIVE_INFINITY;
      else if(lnr < _LogMin + 0.01)
        r = _PosZero;
      else {
        double sy, M, N, a, b;
        sy = y > 0 ? 1 : (y == 0) ? 0 : -1;
        M = Math.floor(x);
        N = Math.floor(ly);
        if(sy == 0) 
          r = 0.0;
        else
          r = sy * Math.exp(M+N) * Math.exp((x-M)+(ly-N));
      }
    }
    return r;
  }
  
  /**
   * Calculate <i>&psi;(x)</i> (the digamma function).
   */
  public static double Psi(double x) {
    double r = 0.0;
    
    double y, s;
    
    y = Math.abs(x);
    
    if(x == 0.0 || x == -1.0 || x == -2.0)
      r = Double.NaN;
    else if (y >= 2.0) {
      r = ChebPsi2.evaluate(8.0 / (y * y) - 1.0);
      if(x < 0.0) {
        s = Math.sin(Math.PI * x);
        
        if(Math.abs(s) < 2.0 * _SqrtMin)
          r = Double.NaN;
        else
          r = Math.log(y) - 0.5 / x + r - Math.PI * Math.cos(Math.PI * x) / s;
      } else r = Math.log(y) - 0.5 / x + r;
    } else {
      final double K, arg;
      if(x < -1.0) {
        K = -(1.0 / x + 1.0 / (x + 1.0) + 1.0 / (x + 2.0));
        arg = 2.0 * (x + 2.0) - 1.0;
      } else if(x < 0.0) {
        K = -(1.0 / x + 1.0 / (x + 1.0));
        arg = 2.0 * (x + 1.0) - 1.0;
      } else if(x < 1.0) {
        K = -(1.0 / x);
        arg = 2.0 * x - 1.0;
      } else {
        K = 0.0;
        arg = 2.0 * (x - 1.0) - 1.0;
      }
      r = K + ChebPsi1.evaluate(arg);
    }
    
    return r;
  }
  
  private static final double LNFACT[] = { //{{{
    0.00000000000000e+000,
    0.00000000000000e+000,
    6.93147180559945e-001,
    1.79175946922805e+000,
    3.17805383034795e+000,
    4.78749174278205e+000,
    6.57925121201010e+000,
    8.52516136106541e+000,
    1.06046029027453e+001,
    1.28018274800815e+001,
    1.51044125730755e+001,
    1.75023078458739e+001,
    1.99872144956619e+001,
    2.25521638531234e+001,
    2.51912211827387e+001,
    2.78992713838409e+001,
    3.06718601060807e+001,
    3.35050734501369e+001,
    3.63954452080331e+001,
    3.93398841871995e+001,
    4.23356164607535e+001,
    4.53801388984769e+001,
    4.84711813518352e+001,
    5.16066755677644e+001,
    5.47847293981123e+001,
    5.80036052229805e+001,
    6.12617017610020e+001,
    6.45575386270063e+001,
    6.78897431371815e+001,
    7.12570389671680e+001,
    7.46582363488302e+001,
    7.80922235533153e+001,
    8.15579594561150e+001,
    8.50544670175815e+001,
    8.85808275421977e+001,
    9.21361756036871e+001,
    9.57196945421432e+001,
    9.93306124547874e+001,
    1.02968198614514e+002,
    1.06631760260643e+002,
    1.10320639714757e+002,
    1.14034211781462e+002,
    1.17771881399745e+002,
    1.21533081515439e+002,
    1.25317271149357e+002,
    1.29123933639127e+002,
    1.32952575035616e+002,
    1.36802722637326e+002,
    1.40673923648234e+002,
    1.44565743946345e+002,
    1.48477766951773e+002,
    1.52409592584497e+002,
    1.56360836303079e+002,
    1.60331128216631e+002,
    1.64320112263195e+002,
    1.68327445448428e+002,
    1.72352797139163e+002,
    1.76395848406997e+002,
    1.80456291417544e+002,
    1.84533828861449e+002,
    1.88628173423672e+002,
    1.92739047287845e+002,
    1.96866181672890e+002,
    2.01009316399282e+002,
    2.05168199482641e+002,
    2.09342586752537e+002,
    2.13532241494563e+002,
    2.17736934113954e+002,
    2.21956441819130e+002,
    2.26190548323728e+002,
    2.30439043565777e+002,
    2.34701723442818e+002,
    2.38978389561834e+002,
    2.43268849002983e+002,
    2.47572914096187e+002,
    2.51890402209723e+002,
    2.56221135550010e+002,
    2.60564940971863e+002,
    2.64921649798553e+002,
    2.69291097651020e+002,
    2.73673124285694e+002,
    2.78067573440366e+002,
    2.82474292687630e+002,
    2.86893133295427e+002,
    2.91323950094270e+002,
    2.95766601350761e+002,
    3.00220948647014e+002,
    3.04686856765669e+002,
    3.09164193580147e+002,
    3.13652829949879e+002,
    3.18152639620209e+002,
    3.22663499126726e+002,
    3.27185287703775e+002,
    3.31717887196928e+002,
    3.36261181979198e+002,
    3.40815058870799e+002,
    3.45379407062267e+002,
    3.49954118040770e+002,
    3.54539085519441e+002,
    3.59134205369575e+002,
    3.63739375555563e+002,
    3.68354496072405e+002,
    3.72979468885689e+002,
    3.77614197873919e+002,
    3.82258588773060e+002,
    3.86912549123218e+002,
    3.91575988217330e+002,
    3.96248817051792e+002,
    4.00930948278916e+002,
    4.05622296161145e+002,
    4.10322776526937e+002,
    4.15032306728250e+002,
    4.19750805599545e+002,
    4.24478193418257e+002,
    4.29214391866652e+002,
    4.33959323995015e+002,
    4.38712914186121e+002,
    4.43475088120919e+002,
    4.48245772745385e+002,
    4.53024896238496e+002,
    4.57812387981278e+002,
    4.62608178526875e+002,
    4.67412199571608e+002,
    4.72224383926981e+002,
    4.77044665492586e+002,
    4.81872979229888e+002,
    4.86709261136839e+002,
    4.91553448223298e+002,
    4.96405478487218e+002,
    5.01265290891579e+002,
    5.06132825342035e+002,
    5.11008022665236e+002,
    5.15890824587822e+002,
    5.20781173716044e+002,
    5.25679013515995e+002,
    5.30584288294433e+002,
    5.35496943180170e+002,
    5.40416924105998e+002,
    5.45344177791155e+002,
    5.50278651724286e+002,
    5.55220294146895e+002,
    5.60169054037273e+002,
    5.65124881094874e+002,
    5.70087725725134e+002,
    5.75057539024710e+002,
    5.80034272767131e+002,
    5.85017879388839e+002,
    5.90008311975618e+002,
    5.95005524249382e+002,
    6.00009470555327e+002,
    6.05020105849424e+002,
    6.10037385686239e+002,
    6.15061266207085e+002,
    6.20091704128477e+002,
    6.25128656730891e+002,
    6.30172081847810e+002,
    6.35221937855060e+002,
    6.40278183660408e+002,
    6.45340778693435e+002,
    6.50409682895655e+002,
    6.55484856710889e+002,
    6.60566261075874e+002,
    6.65653857411106e+002,
    6.70747607611913e+002,
    6.75847474039737e+002,
    6.80953419513637e+002,
    6.86065407301994e+002,
    6.91183401114411e+002,
    6.96307365093814e+002,
    7.01437263808737e+002,
    7.06573062245787e+002    
  }; //}}}
  private static final int MAX_LNFACT = LNFACT.length;
  
  /**
   * Calculate the logarithm of the factorial (log(<i>N!</i>)).
   */
  public static double lnFact(int n) {
    if(n < MAX_LNFACT)
      return LNFACT[n];
    else
      return lnGamma(n + 1.0);
  }
  
  private final static double Shz_MAX_TAIL_BITS = 54.0;
  private final static int Shz_MAX_J = 14;
  private final static int Shz_MAX_K = 10;
  
  private final static double Shz_COEFF[] = { //{{{
     1.000000000000000e+000,
     8.333333333333333e-002,
    -1.388888888888889e-003,
     3.306878306878307e-005,
    -8.267195767195767e-007,
     2.087675698786810e-008,
    -5.284190138687493e-010,
     1.338253653068468e-011,
    -3.389680296322583e-013,
     8.586062056277845e-015,
    -2.174868698558062e-016,
     5.509002828360230e-018,
    -1.395446468581252e-019,
     3.534707039629467e-021,
    -8.953517427037547e-023
  }; //}}}

  /**
   * Calculate the value of the shifted Zeta function (shifted by q).
   */
  public static double ShiftedZeta(double s, double q) {
    double r;
    
    if(s <= 1.0 || q <= 0.0)
      r = Double.NaN;
    else {
      final double ln_term0 = -s * Math.log(q);
      
      if(ln_term0 < _LogMin + 1.0)
        r = _PosZero;
      else if(ln_term0 > _LogMax - 1.0)
        r = Double.POSITIVE_INFINITY;
      else if((s > Shz_MAX_TAIL_BITS && q < 1.0)
              || (s > 0.5 * Shz_MAX_TAIL_BITS && q < 0.25))
        r = Math.pow(q, -s);
      else if(s > 0.5 * Shz_MAX_TAIL_BITS && q < 1.0)
        r = Math.pow(q, -s) 
          * (1.0 + Math.pow(q / (1.0 + q), s) + Math.pow(q / (2.0 + q), s));
      else {
        // Euler-Maclaurin summation formula (Moshier p. 400 + corrections)
        final double maxkpq = (Shz_MAX_K + q);
        final double deninv = 1.0 / (maxkpq * maxkpq);
        final double pmax = Math.pow(maxkpq, -s);
        double scp = s;
        double pcp = pmax / maxkpq;
        double ans = pmax * (maxkpq / (s - 1.0) + 0.5);
        
        for(int k=0; k<Shz_MAX_K; k++) 
          ans += Math.pow(k + q, -s);
        
        for(int j=0; j<Shz_MAX_J; j++) {
          double delta = Shz_COEFF[j+1] * scp * pcp;
          ans += delta;
          if(Math.abs(delta / ans) < 0.5 * _Epsilon) break;
          
          double t = s + ((j<<1) + 1);
          scp *= t * (t + 1.0);
          pcp *= deninv;
        }
        r = ans;        
      }
    }
    
    return r;
  }
  
  public static double PolyGamma(int n, double x) {
    double r;
    
    if(x <= 0.0)
      r = Double.NaN;
    else if(n == 0)
      r = Psi(x);
    else {
      r = YExpX(lnFact(n), ShiftedZeta(n + 1.0, x));
      if(n % 2 == 0) r = -r;
    }
    return r;
  }
  
  /**
   * Calculate the logarithm of the Gamma function with parameter a (a > 0)
   **/
  public static double lnGamma(double a) {
    if(a > 10) {
      final double l2pi = Math.log(2 * Math.PI);
      final double[] p = {12, 360, 1260, 1680, 1188, 360360/691, 156, 
        122400/3617, 244188/43867, 125400/174600};
      double rv = (a - 0.5) * Math.log(a) - a + 0.5 * l2pi;
      double na2 = - a * a;
      double av = a;
      for(int i=0; i<9; i++, av *= na2)
        rv += 1 / (p[i] * av);
      return rv + p[9] / av;
    } else if(a <= 0) {
      return 0;
    } else {
      double av = a;
      double rv = 0;
      while(av <= 10) {
        rv -= Math.log(av);
        av += 1;
      }
      return rv + lnGamma(av);
    }
  }

  public static double gamma(double alpha) {
    return YExpX(lnGamma(alpha), 1.0);
  }

  public static double incGammaC(double alpha, double x) {
    final double tol = 1E-11;
    if(x <= 0 || alpha <= 0)
      return 0.0;
    if(x < alpha - 1)
      return 1.0 - incGamma(alpha, x);
    double lnGa = lnGamma(alpha);
    if(alpha != 1.0) {
      if(x >= (-710.75 + lnGa) / (alpha - 1))
        return 0.0;
    } else {
      if(x >= Math.exp(709.75))
        return 0.0;
    }
    double ax = alpha * Math.log(x) - x - lnGa;
    if(ax < -709.75)
      return 0.0;
    ax = Math.exp(ax);
    
    // calculate continued fraction expansion 
    //  (see Abramowitz & Stegun pg 263 formula 6.5.31)
    //  continued fraction method from Abramowitz & Stegun pg 19 (3.10)
    double Ak, Akm1 = ax, Akm2 = 0.0, Bk, Bkm1 = x, Bkm2 = 1;
    int i=0;
    double t = 1;
    double ans = Akm1 / Bkm1, lans;
    do {
      double a, b;
      i = (i + 1) % 2;
      if(i == 1) {
        b = 1;
        a = t - alpha;
      } else {
        b = x;
        a = t;
        t += 1;
      }
      Ak = b * Akm1 + a * Akm2;
      Akm2 = Akm1;
      Akm1 = Ak;
      Bk = b * Bkm1 + a * Bkm2;
      Bkm2 = Bkm1;
      Bkm1 = Bk;
      lans = ans;
      ans = Ak / Bk;
    } while(Math.abs(ans - lans) > tol);
    return ans;
  }
  
  public static double incGamma(double alpha, double x) {
    if(x <= 0 || alpha <= 0)
      return 0.0;
    if(x > 1.0 && x > alpha)
      return 1.0 - incGammaC(alpha, x);
    double ax = alpha * Math.log(x) - lnGamma(alpha);
    if(ax < -709.75)
      return 0.0;
    ax = Math.exp(ax);
    
    final double tol = 1E-11;
    
    double d = alpha;
    double t = 1.0;
    double rv = 1.0;
    do {
      d += 1;
      t *= x / d;
      rv += t;
    } while (t / rv > tol);
    
    return rv * ax / alpha;
  }
  
}
