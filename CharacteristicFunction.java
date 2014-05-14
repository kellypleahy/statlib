package statlib;

/**
 * Characteristic function interface (and implementing classes)
 * @author KLeahy
 * @version 1.0
 */
public interface CharacteristicFunction {
    /**
    * Compute the value of the characteristic function <i>&psi;</i>(<i>t</i>)
    * at different indices in the array.
    * @param real The destination for the real part of the complex function
    *               value.
    * @param imag The destination for the imaginary part of the complex function
    *               value.
    * @param first The first index (<i>k</i>) to evaluate (normally 0).
    * @param last The last index (<i>k</i>) to evaluate (normally
    *               <code>real.length - 1</code>).
    * @param A The value of <i>t</i> at the first index (<i>k</i> = 0).
    * @param B The value of <i>t</i> at the last index + 1 (<i>k</i> =
    *            <code>real.length</code>).
    **/
    public void compute(double[] real, double[] imag, int first, int last,
                        double A, double B);

    /** Implements the characteristic function for a mixture p.d.f. based on
     *    the characteristic functions and proportions of the component pdf's.
     */
    public static class Mixture implements CharacteristicFunction {

        private double[] m_pProportions;
        private CharacteristicFunction[] m_pCFs;

        /**
        * Constructor of the CFMixturePdf characteristic function object.
        * @param prop The proportion values for the
        *               <i>&psi;<sub>i</sub></i>(<i>t</i>).
        * @param psi The <i>&psi;<sub>i</sub></i>(<i>t</i>) (characteristic
        *              functions for the mixture).
        **/
        public Mixture(double[] prop, CharacteristicFunction[] psi) {
            assert (prop.length == psi.length && prop.length > 0);

            m_pProportions = (double[])prop.clone();
            m_pCFs = (CharacteristicFunction[])psi.clone();
        }

        /**
        * Compute the value of the characteristic function
        *   <i>&psi;</i>(<i>t</i>) at different indices in the array.
        * @param real The destination for the real part of the complex function
        *               value.
        * @param imag The destination for the imaginary part of the complex
        *               function value.
        * @param first The first index (<i>k</i>) to evaluate (normally 0).
        * @param last The last index (<i>k</i>) to evaluate (normally
        *               <code>real.length - 1</code>).
        * @param A The value of <i>t</i> at the first index (<i>k</i> = 0).
        * @param B The value of <i>t</i> at the last index + 1 (<i>k</i> =
        *            <code>real.length</code>).
        **/
        public void compute(double[] real, double[] imag, int first, int last,
                            double A, double B) {
            assert(real.length == imag.length && real.length > last);
            assert(first <= last && A <= B);

            double[] r1 = new double[last - first + 1];
            double[] r2 = new double[last - first + 1];
            double[] i1 = new double[last - first + 1];
            double[] i2 = new double[last - first + 1];

            double min, max, step;

            step = (B - A) / real.length;

            min = A + step * first;
            max = A + step * (last+1);

            // compute the first CF of the mixture.
            m_pCFs[0].compute(r1, i1, 0, last - first, min, max);

            // multiply by the proportion for this element of the mixture.
            for(int j=0; j<r1.length; j++) {
                r1[j] *= m_pProportions[0];
                i1[j] *= m_pProportions[0];
            }

            // compute the other CFs of the mixture.
            for(int i=1; i<m_pProportions.length; i++) {
                // compute the CF
                m_pCFs[i].compute(r2, i2, 0, last - first, min, max);

                // multiply by the applicable proportion and add into the
                //   destination
                for(int j=0; j<r2.length; j++) {
                    r1[j] += r2[j] * m_pProportions[i];
                    i1[j] += i2[j] * m_pProportions[i];
                }
            }

            // finally, move the destination values back to the original arrays
            for(int i=first, j=0; i<=last; i++,j++) {
                real[i] = r1[j];
                imag[i] = i1[j];
            }
        }
    }

    /** Implements the characteristic function for a gaussian kernel mixture
     *    p.d.f. based on the mean and variance supplied and the individual
     *    points supplied.
     */
    public static class GaussianKernel implements CharacteristicFunction {

        private double m_pPoints[];
        private int m_nPoints;
        private double m_dLogScalar;

        /**
         * Constructor, creates a GaussianKernel object with the characteristic
         *   function for the points supplied with the supplied sigma.  The c.f.
         *   represented corresponds to <br>
         *     <image src=gk_formula.bmp></image><br>
         *   where <i>N</i> is the number of points, the <i>x<sub>j</sub></i>
         *   are the points, and <i>&sigma;</i> is the dispersion parameter
         *   for the mixture.
         * @param pts The points at which to compute the mixture.
         * @param sigma The dispersion parameter (<i>&sigma;</i>) for the
         *                mixture components.
         */
        public GaussianKernel(double[] pts, double sigma) {
            m_pPoints = (double[])pts.clone();
            m_nPoints = m_pPoints.length;
            m_dLogScalar = Math.log(1.0/m_nPoints) + sigma * sigma / 2.0;
        }

        /**
        * Compute the value of the characteristic function
        *   <i>&psi;</i>(<i>t</i>) at different indices in the array.
        * @param real The destination for the real part of the complex function
        *               value.
        * @param imag The destination for the imaginary part of the complex
        *               function value.
        * @param first The first index (<i>k</i>) to evaluate (normally 0).
        * @param last The last index (<i>k</i>) to evaluate (normally
        *               <code>real.length - 1</code>).
        * @param A The value of <i>t</i> at the first index (<i>k</i> = 0).
        * @param B The value of <i>t</i> at the last index + 1 (<i>k</i> =
        *            <code>real.length</code>).
        **/
        public void compute(double[] real, double[] imag, int first, int last,
                            double A, double B) {
            assert(real.length == imag.length && real.length > last);
            assert(first <= last && A <= B);

            double step = (B - A) / real.length;
            double t = A + step * first;
            int N = real.length;

            for(int i=first; i<=last; i++) {
                double sreal = 0.0, simag = 0.0;
                for(int j=0; j<m_nPoints; j++) {
                    double xt = m_pPoints[j] * t;
                    sreal += Math.cos(xt);
                    simag += Math.sin(xt);
                }

                double scalar = Math.exp(m_dLogScalar * t * t);

                real[i] = sreal * scalar;
                imag[i] = simag * scalar;

                t += step;
            }
        }
    }

    /** Implements the characteristic function for a convolution of 1 or more
     *    characteristic functions.
     */
    public static class Convolution implements CharacteristicFunction {
        private CharacteristicFunction[] m_pCFs;
        private int m_nCFs;

        /**
         * Create a new Convolution object to convolve the supplied cf's.
         * @param comp the component characteristic function objects
         */
        public Convolution(CharacteristicFunction[] comp) {
            assert(comp.length > 1);
            m_pCFs = (CharacteristicFunction[])comp.clone();
            m_nCFs = comp.length;
        }

        /**
        * Compute the value of the characteristic function
        *   <i>&psi;</i>(<i>t</i>) at different indices in the array.
        * @param real The destination for the real part of the complex function
        *               value.
        * @param imag The destination for the imaginary part of the complex
        *               function value.
        * @param first The first index (<i>k</i>) to evaluate (normally 0).
        * @param last The last index (<i>k</i>) to evaluate (normally
        *               <code>real.length - 1</code>).
        * @param A The value of <i>t</i> at the first index (<i>k</i> = 0).
        * @param B The value of <i>t</i> at the last index + 1 (<i>k</i> =
        *            <code>real.length</code>).
        **/
        public void compute(double[] real, double[] imag, int first, int last,
                            double A, double B) {
            assert(real.length == imag.length && real.length > last);
            assert(first <= last && A <= B);

            double[] r1 = new double[last - first + 1];
            double[] r2 = new double[last - first + 1];
            double[] i1 = new double[last - first + 1];
            double[] i2 = new double[last - first + 1];

            double min, max, step;

            step = (B - A) / real.length;

            min = A + step * first;
            max = A + step * (last + 1);

            m_pCFs[0].compute(r1, i1, 0, last - first, min, max);

            // compute the other CFs of the convolution.
            for(int i=1; i<m_nCFs; i++) {
                // compute the CF
                m_pCFs[i].compute(r2, i2, 0, last - first, min, max);

                // multiply into the destination array
                Fourier.immult(r1, i1, r2, i2, r1, i1);
            }

            // finally, move the destination values back to the original arrays
            for(int i=first, j=0; i<=last; i++,j++) {
                real[i] = r1[j];
                imag[i] = i1[j];
            }
        }
    }

    /** Implements the characteristic function for the Gaussian distribution.
     *    The Gaussian cf is <i>&psi;</i>(<i>t</i>) =
     *       exp(<i>i&mu;t</i> &#43;
     *          <i>&sigma;</i><sup>2</sup><i>t</i><sup>2</sup>/2).
     */
    public static class Gaussian implements CharacteristicFunction {
        private double m_dMu;
        private double m_dSigma;

        /**
         * Create a new Gaussian object to represent the cf.
         * @param mu The mean of the Gaussian distribution
         * @param sigma The stddev of the Gaussian distribution
         */
        public Gaussian(double mu, double sigma) {
            assert(sigma > 0);
            m_dMu = mu;
            m_dSigma = sigma;
        }

        /**
        * Compute the value of the characteristic function
        *   <i>&psi;</i>(<i>t</i>) at different indices in the array.
        * @param real The destination for the real part of the complex function
        *               value.
        * @param imag The destination for the imaginary part of the complex
        *               function value.
        * @param first The first index (<i>k</i>) to evaluate (normally 0).
        * @param last The last index (<i>k</i>) to evaluate (normally
        *               <code>real.length - 1</code>).
        * @param A The value of <i>t</i> at the first index (<i>k</i> = 0).
        * @param B The value of <i>t</i> at the last index + 1 (<i>k</i> =
        *            <code>real.length</code>).
        **/
        public void compute(double[] real, double[] imag, int first, int last,
                            double A, double B) {
            assert(real.length == imag.length && real.length > last);
            assert(first <= last && A <= B);

            double step = (B - A) / real.length;
            double t = A + step * first;

            for(int i=first; i<=last; i++) {
                double scalar = Math.exp(-m_dSigma * m_dSigma * t * t / 2);
                double arg = m_dMu * t;
                real[i] = scalar * Math.cos(arg);
                imag[i] = scalar * Math.sin(arg);
                t += step;
            }
        }
    }

    /** Implements the characteristic function for the Gamma distribution.
     *    The Gamma cf is <i>&psi;</i>(<i>t</i>) =
     *       (1-<i>i&beta;t</i>)<sup><i>-&alpha;</i></sup>.
     */
    public static class Gamma implements CharacteristicFunction {
        private double m_dAlpha;
        private double m_dBeta;

        /**
         * Create a new Gamma object to represent the cf.
         * @param alpha The shape parameter of the Gamma distribution
         * @param beta The inverse scale parameter of the Gamma distribution
         */
        public Gamma(double alpha, double beta) {
            assert(beta > 0);
            m_dAlpha = alpha;
            m_dBeta = beta;
        }

        /**
        * Compute the value of the characteristic function
        *   <i>&psi;</i>(<i>t</i>) at different indices in the array.
        * @param real The destination for the real part of the complex function
        *               value.
        * @param imag The destination for the imaginary part of the complex
        *               function value.
        * @param first The first index (<i>k</i>) to evaluate (normally 0).
        * @param last The last index (<i>k</i>) to evaluate (normally
        *               <code>real.length - 1</code>).
        * @param A The value of <i>t</i> at the first index (<i>k</i> = 0).
        * @param B The value of <i>t</i> at the last index + 1 (<i>k</i> =
        *            <code>real.length</code>).
        **/
        public void compute(double[] real, double[] imag, int first, int last,
                            double A, double B) {
            assert(real.length == imag.length && real.length > last);
            assert(first <= last && A <= B);

            double step = (B - A) / real.length;
            double t = A + step * first;
            
            final double lscalar = m_dAlpha / 2;
            final double bsq = m_dBeta * m_dBeta;

            for(int i=first; i<=last; i++) {
                double larg = 1 + bsq * t * t;
                double scalar = Math.exp(lscalar * Math.log(larg));
                double arg = m_dAlpha * Math.atan(m_dBeta * t);
                real[i] = scalar * Math.cos(arg);
                imag[i] = scalar * Math.sin(arg);
                t += step;
            }
        }
    }
}