/*
 * fourier.java
 *
 * Created on March 12, 2002, 9:23 PM
 */

package statlib;

import java.lang.IllegalArgumentException;

/**
 *
 * @author KLeahy
 * @version 1.0
 */
public class Fourier {

    private static final int numBits(int x)
    {
        return (int)(Math.log(x) / Math.log(2.0));
    }

    private static final int reverseBits(int v, int n)
    {
        int rev = 0, index = v;
        for(int i=0; i<n; i++) {
            rev = (rev << 1) | (index & 0x1);
            index >>= 1;
        }
        return rev;
    }

    private static final void checkArgs(double[] realIn, double[] imagIn,
        double[] realOut, double[] imagOut) throws IllegalArgumentException
    {
        if(realIn == null || realOut == null || imagIn == null
                || imagOut == null)
            throw new IllegalArgumentException("Arrays must not be null");
        if(realIn.length != realOut.length || realIn.length != imagIn.length
                || realIn.length != imagOut.length)
            throw new IllegalArgumentException("Array lengths must all match");
        if(realIn.length != (1<<(int)(Math.log(realIn.length)/Math.log(2.0))))
            throw new IllegalArgumentException(
                "Array length must be power of 2");
    }

    public static final void FFT(double[] realIn, double[] imagIn,
        double[] realOut, double[] imagOut) throws IllegalArgumentException
    {
        checkArgs(realIn, imagIn, realOut, imagOut);
        internalFFT(2.0 * Math.PI, realIn, imagIn, realOut, imagOut);
    }

    public static final void iFFT(double[] realIn, double[] imagIn,
        double[] realOut, double[] imagOut) throws IllegalArgumentException
    {
        checkArgs(realIn, imagIn, realOut, imagOut);
        internalFFT(-2.0 * Math.PI, realIn, imagIn, realOut, imagOut);
        int N = realIn.length;
        for(int i=0; i<N; i++)
        {
            realOut[i] /= N;
            imagOut[i] /= N;
        }
    }

    private static final void internalFFT(double numer, double[] realIn,
        double[] imagIn, double[] realOut, double[] imagOut)
    {
        int nBits = numBits(realIn.length);
        int N = realIn.length;
        for(int i=0; i<N; i++) {
            int j = reverseBits(i, nBits);
            realOut[j] = realIn[i];
            imagOut[j] = imagIn[i];
        }

        int blkEnd = 1, blkSize = 2;
        while(blkSize <= N) {
            double delta = numer / blkSize;
            double alpha = Math.sin(0.5 * delta);
            alpha = 2.0 * alpha * alpha;
            double beta = Math.sin(delta);

            for(int i=0; i<N; i+=blkSize) {
                double ar = 1.0, ai = 0.0;
                for(int m = 0, j=i; m<blkEnd; m++, j++) {
                    int k = j + blkEnd;
                    double tr = ar * realOut[k] - ai * imagOut[k],
                        ti = ar * imagOut[k] + ai * realOut[k];
                    realOut[k] = realOut[j] - tr;
                    imagOut[k] = imagOut[j] - ti;
                    realOut[j] += tr;
                    imagOut[j] += ti;
                    double dar = alpha * ar + beta * ai;
                    ai -= alpha * ai - beta * ar;
                    ar -= dar;
                }
            }
            blkEnd = blkSize;
            blkSize <<= 1;
        }
    }

    public static double immod(double real, double imag) {
      return Math.sqrt(real * real + imag * imag);
    }

    public static double imang(double real, double imag) {
        if(real == 0) {
            if(imag >= 0) {
                return Math.PI * 0.5;
            } else {
                return Math.PI * 1.5;
            }
        } else {
            return Math.atan2(imag, real);
        }
    }

    public static void immult(double[] real1, double[] imag1,
			                  double[] real2, double[] imag2,
			                  double[] realout, double[] imagout) {
        int N = real1.length;
	    assert(N == imag1.length && N == real2.length && N == imag2.length);
        /* old way - using r, theta
	    double r, theta;
	    for(int i=0; i<N; i++) {
	        r = immod(real1[i], imag1[i]) * immod(real2[i], imag2[i]);
	        theta = imang(real1[i], imag1[i]) + imang(real2[i], imag2[i]);
	        realout[i] = r * Math.cos(theta);
	        imagout[i] = r * Math.sin(theta);
	    }
	    */

	    // new way, using R = r1r2 - i1i2, I = r1i2 + i1r2
	    for(int i=0; i<N; i++) {
	        // save inputs in case overlaps of input arrays exist.
	        double r1 = real1[i], r2 = real2[i];
	        double i1 = imag1[i], i2 = imag2[i];
	        realout[i] = r1 * r2 - i1 * i2;
	        imagout[i] = r1 * i2 + i1 * r2;
	    }
    }

    public static void impower(double[] real, double[] imag, int power) {
        int N = real.length;
        assert (N == imag.length);
        for(int i=0; i<N; i++) {
            double r, theta;
            r = Math.pow(immod(real[i], imag[i]), power);
            theta = power * imang(real[i], imag[i]);
            real[i] = r * Math.cos(theta);
            imag[i] = r * Math.sin(theta);
        }
    }


    /**
     * Compute the DFT that corresponds to the analytical cf represented by the
     *   argument f.  Use the domain [A,B) for corresponding pdf with evaluation
     *   at the points <i>A</i>, <i>A</i> + (<i>B</i> - <i>A</i>)/<i>N</i>,
     *   <i>A</i> + 2(<i>B</i> - <i>A</i>)/<i>N</i>, ..., <i>A</i> +
     *   (<i>N</i> - 1)(<i>B</i> - <i>A</i>)/<i>N</i>.
     * @param f The <code>CharacteristicFunction</code> object representing the
     *            distribution's <i>&psi;</i>(<i>t</i>).
     * @param A The minimum value in the domain of the target pdf.
     * @param B The maximum value in the domain of the target pdf.
     * @param N The number of buckets in the domain of the target pdf.
     * @param normalize Normalize the output if true (make it sum to 1).
     */
    public static void CharFnToDFT(CharacteristicFunction f, double A, double B,
                                   int N, double[] real, double[] imag,
                                   boolean normalize) {
        assert (real.length == imag.length && imag.length == N);
        assert (N > 0 && N % 2 == 0);
        assert (B > A);

        double scalar_const = N / (B - A);

        double max = Math.PI * scalar_const;
        double min = -max;

        double[] rscalar = new double[N], iscalar = new double[N];

        double arg_scalar = -2 * Math.PI * A / (B - A);

        if(!normalize) scalar_const = 1;

        double[] tmpreal = new double[N], tmpimag = new double[N];

        // compute the CF at the points k = -N/2 to k = N/2 - 1.
        f.compute(tmpreal, tmpimag, 0, N-1, min, max);

        // compute the scalars at the points k = -N/2 to k = N/2 - 1.
        for(int j=0,k=-(N>>1); j < N; k++,j++) {
            double arg = arg_scalar * k;
            rscalar[j] = scalar_const * Math.cos(arg);
            iscalar[j] = scalar_const * Math.sin(arg);
        }

        // multiply the two arrays.
        Fourier.immult(tmpreal, tmpimag, rscalar, iscalar, tmpreal, tmpimag);

        // now, shift the arrays into the right order (iFFT expects the array
        //   to have k = 0, 1, ..., N/2-1, -N/2, -N/2 + 1, ..., -1 in the
        //   array (in order).
        for(int j=0, k=(N>>1); k<N; j++,k++) {
            real[k] = tmpreal[j];
            imag[k] = tmpimag[j];
            real[j] = tmpreal[k];
            imag[j] = tmpimag[k];
        }
    }

    /**
     * Compute the DFT that corresponds to the analytical cf represented by the
     *   argument f.  Use the domain [A,B) for corresponding pdf with evaluation
     *   at the points <i>A</i>, <i>A</i> + (<i>B</i> - <i>A</i>)/<i>N</i>,
     *   <i>A</i> + 2(<i>B</i> - <i>A</i>)/<i>N</i>, ..., <i>A</i> +
     *   (<i>N</i> - 1)(<i>B</i> - <i>A</i>)/<i>N</i>.
     * @param f The <code>CharacteristicFunction</code> object representing the
     *            distribution's <i>&psi;</i>(<i>t</i>).
     * @param A The minimum value in the domain of the target pdf.
     * @param B The maximum value in the domain of the target pdf.
     * @param N The number of buckets in the domain of the target pdf.
     */
    public static void CharFnToDFT(CharacteristicFunction f, double A, double B,
                                   int N, double[] real, double[] imag) {
        CharFnToDFT(f, A, B, N, real, imag, true);
    }
}

