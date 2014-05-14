package statlib;

/**
 * <p>Title: JavaRT (Washington University)</p>
 * <p>Description: The bucket optimizer for FFT buckets in a freq dist.</p>
 * <p>Copyright: Copyright (c) 2003</p>
 * <p>Company: Washington University</p>
 * @author Kelly Leahy
 * @version 1.0
 */

public class BucketOptimize {
    private static final double zerotol = 0.0000001;

    /**
     * The original right-hand-side vector (b)
     */
    private double rhs[];

    /**
     * The A matrix...
     */
    private double constraints[][];

    /**
     * The current B^{-1} matrix for the problem.
     */
    private double currentBinv[][];

    /**
     * The set of variables representing each row.
     */
    private int basis[];

    /**
     * The most recently leaving basic variable.
     */
    private int lastLeaving;

    /**
     * Compute the dot product of two vectors.
     * @param v1 vector 1
     * @param v2 vector 2
     * @return the dot product.
     */
    private static double dotProduct(double[] v1, double[] v2) {
        double r = 0.0;
        for (int i = 0; i < v1.length; i++) {
            double v = v1[i] * v2[i];
            r += (Math.abs(v) > zerotol ? v : 0);
        }
        return r;
    }

    private static double[][] getIdentity(int n) {
        double[][] ret = new double[n][n];
        for (int i = 0; i < n; i++) {
            ret[i][i] = 1.0;
        }
        return ret;
    }

    private double[] getVarCoeff(int c) {
        double[] ret = new double[8];
        if (c < 8) {
            // this is a slack variable...
            for (int i = 0; i < 8; i++) {
                ret[i] = currentBinv[i][c];
            }
        }
        else {
            // this is a non-slack variable (x(i) or b)
            for (int i = 0; i < 8; i++) {
                // constraints is column, row
                ret[i] = dotProduct(currentBinv[i], constraints[c - 8]);
            }
        }
        return ret;
    }

    private double[] getRHS() {
        double[] ret = new double[8];
        for (int i = 0; i < 8; i++) {
            ret[i] = dotProduct(currentBinv[i], rhs);
        }
        return ret;
    }

    private String varName(int var) {
        if(var <= 7)
            return " [s(" + (var + 1) + ")]";
        else if(var <= 10)
            return " [x(" + (var - 7) + ")]";
        else if(var <= 15)
            return " [l(" + (var - 10) + ")]";
        else if(var == 16)
            return " [a]";
        else
            return " (what??? " + var + " ???)";
    }

    private boolean pivot(int entering) {
        if(entering < 0) return true;

        double[] coeff = getVarCoeff(entering);
        boolean negative = false;
        for(int i=0; i<8 && !negative; i++) {
            negative |= (coeff[i] > -zerotol);
        }
        if(!negative) {
            System.out.println("Entering variable column is non-negative!!!");
            return true;
        }

        double[] thisrhs = getRHS();

        int min = -1;
        double minv;
        if(entering != 16) {
            minv = Double.POSITIVE_INFINITY;
            for (int i = 0; i < 8; i++) {
                double ratio = thisrhs[i] / coeff[i];
                if (coeff[i] > zerotol && minv > ratio) {
                    min = i;
                    minv = ratio;
                }
            }
        } else {
            minv = Double.NEGATIVE_INFINITY;
            for(int i=0; i<8; i++) {
                double ratio = thisrhs[i] / coeff[i];
                if (coeff[i] < -zerotol && minv < ratio) {
                    min = i;
                    minv = ratio;
                }
            }
        }

        // if no max found, no more pivots left.
        if (min == -1) {
            System.out.println("No pivot point found!!! (entering = " + varName(entering) + ")");
            return true;
        }

        lastLeaving = basis[min];
        basis[min] = entering;

        System.out.println("Entering variable: " + entering + varName(entering));
        System.out.println("Leaving variable: " + lastLeaving + varName(lastLeaving));

        // E is row, column.
        double[][] E = getIdentity(8);
        for (int i = 0; i < 8; i++) {
            E[i][min] = -coeff[i] / coeff[min];
        }
        E[min][min] = 1 / coeff[min];

        // now we need nextBinv = E currentBinv.
        nextBinv(E);

        printTableau();

        return (lastLeaving == 16); // 16 is the artificial variable...
    }

    private void nextBinv(double[][] E) {
        double[][] next = new double[8][8];
        // currentBinv is row, column.  E is row, column.
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                for (int k = 0; k < 8; k++) {
                    double v = E[i][k] * currentBinv[k][j];
                    next[i][j] += (Math.abs(v) > zerotol ? v : 0);
                }
            }
        }
        currentBinv = next;
    }

    private void initProblem(double[] x, int[] c,
                             double alpha, double beta,
                             double gamma, double delta,
                             double h0, double h1) {

        // set up the initial Binv = I.
        currentBinv = getIdentity(8);

        // compute sigma values
        int N = x.length;
        double sigma2 = 0.0, sigma3 = 0.0, sigma4 = 0.0, sigma5 = 0.0;
        for (int i = 0; i < N; i++) {
            sigma2 += x[i];
            sigma3 += c[i] * x[i];
            sigma4 += c[i];
            sigma5 += c[i] * c[i];
        }

        // set up the RHS vector.
        rhs = new double[8];
        rhs[0] = rhs[1] = sigma2;
        rhs[2] = sigma3;
        rhs[3] = -alpha;
        rhs[4] = gamma;
        rhs[5] = -h0;
        rhs[6] = h1;

        // set up the A matrix (with the artifical var at the end)
        //  constraints is (column, row) matrix.
        constraints = new double[9][8];
        constraints[0][0] = constraints[1][0] = constraints[0][1]
            = constraints[1][1] = N;
        constraints[0][2] = constraints[1][2] = constraints[2][1]
            = constraints[2][0] = sigma4;
        constraints[0][3] = constraints[3][0] = constraints[2][5]
            = constraints[5][2] = constraints[2][7] = constraints[7][2]
            = -1;
        constraints[0][4] = constraints[4][0] = constraints[1][7]
            = constraints[7][1] = constraints[2][6] = constraints[6][2] = 1;
        constraints[2][2] = sigma5;
        constraints[2][3] = constraints[3][2] = beta;
        constraints[2][4] = constraints[4][2] = -delta;
        for(int i=0; i<8; i++)
            if(rhs[i] < 0)
                constraints[8][i] = -1;
            else
                constraints[8][i] = 0;

        // set up the basis vector
        basis = new int[8];
        basis[0] = 0;
        basis[1] = 1;
        basis[2] = 2;
        basis[3] = 3;
        basis[4] = 4;
        basis[5] = 5;
        basis[6] = 6;
        basis[7] = 7;

        printTableau();

        // pivot a into the basis
        pivot(16); // 16 = the artificial variable...
    }

    private void printTableau() {
        // tableau is { Binv A | Binv | Binv b }
        System.out.println("------------- Current Tableau ----------------");
        double[][] tab = new double[18][];
        for(int c=0; c<17; c++) tab[c] = getVarCoeff(c);
        tab[17] = getRHS();
        for(int i=1; i<=8; i++) System.out.print("s(" + i + ")\t");
        for(int i=1; i<=3; i++) System.out.print("x(" + i + ")\t");
        for(int i=1; i<=5; i++) System.out.print("l(" + i + ")\t");
        System.out.println("a\tRHS");
        for(int i=0; i<8; i++) {
            for (int j = 0; j < 18; j++) {
                System.out.print((Math.round(tab[j][i] * 1000) / 1000.0) + "\t");
            }
            System.out.println();
        }
        System.out.print("Basis: \t");
        for(int i=0; i<8; i++) System.out.print(basis[i] + varName(basis[i]) + "\t");
        System.out.println();
        System.out.println("--------------- End Tableau ------------------");
    }

    public BucketOptimize(double[] points, int[] buckets, double alpha,
                          double beta, double gamma, double delta,
                          double h0, double h1) {
        System.out.println("alpha = " + alpha + ", beta = " + beta
                           + ", gamma = " + gamma + ", delta = " + delta
                           + ", h0 = " + h0 + ", h1 = " + h1);
        initProblem(points, buckets, alpha, beta, gamma, delta, h0, h1);
    }

    public boolean Solve() {
        while(!pivot((lastLeaving + 8) % 16)) ;
        return true;
    }
}