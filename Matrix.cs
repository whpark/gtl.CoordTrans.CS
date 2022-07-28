//=============================================================================
// 
// xMatrix - double
// PWH. 2021.08.06. ==> C#...
//
//

using System;
using System.Collections.Generic;
using System.Text;

namespace gtl.CoordTrans {
    public class xMatrix {
        public double[,] m_mat;
        public int rows => m_mat.GetLength(0);
        public int cols => m_mat.GetLength(1);

        /// <summary>
        /// No Initialize... ?
        /// </summary>
        /// <param name="rows"></param>
        /// <param name="cols"></param>
        public xMatrix(int rows, int cols) {
            m_mat = new double[rows, cols];
        }
        public xMatrix(int rows, int cols, double value) {
            m_mat = new double[rows, cols];
            for (int y = 0; y < rows; y++) {
                for (int x = 0; x < cols; x++)
                    m_mat[y, x] = value;
            }
        }

        public xMatrix(xMatrix m) {
            m_mat = new double[m.rows, m.cols];
            ForAll((y, x) => m.m_mat[y, x]);
        }

        public static xMatrix eye(int nSize) {
            xMatrix M = new xMatrix(nSize, nSize);
            for (int i = 0; i < nSize; i++)
                M.m_mat[i, i] = 1.0;
            return M;
        }

        public ref double At(int row, int col) {
            return ref m_mat[row, col];
        }

        public xMatrix SetAll(double a) {
            for (int y = 0; y < rows; y++) {
                for (int x = 0; x < cols; x++)
                    At(y, x) = a;
            }
            return this;
        }

        public xMatrix ForAll(Func<int, int, double> func) {
            for (int y = 0; y < rows; y++) {
                for (int x = 0; x < cols; x++) {
                    m_mat[y, x] = func(y, x);
                }
            }
            return this;
        }


        // 4칙연산 with double
        public static xMatrix operator +(xMatrix a) => a;
        public static xMatrix operator -(xMatrix a) {
            xMatrix b = new xMatrix(a);
            b.ForAll((y, x) => -a.m_mat[y, x]);
            return b;
        }
        public static xMatrix operator +(xMatrix a, xMatrix b) {
            if ((a.rows != b.rows) || (a.cols != b.cols))
                throw new Exception("Size mismatching");

            xMatrix c = new xMatrix(a.rows, a.cols);
            c.ForAll((y, x) => (a.m_mat[y, x] + b.m_mat[y, x]));
            return c;
        }
        public static xMatrix operator -(xMatrix a, xMatrix b) {
            if ((a.rows != b.rows) || (a.cols != b.cols))
                throw new Exception("Size mismatching");

            xMatrix c = new xMatrix(a.rows, a.cols);
            c.ForAll((y, x) => (a.m_mat[y, x] - b.m_mat[y, x]));
            return c;
        }

        public static xMatrix operator *(xMatrix a, double b) {
            return new xMatrix(a.rows, a.cols).ForAll((y, x) => a.At(y, x) * b);
        }
        public static xMatrix operator *(double b, xMatrix a) {
            return new xMatrix(a.rows, a.cols).ForAll((y, x) => a.At(y, x) * b);
        }
        public static xMatrix operator /(double b, xMatrix a) {
            return new xMatrix(a.rows, a.cols).ForAll((y, x) => b / a.At(y, x));
        }
        public static xMatrix operator /(xMatrix a, double b) {
            return new xMatrix(a.rows, a.cols).ForAll((y, x) => a.At(y, x) / b);
        }

        // xMatrix Multiplication
        public static xMatrix operator *(xMatrix a, xMatrix b) {
            if (a.cols != b.rows)
                throw new Exception("size mismatch");

            xMatrix C = new xMatrix(a.rows, b.cols);

            for (int i = 0; i < a.rows; i++) {
                for (int j = 0; j < b.cols; j++) {
                    for (int k = 0; k < a.cols; k++) {
                        C.At(i, j) += a.At(i, k) * b.At(k, j);
                    }
                }
            }

            return C;
        }

        public xMatrix Trans() {
            return new xMatrix(cols, rows).ForAll((y, x) => m_mat[x, y]);
        }

        public double Determinant() {
            if (rows != cols)
                throw new Exception("rows != cols");

            if (rows == 1) {
                return m_mat[0, 0];
            }

            else if (rows == 2) {
                return m_mat[0, 0] * m_mat[1, 1] - m_mat[0, 1] * m_mat[1, 0];   // ad-bc
            }
            else if (rows == 3) {
                Func<int, int, int, double> M = (i0, i1, i2) => m_mat[i0 / 3, i0 % 3] * m_mat[i1 / 3, i1 % 3] * m_mat[i2 / 3, i2 % 3];

                return M(0, 4, 8) + M(1, 5, 6) + M(2, 3, 7) - M(0, 5, 7) - M(1, 3, 8) - M(2, 4, 6);
            }
            else {
                double d = 0;
                for (int i = 0; i < rows; i++) {
                    if (At(0, i) == 0.0) continue;
                    xMatrix C = new xMatrix(rows - 1, cols - 1);
                    for (int j = 0; j < C.rows; j++) {
                        for (int k = 0; k < C.cols; k++) {
                            C.At(j, k) = At(j + 1, (k + i + 1) % cols);
                        }
                    }
                    d += At(0, i) * C.Determinant();
                }
                return d;
            }
        }

        /* Routine for L-U Decomposition
		* ref> Numerical Recipes
		*/
        public static void ludcmp(ref xMatrix A, ref int[] index) {
            double TINY = 1.0e-20;
            int dim = A.rows;
            int imax = 0;
            double big = 0, dum = 0, sum = 0, temp = 0;

            xMatrix vv = new xMatrix(1, dim, 0.0);
            for (int i = 0; i < dim; i++) {
                big = .0;
                for (int j = 0; j < dim; j++)
                    if ((temp = Math.Abs(A.At(i, j))) > big) big = temp;
                //ASSERT(big != 0.0);
                vv.At(0, i) = 1.0 / big;
            }
            for (int j = 0; j < dim; j++) {
                for (int i = 0; i < j; i++) {
                    sum = A.At(i, j);
                    for (int k = 0; k < i; k++)
                        sum -= A.At(i, k) * A.At(k, j);
                    A.At(i, j) = sum;
                }
                big = .0;
                for (int i = j; i < dim; i++) {
                    sum = A.At(i, j);
                    for (int k = 0; k < j; k++)
                        sum -= A.At(i, k) * A.At(k, j);
                    A.At(i, j) = sum;
                    if ((dum = vv.At(0, i) * Math.Abs(sum)) >= big) {
                        big = dum;
                        imax = i;
                    }
                }
                if (j != imax) {
                    for (int k = 0; k < dim; k++) {
                        dum = A.At(imax, k);
                        A.At(imax, k) = A.At(j, k);
                        A.At(j, k) = dum;
                    }
                    vv.At(0, imax) = vv.At(0, j);
                }
                index[j] = imax;
                if (A.At(j, j) == .0) A.At(j, j) = TINY;
                if (j != (dim - 1)) {
                    dum = 1.0 / A.At(j, j);
                    for (int i = j + 1; i < dim; i++) A.At(i, j) *= dum;
                }
            }
        }

        /* Routine for L-U Back Substitution
		* ref> Numerical Recipes
		*/
        public void lubksb(ref xMatrix A, ref xMatrix Col, ref int[] index) {
            int ii = -1;
            int dim = A.rows;
            double sum = 0;

            for (int i = 0; i < dim; i++) {
                int ip = index[i];
                sum = Col.At(0, ip);
                Col.At(0, ip) = Col.At(0, i);
                if (ii != -1) {
                    for (int j = ii; j < i; j++) {
                        sum -= A.At(i, j) * Col.At(0, j);
                    }
                }
                else {
                    if (sum != 0.0)
                        ii = i;
                }
                Col.At(0, i) = sum;
            }
            for (int i = dim - 1; i >= 0; i--) {
                sum = Col.At(0, i);
                for (int j = i + 1; j < dim; j++) sum -= A.At(i, j) * Col.At(0, j);
                Col.At(0, i) = sum / A.At(i, i);
            }
        }

        public xMatrix Inverse(double tiny, ref bool bOK) {
            if (rows != cols) {
                throw new Exception("rows != cols");
            }
            if (Determinant() <= tiny) {
                bOK = false;
                return eye(rows);
            }
            int dim = rows;
            int[] index = new int[rows];

            //	ASSERT(det() != 0.0);
            xMatrix A = new xMatrix(this);
            xMatrix B = new xMatrix(rows, cols, 0.0);
            xMatrix col = new xMatrix(1, cols, 0.0);

            ludcmp(ref A, ref index);

            for (int j = 0; j < dim; j++) {
                for (int i = 0; i < dim; i++)
                    col.At(0, i) = .0;
                col.At(0, j) = 1.0;
                lubksb(ref A, ref col, ref index);
                for (int i = 0; i < dim; i++)
                    B.At(i, j) = col.At(0, i);
            }

            bOK = true;

            return B;
        }

    }
}
