﻿using System.Text.RegularExpressions;

namespace System.Astronomy;

internal partial class Matrix
{
    public int rows;
    public int cols;
    public double[,] mat;

    public Matrix? L;
    public Matrix? U;
    private int[]? pi;
    private double detOfP = 1;


    public Matrix(int iRows, int iCols)         // Matrix Class constructor
    {
        rows = iRows;
        cols = iCols;
        mat = new double[rows, cols];
    }

    public static Matrix ABS(Matrix C)          // C ABS(C)
    {
        for (int i = 0; i < C.cols; i++)
        {
            C[0, i] = Math.Abs(C[0, i]);
        }

        return C;
    }
    public static Matrix ScalerInverse(Matrix C)          // C ABS(C)
    {
        for (int i = 0; i < C.cols; i++)
        {
            C[0, i] = 1 / C[0, i];
        }

        return C;
    }
    public static Matrix ScalerPower(Matrix C, double n)          // C ABS(C)
    {
        for (int i = 0; i < C.cols; i++)
        {
            C[0, i] = Math.Pow(C[0, i], n);
        }

        return C;
    }
    public bool IsSquare()
    {
        return rows == cols;
    }

    public double this[int iRow, int iCol]      // Access this matrix as a 2D array
    {
        get => mat[iRow, iCol];
        set => mat[iRow, iCol] = value;
    }

    public Matrix GetCol(int k)
    {
        Matrix m = new(rows, 1);
        for (int i = 0; i < rows; i++)
        {
            m[i, 0] = mat[i, k];
        }

        return m;
    }

    public void SetCol(Matrix v, int k)
    {
        for (int i = 0; i < rows; i++)
        {
            mat[i, k] = v[i, 0];
        }
    }

    public void MakeLU()                        // Function for LU decomposition
    {
        if (!IsSquare())
        {
            throw new MatrixException("The matrix is not square!");
        }

        L = IdentityMatrix(rows, cols);
        U = Duplicate();

        pi = new int[rows];
        for (int i = 0; i < rows; i++)
        {
            pi[i] = i;
        }

        double pom2;
        int k0 = 0;

        for (int k = 0; k < cols - 1; k++)
        {
            var p = 0.0;
            for (int i = k; i < rows; i++)      // find the row with the biggest pivot
            {
                if (Math.Abs(U[i, k]) > p)
                {
                    p = Math.Abs(U[i, k]);
                    k0 = i;
                }
            }
            if (p == 0) // samé nuly ve sloupci
            {
                throw new MatrixException("The matrix is singular!");
            }

            // switch two rows in permutation matrix
            (pi[k0], pi[k]) = (pi[k], pi[k0]);

            for (int i = 0; i < k; i++)
            {
                pom2 = L[k, i]; L[k, i] = L[k0, i]; L[k0, i] = pom2;
            }

            if (k != k0)
            {
                detOfP *= -1;
            }

            for (int i = 0; i < cols; i++)                  // Switch rows in U
            {
                pom2 = U[k, i]; U[k, i] = U[k0, i]; U[k0, i] = pom2;
            }

            for (int i = k + 1; i < rows; i++)
            {
                L[i, k] = U[i, k] / U[k, k];
                for (int j = k; j < cols; j++)
                {
                    U[i, j] = U[i, j] - L[i, k] * U[k, j];
                }
            }
        }
    }


    public Matrix SolveWith(Matrix v)                        // Function solves Ax = v in confirmity with solution vector "v"
    {
        if (rows != cols)
        {
            throw new MatrixException("The matrix is not square!");
        }

        if (rows != v.rows)
        {
            throw new MatrixException("Wrong number of results in solution vector!");
        }

        if (L == null)
        {
            MakeLU();
        }

        Matrix b = new(rows, 1);
        for (int i = 0; i < rows; i++)
        {
            b[i, 0] = v[pi![i], 0];   // switch two items in "v" due to permutation matrix
        }

        Matrix z = SubsForth(L!, b);
        Matrix x = SubsBack(U!, z);

        return x;
    }

    public Matrix Invert()                                   // Function returns the inverted matrix
    {
        if (L == null)
        {
            MakeLU();
        }

        Matrix inv = new(rows, cols);

        for (int i = 0; i < rows; i++)
        {
            Matrix Ei = ZeroMatrix(rows, 1);
            Ei[i, 0] = 1;
            Matrix col = SolveWith(Ei);
            inv.SetCol(col, i);
        }
        return inv;
    }


    public double Det()                         // Function for determinant
    {
        if (L == null)
        {
            MakeLU();
        }

        double det = detOfP;
        for (int i = 0; i < rows; i++)
        {
            det *= U![i, i];
        }

        return det;
    }
    public double Size()                         // Function for determinant
    {
        return Math.Max(rows, cols);
    }
    public Matrix GetP()                        // Function returns permutation matrix "P" due to permutation vector "pi"
    {
        if (L == null)
        {
            MakeLU();
        }

        Matrix matrix = ZeroMatrix(rows, cols);
        for (int i = 0; i < rows; i++)
        {
            matrix[pi![i], i] = 1;
        }

        return matrix;
    }

    public Matrix Duplicate()                   // Function returns the copy of this matrix
    {
        Matrix matrix = new(rows, cols);
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                matrix[i, j] = mat[i, j];
            }
        }

        return matrix;
    }

    public static Matrix SubsForth(Matrix A, Matrix b)          // Function solves Ax = b for A as a lower triangular matrix
    {
        if (A.L == null)
        {
            A.MakeLU();
        }

        int n = A.rows;
        Matrix x = new(n, 1);

        for (int i = 0; i < n; i++)
        {
            x[i, 0] = b[i, 0];
            for (int j = 0; j < i; j++)
            {
                x[i, 0] -= A[i, j] * x[j, 0];
            }

            x[i, 0] = x[i, 0] / A[i, i];
        }
        return x;
    }

    public static Matrix SubsBack(Matrix A, Matrix b)           // Function solves Ax = b for A as an upper triangular matrix
    {
        if (A.L == null)
        {
            A.MakeLU();
        }

        int n = A.rows;
        Matrix x = new(n, 1);

        for (int i = n - 1; i > -1; i--)
        {
            x[i, 0] = b[i, 0];
            for (int j = n - 1; j > i; j--)
            {
                x[i, 0] -= A[i, j] * x[j, 0];
            }

            x[i, 0] = x[i, 0] / A[i, i];
        }
        return x;
    }

    public static Matrix ZeroMatrix(int iRows, int iCols)       // Function generates the zero matrix
    {
        Matrix matrix = new(iRows, iCols);
        for (int i = 0; i < iRows; i++)
        {
            for (int j = 0; j < iCols; j++)
            {
                matrix[i, j] = 0;
            }
        }

        return matrix;
    }

    public static Matrix IdentityMatrix(int iRows, int iCols)   // Function generates the identity matrix
    {
        Matrix matrix = ZeroMatrix(iRows, iCols);
        for (int i = 0; i < Math.Min(iRows, iCols); i++)
        {
            matrix[i, i] = 1;
        }

        return matrix;
    }

    public static Matrix RandomMatrix(int iRows, int iCols, int dispersion)       // Function generates the random matrix
    {
        Random random = new();
        Matrix matrix = new(iRows, iCols);
        for (int i = 0; i < iRows; i++)
        {
            for (int j = 0; j < iCols; j++)
            {
                matrix[i, j] = random.Next(-dispersion, dispersion);
            }
        }

        return matrix;
    }

    public static Matrix Parse(string ps)                        // Function parses the matrix from string
    {
        string s = NormalizeMatrixString(ps);
        string[] rows = CRLF().Split(s);
        string[] nums = rows[0].Split(' ');
        Matrix matrix = new(rows.Length, nums.Length);
        try
        {
            for (int i = 0; i < rows.Length; i++)
            {
                nums = rows[i].Split(' ');
                for (int j = 0; j < nums.Length; j++)
                {
                    matrix[i, j] = double.Parse(nums[j]);
                }
            }
        }
        catch (FormatException exc) { throw new MatrixException("Wrong input format!"); }
        return matrix;
    }
    public static Matrix Parse2(string ps)                        // Function parses the matrix from string
    {
        string s = NormalizeMatrixString(ps);
        string[] rows = Semicolon().Split(s);
        string[] nums = rows[0].Split(' ');
        Matrix matrix = new(rows.Length, nums.Length);
        try
        {
            for (int i = 0; i < rows.Length; i++)
            {
                nums = rows[i].Split(' ');
                for (int j = 0; j < nums.Length; j++)
                {
                    matrix[i, j] = double.Parse(nums[j]);
                }
            }
        }
        catch (FormatException exc)
        {
            throw new MatrixException($"Wrong input format! {exc.Message}");
        }

        return matrix;
    }

    public override string ToString()                           // Function returns matrix as a string
    {
        string s = "";
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                s += string.Format("{0,5:0.0000000000000}", mat[i, j]) + " ";
            }

            s += "\r\n";
        }
        return s;
    }

    public static Matrix Transpose(Matrix m)              // Matrix transpose, for any rectangular matrix
    {
        Matrix t = new(m.cols, m.rows);
        for (int i = 0; i < m.rows; i++)
        {
            for (int j = 0; j < m.cols; j++)
            {
                t[j, i] = m[i, j];
            }
        }

        return t;
    }

    public static Matrix Power(Matrix m, int pow)           // Power matrix to exponent
    {
        if (pow == 0)
        {
            return IdentityMatrix(m.rows, m.cols);
        }

        if (pow == 1)
        {
            return m.Duplicate();
        }

        if (pow == -1)
        {
            return m.Invert();
        }

        Matrix x;
        if (pow < 0) { x = m.Invert(); pow *= -1; }
        else
        {
            x = m.Duplicate();
        }

        Matrix ret = IdentityMatrix(m.rows, m.cols);
        while (pow != 0)
        {
            if ((pow & 1) == 1)
            {
                ret *= x;
            }

            x *= x;
            pow >>= 1;
        }
        return ret;
    }

    private static void SafeAplusBintoC(Matrix A, int xa, int ya, Matrix B, int xb, int yb, Matrix C, int size)
    {
        for (int i = 0; i < size; i++)          // rows
        {
            for (int j = 0; j < size; j++)     // cols
            {
                C[i, j] = 0;
                if (xa + j < A.cols && ya + i < A.rows)
                {
                    C[i, j] += A[ya + i, xa + j];
                }

                if (xb + j < B.cols && yb + i < B.rows)
                {
                    C[i, j] += B[yb + i, xb + j];
                }
            }
        }
    }

    private static void SafeAminusBintoC(Matrix A, int xa, int ya, Matrix B, int xb, int yb, Matrix C, int size)
    {
        for (int i = 0; i < size; i++)          // rows
        {
            for (int j = 0; j < size; j++)     // cols
            {
                C[i, j] = 0;
                if (xa + j < A.cols && ya + i < A.rows)
                {
                    C[i, j] += A[ya + i, xa + j];
                }

                if (xb + j < B.cols && yb + i < B.rows)
                {
                    C[i, j] -= B[yb + i, xb + j];
                }
            }
        }
    }

    private static void SafeACopytoC(Matrix A, int xa, int ya, Matrix C, int size)
    {
        for (int i = 0; i < size; i++)          // rows
        {
            for (int j = 0; j < size; j++)     // cols
            {
                C[i, j] = 0;
                if (xa + j < A.cols && ya + i < A.rows)
                {
                    C[i, j] += A[ya + i, xa + j];
                }
            }
        }
    }

    private static void AplusBintoC(Matrix A, int xa, int ya, Matrix B, int xb, int yb, Matrix C, int size)
    {
        for (int i = 0; i < size; i++)          // rows
        {
            for (int j = 0; j < size; j++)
            {
                C[i, j] = A[ya + i, xa + j] + B[yb + i, xb + j];
            }
        }
    }

    private static void AminusBintoC(Matrix A, int xa, int ya, Matrix B, int xb, int yb, Matrix C, int size)
    {
        for (int i = 0; i < size; i++)          // rows
        {
            for (int j = 0; j < size; j++)
            {
                C[i, j] = A[ya + i, xa + j] - B[yb + i, xb + j];
            }
        }
    }

    private static void ACopytoC(Matrix A, int xa, int ya, Matrix C, int size)
    {
        for (int i = 0; i < size; i++)          // rows
        {
            for (int j = 0; j < size; j++)
            {
                C[i, j] = A[ya + i, xa + j];
            }
        }
    }

    private static Matrix StrassenMultiply(Matrix A, Matrix B)                // Smart matrix multiplication
    {
        if (A.cols != B.rows)
        {
            throw new MatrixException("Wrong dimension of matrix!");
        }

        Matrix R;

        int msize = Math.Max(Math.Max(A.rows, A.cols), Math.Max(B.rows, B.cols));

        if (msize < 32)
        {
            R = ZeroMatrix(A.rows, B.cols);
            for (int i = 0; i < R.rows; i++)
            {
                for (int j = 0; j < R.cols; j++)
                {
                    for (int k = 0; k < A.cols; k++)
                    {
                        R[i, j] += A[i, k] * B[k, j];
                    }
                }
            }

            return R;
        }

        int size = 1; int n = 0;
        while (msize > size) { size *= 2; n++; };
        int h = size / 2;


        Matrix[,] mField = new Matrix[n, 9];

        /*
         *  8x8, 8x8, 8x8, ...
         *  4x4, 4x4, 4x4, ...
         *  2x2, 2x2, 2x2, ...
         *  . . .
         */

        int z;
        for (int i = 0; i < n - 4; i++)          // rows
        {
            z = (int)Math.Pow(2, n - i - 1);
            for (int j = 0; j < 9; j++)
            {
                mField[i, j] = new Matrix(z, z);
            }
        }

        SafeAplusBintoC(A, 0, 0, A, h, h, mField[0, 0], h);
        SafeAplusBintoC(B, 0, 0, B, h, h, mField[0, 1], h);
        StrassenMultiplyRun(mField[0, 0], mField[0, 1], mField[0, 1 + 1], 1, mField); // (A11 + A22) * (B11 + B22);

        SafeAplusBintoC(A, 0, h, A, h, h, mField[0, 0], h);
        SafeACopytoC(B, 0, 0, mField[0, 1], h);
        StrassenMultiplyRun(mField[0, 0], mField[0, 1], mField[0, 1 + 2], 1, mField); // (A21 + A22) * B11;

        SafeACopytoC(A, 0, 0, mField[0, 0], h);
        SafeAminusBintoC(B, h, 0, B, h, h, mField[0, 1], h);
        StrassenMultiplyRun(mField[0, 0], mField[0, 1], mField[0, 1 + 3], 1, mField); //A11 * (B12 - B22);

        SafeACopytoC(A, h, h, mField[0, 0], h);
        SafeAminusBintoC(B, 0, h, B, 0, 0, mField[0, 1], h);
        StrassenMultiplyRun(mField[0, 0], mField[0, 1], mField[0, 1 + 4], 1, mField); //A22 * (B21 - B11);

        SafeAplusBintoC(A, 0, 0, A, h, 0, mField[0, 0], h);
        SafeACopytoC(B, h, h, mField[0, 1], h);
        StrassenMultiplyRun(mField[0, 0], mField[0, 1], mField[0, 1 + 5], 1, mField); //(A11 + A12) * B22;

        SafeAminusBintoC(A, 0, h, A, 0, 0, mField[0, 0], h);
        SafeAplusBintoC(B, 0, 0, B, h, 0, mField[0, 1], h);
        StrassenMultiplyRun(mField[0, 0], mField[0, 1], mField[0, 1 + 6], 1, mField); //(A21 - A11) * (B11 + B12);

        SafeAminusBintoC(A, h, 0, A, h, h, mField[0, 0], h);
        SafeAplusBintoC(B, 0, h, B, h, h, mField[0, 1], h);
        StrassenMultiplyRun(mField[0, 0], mField[0, 1], mField[0, 1 + 7], 1, mField); // (A12 - A22) * (B21 + B22);

        R = new Matrix(A.rows, B.cols);                  // result

        /// C11
        for (int i = 0; i < Math.Min(h, R.rows); i++)          // rows
        {
            for (int j = 0; j < Math.Min(h, R.cols); j++)     // cols
            {
                R[i, j] = mField[0, 1 + 1][i, j] + mField[0, 1 + 4][i, j] - mField[0, 1 + 5][i, j] + mField[0, 1 + 7][i, j];
            }
        }

        /// C12
        for (int i = 0; i < Math.Min(h, R.rows); i++)          // rows
        {
            for (int j = h; j < Math.Min(2 * h, R.cols); j++)     // cols
            {
                R[i, j] = mField[0, 1 + 3][i, j - h] + mField[0, 1 + 5][i, j - h];
            }
        }

        /// C21
        for (int i = h; i < Math.Min(2 * h, R.rows); i++)          // rows
        {
            for (int j = 0; j < Math.Min(h, R.cols); j++)     // cols
            {
                R[i, j] = mField[0, 1 + 2][i - h, j] + mField[0, 1 + 4][i - h, j];
            }
        }

        /// C22
        for (int i = h; i < Math.Min(2 * h, R.rows); i++)          // rows
        {
            for (int j = h; j < Math.Min(2 * h, R.cols); j++)     // cols
            {
                R[i, j] = mField[0, 1 + 1][i - h, j - h] - mField[0, 1 + 2][i - h, j - h] + mField[0, 1 + 3][i - h, j - h] + mField[0, 1 + 6][i - h, j - h];
            }
        }

        return R;
    }

    // function for square matrix 2^N x 2^N

    private static void StrassenMultiplyRun(Matrix A, Matrix B, Matrix C, int l, Matrix[,] f)    // A * B into C, level of recursion, matrix field
    {
        int size = A.rows;
        int h = size / 2;

        if (size < 32)
        {
            for (int i = 0; i < C.rows; i++)
            {
                for (int j = 0; j < C.cols; j++)
                {
                    C[i, j] = 0;
                    for (int k = 0; k < A.cols; k++)
                    {
                        C[i, j] += A[i, k] * B[k, j];
                    }
                }
            }

            return;
        }

        AplusBintoC(A, 0, 0, A, h, h, f[l, 0], h);
        AplusBintoC(B, 0, 0, B, h, h, f[l, 1], h);
        StrassenMultiplyRun(f[l, 0], f[l, 1], f[l, 1 + 1], l + 1, f); // (A11 + A22) * (B11 + B22);

        AplusBintoC(A, 0, h, A, h, h, f[l, 0], h);
        ACopytoC(B, 0, 0, f[l, 1], h);
        StrassenMultiplyRun(f[l, 0], f[l, 1], f[l, 1 + 2], l + 1, f); // (A21 + A22) * B11;

        ACopytoC(A, 0, 0, f[l, 0], h);
        AminusBintoC(B, h, 0, B, h, h, f[l, 1], h);
        StrassenMultiplyRun(f[l, 0], f[l, 1], f[l, 1 + 3], l + 1, f); //A11 * (B12 - B22);

        ACopytoC(A, h, h, f[l, 0], h);
        AminusBintoC(B, 0, h, B, 0, 0, f[l, 1], h);
        StrassenMultiplyRun(f[l, 0], f[l, 1], f[l, 1 + 4], l + 1, f); //A22 * (B21 - B11);

        AplusBintoC(A, 0, 0, A, h, 0, f[l, 0], h);
        ACopytoC(B, h, h, f[l, 1], h);
        StrassenMultiplyRun(f[l, 0], f[l, 1], f[l, 1 + 5], l + 1, f); //(A11 + A12) * B22;

        AminusBintoC(A, 0, h, A, 0, 0, f[l, 0], h);
        AplusBintoC(B, 0, 0, B, h, 0, f[l, 1], h);
        StrassenMultiplyRun(f[l, 0], f[l, 1], f[l, 1 + 6], l + 1, f); //(A21 - A11) * (B11 + B12);

        AminusBintoC(A, h, 0, A, h, h, f[l, 0], h);
        AplusBintoC(B, 0, h, B, h, h, f[l, 1], h);
        StrassenMultiplyRun(f[l, 0], f[l, 1], f[l, 1 + 7], l + 1, f); // (A12 - A22) * (B21 + B22);

        /// C11
        for (int i = 0; i < h; i++)          // rows
        {
            for (int j = 0; j < h; j++)     // cols
            {
                C[i, j] = f[l, 1 + 1][i, j] + f[l, 1 + 4][i, j] - f[l, 1 + 5][i, j] + f[l, 1 + 7][i, j];
            }
        }

        /// C12
        for (int i = 0; i < h; i++)          // rows
        {
            for (int j = h; j < size; j++)     // cols
            {
                C[i, j] = f[l, 1 + 3][i, j - h] + f[l, 1 + 5][i, j - h];
            }
        }

        /// C21
        for (int i = h; i < size; i++)          // rows
        {
            for (int j = 0; j < h; j++)     // cols
            {
                C[i, j] = f[l, 1 + 2][i - h, j] + f[l, 1 + 4][i - h, j];
            }
        }

        /// C22
        for (int i = h; i < size; i++)          // rows
        {
            for (int j = h; j < size; j++)     // cols
            {
                C[i, j] = f[l, 1 + 1][i - h, j - h] - f[l, 1 + 2][i - h, j - h] + f[l, 1 + 3][i - h, j - h] + f[l, 1 + 6][i - h, j - h];
            }
        }
    }

    public static Matrix StupidMultiply(Matrix m1, Matrix m2)                  // Stupid matrix multiplication
    {
        if (m1.cols != m2.rows)
        {
            throw new MatrixException("Wrong dimensions of matrix!");
        }

        Matrix result = ZeroMatrix(m1.rows, m2.cols);
        for (int i = 0; i < result.rows; i++)
        {
            for (int j = 0; j < result.cols; j++)
            {
                for (int k = 0; k < m1.cols; k++)
                {
                    result[i, j] += m1[i, k] * m2[k, j];
                }
            }
        }

        return result;
    }
    private static Matrix Multiply(double n, Matrix m)                          // Multiplication by constant n
    {
        Matrix r = new(m.rows, m.cols);
        for (int i = 0; i < m.rows; i++)
        {
            for (int j = 0; j < m.cols; j++)
            {
                r[i, j] = m[i, j] * n;
            }
        }

        return r;
    }
    private static Matrix Add(Matrix m1, Matrix m2)         // Sčítání matic
    {
        if (m1.rows != m2.rows || m1.cols != m2.cols)
        {
            throw new MatrixException("Matrices must have the same dimensions!");
        }

        Matrix r = new(m1.rows, m1.cols);
        for (int i = 0; i < r.rows; i++)
        {
            for (int j = 0; j < r.cols; j++)
            {
                r[i, j] = m1[i, j] + m2[i, j];
            }
        }

        return r;
    }
    public static double Norm(Matrix val)
    {
        //Only For 1*3 
        return val.rows == 1 && val.cols == 3
            ? Math.Sqrt(Math.Pow(val[0, 0], 2) + Math.Pow(val[0, 1], 2) + Math.Pow(val[0, 2], 2))
            : val.rows == 3 && val.cols == 1 ? Math.Sqrt(Math.Pow(val[0, 0], 2) + Math.Pow(val[1, 0], 2) + Math.Pow(val[2, 0], 2)) : -1;
    }
    public static double Dot(Matrix val1, Matrix val2)
    {


        return val1.rows == 1 && val1.cols == 3
            ? val1[0, 0] * val2[0, 0] + val1[0, 1] * val2[0, 1] + val1[0, 2] * val2[0, 2]
            : val1.rows == 3 && val1.cols == 1 ? val1[0, 0] * val2[0, 0] + val1[1, 0] * val2[1, 0] + val1[2, 0] * val2[2, 0] : -1;
    }

    public static Matrix Cross(Matrix val1, Matrix val2)
    {
        Matrix Rc = new(1, 3);
        Rc[0, 0] = val1[0, 1] * val2[0, 2] - val1[0, 2] * val2[0, 1];
        Rc[0, 1] = val1[0, 2] * val2[0, 0] - val1[0, 0] * val2[0, 2];
        Rc[0, 2] = val1[0, 0] * val2[0, 1] - val1[0, 1] * val2[0, 0];
        return Rc;
    }

    public static Matrix Merge(Matrix m1, Matrix m2)
    {
        if (m1.rows != m2.rows)
        {
            throw new MatrixException("Wrong dimensions of matrix!");
        }

        Matrix result = ZeroMatrix(m1.rows, m1.cols + m2.cols);
        for (int i = 0; i < m1.rows; i++)
        {
            for (int j = 0; j < m1.cols; j++)
            {
                result[i, j] = m1[i, j];
            }
        }

        for (int i = 0; i < m2.rows; i++)
        {
            for (int j = 0; j < m2.cols; j++)
            {
                result[i, j + m1.cols] = m2[i, j];
            }
        }

        return result;
    }
    public static string NormalizeMatrixString(string matStr)   // From Andy - thank you! :)
    {
        // Remove any multiple spaces
        while (matStr.IndexOf("  ") != -1)
        {
            matStr = matStr.Replace("  ", " ");
        }

        // Remove any spaces before or after newlines
        matStr = matStr.Replace(" \r\n", "\r\n");
        matStr = matStr.Replace("\r\n ", "\r\n");

        // If the data ends in a newline, remove the trailing newline.
        // Make it easier by first replacing \r\n’s with |’s then
        // restore the |’s with \r\n’s
        matStr = matStr.Replace("\r\n", "|");
        while (matStr.LastIndexOf("|") == matStr.Length - 1)
        {
            matStr = matStr[..^1];
        }

        matStr = matStr.Replace("|", "\r\n");
        return matStr.Trim();
    }

    //   O P E R A T O R S

    public static Matrix operator -(Matrix m)
    { return Multiply(-1, m); }

    public static Matrix operator +(Matrix m1, Matrix m2)
    { return Add(m1, m2); }

    public static Matrix operator -(Matrix m1, Matrix m2)
    { return Add(m1, -m2); }

    public static Matrix operator *(Matrix m1, Matrix m2)
    { return StrassenMultiply(m1, m2); }

    public static Matrix operator *(double n, Matrix m)
    { return Multiply(n, m); }

    [GeneratedRegex("\r\n")]
    private static partial Regex CRLF();
    [GeneratedRegex(";")]
    private static partial Regex Semicolon();
}

internal sealed class MatrixException(string Message) : Exception(Message)
{
}
