using System.Numerics;

namespace System;

public ref struct Vector4<T> where T : INumber<T>, IFloatingPoint<T>
{
    public Vector4(T x, T y, T z, T w)
    {
        X = x;
        Y = y;
        Z = z;
        W = w;
    }

    public readonly T X;
    public readonly T Y;
    public readonly T Z;
    public readonly T W;
}
