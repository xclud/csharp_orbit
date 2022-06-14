using System.Numerics;

namespace System;

public readonly struct Vector4<T> where T : INumber<T>, IFloatingPoint<T>
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

    public static Vector4<T> operator *(Vector4<T> left, T right)
    {
        return new Vector4<T>(left.X * right, left.Y * right, left.Z * right, left.W * right);
    }

    public static Vector4<T> operator *(T left, Vector4<T> right)
    {
        return new Vector4<T>(right.X * left, right.Y * left, right.Z * left, right.W * left);
    }
}
