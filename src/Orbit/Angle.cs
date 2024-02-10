using System.Numerics;

namespace System.Astronomy;

public static class Angle
{
    public static Angle<T> FromDegrees<T>(T degrees) where T : INumber<T>, IFloatingPoint<T>
    {
        return Angle<T>.FromDegrees(degrees);
    }

    public static Angle<T> FromRadians<T>(T radians) where T : INumber<T>, IFloatingPoint<T>
    {
        return Angle<T>.FromRadians(radians);
    }
}

public sealed record class Angle<T> where T : INumber<T>, IFloatingPoint<T>
{
    public static readonly Angle<T> Zero = new Angle<T>(T.Zero, T.Zero);
    public static readonly Angle<T> Pi = Angle<T>.FromRadians(T.Pi);
    public static readonly Angle<T> TwoPi = Angle<T>.FromRadians(T.Pi * (T.One + T.One));

    private static readonly T _180 = T.Parse("180", null);

    public readonly T Degrees;
    public readonly T Radians;

    private Angle(T degrees, T radians)
    {
        Degrees = degrees;
        Radians = radians;
    }

    public static Angle<T> operator *(Angle<T> angle, T scalar)
    {
        return new Angle<T>(angle.Degrees * scalar, angle.Radians * scalar);
    }

    public static Angle<T> operator *(T scalar, Angle<T> angle)
    {
        return new Angle<T>(angle.Degrees * scalar, angle.Radians * scalar);
    }

    public static Angle<T> operator /(Angle<T> angle, T scalar)
    {
        return new Angle<T>(angle.Degrees / scalar, angle.Radians / scalar);
    }

    public static Angle<T> operator +(Angle<T> left, Angle<T> right)
    {
        return new Angle<T>(right.Degrees + left.Degrees, right.Radians + left.Radians);
    }
    public static Angle<T> operator -(Angle<T> left, Angle<T> right)
    {
        return new Angle<T>(right.Degrees - left.Degrees, right.Radians - left.Radians);
    }

    public static Angle<T> operator -(Angle<T> angle)
    {
        return new Angle<T>(-angle.Degrees, -angle.Radians);
    }

    public static Angle<T> FromDegrees(T degrees)
    {
        return new Angle<T>(degrees, degrees / _180 * T.Pi);
    }

    public static Angle<T> FromRadians(T radians)
    {
        return new Angle<T>(radians * _180 / T.Pi, radians);
    }

    public override string ToString()
    {
        return $"{Degrees}° ({Radians} rad)";
    }
}