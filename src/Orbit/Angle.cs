using System.Numerics;

namespace System.Astronomy;

public sealed record class Angle<T> where T : INumber<T>, IFloatingPoint<T>
{
    public readonly T Degrees;
    public readonly T Radians;

    private static readonly T _180 = T.Parse("180", null);

    private Angle(T degrees, T radians)
    {
        Degrees = degrees;
        Radians = radians;
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