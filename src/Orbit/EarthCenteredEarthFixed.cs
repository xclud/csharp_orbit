using System.Numerics;

namespace System.Astronomy;

public sealed class EarthCenteredEarthFixed<T>(T x, T y, T z) where T : INumber<T>, IFloatingPoint<T>
{
    public readonly T X = x;
    public readonly T Y = y;
    public readonly T Z = z;
}