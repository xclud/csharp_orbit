using System.Numerics;

namespace System.Astronomy;

public sealed class Topocentric<T>(T south, T east, T normal) where T : INumber<T>, IFloatingPoint<T>
{
    public readonly T South = south;
    public readonly T East = east;
    public readonly T Normal = normal;
}
