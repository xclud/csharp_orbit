using System.Numerics;

namespace System;


public sealed class Topocentric<T> where T : INumber<T>, IFloatingPoint<T>
{
    public Topocentric(T south, T east, T normal)
    {
        this.South = south;
        this.East = east;
        this.Normal = normal;
    }

    public readonly T South;
    public readonly T East;
    public readonly T Normal;
}
