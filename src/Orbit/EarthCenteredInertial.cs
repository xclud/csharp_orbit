using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace System;

public sealed class EarthCenteredInertial<T> where T : INumber<T>, IFloatingPoint<T>
{
    public EarthCenteredInertial(T x, T y, T z)
    {
        X = x;
        Y = y;
        Z = z;
    }

    public readonly T X;
    public readonly T Y;
    public readonly T Z;
}

public sealed class EarthCenteredEarthFixed<T> where T : INumber<T>, IFloatingPoint<T>
{
    public EarthCenteredEarthFixed(T x, T y, T z)
    {
        X = x;
        Y = y;
        Z = z;
    }

    public readonly T X;
    public readonly T Y;
    public readonly T Z;
}