using System.Numerics;

namespace System;

public sealed class KeplerianElements<T> where T : INumber<T>, IFloatingPoint<T>
{
    public KeplerianElements(T epoch, T eccentricity, T inclination, T rightAscensionOfAscendingNode, T argumentOfPeriapsis, T meanMotion, T meanAnomaly, T drag, int revolutionNumberAtEpoch)
    {
        Epoch = epoch;
        Eccentricity = eccentricity;
        Inclination = inclination;
        RightAscensionOfAscendingNode = rightAscensionOfAscendingNode;
        ArgumentOfPeriapsis = argumentOfPeriapsis;
        MeanMotion = meanMotion;
        MeanAnomaly = meanAnomaly;
        Drag = drag;
        RevolutionNumberAtEpoch = revolutionNumberAtEpoch;
    }

    public readonly T Epoch;

    public readonly T Eccentricity;

    public readonly T Inclination;

    public readonly T RightAscensionOfAscendingNode;

    public readonly T ArgumentOfPeriapsis;

    public readonly T MeanMotion;

    public readonly T MeanAnomaly;

    public readonly T Drag;

    public readonly int RevolutionNumberAtEpoch;
}