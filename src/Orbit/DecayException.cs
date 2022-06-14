using System.Numerics;

namespace System;

public sealed class DecayException<T> : PropagationException where T : INumber<T>, IFloatingPoint<T>
{
    /// <summary>
    /// The GMT when the satellite orbit decays.
    /// </summary>
    public DateTime DecayTime { get; private set; }

    /// <summary>
    /// The name of the satellite whose orbit decayed.
    /// </summary>
    public IKeplerianElements<T> KeplerianElements { get; private set; }

    internal DecayException(Julian decayTime, IKeplerianElements<T> keplerianElements)
       : this(decayTime.ToTime(), keplerianElements)
    {

    }

    public DecayException(DateTime decayTime, IKeplerianElements<T> keplerianElements)
    {
        DecayTime = decayTime;
        KeplerianElements = keplerianElements;
    }
}