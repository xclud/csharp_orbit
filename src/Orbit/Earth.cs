namespace System;

public static class Earth
{
    public static readonly IPlanet WGS72 = new _WGS72();
    public static readonly IPlanet WGS84 = new _WGS84();

    private sealed class _WGS72 : IPlanet
    {
        public double Radius => 6378.135;
        public double Mu => 398600.8;
        public double J2 => 0.001082616;
        public double J3 => -0.00000253881;
        public double J4 => -0.00000165597;
        public double Flattening => 1 / 298.26;
    }

    private sealed class _WGS84 : IPlanet
    {
        public double Radius => 6378.137;
        public double Mu => 398600.5;
        public double J2 => 0.00108262998905;
        public double J3 => -0.00000253215306;
        public double J4 => -0.00000161098761;
        public double Flattening => 1 / 298.257223563;
    }
}
