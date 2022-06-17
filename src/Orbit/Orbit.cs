namespace System;


/// <summary>
/// This class accepts a single satellite's NORAD two-line element
/// set and provides information regarding the satellite's orbit 
/// such as period, axis length, ECI coordinates, velocity, etc.
/// </summary>
public sealed class Orbit
{
    public readonly TimeSpan Period;

    private readonly double m_Inclination;
    private readonly double m_Eccentricity;
    private readonly double m_RAAN;
    private readonly double m_ArgPerigee;
    private readonly double m_BStar;
    private readonly double m_MeanAnomaly;

    // Caching variables recovered from the input TLE elements
    private readonly double m_aeAxisSemiMajorRec;  // semi-major axis, in AE units
    private readonly double m_aeAxisSemiMinorRec;  // semi-minor axis, in AE units
    private readonly double m_rmMeanMotionRec;     // radians per minute

    /// <summary>
    /// In Kilometers.
    /// </summary>
    public readonly double Periapsis;

    /// <summary>
    /// In Kilometers.
    /// </summary>
    public readonly double Apoapsis;


    internal readonly double Xke;

    #region Properties

    internal Julian EpochJ { get; private set; }
    public DateTime Epoch => EpochJ.ToTime();

    private ICartesianElements NoradModel { get; set; }

    // "Recovered" from the input elements
    public double SemiMajor => m_aeAxisSemiMajorRec;
    public double SemiMinor => m_aeAxisSemiMinorRec;
    public double MeanMotion => m_rmMeanMotionRec;
    public double Major => 2.0 * SemiMajor;
    public double Minor => 2.0 * SemiMinor;

    public double Inclination => m_Inclination;
    public double Eccentricity => m_Eccentricity;
    public double RAAN => m_RAAN;
    public double ArgPerigee => m_ArgPerigee;
    public double BStar => m_BStar;
    public double MeanAnomaly => m_MeanAnomaly;

    public readonly IPlanet Planet;
    public readonly IKeplerianElements<double> KeplerianElements;

    #endregion


    /// <summary>
    /// Standard constructor.
    /// </summary>
    /// <param name="keplerianElements">Two-line element orbital parameters.</param>
    public Orbit(IKeplerianElements<double> keplerianElements, IPlanet planet)
    {
        this.Planet = planet;
        this.KeplerianElements = keplerianElements;

        var year = (int)(this.KeplerianElements.Epoch / 1000.0);
        var doy = this.KeplerianElements.Epoch - (year * 1000.0);

        year += year > 57 ? 1900 : 2000;
        EpochJ = new Julian(year, doy);


        m_Inclination = GetRad(this.KeplerianElements.Inclination);
        m_Eccentricity = this.KeplerianElements.Eccentricity;
        m_RAAN = GetRad(this.KeplerianElements.RightAscensionOfAscendingNode);
        m_ArgPerigee = GetRad(this.KeplerianElements.ArgumentOfPeriapsis);
        m_BStar = this.KeplerianElements.Drag;
        m_MeanAnomaly = GetRad(this.KeplerianElements.MeanAnomaly);

        // Recover the original mean motion and semi-major axis from the
        // input elements.
        double mm = keplerianElements.MeanMotion;
        double rpmin = mm * Globals.TwoPi / Globals.MinPerDay;   // rads per minute

        var mu = planet.Mu;
        var Xkmper = planet.Radius;

        this.Xke = Math.Sqrt(3600.0 * mu / (Xkmper * Xkmper * Xkmper));
        var Ck2 = planet.J2 / 2;

        double a1 = Math.Pow(Xke / rpmin, 2.0 / 3.0);
        double e = Eccentricity;
        double i = Inclination;
        double temp = 1.5 * Ck2 * ((3.0 * Globals.Sqr(Math.Cos(i))) - 1.0) /
                        Math.Pow(1.0 - (e * e), 1.5);
        double delta1 = temp / (a1 * a1);
        double a0 = a1 *
                       (1.0 - (delta1 *
                       ((1.0 / 3.0) + (delta1 *
                       (1.0 + (134.0 / 81.0 * delta1))))));

        double delta0 = temp / (a0 * a0);

        m_rmMeanMotionRec = rpmin / (1.0 + delta0);
        m_aeAxisSemiMajorRec = a0 / (1.0 - delta0);
        m_aeAxisSemiMinorRec = m_aeAxisSemiMajorRec * Math.Sqrt(1.0 - (e * e));
        Periapsis = Xkmper * ((m_aeAxisSemiMajorRec * (1.0 - e)) - Globals.Ae);
        Apoapsis = Xkmper * ((m_aeAxisSemiMajorRec * (1.0 + e)) - Globals.Ae);

        // Calculate the period using the recovered mean motion.
        if (m_rmMeanMotionRec == 0)
        {
            Period = new TimeSpan(0, 0, 0);
        }
        else
        {
            double secs = Globals.TwoPi / m_rmMeanMotionRec * 60.0;
            int msecs = (int)((secs - (int)secs) * 1000);

            Period = new TimeSpan(0, 0, 0, (int)secs, msecs);
        }

        if (Period.TotalMinutes >= 225.0)
        {
            // SDP4 - period >= 225 minutes.
            NoradModel = new NoradSDP4(this);
        }
        else
        {
            // SGP4 - period < 225 minutes
            NoradModel = new SGP4(this.KeplerianElements, this.Planet);
        }
    }


    /// <summary>
    /// Calculate satellite ECI position/velocity for a given time.
    /// </summary>
    /// <param name="minutes">Target time, in minutes past the Keplerian epoch.</param>
    /// <returns>Kilometer-based position/velocity ECI coordinates.</returns>
    public OrbitalState<double> GetPosition(double minutes)
    {
        var eci = NoradModel.GetPosition(minutes);

        if (NoradModel is not SGP4)
        {
            // Convert ECI vector units from AU to kilometers
            double radiusAe = Planet.Radius / Globals.Ae;  // km
            var velocityScale = radiusAe * (Globals.MinPerDay / 86400.0);// km/sec

            eci = new OrbitalState<double>(eci.Position * radiusAe, eci.Velocity * velocityScale);
        }
        return eci;
    }

    /// <summary>
    /// Calculate ECI position/velocity for a given time.
    /// </summary>
    /// <param name="utc">Target time (UTC).</param>
    /// <returns>Kilometer-based position/velocity ECI coordinates.</returns>
    public OrbitalState<double> GetOrbitalState(DateTime utc)
    {
        return GetPosition(TPlusEpoch(utc).TotalMinutes);
    }



    // ///////////////////////////////////////////////////////////////////////////
    // Returns elapsed time from epoch to given time.
    // Note: "Predicted" TLEs can have epochs in the future.
    private TimeSpan TPlusEpoch(DateTime utc)
    {
        return utc - Epoch;
    }

    #region Utility

    // ///////////////////////////////////////////////////////////////////
    private static double GetRad(double v)
    {
        return Globals.ToRadians(v);
    }

    #endregion
}
