using System.Numerics;

namespace System.Astronomy;

/// <summary>
/// This class encapsulates a single set of standard NORAD two-line element.
/// 
/// <para>Seven numbers are required to define a satellite orbit. This set of seven numbers is called the satellite orbital elements, or sometimes “Keplerian” elements (after Johann Kepler [1571-1630]), or just elements. These numbers define an ellipse, orient it about the earth, and place the satellite on the ellipse at a particular time. In the Keplerian model, satellites orbit in an ellipse of constant shape and orientation. The Earth is at one focus of the ellipse, not the center (unless the orbit ellipse is actually a perfect circle).</para>
/// <para>The real world is slightly more complex than the Keplerian model, and tracking programs compensate for this by introducing minor corrections to the Keplerian model. These corrections are known as perturbations.The perturbations that amateur tracking programs know about are due to the lumpiness of the earth’s gravitational field (which luckily you don’t have to specify), and the “drag” on the satellite due to atmosphere. Drag becomes an optional eighth orbital element.</para>
/// </summary>
public sealed class TwoLineElement<T> : IKeplerianElements<T> where T : INumber<T>, IFloatingPoint<T>
{
    public sealed class SatelliteCatalogNumber
    {
        public readonly int Year;
        public readonly int LaunchNumber;
        public readonly string Piece;

        internal SatelliteCatalogNumber(int year, int launchNumber, string piece)
        {
            Year = year;
            LaunchNumber = launchNumber;
            Piece = piece;
        }
    }

    public readonly string Name;
    public readonly string Id;
    public readonly string Classification;
    public readonly T SetNumber;
    public readonly SatelliteCatalogNumber InternationalDesignator;

    public readonly KeplerianElements<T> KeplerianElements;


    public readonly T MeanMotionFirstDerivation;
    public readonly T MeanMotionSecondDerivation;

    /// <summary>
    /// This tells the tracking program how many times the satellite has orbited from the time
    /// it was launched until the time specified by “Epoch”.
    /// 
    /// <para>Epoch Rev is used to calculate the revolution number displayed by the tracking program.</para>
    /// </summary>
    public readonly int RevolutionNumberAtEpoch;


    #region Column Offsets

    // Note: The column offsets are zero-based.

    // Name
    private const int TLE_LEN_LINE_DATA = 69; private const int TLE_LEN_LINE_NAME = 24;

    // Line 1
    private const int TLE1_COL_SATNUM = 2; private const int TLE1_LEN_SATNUM = 5;
    private const int TLE1_COL_INTLDESC_A = 9; private const int TLE1_LEN_INTLDESC_A = 2;
    private const int TLE1_COL_INTLDESC_B = 11; private const int TLE1_LEN_INTLDESC_B = 3;
    private const int TLE1_COL_INTLDESC_C = 14; private const int TLE1_LEN_INTLDESC_C = 3;
    private const int TLE1_COL_EPOCH_A = 18; private const int TLE1_LEN_EPOCH_A = 2;
    private const int TLE1_COL_EPOCH_B = 20; private const int TLE1_LEN_EPOCH_B = 12;
    private const int TLE1_COL_MEANMOTIONDT = 33; private const int TLE1_LEN_MEANMOTIONDT = 10;
    private const int TLE1_COL_MEANMOTIONDT2 = 44; private const int TLE1_LEN_MEANMOTIONDT2 = 8;
    private const int TLE1_COL_BSTAR = 53; private const int TLE1_LEN_BSTAR = 8;
    private const int TLE1_COL_EPHEMTYPE = 62; private const int TLE1_LEN_EPHEMTYPE = 1;
    private const int TLE1_COL_ELNUM = 64; private const int TLE1_LEN_ELNUM = 4;

    // Line 2
    private const int TLE2_COL_SATNUM = 2; private const int TLE2_LEN_SATNUM = 5;
    private const int TLE2_COL_INCLINATION = 8; private const int TLE2_LEN_INCLINATION = 8;
    private const int TLE2_COL_RAASCENDNODE = 17; private const int TLE2_LEN_RAASCENDNODE = 8;
    private const int TLE2_COL_ECCENTRICITY = 26; private const int TLE2_LEN_ECCENTRICITY = 7;
    private const int TLE2_COL_ARGPERIGEE = 34; private const int TLE2_LEN_ARGPERIGEE = 8;
    private const int TLE2_COL_MEANANOMALY = 43; private const int TLE2_LEN_MEANANOMALY = 8;
    private const int TLE2_COL_MEANMOTION = 52; private const int TLE2_LEN_MEANMOTION = 11;
    private const int TLE2_COL_REVATEPOCH = 63; private const int TLE2_LEN_REVATEPOCH = 5;

    T IKeplerianElements<T>.Epoch => KeplerianElements.Epoch;

    T IKeplerianElements<T>.Eccentricity => KeplerianElements.Eccentricity;

    T IKeplerianElements<T>.MeanMotion => KeplerianElements.MeanMotion;

    T IKeplerianElements<T>.Inclination => KeplerianElements.Inclination;

    T IKeplerianElements<T>.RightAscensionOfAscendingNode => KeplerianElements.RightAscensionOfAscendingNode;

    T IKeplerianElements<T>.MeanAnomaly => KeplerianElements.MeanAnomaly;

    T IKeplerianElements<T>.ArgumentOfPeriapsis => KeplerianElements.ArgumentOfPeriapsis;

    T IKeplerianElements<T>.Drag => KeplerianElements.Drag;

    #endregion

    public TwoLineElement(string name, string classification, string noradNumber, SatelliteCatalogNumber internationalDesignator, T epoch, T meanMotionFirstDerivation, T meanMotionSecondDerivation, T drag, T setNumber, T eccentricity, T inclination, T rightAscensionOfAscendingNode, T argumentOfPeriapsis, T meanAnomaly, T meanMotion, int revolutionNumberAtEpoch)
    {
        Name = name;
        Id = noradNumber;
        Classification = classification;
        SetNumber = setNumber;
        InternationalDesignator = internationalDesignator;
        MeanMotionFirstDerivation = meanMotionFirstDerivation;
        MeanMotionSecondDerivation = meanMotionSecondDerivation;
        RevolutionNumberAtEpoch = revolutionNumberAtEpoch;

        KeplerianElements = new KeplerianElements<T>(epoch, eccentricity, inclination, rightAscensionOfAscendingNode, argumentOfPeriapsis, meanMotion, meanAnomaly, drag);
    }

    public static TwoLineElement<T> Parse(string text)
    {
        var lines = text.Split('\n');
        return Parse(lines[0], lines[1], lines[2]);
    }

    public static TwoLineElement<T> Parse(string[] lines)
    {
        return Parse(lines[0], lines[1], lines[2]);
    }

    public static TwoLineElement<T> Parse(string name, string line1, string line2)
    {
        var noradNumber = line1[TLE1_COL_SATNUM..(TLE1_LEN_SATNUM + TLE1_COL_SATNUM)];
        var classification = line1.Substring(7, 1);

        var intlYY = line1.Substring(TLE1_COL_INTLDESC_A, TLE1_LEN_INTLDESC_A);
        var intlLN = line1.Substring(TLE1_COL_INTLDESC_A + TLE1_LEN_INTLDESC_A, TLE1_LEN_INTLDESC_B);
        var intlPL = line1.Substring(TLE1_COL_INTLDESC_A + TLE1_LEN_INTLDESC_A + TLE1_LEN_INTLDESC_B, TLE1_LEN_INTLDESC_C);

        var idYY = 0;
        if (!string.IsNullOrWhiteSpace(intlYY))
        {
            idYY = int.Parse(intlYY.Trim(), null);
        }
        var idLN = 0;
        if (!string.IsNullOrWhiteSpace(intlLN))
        {
            idLN = int.Parse(intlLN.Trim(), null);
        }

        var idPL = intlPL.Trim();

        var intl = new SatelliteCatalogNumber(idYY, idLN, idPL);


        var epoch = ParseDecimal(line1, TLE1_COL_EPOCH_A, TLE1_LEN_EPOCH_A + TLE1_LEN_EPOCH_B);

        var setNumber = line1.Substring(TLE1_COL_ELNUM, TLE1_LEN_ELNUM).TrimStart();
        var SetNumber = string.IsNullOrWhiteSpace(setNumber) ? T.Zero : T.Parse(setNumber, null);


        var meanMotionDt = "0";
        if (line1[TLE1_COL_MEANMOTIONDT] == '-')
        {
            // value is negative
            meanMotionDt = "-0";
        }

        meanMotionDt += line1.Substring(TLE1_COL_MEANMOTIONDT + 1, TLE1_LEN_MEANMOTIONDT);
        var MeanMotionDt = T.Parse(meanMotionDt, null);

        // T point assumed; exponential notation
        var MeanMotionDt2 = ExpToDecimal(line1.Substring(TLE1_COL_MEANMOTIONDT2, TLE1_LEN_MEANMOTIONDT2));


        // T point assumed; exponential notation
        var BStarDrag = ExpToDecimal(line1.Substring(TLE1_COL_BSTAR, TLE1_LEN_BSTAR));

        var Eccentricity = T.Parse("0." + line2.Substring(TLE2_COL_ECCENTRICITY, TLE2_LEN_ECCENTRICITY), null);


        var Inclination = T.Parse(line2.Substring(TLE2_COL_INCLINATION, TLE2_LEN_INCLINATION).TrimStart(), null);



        var raan = line2.Substring(TLE2_COL_RAASCENDNODE, TLE2_LEN_RAASCENDNODE).TrimStart();
        var Raan = T.Parse(raan, null);

        var argPerigee = line2.Substring(TLE2_COL_ARGPERIGEE, TLE2_LEN_ARGPERIGEE).TrimStart();
        var ArgPerigee = T.Parse(argPerigee, null);

        var meanAnomaly = line2.Substring(TLE2_COL_MEANANOMALY, TLE2_LEN_MEANANOMALY).TrimStart();
        var MeanAnomaly = T.Parse(meanAnomaly, null);

        var meanMotion = line2.Substring(TLE2_COL_MEANMOTION, TLE2_LEN_MEANMOTION).TrimStart();
        var MeanMotion = T.Parse(meanMotion, null);

        var revAtEpoch = line2.Substring(TLE2_COL_REVATEPOCH, TLE2_LEN_REVATEPOCH).TrimStart();

        var RevAtEpoch = string.IsNullOrEmpty(revAtEpoch) ? 0 : int.Parse(revAtEpoch, null);

        return new TwoLineElement<T>(name, classification, noradNumber, intl, epoch, MeanMotionDt, MeanMotionDt2, BStarDrag, SetNumber, Eccentricity, Inclination, Raan, ArgPerigee, MeanAnomaly, MeanMotion, RevAtEpoch);
    }

    private static T ParseDecimal(string str, int start, int end)
    {
        var sub = str.Substring(start, end);

        return T.Parse(sub, null);
    }

    private static int ParseInt(string str, int start, int end)
    {
        var sub = str.Substring(start, end);

        return int.Parse(sub, null);
    }

    /// <summary>
    /// Converts TLE-style exponential notation of the form 
    ///       [ |+|-]00000[ |+|-]0
    /// to T notation. Assumes implied T point to the left
    /// of the first number in the string, i.e., 
    ///       " 12345-3" =  0.00012345
    ///       "-23429-5" = -0.0000023429   
    ///       " 40436+1" =  4.0436
    /// No sign character implies a positive value, i.e.,
    ///       " 00000 0" =  0.00000
    ///       " 31415 1" =  3.1415
    /// </summary>
    private static T ExpToDecimal(string str)
    {
        const int LEN_SIGN = 1;

        const int COL_MANTISSA = 1;
        const int LEN_MANTISSA = 5;

        const int COL_EXPONENT = 6;
        const int LEN_EXPONENT = 2;

        string sign = str[..LEN_SIGN];
        string mantissa = str.Substring(COL_MANTISSA, LEN_MANTISSA);
        string exponent = str.Substring(COL_EXPONENT, LEN_EXPONENT).TrimStart();

        var val = T.Parse(sign + "0." + mantissa + "e" + exponent, style: System.Globalization.NumberStyles.AllowExponent | System.Globalization.NumberStyles.Float, null);


        return val;
    }
}
