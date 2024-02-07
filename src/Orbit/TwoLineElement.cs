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

    T IKeplerianElements<T>.Epoch => KeplerianElements.Epoch;

    T IKeplerianElements<T>.Eccentricity => KeplerianElements.Eccentricity;

    T IKeplerianElements<T>.MeanMotion => KeplerianElements.MeanMotion;

    Angle<T> IKeplerianElements<T>.Inclination => KeplerianElements.Inclination;

    Angle<T> IKeplerianElements<T>.RightAscensionOfAscendingNode => KeplerianElements.RightAscensionOfAscendingNode;

    Angle<T> IKeplerianElements<T>.MeanAnomaly => KeplerianElements.MeanAnomaly;

    Angle<T> IKeplerianElements<T>.ArgumentOfPeriapsis => KeplerianElements.ArgumentOfPeriapsis;

    T IKeplerianElements<T>.Drag => KeplerianElements.Drag;


    public TwoLineElement(string name, string classification, string noradNumber, SatelliteCatalogNumber internationalDesignator, T epoch, T meanMotionFirstDerivation, T meanMotionSecondDerivation, T drag, T setNumber, T eccentricity, Angle<T> inclination, Angle<T> rightAscensionOfAscendingNode, Angle<T> argumentOfPeriapsis, Angle<T> meanAnomaly, T meanMotion, int revolutionNumberAtEpoch)
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
        var noradNumber = line1[2..7];
        var classification = line1.Substring(7, 1);

        var intlYY = line1.Substring(9, 2);
        var intlLN = line1.Substring(11, 3);
        var intlPL = line1.Substring(14, 3);

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


        var epoch = ParseDecimal(line1, 18, 2 + 12);

        var setNumber = line1.Substring(64, 4).TrimStart();
        var SetNumber = string.IsNullOrWhiteSpace(setNumber) ? T.Zero : T.Parse(setNumber, null);


        var meanMotionDt = "0";
        if (line1[33] == '-')
        {
            // value is negative
            meanMotionDt = "-0";
        }

        meanMotionDt += line1.Substring(33 + 1, 10);
        var MeanMotionDt = T.Parse(meanMotionDt, null);

        // T point assumed; exponential notation
        var MeanMotionDt2 = ExpToDecimal(line1.Substring(44, 8));


        // T point assumed; exponential notation
        var BStarDrag = ExpToDecimal(line1.Substring(53, 8));

        var Eccentricity = T.Parse("0." + line2.Substring(26, 7), null);


        var Inclination = T.Parse(line2.Substring(8, 8).TrimStart(), null);



        var raan = line2.Substring(17, 8).TrimStart();
        var Raan = T.Parse(raan, null);

        var argPerigee = line2.Substring(34, 8).TrimStart();
        var ArgPerigee = T.Parse(argPerigee, null);

        var meanAnomaly = line2.Substring(43, 8).TrimStart();
        var MeanAnomaly = T.Parse(meanAnomaly, null);

        var meanMotion = line2.Substring(52, 11).TrimStart();
        var MeanMotion = T.Parse(meanMotion, null);

        var revAtEpoch = line2.Substring(63, 5).TrimStart();

        var RevAtEpoch = string.IsNullOrEmpty(revAtEpoch) ? 0 : int.Parse(revAtEpoch, null);

        return new TwoLineElement<T>(name, classification, noradNumber, intl, epoch, MeanMotionDt, MeanMotionDt2, BStarDrag, SetNumber, Eccentricity, Angle<T>.FromDegrees(Inclination), Angle<T>.FromDegrees(Raan), Angle<T>.FromDegrees(ArgPerigee), Angle<T>.FromDegrees(MeanAnomaly), MeanMotion, RevAtEpoch);
    }

    private static T ParseDecimal(string str, int start, int length)
    {
        var sub = str.Substring(start, length);

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
