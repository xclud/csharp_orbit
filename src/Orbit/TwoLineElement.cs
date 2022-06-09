using System.Numerics;

namespace System;

/// <summary>
/// This class encapsulates a single set of standard NORAD two-line element.
/// 
/// <para>Seven numbers are required to define a satellite orbit. This set of seven numbers is called the satellite orbital elements, or sometimes “Keplerian” elements (after Johann Kepler [1571-1630]), or just elements. These numbers define an ellipse, orient it about the earth, and place the satellite on the ellipse at a particular time. In the Keplerian model, satellites orbit in an ellipse of constant shape and orientation. The Earth is at one focus of the ellipse, not the center (unless the orbit ellipse is actually a perfect circle).</para>
/// <para>The real world is slightly more complex than the Keplerian model, and tracking programs compensate for this by introducing minor corrections to the Keplerian model. These corrections are known as perturbations.The perturbations that amateur tracking programs know about are due to the lumpiness of the earth’s gravitational field (which luckily you don’t have to specify), and the “drag” on the satellite due to atmosphere. Drag becomes an optional eighth orbital element.</para>
/// </summary>
public readonly ref struct TwoLineElement<T> where T : INumber<T>
{

    public class SatelliteCatalogNumber
    {
        public int Year;
        public int LaunchNumber;
        public string Piece;

        internal SatelliteCatalogNumber(int year, int launchNumber, string piece)
        {
            Year = year;
            LaunchNumber = launchNumber;
            Piece = piece;
        }
    }

    public readonly string Name;
    public readonly T NoradNumber;
    public readonly string Classification;
    public readonly T SetNumber;
    public readonly SatelliteCatalogNumber InternationalDesignator;
    public readonly T Epoch;

    public readonly T MeanMotionFirstDerivation;
    public readonly T MeanMotionSecondDerivation;

    public readonly T BStarDrag;

    /// <summary>
    /// Eccentricity in the range 0 <= e < 1.
    /// 
    /// <para>In the Keplerian orbit model, the satellite orbit is an ellipse.
    /// Eccentricity tells us the "shape" of the ellipse. When e = 0, the ellipse is a circle.
    /// When e is very near 1, the ellipse is very long and skinny.
    /// </para>
    /// <para>To be precise, the Keplerian orbit is a conic section, which can be either an ellipse, which includes circles, a parabola, a hyperbola, or a straight line.
    /// But here, we are only interested in elliptical orbits.
    /// The other kinds of orbits are not used for satellites, at least not on purpose, and tracking programs typically aren't programmed to handle them.
    /// </para>
    /// <para>For our purposes, eccentricity must be in the range 0 <= e < 1.</para>
    /// </summary>
    public readonly T Eccentricity;

    public readonly T MeanMotion;

    /// <summary>
    /// This tells the tracking program how many times the satellite has orbited from the time
    /// it was launched until the time specified by “Epoch”.
    /// 
    /// <para>Epoch Rev is used to calculate the revolution number displayed by the tracking program.</para>
    /// </summary>
    public readonly T RevolutionNumberAtEpoch;


    /// <summary>
    /// The orbit ellipse lies in a plane known as the orbital plane. The orbital plane always goes through the center of the earth,
    /// but may be tilted any angle relative to the equator.
    /// Inclination is the angle between the orbital plane and the equatorial plane.
    /// By convention, inclination is a number between 0 and 180 degrees.
    /// <para>Orbits with inclination near 0 degrees are called equatorial orbits (because the satellite stays nearly over the equator).</para>
    /// <para>Orbits with inclination near 90 degrees are called polar (because the satellite crosses over the north and south poles).</para>
    /// <para>The intersection of the equatorial plane and the orbital plane is a line which is called the line of nodes.</para>
    /// </summary>
    public readonly T Inclination;

    /// <summary>
    /// Right Ascension of the Ascending Node (RAAN) In Degrees.
    /// <para>Two numbers orient the orbital plane in space. The first number was Inclination. This is the second. After we’ve specified inclination, there are still an infinite number of orbital planes possible. The line of nodes can poke out the anywhere along the equator. If we specify where along the equator the line of nodes pokes out, we will have the orbital plane fully specified. The line of nodes pokes out two places, of course. We only need to specify one of them. One is called the ascending node (where the satellite crosses the equator going from south to north). The other is called the descending node (where the satellite crosses the equator going from north to south). By convention, we specify the location of the ascending node.</para>
    /// <para>Now, the earth is spinning.This means that we can’t use the common latitude/longitude coordinate system to specify where the line of nodes points.Instead, we use an astronomical coordinate system, known as the right ascension / declination coordinate system, which does not spin with the earth.Right ascension is another fancy word for an angle, in this case, an angle measured in the equatorial plane from a reference point in the sky where right ascension is defined to be zero.Astronomers call this point the vernal equinox..</para>
    /// <para>Finally, “right ascension of ascending node” is an angle, measured at the center of the earth, from the vernal equinox to the ascending node.</para>
    /// <para>I know this is getting complicated. Here’s an example.Draw a line from the center of the earth to the point where our satellite crosses the equator (going from south to north). If this line points directly at the vernal equinox, then RAAN = 0 degrees..</para>
    /// <para>By convention, RAAN is a number in the range 0 to 360 degrees.</para>
    /// <para>I used the term “vernal equinox” above without really defining it.If you can tolerate a minor digression, I’ll do that now. Teachers have told children for years that the vernal equinox is “the place in the sky where the sun rises on the first day of Spring”. This is a horrible definition.Most teachers, and students, have no idea what the first day of spring is (except a date on a calendar), and no idea why the sun should be in the same place in the sky on that date every year.</para>
    /// <para>You now have enough astronomy vocabulary to get a better definition. Consider the orbit of the sun around the earth. I know in school they told you the earth orbits around the sun, but the math is equally valid either way, and it suits our needs at this instant to think of the sun orbiting the earth. The orbit of the sun has an inclination of about 23.5 degrees. (Astronomers don’t usually call this 23.5 degree angle an ‘inclination’, by the way.They use an infinitely more obscure name: The Obliquity of The Ecliptic.) The orbit of the sun is divided (by humans) into four equally sized portions called seasons.The one called Spring begins when the sun pops up past the equator. In other words, the first day of Spring is the day that the sun crosses through the equatorial plane going from South to North.We have a name for that! It’s the ascending node of the Sun’s orbit. So finally, the vernal equinox is nothing more than the ascending node of the Sun’s orbit. The Sun’s orbit has RAAN = 0 simply because we’ve defined the Sun’s ascending node as the place from which all ascending nodes are measured.The RAAN of your satellite’s orbit is just the angle (measured at the center of the earth) between the place the Sun’s orbit pops up past the equator, and the place your satellite’s orbit pops up past the equator.</para>
    /// </summary>
    public readonly T RightAscensionOfAscendingNode;


    /// <summary>
    /// Mean Anomaly in Degrees.
    /// Now that we have the size, shape, and orientation of the orbit firmly established, the only thing left to do is specify where exactly the satellite is on this orbit ellipse at some particular time. Our very first orbital element (Epoch) specified a particular time, so all we need to do now is specify where, on the ellipse, our satellite was exactly at the Epoch time.
    /// Anomaly is yet another astronomer-word for angle.Mean anomaly is simply an angle that marches uniformly in time from 0 to 360 degrees during one revolution.It is defined to be 0 degrees at perigee, and therefore is 180 degrees at apogee.
    /// If you had a satellite in a circular orbit (therefore moving at constant speed) and you stood in the center of the earth and measured this angle from perigee, you would point directly at the satellite.Satellites in non-circular orbits move at a non-constant speed, so this simple relation doesn’t hold. This relation does hold for two important points on the orbit, however, no matter what the eccentricity.Perigee always occurs at MA = 0, and apogee always occurs at MA = 180 degrees.
    /// It has become common practice with radio amateur satellites to use Mean Anomaly to schedule satellite operations.Satellites commonly change modes or turn on or off at specific places in their orbits, specified by Mean Anomaly. Unfortunately, when used this way, it is common to specify MA in units of 256ths of a circle instead of degrees! Some tracking programs use the term “phase” when they display MA in these units. It is still specified in degrees, between 0 and 360, when entered as an orbital element.
    /// Example: Suppose Oscar-99 has a period of 12 hours, and is turned off from Phase 240 to 16. That means it’s off for 32 ticks of phase.There are 256 of these ticks in the entire 12 hour orbit, so it’s off for (32/256)x12hrs = 1.5 hours. Note that the off time is centered on perigee. Satellites in highly eccentric orbits are often turned off near perigee when they’re moving the fastest, and therefore difficult to use.
    /// </summary>
    public readonly T MeanAnomaly;

    /// <summary>
    /// Argument of Periapsis In Degrees.
    /// 
    /// <para>Now that we’ve oriented the orbital plane in space, we need to orient the orbit ellipse in the orbital plane. We do this by specifying a single angle known as argument of perigee.</para>
    /// <para>A few words about elliptical orbits… The point where the satellite is closest to the earth is called perigee, although it’s sometimes called periapsis or perifocus.We’ll call it perigee. The point where the satellite is farthest from earth is called apogee (aka apoapsis, or apifocus). If we draw a line from perigee to apogee, this line is called the line-of-apsides. (Apsides is, of course, the plural of apsis.) I know, this is getting complicated again.Sometimes the line-of-apsides is called the major-axis of the ellipse.It’s just a line drawn through the ellipse the “long way”.</para>
    /// <para>The line-of-apsides passes through the center of the earth. We’ve already identified another line passing through the center of the earth: the line of nodes. The angle between these two lines is called the argument of perigee.Where any two lines intersect, they form two supplementary angles, so to be specific, we say that argument of perigee is the angle (measured at the center of the earth) from the ascending node to perigee.</para>
    /// <para>Example: When ARGP = 0, the perigee occurs at the same place as the ascending node.That means that the satellite would be closest to earth just as it rises up over the equator. When ARGP = 180 degrees, apogee would occur at the same place as the ascending node.That means that the satellite would be farthest from earth just as it rises up over the equator.</para>
    /// <para>By convention, ARGP is an angle between 0 and 360 degrees.</para>
    /// </summary>
    public readonly T ArgumentOfPeriapsis;

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

    #endregion

    public TwoLineElement(string name, string classification, T noradNumber, SatelliteCatalogNumber internationalDesignator, T epoch, T meanMotionFirstDerivation, T meanMotionSecondDerivation, T bStarDrag, T setNumber, T eccentricity, T inclination, T rightAscensionOfAscendingNode, T argumentOfPeriapsis, T meanAnomaly, T meanMotion, T revolutionNumberAtEpoch)
    {
        Name = name;
        NoradNumber = noradNumber;
        Classification = classification;
        SetNumber = setNumber;
        InternationalDesignator = internationalDesignator;
        Epoch = epoch;
        MeanMotionFirstDerivation = meanMotionFirstDerivation;
        MeanMotionSecondDerivation = meanMotionSecondDerivation;
        BStarDrag = bStarDrag;
        Eccentricity = eccentricity;
        RightAscensionOfAscendingNode = rightAscensionOfAscendingNode;
        ArgumentOfPeriapsis = argumentOfPeriapsis;
        MeanAnomaly = meanAnomaly;
        MeanMotion = meanMotion;

        RevolutionNumberAtEpoch = revolutionNumberAtEpoch;
        Inclination = inclination;
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
        var noradNumber = ParseDecimal(line1, TLE1_COL_SATNUM, TLE1_LEN_SATNUM);
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

        var RevAtEpoch = string.IsNullOrEmpty(revAtEpoch) ? T.Zero : T.Parse(revAtEpoch, null);

        return new TwoLineElement<T>(name, classification, noradNumber, intl, epoch, MeanMotionDt, MeanMotionDt2, BStarDrag, SetNumber, Eccentricity, Inclination, Raan, ArgPerigee, MeanAnomaly, MeanMotion, RevAtEpoch);
    }

    private static T ParseDecimal(string str, int start, int end)
    {
        var sub = str.Substring(start, end);

        return T.Parse(sub, null);
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
        const int COL_SIGN = 0;
        const int LEN_SIGN = 1;

        const int COL_MANTISSA = 1;
        const int LEN_MANTISSA = 5;

        const int COL_EXPONENT = 6;
        const int LEN_EXPONENT = 2;

        string sign = str.Substring(COL_SIGN, LEN_SIGN);
        string mantissa = str.Substring(COL_MANTISSA, LEN_MANTISSA);
        string exponent = str.Substring(COL_EXPONENT, LEN_EXPONENT).TrimStart();

        var val = T.Parse(sign + "0." + mantissa + "e" + exponent, style: System.Globalization.NumberStyles.AllowExponent | System.Globalization.NumberStyles.Float, null);


        return val;
    }
}
