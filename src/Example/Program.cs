using System.Astronomy;

var a1 = @"STARLINK-30591
1 56003U 23042T   24037.86834427  .00001215  00000-0  10499-3 0  9998
2 56003  43.0022 235.5557 0001246 278.4576  81.6125 15.02545537 48317";

//var a1 = "STARLINK-31155";
//var a2 = "1 58728U 24005A   24036.91667824  .00091608  00000+0  36484-2 0  9999";
//var a3 = "2 58728  43.0002  49.7284 0001969 257.5625   0.2160 15.24223380  5747";

var s2763 = @"STARLINK-2763
1 48670U 21044AJ  24036.44046090  .00001946  00000+0  14948-3 0  9996
2 48670  53.0543 297.9783 0001510  99.0250 261.0910 15.06413955148803";

/*
STARLINK-31155          
1 58728U 24005A   24036.91667824  .00091608  00000+0  36484-2 0  9999
2 58728  43.0002  49.7284 0001969 257.5625   0.2160 15.24223380  5747
 */

var tlex = TwoLineElement<double>.Parse(a1);

var sgp = new SGP4(tlex, Earth.WGS84);

//var now = new DateTime(2024, 2, 6, 20, 54, 45, DateTimeKind.Utc);
var now = new DateTime(2024, 2, 7, 10, 3, 0, DateTimeKind.Utc);
var observer = new LatLongAlt<double>(Angle<double>.FromDegrees(35.764472), Angle<double>.FromDegrees(50.786492), 1185.9);

var rv = sgp.GetPosition(now);
var ecf = rv.Position.ToEcf(now);
var tp = Earth.WGS84.Topocentric(observer, ecf);
var la = tp.ToLookAngle();
var el = la.Elevation.Degrees;
var az = la.Azimuth.Degrees;

string str1 = "SGP4 Test";
string str2 = "1 88888U          80275.98708465  .00073094  13844-3  66816-4 0     8";
string str3 = "2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518   105";

var tle1 = TwoLineElement<double>.Parse(str1, str2, str3);

PrintPosVel(tle1);

Console.WriteLine();

// Test SDP4
str1 = "SDP4 Test";
str2 = "1 11801U          80230.29629788  .01431103  00000-0  14311-1       8";
str3 = "2 11801  46.7916 230.4354 7318036  47.4722  10.4117  2.28537848     6";

var tle2 = TwoLineElement<double>.Parse(str1, str2, str3);

PrintPosVel(tle2);

Console.WriteLine("\nExample output:");

// Example: Define a location on the earth, then determine the look-angle
// to the SDP4 satellite defined above.

// Get the location of the satellite from the Orbit object. The 
// earth-centered inertial information is placed into eciSDP4.
// Here we ask for the location of the satellite 90 minutes after
// the TLE epoch.


// Now create a site object. Site objects represent a location on the 
// surface of the earth. Here we arbitrarily select a point on the
// equator.
// 0.00 N, 100.00 W, 0 km altitude

var eci = sgp.GetPosition(now);

var zone = CalculZoneVisibilite(51, 35, eci.Position);

var ll = Earth.WGS84.GetLookAngle(observer, eci.Position, now);

var c = zone;

//var when = satellite.Epoch.AddDays(900);
//// Now get the "look angle" from the site to the satellite. 
//var nn = Earth.WGS84.GetLookAngle(observer, satellite, now);
//var topoLook = Earth.WGS84.GetLookAngle(observer, satellite, when);
//var topoLook1 = Earth.WGS84.GetLookAngle(observer, satellite, when.AddMinutes(-5));
//var topoLook2 = Earth.WGS84.GetLookAngle(observer, satellite, when.AddMinutes(5));

//// Print out the results. Note that the Azimuth and Elevation are
//// stored in the Topocentric object as radians. Here we convert
//// to degrees using 180 / Math.PI.
//Console.WriteLine($"AZ: {topoLook1.Azimuth * 180 / Math.PI:f3}  EL: {topoLook1.Elevation * 180 / Math.PI:f3}");
//Console.WriteLine($"AZ: {topoLook.Azimuth * 180 / Math.PI:f3}  EL: {topoLook.Elevation * 180 / Math.PI:f3}");
//Console.WriteLine($"AZ: {topoLook2.Azimuth * 180 / Math.PI:f3}  EL: {topoLook2.Elevation * 180 / Math.PI:f3}");


void PrintPosVel(IKeplerianElements<double> tle)
{
    const int Step = 360;

    var sat = new Orbit(tle, Earth.WGS72);
    var coords = new List<OrbitalState<double>>();

    // Calculate position, velocity
    // mpe = "minutes past epoch"
    for (int mpe = 0; mpe <= (Step * 4); mpe += Step)
    {
        // Get the position of the satellite at time "mpe".
        // The coordinates are placed into the variable "eci".
        var eci = sat.GetPosition(mpe);

        // Add the coordinate object to the list
        coords.Add(eci);
    }

    // Print TLE data
    //Console.Write("{0}\n", tle.Name);
    //Console.Write("{0}\n", tle.Line1);
    //Console.Write("{0}\n", tle.Line2);

    // Header
    Console.Write("\n  TSINCE            X                Y                Z\n\n");

    // Iterate over each of the ECI position objects pushed onto the
    // coordinate list, above, printing the ECI position information
    // as we go.
    for (int i = 0; i < coords.Count; i++)
    {
        var e = coords[i];

        Console.Write("{0,8}.00 {1,16:f8} {2,16:f8} {3,16:f8}\n",
                      i * Step,
                      e.Position.X,
                      e.Position.Y,
                      e.Position.Z);
    }

    Console.Write("\n                  XDOT             YDOT             ZDOT\n\n");

    // Iterate over each of the ECI position objects in the coordinate
    // list again, but this time print the velocity information.
    for (int i = 0; i < coords.Count; i++)
    {
        var e = coords[i];

        Console.Write("{0,24:f8} {1,16:f8} {2,16:f8}\n",
                      e.Velocity.X,
                      e.Velocity.Y,
                      e.Velocity.Z);
    }
}

LatLongAlt<double>[] CalculZoneVisibilite(double longitude, double latitude, EarthCenteredInertial<double> position)
{
    var _zone = new LatLongAlt<double>[361];

    double num = longitude;
    if (num > 0.0)
    {
        num -= Math.PI * 2.0;
    }
    double num2 = Math.Cos(latitude);
    double num3 = Math.Sin(latitude);
    double num4 = Math.Acos(6363.136658 / position.Length());
    if (double.IsNaN(num4))
    {
        num4 = 0.0;
    }
    double num5 = Math.Cos(num4);
    double num6 = Math.Sin(num4);
    int i = 0;
    do
    {
        double d = Math.PI / 180.0 * i;
        double lat = Math.Asin(num3 * num5 + Math.Cos(d) * num6 * num2);
        double num9 = (num5 - num3 * Math.Sin(lat)) / (num2 * Math.Cos(lat));
        double lng = (((i != 0 || !(num4 > Math.PI / 2.0 - latitude)) && 0 == 0) ? (((i == 180 && num4 > Math.PI / 2.0 + latitude) ? true : false) ? (num + Math.PI) : ((Math.Abs(num9) > 1.0) ? num : ((i > 180) ? (num - Math.Acos(num9)) : (num + Math.Acos(num9))))) : (num + Math.PI));

        const double twoPi = Math.PI * 2.0;

        var finalLng = (Math.PI - lng) % twoPi;
        var finalLat = (Math.PI / 2.0 - lat);

        var z = new LatLongAlt<double>(Angle<double>.FromRadians(finalLat), Angle<double>.FromRadians(finalLng), 0);

        _zone[i] = z;

        i++;
    }


    while (i <= 359);

    _zone[360] = new LatLongAlt<double>(_zone[0].Latitude, _zone[0].Longitude, _zone[0].Altitude);

    return _zone;
}