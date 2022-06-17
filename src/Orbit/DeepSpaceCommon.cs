using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace System;


/*-----------------------------------------------------------------------------
*
*                           procedure dscom
*
*  this procedure provides deep space common items used by both the secular
*    and periodics subroutines. input is provided as shown. this routine
*    used to be called dpper, but the functions inside weren't well organized.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    epoch       -
*    ep          - eccentricity
*    argpp       - argument of perigee
*    tc          -
*    inclp       - inclination
*    nodep       - right ascension of ascending node
*    np          - mean motion
*
*  outputs       :
*    sinim  , cosim  , sinomm , cosomm , snodm  , cnodm
*    day         -
*    e3          -
*    ee2         -
*    em          - eccentricity
*    emsq        - eccentricity squared
*    gam         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    rtemsq      -
*    se2, se3         -
*    sgh2, sgh3, sgh4        -
*    sh2, sh3, si2, si3, sl2, sl3, sl4         -
*    s1, s2, s3, s4, s5, s6, s7          -
*    ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3         -
*    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
*    xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4         -
*    nm          - mean motion
*    z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33         -
*    zmol        -
*    zmos        -
*
*  locals        :
*    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10         -
*    betasq      -
*    cc          -
*    ctem, stem        -
*    x1, x2, x3, x4, x5, x6, x7, x8          -
*    xnodce      -
*    xnoi        -
*    zcosg  , zsing  , zcosgl , zsingl , zcosh  , zsinh  , zcoshl , zsinhl ,
*    zcosi  , zsini  , zcosil , zsinil ,
*    zx          -
*    zy          -
*
*  coupling      :
*    none.
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/


/// <summary>
/// Provides deep space common items used by both the secular and periodics subroutines.
/// </summary>
internal sealed class DeepSpaceCommon
{
    private const double twoPi = Math.PI * 2;
    private const double pi = Math.PI;
    private const double deg2rad = pi / 180.0;
    private const double rad2deg = 180 / pi;
    private const double x2o3 = 2.0 / 3.0;

    internal DeepSpaceCommon(double epoch, double ep, double argpp, double tc, double inclp, double nodep, double np)
    {
        var ss1 = 0.0;
        var ss2 = 0.0;
        var ss3 = 0.0;
        var ss4 = 0.0;
        var ss5 = 0.0;
        var ss6 = 0.0;
        var ss7 = 0.0;
        var sz1 = 0.0;
        var sz2 = 0.0;
        var sz3 = 0.0;

        var sz11 = 0.0;
        var sz12 = 0.0;
        var sz13 = 0.0;
        var sz21 = 0.0;
        var sz22 = 0.0;
        var sz23 = 0.0;
        var sz31 = 0.0;
        var sz32 = 0.0;
        var sz33 = 0.0;

        var s1 = 0.0;
        var s2 = 0.0;
        var s3 = 0.0;
        var s4 = 0.0;
        var s5 = 0.0;
        var s6 = 0.0;
        var s7 = 0.0;
        var z1 = 0.0;
        var z2 = 0.0;
        var z3 = 0.0;

        var z11 = 0.0;
        var z12 = 0.0;
        var z13 = 0.0;
        var z21 = 0.0;
        var z22 = 0.0;
        var z23 = 0.0;
        var z31 = 0.0;
        var z32 = 0.0;
        var z33 = 0.0;

        // -------------------------- constants -------------------------
        const double zes = 0.01675;
        const double zel = 0.05490;
        const double c1ss = 2.9864797e-6;
        const double c1l = 4.7968065e-7;
        const double zsinis = 0.39785416;
        const double zcosis = 0.91744867;
        const double zcosgs = 0.1945905;
        const double zsings = -0.98088458;

        //  --------------------- local variables ------------------------
        double nm = np;
        double em = ep;
        double snodm = Math.Sin(nodep);
        double cnodm = Math.Cos(nodep);
        double sinomm = Math.Sin(argpp);
        double cosomm = Math.Cos(argpp);
        double sinim = Math.Sin(inclp);
        double cosim = Math.Cos(inclp);
        double emsq = em * em;
        double betasq = 1.0 - emsq;
        double rtemsq = Math.Sqrt(betasq);

        //  ----------------- initialize lunar solar terms ---------------
        const double peo = 0.0;
        const double pinco = 0.0;
        const double plo = 0.0;
        const double pgho = 0.0;
        const double pho = 0.0;
        var day = epoch + 18261.5 + (tc / 1440.0);
        var xnodce = (4.5236020 - (9.2422029e-4 * day)) % twoPi;
        var stem = Math.Sin(xnodce);
        var ctem = Math.Cos(xnodce);
        var zcosil = 0.91375164 - (0.03568096 * ctem);
        var zsinil = Math.Sqrt(1.0 - (zcosil * zcosil));
        var zsinhl = 0.089683511 * stem / zsinil;
        var zcoshl = Math.Sqrt(1.0 - (zsinhl * zsinhl));
        var gam = 5.8351514 + (0.0019443680 * day);
        var zx = 0.39785416 * stem / zsinil;
        var zy = (zcoshl * ctem) + (0.91744867 * zsinhl * stem);
        zx = Math.Atan2(zx, zy);
        zx += gam - xnodce;
        var zcosgl = Math.Cos(zx);
        var zsingl = Math.Sin(zx);

        //  ------------------------- do solar terms ---------------------
        var zcosg = zcosgs;
        var zsing = zsings;
        var zcosi = zcosis;
        var zsini = zsinis;
        var zcosh = cnodm;
        var zsinh = snodm;
        var cc = c1ss;
        double xnoi = 1.0 / nm;

        var lsflg = 0;
        while (lsflg < 2)
        {
            lsflg += 1;
            var a1 = (zcosg * zcosh) + (zsing * zcosi * zsinh);
            var a3 = (-zsing * zcosh) + (zcosg * zcosi * zsinh);
            var a7 = (-zcosg * zsinh) + (zsing * zcosi * zcosh);
            var a8 = zsing * zsini;
            var a9 = (zsing * zsinh) + (zcosg * zcosi * zcosh);
            var a10 = zcosg * zsini;
            var a2 = (cosim * a7) + (sinim * a8);
            var a4 = (cosim * a9) + (sinim * a10);
            var a5 = (-sinim * a7) + (cosim * a8);
            var a6 = (-sinim * a9) + (cosim * a10);

            var x1 = (a1 * cosomm) + (a2 * sinomm);
            var x2 = (a3 * cosomm) + (a4 * sinomm);
            var x3 = (-a1 * sinomm) + (a2 * cosomm);
            var x4 = (-a3 * sinomm) + (a4 * cosomm);
            var x5 = a5 * sinomm;
            var x6 = a6 * sinomm;
            var x7 = a5 * cosomm;
            var x8 = a6 * cosomm;

            z31 = (12.0 * x1 * x1) - (3.0 * x3 * x3);
            z32 = (24.0 * x1 * x2) - (6.0 * x3 * x4);
            z33 = (12.0 * x2 * x2) - (3.0 * x4 * x4);

            z1 = (3.0 * ((a1 * a1) + (a2 * a2))) + (z31 * emsq);
            z2 = (6.0 * ((a1 * a3) + (a2 * a4))) + (z32 * emsq);
            z3 = (3.0 * ((a3 * a3) + (a4 * a4))) + (z33 * emsq);

            z11 = (-6.0 * a1 * a5)
              + (emsq * ((-24.0 * x1 * x7) - (6.0 * x3 * x5)));
            z12 = (-6.0 * ((a1 * a6) + (a3 * a5)))
              + (emsq * ((-24.0 * ((x2 * x7) + (x1 * x8))) + (-6.0 * ((x3 * x6) + (x4 * x5)))));

            z13 = (-6.0 * a3 * a6)
              + (emsq * ((-24.0 * x2 * x8) - (6.0 * x4 * x6)));

            z21 = (6.0 * a2 * a5)
              + (emsq * ((24.0 * x1 * x5) - (6.0 * x3 * x7)));
            z22 = (6.0 * ((a4 * a5) + (a2 * a6)))
              + (emsq * ((24.0 * ((x2 * x5) + (x1 * x6))) - (6.0 * ((x4 * x7) + (x3 * x8)))));
            z23 = (6.0 * a4 * a6)
              + (emsq * ((24.0 * x2 * x6) - (6.0 * x4 * x8)));

            z1 = z1 + z1 + (betasq * z31);
            z2 = z2 + z2 + (betasq * z32);
            z3 = z3 + z3 + (betasq * z33);
            s3 = cc * xnoi;
            s2 = -0.5 * s3 / rtemsq;
            s4 = s3 * rtemsq;
            s1 = -15.0 * em * s4;
            s5 = (x1 * x3) + (x2 * x4);
            s6 = (x2 * x3) + (x1 * x4);
            s7 = (x2 * x4) - (x1 * x3);

            //  ----------------------- do lunar terms -------------------
            if (lsflg == 1)
            {
                ss1 = s1;
                ss2 = s2;
                ss3 = s3;
                ss4 = s4;
                ss5 = s5;
                ss6 = s6;
                ss7 = s7;
                sz1 = z1;
                sz2 = z2;
                sz3 = z3;
                sz11 = z11;
                sz12 = z12;
                sz13 = z13;
                sz21 = z21;
                sz22 = z22;
                sz23 = z23;
                sz31 = z31;
                sz32 = z32;
                sz33 = z33;
                zcosg = zcosgl;
                zsing = zsingl;
                zcosi = zcosil;
                zsini = zsinil;
                zcosh = (zcoshl * cnodm) + (zsinhl * snodm);
                zsinh = (snodm * zcoshl) - (cnodm * zsinhl);
                cc = c1l;
            }
        }

        double zmol = (4.7199672 + ((0.22997150 * day) - gam)) % twoPi;
        double zmos = (6.2565837 + (0.017201977 * day)) % twoPi;

        //  ------------------------ do solar terms ----------------------
        double se2 = 2.0 * ss1 * ss6;
        double se3 = 2.0 * ss1 * ss7;
        double si2 = 2.0 * ss2 * sz12;
        double si3 = 2.0 * ss2 * (sz13 - sz11);
        double sl2 = -2.0 * ss3 * sz2;
        double sl3 = -2.0 * ss3 * (sz3 - sz1);
        double sl4 = -2.0 * ss3 * (-21.0 - (9.0 * emsq)) * zes;
        double sgh2 = 2.0 * ss4 * sz32;
        double sgh3 = 2.0 * ss4 * (sz33 - sz31);
        double sgh4 = -18.0 * ss4 * zes;
        double sh2 = -2.0 * ss2 * sz22;
        double sh3 = -2.0 * ss2 * (sz23 - sz21);

        //  ------------------------ do lunar terms ----------------------
        double ee2 = 2.0 * s1 * s6;
        double e3 = 2.0 * s1 * s7;
        double xi2 = 2.0 * s2 * z12;
        double xi3 = 2.0 * s2 * (z13 - z11);
        double xl2 = -2.0 * s3 * z2;
        double xl3 = -2.0 * s3 * (z3 - z1);
        double xl4 = -2.0 * s3 * (-21.0 - (9.0 * emsq)) * zel;
        double xgh2 = 2.0 * s4 * z32;
        double xgh3 = 2.0 * s4 * (z33 - z31);
        double xgh4 = -18.0 * s4 * zel;
        double xh2 = -2.0 * s2 * z22;
        double xh3 = -2.0 * s2 * (z23 - z21);

        this.snodm = snodm;
        this.cnodm = cnodm;
        this.sinim = sinim;
        this.cosim = cosim;
        this.sinomm = sinomm;
        this.cosomm = cosomm;
        this.day = day;
        this.e3 = e3;
        this.ee2 = ee2;
        this.em = em;
        this.emsq = emsq;
        this.gam = gam;
        this.peo = peo;
        this.pgho = pgho;
        this.pho = pho;
        this.pinco = pinco;
        this.plo = plo;
        this.rtemsq = rtemsq;
        this.se2 = se2;
        this.se3 = se3;
        this.sgh2 = sgh2;
        this.sgh3 = sgh3;
        this.sgh4 = sgh4;
        this.sh2 = sh2;
        this.sh3 = sh3;
        this.si2 = si2;
        this.si3 = si3;
        this.sl2 = sl2;
        this.sl3 = sl3;
        this.sl4 = sl4;
        this.s1 = s1;
        this.s2 = s2;
        this.s3 = s3;
        this.s4 = s4;
        this.s5 = s5;
        this.s6 = s6;
        this.s7 = s7;
        this.ss1 = ss1;
        this.ss2 = ss2;
        this.ss3 = ss3;
        this.ss4 = ss4;
        this.ss5 = ss5;
        this.ss6 = ss6;
        this.ss7 = ss7;
        this.sz1 = sz1;
        this.sz2 = sz2;
        this.sz3 = sz3;
        this.sz11 = sz11;
        this.sz12 = sz12;
        this.sz13 = sz13;
        this.sz21 = sz21;
        this.sz22 = sz22;
        this.sz23 = sz23;
        this.sz31 = sz31;
        this.sz32 = sz32;
        this.sz33 = sz33;
        this.xgh2 = xgh2;
        this.xgh3 = xgh3;
        this.xgh4 = xgh4;
        this.xh2 = xh2;
        this.xh3 = xh3;
        this.xi2 = xi2;
        this.xi3 = xi3;
        this.xl2 = xl2;
        this.xl3 = xl3;
        this.xl4 = xl4;
        this.nm = nm;
        this.z1 = z1;
        this.z2 = z2;
        this.z3 = z3;
        this.z11 = z11;
        this.z12 = z12;
        this.z13 = z13;
        this.z21 = z21;
        this.z22 = z22;
        this.z23 = z23;
        this.z31 = z31;
        this.z32 = z32;
        this.z33 = z33;
        this.zmol = zmol;
        this.zmos = zmos;
    }

    public DeepSpaceCommon(double snodm, double cnodm, double sinim, double cosim, double sinomm, double cosomm, double day, double e3, double ee2, double em, double emsq, double gam, double peo, double pgho, double pho, double pinco, double plo, double rtemsq, double se2, double se3, double sgh2, double sgh3, double sgh4, double sh2, double sh3, double si2, double si3, double sl2, double sl3, double sl4, double s1, double s2, double s3, double s4, double s5, double s6, double s7, double ss1, double ss2, double ss3, double ss4, double ss5, double ss6, double ss7, double sz1, double sz2, double sz3, double sz11, double sz12, double sz13, double sz21, double sz22, double sz23, double sz31, double sz32, double sz33, double xgh2, double xgh3, double xgh4, double xh2, double xh3, double xi2, double xi3, double xl2, double xl3, double xl4, double nm, double z1, double z2, double z3, double z11, double z12, double z13, double z21, double z22, double z23, double z31, double z32, double z33, double zmol, double zmos)
    {
        this.snodm = snodm;
        this.cnodm = cnodm;
        this.sinim = sinim;
        this.cosim = cosim;
        this.sinomm = sinomm;
        this.cosomm = cosomm;
        this.day = day;
        this.e3 = e3;
        this.ee2 = ee2;
        this.em = em;
        this.emsq = emsq;
        this.gam = gam;
        this.peo = peo;
        this.pgho = pgho;
        this.pho = pho;
        this.pinco = pinco;
        this.plo = plo;
        this.rtemsq = rtemsq;
        this.se2 = se2;
        this.se3 = se3;
        this.sgh2 = sgh2;
        this.sgh3 = sgh3;
        this.sgh4 = sgh4;
        this.sh2 = sh2;
        this.sh3 = sh3;
        this.si2 = si2;
        this.si3 = si3;
        this.sl2 = sl2;
        this.sl3 = sl3;
        this.sl4 = sl4;
        this.s1 = s1;
        this.s2 = s2;
        this.s3 = s3;
        this.s4 = s4;
        this.s5 = s5;
        this.s6 = s6;
        this.s7 = s7;
        this.ss1 = ss1;
        this.ss2 = ss2;
        this.ss3 = ss3;
        this.ss4 = ss4;
        this.ss5 = ss5;
        this.ss6 = ss6;
        this.ss7 = ss7;
        this.sz1 = sz1;
        this.sz2 = sz2;
        this.sz3 = sz3;
        this.sz11 = sz11;
        this.sz12 = sz12;
        this.sz13 = sz13;
        this.sz21 = sz21;
        this.sz22 = sz22;
        this.sz23 = sz23;
        this.sz31 = sz31;
        this.sz32 = sz32;
        this.sz33 = sz33;
        this.xgh2 = xgh2;
        this.xgh3 = xgh3;
        this.xgh4 = xgh4;
        this.xh2 = xh2;
        this.xh3 = xh3;
        this.xi2 = xi2;
        this.xi3 = xi3;
        this.xl2 = xl2;
        this.xl3 = xl3;
        this.xl4 = xl4;
        this.nm = nm;
        this.z1 = z1;
        this.z2 = z2;
        this.z3 = z3;
        this.z11 = z11;
        this.z12 = z12;
        this.z13 = z13;
        this.z21 = z21;
        this.z22 = z22;
        this.z23 = z23;
        this.z31 = z31;
        this.z32 = z32;
        this.z33 = z33;
        this.zmol = zmol;
        this.zmos = zmos;
    }

    public readonly double snodm;
    public readonly double cnodm;
    public readonly double sinim;
    public readonly double cosim;
    public readonly double sinomm;
    public readonly double cosomm;
    public readonly double day;
    public readonly double e3;
    public readonly double ee2;
    public readonly double em;
    public readonly double emsq;
    public readonly double gam;
    public readonly double peo;
    public readonly double pgho;
    public readonly double pho;
    public readonly double pinco;
    public readonly double plo;
    public readonly double rtemsq;
    public readonly double se2;
    public readonly double se3;
    public readonly double sgh2;
    public readonly double sgh3;
    public readonly double sgh4;
    public readonly double sh2;
    public readonly double sh3;
    public readonly double si2;
    public readonly double si3;
    public readonly double sl2;
    public readonly double sl3;
    public readonly double sl4;
    public readonly double s1;
    public readonly double s2;
    public readonly double s3;
    public readonly double s4;
    public readonly double s5;
    public readonly double s6;
    public readonly double s7;
    public readonly double ss1;
    public readonly double ss2;
    public readonly double ss3;
    public readonly double ss4;
    public readonly double ss5;
    public readonly double ss6;
    public readonly double ss7;
    public readonly double sz1;
    public readonly double sz2;
    public readonly double sz3;
    public readonly double sz11;
    public readonly double sz12;
    public readonly double sz13;
    public readonly double sz21;
    public readonly double sz22;
    public readonly double sz23;
    public readonly double sz31;
    public readonly double sz32;
    public readonly double sz33;
    public readonly double xgh2;
    public readonly double xgh3;
    public readonly double xgh4;
    public readonly double xh2;
    public readonly double xh3;
    public readonly double xi2;
    public readonly double xi3;
    public readonly double xl2;
    public readonly double xl3;
    public readonly double xl4;
    public readonly double nm;
    public readonly double z1;
    public readonly double z2;
    public readonly double z3;
    public readonly double z11;
    public readonly double z12;
    public readonly double z13;
    public readonly double z21;
    public readonly double z22;
    public readonly double z23;
    public readonly double z31;
    public readonly double z32;
    public readonly double z33;
    public readonly double zmol;
    public readonly double zmos;
}

