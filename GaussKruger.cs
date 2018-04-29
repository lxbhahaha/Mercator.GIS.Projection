using Mercator.Mathematics.Calculus;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Mercator.GIS.Projection
{
    public class GaussKrugerProjection
    {
        /// <summary>
        /// 比例因子
        /// </summary>
        public double ScaleFactor
        {
            get
            {
                return k0;
            }
        }

        /// <summary>
        /// 起始纬度（赤道）
        /// </summary>
        public double LatitudeOfOrigin
        {
            get
            {
                return φ0;
            }
        }

        /// <summary>
        /// 起始经度（中央子午线）
        /// </summary>
        public double LongitudeOfOrigin
        {
            get
            {
                return λ0;
            }
            set
            {
                λ0 = value;
            }
        }

        /// <summary>
        /// 东伪偏移
        /// </summary>
        public double FalseEasting
        {
            get
            {
                return FE;
            }
            set
            {
                FE = value;
            }
        }

        /// <summary>
        /// 北伪偏移
        /// </summary>
        public double FalseNorthing
        {
            get
            {
                return FN;
            }
            set
            {
                FN = value;
            }
        }

        /// <summary>
        /// Scale factor at natural origin
        /// </summary>
        private double k0 = 1.0;

        /// <summary>
        /// Latitude of natural origin,Usually 0°
        /// </summary>
        private double φ0 = 0.0;

        /// <summary>
        /// Longitude of natural origin,central meridian (CM) 
        /// </summary>
        private double λ0 = 0.0;

        /// <summary>
        /// the meridional arc distance from equator to the projection origin
        /// </summary>
        private double M0
        {
            get
            {
                double[] m0 = new double[4];
                m0[0] = (1 - Math.Pow(e, 2) / 4.0 - 3.0 * Math.Pow(e, 4) / 64.0 - 5.0 * Math.Pow(e, 6) / 256.0) * φ0;
                m0[1] = (3.0 * Math.Pow(e, 2) / 8.0 + 3.0 * Math.Pow(e, 4) / 32.0 + 45.0 * Math.Pow(e, 6) / 1024.0) * Math.Sin(2 * φ0);
                m0[2] = (15.0 * Math.Pow(e, 4) / 256.0 + 45.0 * Math.Pow(e, 6) / 1024.0) * Math.Sin(4 * φ0);
                m0[3] = (35.0 * Math.Pow(e, 6) / 3072.0) * Math.Sin(6 * φ0);

                return a * (m0[0] - m0[1] + m0[2] - m0[3]);
            }
        }

        /// <summary>
        /// False easting 
        /// </summary>
        private double FE;

        /// <summary>
        /// False northing 
        /// </summary>
        private double FN;

        private double n
        {
            get
            {
                double f = Ellipsoid.f;
                return f / (2 - f);
            }
        }

        private double B
        {
            get
            {
                double a = Ellipsoid.a;
                return (a / (1 + n)) * (1 + Math.Pow(n, 2) / 4.0 + Math.Pow(n, 4) / 64.0);
            }
        }

        private double a
        {
            get
            {
                return Ellipsoid.a;
            }
        }

        private double e
        {
            get
            {
                return Ellipsoid.e;
            }
        }

        public ReferenceEllipsoid Ellipsoid;

        public void Forward(double φ, double λ, out double E, out double N)
        {
            double h1 = n / 2 - 2.0 / 3.0 * Math.Pow(n, 2) + 5.0 / 16.0 * Math.Pow(n, 3) + 41.0 / 180.0 * Math.Pow(n, 4);
            double h2 = 13.0 / 48.0 * Math.Pow(n, 2) - 3.0 / 4.0 * Math.Pow(n, 3) + 557.0 / 1440.0 * Math.Pow(n, 4);
            double h3 = 61.0 / 240.0 * Math.Pow(n, 3) - 103.0 / 140.0 * Math.Pow(n, 4);
            double h4 = 49561.0 / 161280.0 * Math.Pow(n, 4);

            double Q = InverseHyperbolic.Asinh(Math.Tan(φ)) - (Ellipsoid.e * InverseHyperbolic.Atanh(Ellipsoid.e * Math.Sin(φ)));
            double β = Math.Atan(Math.Sinh(Q));
            double η0 = InverseHyperbolic.Atanh(Math.Cos(β) * Math.Sin(λ - λ0));
            double ξ0 = Math.Asin(Math.Sin(β) * Math.Cosh(η0));

            double ξ1 = h1 * Math.Sin(2 * ξ0) * Math.Cosh(2 * η0);
            double ξ2 = h2 * Math.Sin(4 * ξ0) * Math.Cosh(4 * η0);
            double ξ3 = h3 * Math.Sin(6 * ξ0) * Math.Cosh(6 * η0);
            double ξ4 = h4 * Math.Sin(8 * ξ0) * Math.Cosh(8 * η0);
            double ξ = ξ0 + ξ1 + ξ2 + ξ3 + ξ4;

            double η1 = h1 * Math.Cos(2 * ξ0) * Math.Sinh(2 * η0);
            double η2 = h2 * Math.Cos(4 * ξ0) * Math.Sinh(4 * η0);
            double η3 = h3 * Math.Cos(6 * ξ0) * Math.Sinh(6 * η0);
            double η4 = h4 * Math.Cos(8 * ξ0) * Math.Sinh(8 * η0);
            double η = η0 + η1 + η2 + η3 + η4;

            E = FE + k0 * B * η;
            N = FN + k0 * (B * ξ - M0);
        }

        public void Reverse(double E, double N, out double φ, out double λ)
        {
            double h1ˊ = n / 2.0 - 2.0 / 3.0 * Math.Pow(n, 2) + 37.0 / 96.0 * Math.Pow(n, 3) - 1 / 360.0 * Math.Pow(n, 4);
            double h2ˊ = 1 / 48.0 * Math.Pow(n, 2) + 1 / 15.0 * Math.Pow(n, 3) - 437.0 / 1440.0 * Math.Pow(n, 4);
            double h3ˊ = 17.0 / 480.0 * Math.Pow(n, 3) - 37.0 / 840.0 * Math.Pow(n, 4);
            double h4ˊ = 4397.0 / 161280.0 * Math.Pow(n, 4);

            double ηˊ = (E - FE) / (B * k0);
            double ξˊ = ((N - FN) + k0 * M0) / (B * k0);

            double ξ1ˊ = h1ˊ * Math.Sin(2 * ξˊ) * Math.Cosh(2 * ηˊ);
            double ξ2ˊ = h2ˊ * Math.Sin(4 * ξˊ) * Math.Cosh(4 * ηˊ);
            double ξ3ˊ = h3ˊ * Math.Sin(6 * ξˊ) * Math.Cosh(6 * ηˊ);
            double ξ4ˊ = h4ˊ * Math.Sin(8 * ξˊ) * Math.Cosh(8 * ηˊ);
            double ξ0ˊ = ξˊ - (ξ1ˊ + ξ2ˊ + ξ3ˊ + ξ4ˊ);

            double η1ˊ = h1ˊ * Math.Cos(2 * ξˊ) * Math.Sinh(2 * ηˊ);
            double η2ˊ = h2ˊ * Math.Cos(4 * ξˊ) * Math.Sinh(4 * ηˊ);
            double η3ˊ = h3ˊ * Math.Cos(6 * ξˊ) * Math.Sinh(6 * ηˊ);
            double η4ˊ = h4ˊ * Math.Cos(8 * ξˊ) * Math.Sinh(8 * ηˊ);
            double η0ˊ = ηˊ - (η1ˊ + η2ˊ + η3ˊ + η4ˊ);

            double βˊ = Math.Asin(Math.Sin(ξ0ˊ) / Math.Cosh(η0ˊ));
            double Qˊ = InverseHyperbolic.Asinh(Math.Tan(βˊ));
            double Qˊˊ = Qˊ + (e * InverseHyperbolic.Atanh(e * Math.Tanh(Qˊ)));
            double q = Qˊ + (e * InverseHyperbolic.Atanh(e * Math.Tanh(Qˊˊ)));
            while (Math.Abs(Qˊˊ - q) > 0.0000000001)
            {
                Qˊˊ = q;
                q = Qˊ + (e * InverseHyperbolic.Atanh(e * Math.Tanh(Qˊˊ)));
            }

            φ = Math.Atan(Math.Sinh(Qˊˊ));
            λ = λ0 + Math.Asin(Math.Tanh(η0ˊ) / Math.Cos(βˊ));
        }
    }
}
