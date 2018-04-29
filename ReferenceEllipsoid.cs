using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Mercator.GIS.Projection
{
    /// <summary>
    /// 参考椭球
    /// </summary>
    public class ReferenceEllipsoid
    {
        /// <summary>
        /// 长半轴
        /// </summary>
        public double SemimajorAxis
        {
            get
            {
                return a;
            }
            set
            {
                a = value;
            }
        }

        /// <summary>
        /// 短半轴
        /// </summary>
        public double SemiminorAxis
        {
            get
            {
                return b;
            }
            set
            {
                b = value;
            }
        }

        /// <summary>
        /// 反扁率
        /// </summary>
        public double InverseFlattening
        {
            set
            {
                b = a - a / value;
            }
        }

        /// <summary>
        /// semi-major axis
        /// </summary>
        internal double a;

        /// <summary>
        /// semi-minor axis
        /// </summary>
        internal double b;

        /// <summary>
        /// flattening
        /// </summary>
        internal double f
        {
            get
            {
                return (a - b) / a;
            }
        }

        /// <summary>
        /// eccentricity
        /// </summary>
        internal double e
        {
            get
            {
                return Math.Sqrt(2 * f - Math.Pow(f, 2));
            }
        }

        /// <summary>
        /// second eccentricity
        /// </summary>
        internal double eˊ
        {
            get
            {
                return Math.Sqrt(Math.Pow(e, 2) / (1 - Math.Pow(e, 2)));
            }
        }

        /// <summary>
        /// 1975年国际会议推荐的参考椭球
        /// </summary>
        public static ReferenceEllipsoid International1975
        {
            get
            {
                var ellipsoid = new ReferenceEllipsoid();
                ellipsoid.a = 6378140.0000000000;
                ellipsoid.b = 6356755.2881575287;
                return ellipsoid;
            }
        }

        public static ReferenceEllipsoid OSGB1936
        {
            get
            {
                var ellipsoid = new ReferenceEllipsoid();
                ellipsoid.a = 6377563.396;
                ellipsoid.InverseFlattening = 299.3249646;
                return ellipsoid;
            }
        }

        /// <summary>
        /// WGS 1984 椭球
        /// </summary>
        public static ReferenceEllipsoid WGS84
        {
            get
            {
                var ellipsoid = new ReferenceEllipsoid();
                ellipsoid.a = 6378137.0000000000;
                ellipsoid.b = 6356752.3142451795;
                return ellipsoid;
            }
        }

        /// <summary>
        /// 克拉索夫斯基（Krasovsky）1940
        /// </summary>
        public static ReferenceEllipsoid Krasovsky1940
        {
            get
            {
                var ellipsoid = new ReferenceEllipsoid();
                ellipsoid.a = 6378245.0000000000;
                ellipsoid.b = 6356863.0187730473;
                return ellipsoid;
            }
        }

        /// <summary>
        /// ITRF97
        /// </summary>
        public static ReferenceEllipsoid ITRF97
        {
            get
            {
                var ellipsoid = new ReferenceEllipsoid();
                ellipsoid.a = 6378137.0000000000;
                ellipsoid.b = 6356752.3141403358;
                return ellipsoid;
            }
        }
    }
}
