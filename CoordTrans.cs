//=============================================================================
// 
// CoordTrans (Transforms points to points)
//      Affine Transform, Mesh Transform ...
// PWH. 2021.08.06. ==> C#...
//
//

using System;
using System.Collections.Generic;
using System.Text;

namespace gtl.CoordTrans {

    public abstract class ICoordTrans {
        /// <summary>
        /// Transforms point to new coordinate
        /// </summary>
        /// <param name="pt"></param>
        /// <returns></returns>
        public abstract xPoint2d Trans(xPoint2d pt);

        /// <summary>
        /// Inverse Transforms point
        /// </summary>
        /// <param name="pt"></param>
        /// <returns></returns>
        public abstract xPoint2d TransI(xPoint2d pt);


        public bool PtInRect(xPoint2d ptLeftBottom, xPoint2d ptRightTop, xPoint2d pt) {
            return (pt.x >= ptLeftBottom.x) && (pt.x <= ptRightTop.x)
                && (pt.y >= ptLeftBottom.y) && (pt.y <= ptRightTop.y);
        }

        /// <summary>
        /// Transforms pt
        /// </summary>
        /// <param name="ptsSrc">source [2,2] enclosing Rectangle.</param>
        /// <param name="ptsDest">dest [2,2] enclosing Rectangle.</param>
        /// <param name="pt">Point</param>
        /// <returns>Transformed pt</returns>
        protected static xPoint2d Transform2dLinear(xPoint2d[,] ptsSrc, xPoint2d[,] ptsDest, xPoint2d pt) {
            Func<xPoint2d, double> TERM1 = p => p.x;
            Func<xPoint2d, double> TERM2 = p => p.y;
            Func<xPoint2d, double> TERM3 = p => p.x * p.y;
            Func<xPoint2d, double> TERM4 = p => 1.0;

            xPoint2d[] pts = new xPoint2d[4];
            pts[0] = ptsSrc[0, 0];
            pts[1] = ptsSrc[1, 0];
            pts[2] = ptsSrc[1, 1];
            pts[3] = ptsSrc[0, 1];

            xMatrix m2, m2i;

            m2 = new xMatrix(4, 4, 0);
            for (int i = 0; i < 4; i++) {
                m2.At(i, 0) = TERM1(pts[i]);
                m2.At(i, 1) = TERM2(pts[i]);
                m2.At(i, 2) = TERM3(pts[i]);
                m2.At(i, 3) = TERM4(pts[i]);
            }
            bool bOK = false;
            m2i = m2.Inverse(0.0, ref bOK);
            if (!bOK)
                throw new Exception("Wrong xMatrix");

            xPoint2d ptTrans = new xPoint2d(pt);
            xMatrix m1 = new xMatrix(4, 1, 0.0);
            xMatrix m3;

            // x
            m1.At(0, 0) = ptsDest[0, 0].x;
            m1.At(1, 0) = ptsDest[1, 0].x;
            m1.At(2, 0) = ptsDest[1, 1].x;
            m1.At(3, 0) = ptsDest[0, 1].x;

            m3 = m2i * m1;
            ptTrans.x = m3.At(0, 0) * TERM1(pt) + m3.At(1, 0) * TERM2(pt) + m3.At(2, 0) * TERM3(pt) + m3.At(3, 0) * TERM4(pt);

            m1.At(0, 0) = ptsDest[0, 0].y;
            m1.At(1, 0) = ptsDest[1, 0].y;
            m1.At(2, 0) = ptsDest[1, 1].y;
            m1.At(3, 0) = ptsDest[0, 1].y;

            m3 = m2i * m1;
            ptTrans.y = m3.At(0, 0) * TERM1(pt) + m3.At(1, 0) * TERM2(pt) + m3.At(2, 0) * TERM3(pt) + m3.At(3, 0) * TERM4(pt);

            return ptTrans;
        }
    }

    public class xCoordTrans2d : ICoordTrans {
        // ptNew = m_mat(pt - ptShift) + ptOffset;

        public xMatrix m_mat;
        public xPoint2d m_ptShift;
        public xPoint2d m_ptOffset;

        public xCoordTrans2d() {
            m_mat = xMatrix.eye(2);
            m_ptShift = new xPoint2d(0, 0);
            m_ptOffset = new xPoint2d(0, 0);
        }
        public xCoordTrans2d(double dAngleRad) {
            m_mat = GetRotationMatrix(dAngleRad);
            m_ptShift = new xPoint2d(0, 0);
            m_ptOffset = new xPoint2d(0, 0);
        }

        public static xMatrix GetRotationMatrix(double dAngleRad) {
            xMatrix mat = new xMatrix(2, 2);
            double c = Math.Cos(dAngleRad);
            double s = Math.Sin(dAngleRad);
            mat.At(0, 0) = c;
            mat.At(0, 1) = -s;
            mat.At(1, 0) = s;
            mat.At(1, 1) = c;
            return mat;
        }

        public override xPoint2d Trans(xPoint2d pt) {
            xPoint2d pt2 = new xPoint2d(pt - m_ptShift);
            return new xPoint2d(m_mat.At(0, 0) * pt2.x + m_mat.At(0, 1) * pt2.y + m_ptOffset.x,
                m_mat.At(1, 0) * pt2.x + m_mat.At(1, 1) * pt2.y + m_ptOffset.y);
        }
        public override xPoint2d TransI(xPoint2d pt) {
            bool bOK = false;
            xMatrix mat = m_mat.Inverse(0.0, ref bOK);
            if (!bOK)
                throw new Exception("No Inverse xMatrix.");

            xPoint2d pt2 = new xPoint2d(pt - m_ptOffset);
            return new xPoint2d(mat.At(0, 0) * pt2.x + mat.At(0, 1) * pt2.y + m_ptShift.x,
                mat.At(1, 0) * pt2.x + mat.At(1, 1) * pt2.y + m_ptShift.y);
        }
    }

    public class xCoordTrans1dXY : ICoordTrans {
        private double[] m_xs;
        private xPoint2d[] m_vxs;
        private double[] m_ys;
        private xPoint2d[] m_vys;

        protected bool CheckAxisValue(double[] series, xPoint2d[] vdisp) {
            if (series.Length != vdisp.Length)
                return false;

            for (int i = 1; i < series.Length; i++) {
                if (series[i] <= series[i - 1])
                    return false;
            }

            return true;
        }

        public bool SetAxisX(double[] xs, xPoint2d[] vxs) {
            if (!CheckAxisValue(xs, vxs))
                return false;

            m_xs = xs.Clone() as double[];
            m_vxs = vxs.Clone() as xPoint2d[];
            return true;
        }
        public bool SetAxisY(double[] ys, xPoint2d[] vys) {
            if (!CheckAxisValue(ys, vys))
                return false;

            m_ys = ys.Clone() as double[];
            m_vys = vys.Clone() as xPoint2d[];
            return true;
        }

        protected int FindRange(double[] series, Func<int, bool> funcCompare) {
            if (series == null)
                return -1;
            for (int i = 0; i < series.Length; i++) {
                if (funcCompare(i))
                    return i;
            }
            return -1;
        }


        public override xPoint2d Trans(xPoint2d pt) {
            if ((m_xs == null) && (m_ys == null))
                return pt;

            int ix = FindRange(m_xs, i => pt.x <= m_xs[i]);
            int iy = FindRange(m_ys, i => pt.y <= m_ys[i]);

            xPoint2d ptTrans = new xPoint2d(pt);
            if (ix > 0) {
                double t = (m_xs[ix] - pt.x) / (m_xs[ix] - m_xs[ix - 1]);
                ptTrans += m_vxs[ix - 1] * t + m_vxs[ix] * (1 - t);
            }
            if (iy > 0) {
                double t = (m_ys[iy] - pt.y) / (m_ys[iy] - m_ys[iy - 1]);
                ptTrans += m_vys[iy - 1] * t + m_vys[iy] * (1 - t);
            }

            return ptTrans;
        }

        public override xPoint2d TransI(xPoint2d pt) {
            if ((m_xs == null) && (m_ys == null))
                return pt;

            int ix = FindRange(m_xs, i => pt.x <= m_xs[i] + m_vxs[i].x);
            int iy = FindRange(m_ys, i => pt.y <= m_ys[i] + m_vys[i].y);

            xPoint2d ptTrans = new xPoint2d(pt);
            if (ix > 0) {
                double t = (m_xs[ix] + m_vxs[ix].x - pt.x) / (m_xs[ix] + m_vxs[ix].x - m_xs[ix - 1] - m_vxs[ix - 1].x);
                //double t = (m_xs[ix] - pt.x) / (m_xs[ix] - m_xs[ix - 1] );
                ptTrans -= m_vxs[ix - 1] * t + m_vxs[ix] * (1 - t);
            }
            if (iy > 0) {
                double t = (m_ys[iy] + m_vys[iy].y - pt.y) / (m_ys[iy] + m_vys[iy].y - m_ys[iy - 1] - m_vys[iy - 1].y);
                //double t = (m_ys[iy] - pt.y) / (m_ys[iy] - m_ys[iy - 1]);
                ptTrans -= m_vys[iy - 1] * t + m_vys[iy] * (1 - t);
            }

            return ptTrans;
        }
    }


    public class xCoordTrans2dMesh : ICoordTrans {
        private xMeshTable m_src;
        private xMeshTable m_dest;

        public xCoordTrans2dMesh(xMeshTable src, xMeshTable dest) {
            if (src.rows != dest.rows || src.cols != dest.cols)
                throw new Exception("src.size != dest.size");
            m_src = new xMeshTable(src);
            m_dest = new xMeshTable(dest);
        }

        public override xPoint2d Trans(xPoint2d pt) {
			try
			{
                pt = Trans(m_src, m_dest, pt);
			}
			catch {}
            return pt;
        }

        public override xPoint2d TransI(xPoint2d pt) {
            return Trans(m_dest, m_src, pt);
        }


        protected static xPoint2d Trans(xMeshTable src, xMeshTable dest, xPoint2d pt) {
            if ((src == null) || (dest == null))
                return pt;

            int iy = 0, ix = 0;
            if (!src.FindEnclosingPTS(pt, ref iy, ref ix))
			{
                throw new Exception("OutBound");
                //return pt;
			}

            xPoint2d[,] ptsSrc = new xPoint2d[2, 2];
            xPoint2d[,] ptsDest = new xPoint2d[2, 2];
            ptsSrc[0, 0] = src.At(iy-1, ix-1);
            ptsSrc[0, 1] = src.At(iy-1, ix);
            ptsSrc[1, 0] = src.At(iy, ix-1);
            ptsSrc[1, 1] = src.At(iy, ix);
            ptsDest[0, 0] = dest.At(iy-1, ix-1);
            ptsDest[0, 1] = dest.At(iy-1, ix);
            ptsDest[1, 0] = dest.At(iy, ix-1);
            ptsDest[1, 1] = dest.At(iy, ix);
            return Transform2dLinear(ptsSrc, ptsDest, pt);
        }
    }

}
