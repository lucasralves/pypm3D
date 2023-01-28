#include <math.h>


/*
    Struct
*/
typedef struct { double x, y, z; } Vec3D;
typedef struct { double x, y; } Vec2D;

/*
    Custom math
*/
double divide(double a, double b)
{
    if ((0.0 <= b) && (b < 1e-12)) {
        return a / (b + 1e-12);
    } else if ((-1e12 <= b) && (b < 0.0)) {
        return a / (b - 1e-12);
    } else {
        return a / b;
    }
}

double dot(Vec3D a, Vec3D b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/*
    Potential side
*/
void phi_panel_side(Vec3D p, Vec2D pa, Vec2D pb, double *source, double *doublet)
{
    double dab = sqrt((pb.x - pa.x) * (pb.x - pa.x) + (pb.y - pa.y) * (pb.y - pa.y));
    double ra = sqrt((p.x - pa.x) * (p.x - pa.x) + (p.y - pa.y) * (p.y - pa.y) + p.z * p.z);
    double rb = sqrt((p.x - pb.x) * (p.x - pb.x) + (p.y - pb.y) * (p.y - pb.y) + p.z * p.z);
    double mab = divide(pb.y - pa.y, pb.x - pa.x);
    double ea = (p.x - pa.x) * (p.x - pa.x) + p.z * p.z;
    double eb = (p.x - pb.x) * (p.x - pb.x) + p.z * p.z;
    double ha = (p.x - pa.x) * (p.y - pa.y);
    double hb = (p.x - pb.x) * (p.y - pb.y);
    double l1 = (mab * ea - ha) / (p.z * ra);
    double l2 = (mab * eb - hb) / (p.z * rb);

    double f1 = ((p.x - pa.x) * (pb.y - pa.y) - (p.y - pa.y) * (pb.x - pa.x)) / dab;
    double f2 = log((ra + rb + dab) / (ra + rb - dab));
    double f3 = atan2(l1 - l2, 1 + l1 * l2);

    *source = - (1 / (4.0 * M_PI)) * (f1 * f2 - p.z * f3);
    *doublet = (1 / (4.0 * M_PI)) * f3;
}

/*
    Panel potential
*/
void panel_potential(int n,
                     int sides,
                     Vec2D p1, Vec2D p2, Vec2D p3, Vec2D p4,
                     Vec3D e1, Vec3D e2, Vec3D e3,
                     Vec3D p_avg,
                     double p_ctrls[],
                     double source[],
                     double doublet[])
{

    Vec3D vec;
    Vec3D p;

    double s1, s2, s3, s4;
    double d1, d2, d3, d4;

    for (int i = 0; i < n; i++)
    {

        vec.x = p_ctrls[3 * i] - p_avg.x;
        vec.y = p_ctrls[3 * i + 1] - p_avg.y;
        vec.z = p_ctrls[3 * i + 2] - p_avg.z;

        p.x = dot(vec, e1);
        p.y = dot(vec, e2);
        p.z = dot(vec, e3);
        
        phi_panel_side(p, p1, p2, &s1, &d1);
        phi_panel_side(p, p2, p3, &s2, &d2);

        if (sides == 4) {
            phi_panel_side(p, p3, p4, &s3, &d3);
            phi_panel_side(p, p4, p1, &s4, &d4);
        } else if (sides == 3)
        {
            phi_panel_side(p, p3, p1, &s3, &d3);
            s4 = 0.0;
            d4 = 0.0;
        }

        source[i] = s1 + s2 + s3 + s4;
        doublet[i] = d1 + d2 + d3 + d4;
    }

}
