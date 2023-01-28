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

Vec3D cross(Vec3D a, Vec3D b)
{
    Vec3D out;
    out.x = a.y * b.z - a.z * b.y;
    out.y = a.z * b.x - a.x * b.z;
    out.z = a.x * b.y - a.y * b.x;
    return out;
}


double norm(Vec3D a)
{
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

/*
    Potential side
*/
void vel_panel_side(Vec3D p, Vec2D pa, Vec2D pb, Vec3D *source, Vec3D *doublet)
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

    source->x = (1 / (4.0 * M_PI)) * ((pb.y - pa.y) / dab) * log((ra + rb - dab) / (ra + rb + dab));
    source->y = (1 / (4.0 * M_PI)) * ((pa.x - pb.x) / dab) * log((ra + rb - dab) / (ra + rb + dab));
    source->z = (1 / (4.0 * M_PI)) * atan2(l1 - l2, 1 + l1 * l2);

    Vec3D r1 = {p.x - pa.x, p.y - pa.y, p.z};
    Vec3D r2 = {p.x - pb.x, p.y - pb.y, p.z};
    Vec3D r0 = {pa.x - pb.x, pa.y - pb.y, 0.0};
    Vec3D r1_x_r2 = cross(r1, r2);

    double r1_norm = norm(r1);
    double r2_norm = norm(r2);
    double r1_x_r2_norm = norm(r1_x_r2);
    double k = (1 / (4.0 * M_PI)) * (dot(r0, r1) / r1_norm - dot(r0, r2) / r2_norm) / r1_x_r2_norm;

    doublet->x = k * r1_x_r2.x;
    doublet->y = k * r1_x_r2.y;
    doublet->z = k * r1_x_r2.z;

}

/*
    Panel potential
*/
void panel_velocity(int n,
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

    Vec3D s1, s2, s3, s4;
    Vec3D d1, d2, d3, d4;

    double u, v, w;

    for (int i = 0; i < n; i++)
    {

        vec.x = p_ctrls[3 * i] - p_avg.x;
        vec.y = p_ctrls[3 * i + 1] - p_avg.y;
        vec.z = p_ctrls[3 * i + 2] - p_avg.z;

        p.x = dot(vec, e1);
        p.y = dot(vec, e2);
        p.z = dot(vec, e3);
        
        vel_panel_side(p, p1, p2, &s1, &d1);
        vel_panel_side(p, p2, p3, &s2, &d2);

        if (sides == 4) {
            vel_panel_side(p, p3, p4, &s3, &d3);
            vel_panel_side(p, p4, p1, &s4, &d4);
        } else if (sides == 3)
        {
            vel_panel_side(p, p3, p1, &s3, &d3);
            s4.x = 0.0; s4.y = 0.0; s4.z = 0.0;
            d4.x = 0.0; d4.y = 0.0; d4.z = 0.0;
        }

        u = s1.x + s2.x + s3.x + s4.x;
        v = s1.y + s2.y + s3.y + s4.y;
        w = s1.z + s2.z + s3.z + s4.z;

        source[3 * i] = u * e1.x + v * e2.x + w * e3.x;
        source[3 * i + 1] = u * e1.y + v * e2.y + w * e3.y;
        source[3 * i + 2] = u * e1.z + v * e2.z + w * e3.z;

        u = d1.x + d2.x + d3.x + d4.x;
        v = d1.y + d2.y + d3.y + d4.y;
        w = d1.z + d2.z + d3.z + d4.z;

        doublet[3 * i] = u * e1.x + v * e2.x + w * e3.x;
        doublet[3 * i + 1] = u * e1.y + v * e2.y + w * e3.y;
        doublet[3 * i + 2] = u * e1.z + v * e2.z + w * e3.z;
    }

}
