#include <stdio.h>
#include <math.h>

/*-------------------------------------------
    CONSTANTS
--------------------------------------------*/
const double ZERO_ERROR = 1e-12;

/*-------------------------------------------
    STRUCTS
--------------------------------------------*/
typedef struct {
    double x;
    double y;
    double z;
} Vec3D;

typedef struct {
    double x;
    double y;
} Vec2D;

/*-------------------------------------------
    MATH
--------------------------------------------*/
Vec3D cross_func(Vec3D a, Vec3D b)
{
    Vec3D out;
    out.x = a.y * b.z - a.z * b.y;
    out.y = a.z * b.x - a.x * b.z;
    out.z = a.x * b.y - a.y * b.x;
    return out;
}

double norm_func(Vec3D a)
{
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

double dot_func(Vec3D a, Vec3D b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double tri_plane_area(Vec2D a, Vec2D b, Vec2D c)
{
    return 0.5 * fabs(a.x * b.y + b.x * c.y + c.x * a.y - (a.y * b.x + b.y * c.x + c.y * a.x));
}

/*-------------------------------------------
    INDUCED VELOCITY
--------------------------------------------*/
void quad_planar_panel_velocity(Vec3D p,
                                Vec2D p1, Vec2D p2, Vec2D p3, Vec2D p4,
                                Vec3D e1, Vec3D e2, Vec3D e3,
                                Vec3D *source_vel, Vec3D *doublet_vel)
{
    double d12 = sqrt( (p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y) );
    double d23 = sqrt( (p3.x - p2.x) * (p3.x - p2.x) + (p3.y - p2.y) * (p3.y - p2.y) );
    double d34 = sqrt( (p4.x - p3.x) * (p4.x - p3.x) + (p4.y - p3.y) * (p4.y - p3.y) );
    double d41 = sqrt( (p1.x - p4.x) * (p1.x - p4.x) + (p1.y - p4.y) * (p1.y - p4.y) );

    double C12 = (p2.x - p1.x) / d12;
    double C23 = (p3.x - p2.x) / d23;
    double C34 = (p4.x - p3.x) / d34;
    double C41 = (p1.x - p4.x) / d41;

    double S12 = (p2.y - p1.y) / d12;
    double S23 = (p3.y - p2.y) / d23;
    double S34 = (p4.y - p3.y) / d34;
    double S41 = (p1.y - p4.y) / d41;

    double s12_1 = (p1.x - p.x) * C12 + (p1.y - p.y) * S12;
    double s12_2 = (p2.x - p.x) * C12 + (p2.y - p.y) * S12;
    double s23_2 = (p2.x - p.x) * C23 + (p2.y - p.y) * S23;
    double s23_3 = (p3.x - p.x) * C23 + (p3.y - p.y) * S23;
    double s34_3 = (p3.x - p.x) * C34 + (p3.y - p.y) * S34;
    double s34_4 = (p4.x - p.x) * C34 + (p4.y - p.y) * S34;
    double s41_4 = (p4.x - p.x) * C41 + (p4.y - p.y) * S41;
    double s41_1 = (p1.x - p.x) * C41 + (p1.y - p.y) * S41;

    double R12 = (p.x - p1.x) * S12 - (p.y - p1.y) * C12;
    double R23 = (p.x - p2.x) * S23 - (p.y - p2.y) * C23;
    double R34 = (p.x - p3.x) * S34 - (p.y - p3.y) * C34;
    double R41 = (p.x - p4.x) * S41 - (p.y - p4.y) * C41;

    double r1 = sqrt((p.x - p1.x) * (p.x - p1.x) + (p.y - p1.y) * (p.y - p1.y) + p.z * p.z);
    double r2 = sqrt((p.x - p2.x) * (p.x - p2.x) + (p.y - p2.y) * (p.y - p2.y) + p.z * p.z);
    double r3 = sqrt((p.x - p3.x) * (p.x - p3.x) + (p.y - p3.y) * (p.y - p3.y) + p.z * p.z);
    double r4 = sqrt((p.x - p4.x) * (p.x - p4.x) + (p.y - p4.y) * (p.y - p4.y) + p.z * p.z);

    double Q12 = log((r1 + r2 + d12) / (r1 + r2 - d12));
    double Q23 = log((r2 + r3 + d23) / (r2 + r3 - d23));
    double Q34 = log((r3 + r4 + d34) / (r3 + r4 - d34));
    double Q41 = log((r4 + r1 + d41) / (r4 + r1 - d41));

    double F12 = R12 * fabs(p.z) * (r1 * s12_2 - r2 * s12_1);
    double F23 = R23 * fabs(p.z) * (r2 * s23_3 - r3 * s23_2);
    double F34 = R34 * fabs(p.z) * (r3 * s34_4 - r4 * s34_3);
    double F41 = R41 * fabs(p.z) * (r4 * s41_1 - r1 * s41_4);

    double G12 = r1 * r2 * R12 * R12 + p.z * p.z * s12_1 * s12_2;
    double G23 = r2 * r3 * R23 * R23 + p.z * p.z * s23_2 * s23_3;
    double G34 = r3 * r4 * R34 * R34 + p.z * p.z * s34_3 * s34_4;
    double G41 = r4 * r1 * R41 * R41 + p.z * p.z * s41_4 * s41_1;

    double J12 = atan2(F12, G12);
    double J23 = atan2(F23, G23);
    double J34 = atan2(F34, G34);
    double J41 = atan2(F41, G41);

    double R12_x = S12; double R12_y = - C12;
    double R23_x = S23; double R23_y = - C23;
    double R34_x = S34; double R34_y = - C34;
    double R41_x = S41; double R41_y = - C41;

    double r1_x = (p.x - p1.x) / r1; double r1_y = (p.y - p1.y) / r1; double r1_z = p.z / r1;
    double r2_x = (p.x - p2.x) / r2; double r2_y = (p.y - p2.y) / r2; double r2_z = p.z / r2;
    double r3_x = (p.x - p3.x) / r3; double r3_y = (p.y - p3.y) / r3; double r3_z = p.z / r3;
    double r4_x = (p.x - p4.x) / r4; double r4_y = (p.y - p4.y) / r4; double r4_z = p.z / r4;

    double s12_1_x = - C12; double s12_2_x = - C12; double s12_1_y = - S12; double s12_2_y = - S12;
    double s23_2_x = - C23; double s23_3_x = - C23; double s23_2_y = - S23; double s23_3_y = - S23;
    double s34_3_x = - C34; double s34_4_x = - C34; double s34_3_y = - S34; double s34_4_y = - S34;
    double s41_4_x = - C41; double s41_1_x = - C41; double s41_4_y = - S41; double s41_1_y = - S41;

    double F12_x = fabs(p.z) * (R12_x * (r1 * s12_2 - r2 * s12_1) + R12 * (r1_x * s12_2 + r1 * s12_2_x - (r2_x * s12_1 + r2 * s12_1_x)));
    double F23_x = fabs(p.z) * (R23_x * (r2 * s23_3 - r3 * s23_2) + R23 * (r2_x * s23_3 + r2 * s23_3_x - (r3_x * s23_2 + r3 * s23_2_x)));
    double F34_x = fabs(p.z) * (R34_x * (r3 * s34_4 - r4 * s34_3) + R34 * (r3_x * s34_4 + r3 * s34_4_x - (r4_x * s34_3 + r4 * s34_3_x)));
    double F41_x = fabs(p.z) * (R41_x * (r4 * s41_1 - r1 * s41_4) + R41 * (r4_x * s41_1 + r4 * s41_1_x - (r1_x * s41_4 + r1 * s41_4_x)));

    double F12_y = fabs(p.z) * (R12_y * (r1 * s12_2 - r2 * s12_1) + R12 * (r1_y * s12_2 + r1 * s12_2_y - (r2_y * s12_1 + r2 * s12_1_y)));
    double F23_y = fabs(p.z) * (R23_y * (r2 * s23_3 - r3 * s23_2) + R23 * (r2_y * s23_3 + r2 * s23_3_y - (r3_y * s23_2 + r3 * s23_2_y)));
    double F34_y = fabs(p.z) * (R34_y * (r3 * s34_4 - r4 * s34_3) + R34 * (r3_y * s34_4 + r3 * s34_4_y - (r4_y * s34_3 + r4 * s34_3_y)));
    double F41_y = fabs(p.z) * (R41_y * (r4 * s41_1 - r1 * s41_4) + R41 * (r4_y * s41_1 + r4 * s41_1_y - (r1_y * s41_4 + r1 * s41_4_y)));

    double F12_z = (p.z / fabs(p.z)) * R12 * (r1 * s12_2 - r2 * s12_1);
    double F23_z = (p.z / fabs(p.z)) * R23 * (r2 * s23_3 - r3 * s23_2);
    double F34_z = (p.z / fabs(p.z)) * R34 * (r3 * s34_4 - r4 * s34_3);
    double F41_z = (p.z / fabs(p.z)) * R41 * (r4 * s41_1 - r1 * s41_4);

    double G12_x = R12 * R12 * (r1_x * r2 + r1 * r2_x) + 2 * r1 * r2 * R12_x * R12 + p.z * p.z * (s12_1_x * s12_2 + s12_1 * s12_1_x);
    double G23_x = R23 * R23 * (r2_x * r3 + r2 * r3_x) + 2 * r2 * r3 * R23_x * R23 + p.z * p.z * (s23_2_x * s23_3 + s23_2 * s23_2_x);
    double G34_x = R34 * R34 * (r3_x * r4 + r3 * r4_x) + 2 * r3 * r4 * R34_x * R34 + p.z * p.z * (s34_3_x * s34_4 + s34_3 * s34_3_x);
    double G41_x = R41 * R41 * (r4_x * r1 + r4 * r1_x) + 2 * r4 * r1 * R41_x * R41 + p.z * p.z * (s41_4_x * s41_1 + s41_4 * s41_4_x);

    double G12_y = R12 * R12 * (r1_y * r2 + r1 * r2_y) + 2 * r1 * r2 * R12_y * R12 + p.z * p.z * (s12_1_y * s12_2 + s12_1 * s12_1_y);
    double G23_y = R23 * R23 * (r2_y * r3 + r2 * r3_y) + 2 * r2 * r3 * R23_y * R23 + p.z * p.z * (s23_2_y * s23_3 + s23_2 * s23_2_y);
    double G34_y = R34 * R34 * (r3_y * r4 + r3 * r4_y) + 2 * r3 * r4 * R34_y * R34 + p.z * p.z * (s34_3_y * s34_4 + s34_3 * s34_3_y);
    double G41_y = R41 * R41 * (r4_y * r1 + r4 * r1_y) + 2 * r4 * r1 * R41_y * R41 + p.z * p.z * (s41_4_y * s41_1 + s41_4 * s41_4_y);
    
    double G12_z = R12 * R12 * (r1_z * r2 + r1 * r1_z) + 2 * p.z * s12_1 * s12_2;
    double G23_z = R23 * R23 * (r2_z * r3 + r2 * r2_z) + 2 * p.z * s23_2 * s23_3;
    double G34_z = R34 * R34 * (r3_z * r4 + r3 * r3_z) + 2 * p.z * s34_3 * s34_4;
    double G41_z = R41 * R41 * (r4_z * r1 + r4 * r4_z) + 2 * p.z * s41_4 * s41_1;

    double J12_x = ((F12_x * G12 - F12 * G12_x) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12));
    double J23_x = ((F23_x * G23 - F23 * G23_x) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23));
    double J34_x = ((F34_x * G34 - F34 * G34_x) / (G34 * G34)) / (1 + (F34 / G34) * (F34 / G34));
    double J41_x = ((F41_x * G41 - F41 * G41_x) / (G41 * G41)) / (1 + (F41 / G41) * (F41 / G41));

    double J12_y = ((F12_y * G12 - F12 * G12_y) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12));
    double J23_y = ((F23_y * G23 - F23 * G23_y) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23));
    double J34_y = ((F34_y * G34 - F34 * G34_y) / (G34 * G34)) / (1 + (F34 / G34) * (F34 / G34));
    double J41_y = ((F41_y * G41 - F41 * G41_y) / (G41 * G41)) / (1 + (F41 / G41) * (F41 / G41));

    double J12_z = ((F12_z * G12 - F12 * G12_z) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12));
    double J23_z = ((F23_z * G23 - F23 * G23_z) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23));
    double J34_z = ((F34_z * G34 - F34 * G34_z) / (G34 * G34)) / (1 + (F34 / G34) * (F34 / G34));
    double J41_z = ((F41_z * G41 - F41 * G41_z) / (G41 * G41)) / (1 + (F41 / G41) * (F41 / G41));

    double delta;
    if ((R12 < 0.0) && (R23 < 0.0) && (R34 < 0.0) && (R41 < 0.0)) {
        delta = 2 * M_PI;
    } else {
        delta = 0.0;
    }

    // source
    double source_u = (S12 * Q12 + S23 * Q23 + S34 * Q34 + S41 * Q41);
    double source_v = - (C12 * Q12 + C23 * Q23 + C34 * Q34 + C41 * Q41);
    double source_w = (p.z / fabs(p.z)) * (delta + J12 + J23 + J34 + J41);

    source_vel->x = source_u * e1.x + source_v * e2.x + source_w * e3.x;
    source_vel->y = source_u * e1.y + source_v * e2.y + source_w * e3.y;
    source_vel->z = source_u * e1.z + source_v * e2.z + source_w * e3.z;

    // doublet
    double doublet_u = - (p.z / fabs(p.z)) * (J12_x + J23_x + J34_x + J41_x);
    double doublet_v = - (p.z / fabs(p.z)) * (J12_y + J23_y + J34_y + J41_y);
    double doublet_w = - (p.z / fabs(p.z)) * (J12_z + J23_z + J34_z + J41_z);

    doublet_vel->x = doublet_u * e1.x + doublet_v * e2.x + doublet_w * e3.x;
    doublet_vel->y = doublet_u * e1.y + doublet_v * e2.y + doublet_w * e3.y;
    doublet_vel->z = doublet_u * e1.z + doublet_v * e2.z + doublet_w * e3.z;

}

void tri_planar_panel_velocity(Vec3D p,
                               Vec2D p1, Vec2D p2, Vec2D p3,
                               Vec3D e1, Vec3D e2, Vec3D e3,
                               Vec3D *source_vel, Vec3D *doublet_vel)
{
    double d12 = sqrt( (p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y) );
    double d23 = sqrt( (p3.x - p2.x) * (p3.x - p2.x) + (p3.y - p2.y) * (p3.y - p2.y) );
    double d31 = sqrt( (p1.x - p3.x) * (p1.x - p3.x) + (p1.y - p3.y) * (p1.y - p3.y) );

    double C12 = (p2.x - p1.x) / d12;
    double C23 = (p3.x - p2.x) / d23;
    double C31 = (p1.x - p3.x) / d31;

    double S12 = (p2.y - p1.y) / d12;
    double S23 = (p3.y - p2.y) / d23;
    double S31 = (p1.y - p3.y) / d31;

    double s12_1 = (p1.x - p.x) * C12 + (p1.y - p.y) * S12;
    double s12_2 = (p2.x - p.x) * C12 + (p2.y - p.y) * S12;
    double s23_2 = (p2.x - p.x) * C23 + (p2.y - p.y) * S23;
    double s23_3 = (p3.x - p.x) * C23 + (p3.y - p.y) * S23;
    double s31_3 = (p3.x - p.x) * C31 + (p3.y - p.y) * S31;
    double s31_1 = (p1.x - p.x) * C31 + (p1.y - p.y) * S31;

    double R12 = (p.x - p1.x) * S12 - (p.y - p1.y) * C12;
    double R23 = (p.x - p2.x) * S23 - (p.y - p2.y) * C23;
    double R31 = (p.x - p3.x) * S31 - (p.y - p3.y) * C31;

    double r1 = sqrt((p.x - p1.x) * (p.x - p1.x) + (p.y - p1.y) * (p.y - p1.y) + p.z * p.z);
    double r2 = sqrt((p.x - p2.x) * (p.x - p2.x) + (p.y - p2.y) * (p.y - p2.y) + p.z * p.z);
    double r3 = sqrt((p.x - p3.x) * (p.x - p3.x) + (p.y - p3.y) * (p.y - p3.y) + p.z * p.z);

    double Q12 = log((r1 + r2 + d12) / (r1 + r2 - d12));
    double Q23 = log((r2 + r3 + d23) / (r2 + r3 - d23));
    double Q31 = log((r3 + r1 + d31) / (r3 + r1 - d31));

    double F12 = R12 * fabs(p.z) * (r1 * s12_2 - r2 * s12_1);
    double F23 = R23 * fabs(p.z) * (r2 * s23_3 - r3 * s23_2);
    double F31 = R31 * fabs(p.z) * (r3 * s31_1 - r1 * s31_3);

    double G12 = r1 * r2 * R12 * R12 + p.z * p.z * s12_1 * s12_2;
    double G23 = r2 * r3 * R23 * R23 + p.z * p.z * s23_2 * s23_3;
    double G31 = r3 * r1 * R31 * R31 + p.z * p.z * s31_3 * s31_1;

    double J12 = atan2(F12, G12);
    double J23 = atan2(F23, G23);
    double J31 = atan2(F31, G31);

    double R12_x = S12; double R12_y = - C12;
    double R23_x = S23; double R23_y = - C23;
    double R31_x = S31; double R31_y = - C31;

    double r1_x = (p.x - p1.x) / r1; double r1_y = (p.y - p1.y) / r1; double r1_z = p.z / r1;
    double r2_x = (p.x - p2.x) / r2; double r2_y = (p.y - p2.y) / r2; double r2_z = p.z / r2;
    double r3_x = (p.x - p3.x) / r3; double r3_y = (p.y - p3.y) / r3; double r3_z = p.z / r3;

    double s12_1_x = - C12; double s12_2_x = - C12; double s12_1_y = - S12; double s12_2_y = - S12;
    double s23_2_x = - C23; double s23_3_x = - C23; double s23_2_y = - S23; double s23_3_y = - S23;
    double s31_3_x = - C31; double s31_1_x = - C31; double s31_3_y = - S31; double s31_1_y = - S31;

    double F12_x = fabs(p.z) * (R12_x * (r1 * s12_2 - r2 * s12_1) + R12 * (r1_x * s12_2 + r1 * s12_2_x - (r2_x * s12_1 + r2 * s12_1_x)));
    double F23_x = fabs(p.z) * (R23_x * (r2 * s23_3 - r3 * s23_2) + R23 * (r2_x * s23_3 + r2 * s23_3_x - (r3_x * s23_2 + r3 * s23_2_x)));
    double F31_x = fabs(p.z) * (R31_x * (r3 * s31_1 - r1 * s31_3) + R31 * (r3_x * s31_1 + r3 * s31_1_x - (r1_x * s31_3 + r1 * s31_3_x)));

    double F12_y = fabs(p.z) * (R12_y * (r1 * s12_2 - r2 * s12_1) + R12 * (r1_y * s12_2 + r1 * s12_2_y - (r2_y * s12_1 + r2 * s12_1_y)));
    double F23_y = fabs(p.z) * (R23_y * (r2 * s23_3 - r3 * s23_2) + R23 * (r2_y * s23_3 + r2 * s23_3_y - (r3_y * s23_2 + r3 * s23_2_y)));
    double F31_y = fabs(p.z) * (R31_y * (r3 * s31_1 - r1 * s31_3) + R31 * (r3_y * s31_1 + r3 * s31_1_y - (r1_y * s31_3 + r1 * s31_3_y)));

    double F12_z = (p.z / fabs(p.z)) * R12 * (r1 * s12_2 - r2 * s12_1);
    double F23_z = (p.z / fabs(p.z)) * R23 * (r2 * s23_3 - r3 * s23_2);
    double F31_z = (p.z / fabs(p.z)) * R31 * (r3 * s31_1 - r1 * s31_3);

    double G12_x = R12 * R12 * (r1_x * r2 + r1 * r2_x) + 2 * r1 * r2 * R12_x * R12 + p.z * p.z * (s12_1_x * s12_2 + s12_1 * s12_1_x);
    double G23_x = R23 * R23 * (r2_x * r3 + r2 * r3_x) + 2 * r2 * r3 * R23_x * R23 + p.z * p.z * (s23_2_x * s23_3 + s23_2 * s23_2_x);
    double G31_x = R31 * R31 * (r3_x * r1 + r3 * r1_x) + 2 * r3 * r1 * R31_x * R31 + p.z * p.z * (s31_3_x * s31_1 + s31_3 * s31_3_x);

    double G12_y = R12 * R12 * (r1_y * r2 + r1 * r2_y) + 2 * r1 * r2 * R12_y * R12 + p.z * p.z * (s12_1_y * s12_2 + s12_1 * s12_1_y);
    double G23_y = R23 * R23 * (r2_y * r3 + r2 * r3_y) + 2 * r2 * r3 * R23_y * R23 + p.z * p.z * (s23_2_y * s23_3 + s23_2 * s23_2_y);
    double G31_y = R31 * R31 * (r3_y * r1 + r3 * r1_y) + 2 * r3 * r1 * R31_y * R31 + p.z * p.z * (s31_3_y * s31_1 + s31_3 * s31_3_y);

    double G12_z = R12 * R12 * (r1_z * r2 + r1 * r1_z) + 2 * p.z * s12_1 * s12_2;
    double G23_z = R23 * R23 * (r2_z * r3 + r2 * r2_z) + 2 * p.z * s23_2 * s23_3;
    double G31_z = R31 * R31 * (r3_z * r1 + r3 * r3_z) + 2 * p.z * s31_3 * s31_1;

    double J12_x = ((F12_x * G12 - F12 * G12_x) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12));
    double J23_x = ((F23_x * G23 - F23 * G23_x) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23));
    double J31_x = ((F31_x * G31 - F31 * G31_x) / (G31 * G31)) / (1 + (F31 / G31) * (F31 / G31));

    double J12_y = ((F12_y * G12 - F12 * G12_y) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12));
    double J23_y = ((F23_y * G23 - F23 * G23_y) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23));
    double J31_y = ((F31_y * G31 - F31 * G31_y) / (G31 * G31)) / (1 + (F31 / G31) * (F31 / G31));

    double J12_z = ((F12_z * G12 - F12 * G12_z) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12));
    double J23_z = ((F23_z * G23 - F23 * G23_z) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23));
    double J31_z = ((F31_z * G31 - F31 * G31_z) / (G31 * G31)) / (1 + (F31 / G31) * (F31 / G31));

    double delta;
    if ((R12 < 0.0) && (R23 < 0.0) && (R31 < 0.0)) {
        delta = 2 * M_PI;
    } else {
        delta = 0.0;
    }

    // source
    double source_u = (S12 * Q12 + S23 * Q23 + S31 * Q31);
    double source_v = - (C12 * Q12 + C23 * Q23 + C31 * Q31);
    double source_w = (p.z / fabs(p.z)) * (delta + J12 + J23 + J31);

    source_vel->x = source_u * e1.x + source_v * e2.x + source_w * e3.x;
    source_vel->y = source_u * e1.y + source_v * e2.y + source_w * e3.y;
    source_vel->z = source_u * e1.z + source_v * e2.z + source_w * e3.z;

    // doublet
    double doublet_u = - (p.z / fabs(p.z)) * (J12_x + J23_x + J31_x);
    double doublet_v = - (p.z / fabs(p.z)) * (J12_y + J23_y + J31_y);
    double doublet_w = - (p.z / fabs(p.z)) * (J12_z + J23_z + J31_z);

    doublet_vel->x = doublet_u * e1.x + doublet_v * e2.x + doublet_w * e3.x;
    doublet_vel->y = doublet_u * e1.y + doublet_v * e2.y + doublet_w * e3.y;
    doublet_vel->z = doublet_u * e1.z + doublet_v * e2.z + doublet_w * e3.z;

}


void line_vortex(Vec3D p,
                 Vec3D p1, Vec3D p2,
                 Vec3D *vel)
{

    Vec3D r1 = {p.x - p1.x, p.y - p1.y, p.z - p1.z};
    Vec3D r2 = {p.x - p2.x, p.y - p2.y, p.z - p2.z};
    Vec3D r0 = {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};

    Vec3D r1xr2 = {r1.y * r2.z - r1.z * r2.y, r1.z * r2.x - r1.x * r2.z, r1.x * r2.y - r1.y * r2.x};

    double r1xr2_square = r1xr2.x * r1xr2.x + r1xr2.y * r1xr2.y + r1xr2.z * r1xr2.z;

    double r1_norm = sqrt(r1.x * r1.x + r1.y * r1.y + r1.z * r1.z);
    double r2_norm = sqrt(r2.x * r2.x + r2.y * r2.y + r2.z * r2.z);

    if ((r1_norm < ZERO_ERROR) || (r2_norm < ZERO_ERROR) || (r1xr2_square < ZERO_ERROR))
    {
        vel->x = 0.0;
        vel->y = 0.0;
        vel->z = 0.0;
    }
    else
    {
        double r0r1 = r0.x * r1.x + r0.y * r1.y + r0.z * r1.z;
        double r0r2 = r0.x * r2.x + r0.y * r2.y + r0.z * r2.z;
        
        double k = (1 / (4 * M_PI * r1xr2_square)) * (r0r1 / r1_norm - r0r2 / r2_norm);

        vel->x = k * r1xr2.x;
        vel->y = k * r1xr2.y;
        vel->z = k * r1xr2.z;
    }

}

/*-------------------------------------------
    INDUCED POTENTIAL
--------------------------------------------*/
void quad_planar_panel_potential(Vec3D p,
                                Vec2D p1, Vec2D p2, Vec2D p3, Vec2D p4,
                                double *source_pot, double *doublet_pot)
{
    double d12 = sqrt( (p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y) );
    double d23 = sqrt( (p3.x - p2.x) * (p3.x - p2.x) + (p3.y - p2.y) * (p3.y - p2.y) );
    double d34 = sqrt( (p4.x - p3.x) * (p4.x - p3.x) + (p4.y - p3.y) * (p4.y - p3.y) );
    double d41 = sqrt( (p1.x - p4.x) * (p1.x - p4.x) + (p1.y - p4.y) * (p1.y - p4.y) );

    double C12 = (p2.x - p1.x) / d12;
    double C23 = (p3.x - p2.x) / d23;
    double C34 = (p4.x - p3.x) / d34;
    double C41 = (p1.x - p4.x) / d41;

    double S12 = (p2.y - p1.y) / d12;
    double S23 = (p3.y - p2.y) / d23;
    double S34 = (p4.y - p3.y) / d34;
    double S41 = (p1.y - p4.y) / d41;

    double s12_1 = (p1.x - p.x) * C12 + (p1.y - p.y) * S12;
    double s12_2 = (p2.x - p.x) * C12 + (p2.y - p.y) * S12;
    double s23_2 = (p2.x - p.x) * C23 + (p2.y - p.y) * S23;
    double s23_3 = (p3.x - p.x) * C23 + (p3.y - p.y) * S23;
    double s34_3 = (p3.x - p.x) * C34 + (p3.y - p.y) * S34;
    double s34_4 = (p4.x - p.x) * C34 + (p4.y - p.y) * S34;
    double s41_4 = (p4.x - p.x) * C41 + (p4.y - p.y) * S41;
    double s41_1 = (p1.x - p.x) * C41 + (p1.y - p.y) * S41;

    double R12 = (p.x - p1.x) * S12 - (p.y - p1.y) * C12;
    double R23 = (p.x - p2.x) * S23 - (p.y - p2.y) * C23;
    double R34 = (p.x - p3.x) * S34 - (p.y - p3.y) * C34;
    double R41 = (p.x - p4.x) * S41 - (p.y - p4.y) * C41;

    double r1 = sqrt((p.x - p1.x) * (p.x - p1.x) + (p.y - p1.y) * (p.y - p1.y) + p.z * p.z);
    double r2 = sqrt((p.x - p2.x) * (p.x - p2.x) + (p.y - p2.y) * (p.y - p2.y) + p.z * p.z);
    double r3 = sqrt((p.x - p3.x) * (p.x - p3.x) + (p.y - p3.y) * (p.y - p3.y) + p.z * p.z);
    double r4 = sqrt((p.x - p4.x) * (p.x - p4.x) + (p.y - p4.y) * (p.y - p4.y) + p.z * p.z);

    double Q12 = log((r1 + r2 + d12) / (r1 + r2 - d12));
    double Q23 = log((r2 + r3 + d23) / (r2 + r3 - d23));
    double Q34 = log((r3 + r4 + d34) / (r3 + r4 - d34));
    double Q41 = log((r4 + r1 + d41) / (r4 + r1 - d41));

    double F12 = R12 * fabs(p.z) * (r1 * s12_2 - r2 * s12_1);
    double F23 = R23 * fabs(p.z) * (r2 * s23_3 - r3 * s23_2);
    double F34 = R34 * fabs(p.z) * (r3 * s34_4 - r4 * s34_3);
    double F41 = R41 * fabs(p.z) * (r4 * s41_1 - r1 * s41_4);

    double G12 = r1 * r2 * R12 * R12 + p.z * p.z * s12_1 * s12_2;
    double G23 = r2 * r3 * R23 * R23 + p.z * p.z * s23_2 * s23_3;
    double G34 = r3 * r4 * R34 * R34 + p.z * p.z * s34_3 * s34_4;
    double G41 = r4 * r1 * R41 * R41 + p.z * p.z * s41_4 * s41_1;

    double J12 = atan2(F12, G12);
    double J23 = atan2(F23, G23);
    double J34 = atan2(F34, G34);
    double J41 = atan2(F41, G41);

    double R12_x = S12; double R12_y = - C12;
    double R23_x = S23; double R23_y = - C23;
    double R34_x = S34; double R34_y = - C34;
    double R41_x = S41; double R41_y = - C41;

    double r1_x = (p.x - p1.x) / r1; double r1_y = (p.y - p1.y) / r1; double r1_z = p.z / r1;
    double r2_x = (p.x - p2.x) / r2; double r2_y = (p.y - p2.y) / r2; double r2_z = p.z / r2;
    double r3_x = (p.x - p3.x) / r3; double r3_y = (p.y - p3.y) / r3; double r3_z = p.z / r3;
    double r4_x = (p.x - p4.x) / r4; double r4_y = (p.y - p4.y) / r4; double r4_z = p.z / r4;

    double s12_1_x = - C12; double s12_2_x = - C12; double s12_1_y = - S12; double s12_2_y = - S12;
    double s23_2_x = - C23; double s23_3_x = - C23; double s23_2_y = - S23; double s23_3_y = - S23;
    double s34_3_x = - C34; double s34_4_x = - C34; double s34_3_y = - S34; double s34_4_y = - S34;
    double s41_4_x = - C41; double s41_1_x = - C41; double s41_4_y = - S41; double s41_1_y = - S41;

    double F12_x = fabs(p.z) * (R12_x * (r1 * s12_2 - r2 * s12_1) + R12 * (r1_x * s12_2 + r1 * s12_2_x - (r2_x * s12_1 + r2 * s12_1_x)));
    double F23_x = fabs(p.z) * (R23_x * (r2 * s23_3 - r3 * s23_2) + R23 * (r2_x * s23_3 + r2 * s23_3_x - (r3_x * s23_2 + r3 * s23_2_x)));
    double F34_x = fabs(p.z) * (R34_x * (r3 * s34_4 - r4 * s34_3) + R34 * (r3_x * s34_4 + r3 * s34_4_x - (r4_x * s34_3 + r4 * s34_3_x)));
    double F41_x = fabs(p.z) * (R41_x * (r4 * s41_1 - r1 * s41_4) + R41 * (r4_x * s41_1 + r4 * s41_1_x - (r1_x * s41_4 + r1 * s41_4_x)));

    double F12_y = fabs(p.z) * (R12_y * (r1 * s12_2 - r2 * s12_1) + R12 * (r1_y * s12_2 + r1 * s12_2_y - (r2_y * s12_1 + r2 * s12_1_y)));
    double F23_y = fabs(p.z) * (R23_y * (r2 * s23_3 - r3 * s23_2) + R23 * (r2_y * s23_3 + r2 * s23_3_y - (r3_y * s23_2 + r3 * s23_2_y)));
    double F34_y = fabs(p.z) * (R34_y * (r3 * s34_4 - r4 * s34_3) + R34 * (r3_y * s34_4 + r3 * s34_4_y - (r4_y * s34_3 + r4 * s34_3_y)));
    double F41_y = fabs(p.z) * (R41_y * (r4 * s41_1 - r1 * s41_4) + R41 * (r4_y * s41_1 + r4 * s41_1_y - (r1_y * s41_4 + r1 * s41_4_y)));

    double F12_z = (p.z / fabs(p.z)) * R12 * (r1 * s12_2 - r2 * s12_1);
    double F23_z = (p.z / fabs(p.z)) * R23 * (r2 * s23_3 - r3 * s23_2);
    double F34_z = (p.z / fabs(p.z)) * R34 * (r3 * s34_4 - r4 * s34_3);
    double F41_z = (p.z / fabs(p.z)) * R41 * (r4 * s41_1 - r1 * s41_4);

    double G12_x = R12 * R12 * (r1_x * r2 + r1 * r2_x) + 2 * r1 * r2 * R12_x * R12 + p.z * p.z * (s12_1_x * s12_2 + s12_1 * s12_1_x);
    double G23_x = R23 * R23 * (r2_x * r3 + r2 * r3_x) + 2 * r2 * r3 * R23_x * R23 + p.z * p.z * (s23_2_x * s23_3 + s23_2 * s23_2_x);
    double G34_x = R34 * R34 * (r3_x * r4 + r3 * r4_x) + 2 * r3 * r4 * R34_x * R34 + p.z * p.z * (s34_3_x * s34_4 + s34_3 * s34_3_x);
    double G41_x = R41 * R41 * (r4_x * r1 + r4 * r1_x) + 2 * r4 * r1 * R41_x * R41 + p.z * p.z * (s41_4_x * s41_1 + s41_4 * s41_4_x);

    double G12_y = R12 * R12 * (r1_y * r2 + r1 * r2_y) + 2 * r1 * r2 * R12_y * R12 + p.z * p.z * (s12_1_y * s12_2 + s12_1 * s12_1_y);
    double G23_y = R23 * R23 * (r2_y * r3 + r2 * r3_y) + 2 * r2 * r3 * R23_y * R23 + p.z * p.z * (s23_2_y * s23_3 + s23_2 * s23_2_y);
    double G34_y = R34 * R34 * (r3_y * r4 + r3 * r4_y) + 2 * r3 * r4 * R34_y * R34 + p.z * p.z * (s34_3_y * s34_4 + s34_3 * s34_3_y);
    double G41_y = R41 * R41 * (r4_y * r1 + r4 * r1_y) + 2 * r4 * r1 * R41_y * R41 + p.z * p.z * (s41_4_y * s41_1 + s41_4 * s41_4_y);
    
    double G12_z = R12 * R12 * (r1_z * r2 + r1 * r1_z) + 2 * p.z * s12_1 * s12_2;
    double G23_z = R23 * R23 * (r2_z * r3 + r2 * r2_z) + 2 * p.z * s23_2 * s23_3;
    double G34_z = R34 * R34 * (r3_z * r4 + r3 * r3_z) + 2 * p.z * s34_3 * s34_4;
    double G41_z = R41 * R41 * (r4_z * r1 + r4 * r4_z) + 2 * p.z * s41_4 * s41_1;

    double J12_x = ((F12_x * G12 - F12 * G12_x) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12));
    double J23_x = ((F23_x * G23 - F23 * G23_x) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23));
    double J34_x = ((F34_x * G34 - F34 * G34_x) / (G34 * G34)) / (1 + (F34 / G34) * (F34 / G34));
    double J41_x = ((F41_x * G41 - F41 * G41_x) / (G41 * G41)) / (1 + (F41 / G41) * (F41 / G41));

    double J12_y = ((F12_y * G12 - F12 * G12_y) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12));
    double J23_y = ((F23_y * G23 - F23 * G23_y) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23));
    double J34_y = ((F34_y * G34 - F34 * G34_y) / (G34 * G34)) / (1 + (F34 / G34) * (F34 / G34));
    double J41_y = ((F41_y * G41 - F41 * G41_y) / (G41 * G41)) / (1 + (F41 / G41) * (F41 / G41));

    double J12_z = ((F12_z * G12 - F12 * G12_z) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12));
    double J23_z = ((F23_z * G23 - F23 * G23_z) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23));
    double J34_z = ((F34_z * G34 - F34 * G34_z) / (G34 * G34)) / (1 + (F34 / G34) * (F34 / G34));
    double J41_z = ((F41_z * G41 - F41 * G41_z) / (G41 * G41)) / (1 + (F41 / G41) * (F41 / G41));

    double delta;
    if ((R12 < 0.0) && (R23 < 0.0) && (R34 < 0.0) && (R41 < 0.0)) {
        delta = 2 * M_PI;
    } else {
        delta = 0.0;
    }

    // source
    *source_pot = R12 * Q12 + R23 * Q23 + R34 * Q34 + R41 * Q41 + fabs(p.z) * (delta + J12 + J23 + J34 + J41);

    // doublet
    *doublet_pot = - (p.z / fabs(p.z)) * (delta + J12 + J23 + J34 + J41);

}

void tri_planar_panel_potential(Vec3D p,
                               Vec2D p1, Vec2D p2, Vec2D p3,
                               double *source_pot, double *doublet_pot)
{
    double d12 = sqrt( (p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y) );
    double d23 = sqrt( (p3.x - p2.x) * (p3.x - p2.x) + (p3.y - p2.y) * (p3.y - p2.y) );
    double d31 = sqrt( (p1.x - p3.x) * (p1.x - p3.x) + (p1.y - p3.y) * (p1.y - p3.y) );

    double C12 = (p2.x - p1.x) / d12;
    double C23 = (p3.x - p2.x) / d23;
    double C31 = (p1.x - p3.x) / d31;

    double S12 = (p2.y - p1.y) / d12;
    double S23 = (p3.y - p2.y) / d23;
    double S31 = (p1.y - p3.y) / d31;

    double s12_1 = (p1.x - p.x) * C12 + (p1.y - p.y) * S12;
    double s12_2 = (p2.x - p.x) * C12 + (p2.y - p.y) * S12;
    double s23_2 = (p2.x - p.x) * C23 + (p2.y - p.y) * S23;
    double s23_3 = (p3.x - p.x) * C23 + (p3.y - p.y) * S23;
    double s31_3 = (p3.x - p.x) * C31 + (p3.y - p.y) * S31;
    double s31_1 = (p1.x - p.x) * C31 + (p1.y - p.y) * S31;

    double R12 = (p.x - p1.x) * S12 - (p.y - p1.y) * C12;
    double R23 = (p.x - p2.x) * S23 - (p.y - p2.y) * C23;
    double R31 = (p.x - p3.x) * S31 - (p.y - p3.y) * C31;

    double r1 = sqrt((p.x - p1.x) * (p.x - p1.x) + (p.y - p1.y) * (p.y - p1.y) + p.z * p.z);
    double r2 = sqrt((p.x - p2.x) * (p.x - p2.x) + (p.y - p2.y) * (p.y - p2.y) + p.z * p.z);
    double r3 = sqrt((p.x - p3.x) * (p.x - p3.x) + (p.y - p3.y) * (p.y - p3.y) + p.z * p.z);

    double Q12 = log((r1 + r2 + d12) / (r1 + r2 - d12));
    double Q23 = log((r2 + r3 + d23) / (r2 + r3 - d23));
    double Q31 = log((r3 + r1 + d31) / (r3 + r1 - d31));

    double F12 = R12 * fabs(p.z) * (r1 * s12_2 - r2 * s12_1);
    double F23 = R23 * fabs(p.z) * (r2 * s23_3 - r3 * s23_2);
    double F31 = R31 * fabs(p.z) * (r3 * s31_1 - r1 * s31_3);

    double G12 = r1 * r2 * R12 * R12 + p.z * p.z * s12_1 * s12_2;
    double G23 = r2 * r3 * R23 * R23 + p.z * p.z * s23_2 * s23_3;
    double G31 = r3 * r1 * R31 * R31 + p.z * p.z * s31_3 * s31_1;

    double J12 = atan2(F12, G12);
    double J23 = atan2(F23, G23);
    double J31 = atan2(F31, G31);

    double R12_x = S12; double R12_y = - C12;
    double R23_x = S23; double R23_y = - C23;
    double R31_x = S31; double R31_y = - C31;

    double r1_x = (p.x - p1.x) / r1; double r1_y = (p.y - p1.y) / r1; double r1_z = p.z / r1;
    double r2_x = (p.x - p2.x) / r2; double r2_y = (p.y - p2.y) / r2; double r2_z = p.z / r2;
    double r3_x = (p.x - p3.x) / r3; double r3_y = (p.y - p3.y) / r3; double r3_z = p.z / r3;

    double s12_1_x = - C12; double s12_2_x = - C12; double s12_1_y = - S12; double s12_2_y = - S12;
    double s23_2_x = - C23; double s23_3_x = - C23; double s23_2_y = - S23; double s23_3_y = - S23;
    double s31_3_x = - C31; double s31_1_x = - C31; double s31_3_y = - S31; double s31_1_y = - S31;

    double F12_x = fabs(p.z) * (R12_x * (r1 * s12_2 - r2 * s12_1) + R12 * (r1_x * s12_2 + r1 * s12_2_x - (r2_x * s12_1 + r2 * s12_1_x)));
    double F23_x = fabs(p.z) * (R23_x * (r2 * s23_3 - r3 * s23_2) + R23 * (r2_x * s23_3 + r2 * s23_3_x - (r3_x * s23_2 + r3 * s23_2_x)));
    double F31_x = fabs(p.z) * (R31_x * (r3 * s31_1 - r1 * s31_3) + R31 * (r3_x * s31_1 + r3 * s31_1_x - (r1_x * s31_3 + r1 * s31_3_x)));

    double F12_y = fabs(p.z) * (R12_y * (r1 * s12_2 - r2 * s12_1) + R12 * (r1_y * s12_2 + r1 * s12_2_y - (r2_y * s12_1 + r2 * s12_1_y)));
    double F23_y = fabs(p.z) * (R23_y * (r2 * s23_3 - r3 * s23_2) + R23 * (r2_y * s23_3 + r2 * s23_3_y - (r3_y * s23_2 + r3 * s23_2_y)));
    double F31_y = fabs(p.z) * (R31_y * (r3 * s31_1 - r1 * s31_3) + R31 * (r3_y * s31_1 + r3 * s31_1_y - (r1_y * s31_3 + r1 * s31_3_y)));

    double F12_z = (p.z / fabs(p.z)) * R12 * (r1 * s12_2 - r2 * s12_1);
    double F23_z = (p.z / fabs(p.z)) * R23 * (r2 * s23_3 - r3 * s23_2);
    double F31_z = (p.z / fabs(p.z)) * R31 * (r3 * s31_1 - r1 * s31_3);

    double G12_x = R12 * R12 * (r1_x * r2 + r1 * r2_x) + 2 * r1 * r2 * R12_x * R12 + p.z * p.z * (s12_1_x * s12_2 + s12_1 * s12_1_x);
    double G23_x = R23 * R23 * (r2_x * r3 + r2 * r3_x) + 2 * r2 * r3 * R23_x * R23 + p.z * p.z * (s23_2_x * s23_3 + s23_2 * s23_2_x);
    double G31_x = R31 * R31 * (r3_x * r1 + r3 * r1_x) + 2 * r3 * r1 * R31_x * R31 + p.z * p.z * (s31_3_x * s31_1 + s31_3 * s31_3_x);

    double G12_y = R12 * R12 * (r1_y * r2 + r1 * r2_y) + 2 * r1 * r2 * R12_y * R12 + p.z * p.z * (s12_1_y * s12_2 + s12_1 * s12_1_y);
    double G23_y = R23 * R23 * (r2_y * r3 + r2 * r3_y) + 2 * r2 * r3 * R23_y * R23 + p.z * p.z * (s23_2_y * s23_3 + s23_2 * s23_2_y);
    double G31_y = R31 * R31 * (r3_y * r1 + r3 * r1_y) + 2 * r3 * r1 * R31_y * R31 + p.z * p.z * (s31_3_y * s31_1 + s31_3 * s31_3_y);

    double G12_z = R12 * R12 * (r1_z * r2 + r1 * r1_z) + 2 * p.z * s12_1 * s12_2;
    double G23_z = R23 * R23 * (r2_z * r3 + r2 * r2_z) + 2 * p.z * s23_2 * s23_3;
    double G31_z = R31 * R31 * (r3_z * r1 + r3 * r3_z) + 2 * p.z * s31_3 * s31_1;

    double J12_x = ((F12_x * G12 - F12 * G12_x) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12));
    double J23_x = ((F23_x * G23 - F23 * G23_x) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23));
    double J31_x = ((F31_x * G31 - F31 * G31_x) / (G31 * G31)) / (1 + (F31 / G31) * (F31 / G31));

    double J12_y = ((F12_y * G12 - F12 * G12_y) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12));
    double J23_y = ((F23_y * G23 - F23 * G23_y) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23));
    double J31_y = ((F31_y * G31 - F31 * G31_y) / (G31 * G31)) / (1 + (F31 / G31) * (F31 / G31));

    double J12_z = ((F12_z * G12 - F12 * G12_z) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12));
    double J23_z = ((F23_z * G23 - F23 * G23_z) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23));
    double J31_z = ((F31_z * G31 - F31 * G31_z) / (G31 * G31)) / (1 + (F31 / G31) * (F31 / G31));

    double delta;
    if ((R12 < 0.0) && (R23 < 0.0) && (R31 < 0.0)) {
        delta = 2 * M_PI;
    } else {
        delta = 0.0;
    }

    // source
    *source_pot = R12 * Q12 + R23 * Q23 + R31 * Q31 + fabs(p.z) * (delta + J12 + J23 + J31);

    // doublet
    *doublet_pot = - ((p.z / fabs(p.z)) * (delta + J12 + J23 + J31));

}

/*-------------------------------------------
    PANELS
--------------------------------------------*/
void planar_potential_panel(int n_sides,
                            Vec3D p_ctrl,
                            Vec2D p1, Vec2D p2, Vec2D p3, Vec2D p4,
                            double *source_pot,
                            double *doublet_pot)
{
    if (n_sides == 4) {
        quad_planar_panel_potential(p_ctrl, p1, p2, p3, p4, source_pot, doublet_pot);
    } else {
        tri_planar_panel_potential(p_ctrl, p1, p2, p3, source_pot, doublet_pot);
    }
}

void planar_velocity_panel(int n_sides,
                           Vec3D p_ctrl,
                           Vec3D e1, Vec3D e2, Vec3D e3,
                           Vec2D p1, Vec2D p2, Vec2D p3, Vec2D p4,
                           Vec3D *source_vel,
                           Vec3D *doublet_vel)
{
    if (n_sides == 4) {
        quad_planar_panel_velocity(p_ctrl, p1, p2, p3, p4, e1, e2, e3, source_vel, doublet_vel);
    } else {
        tri_planar_panel_velocity(p_ctrl, p1, p2, p3, e1, e2, e3, source_vel, doublet_vel);
    }
}

void non_planar_quad_doublet_panel(Vec3D p_ctrl,
                                   Vec3D p1, Vec3D p2, Vec3D p3, Vec3D p4,
                                   Vec3D *doublet_vel)
{
    
    Vec3D vel_1, vel_2, vel_3, vel_4;

    line_vortex(p_ctrl, p1, p2, &vel_1);
    line_vortex(p_ctrl, p2, p3, &vel_2);
    line_vortex(p_ctrl, p3, p4, &vel_3);
    line_vortex(p_ctrl, p4, p1, &vel_4);

    doublet_vel->x = vel_1.x + vel_2.x + vel_3.x + vel_4.x;
    doublet_vel->y = vel_1.y + vel_2.y + vel_3.y + vel_4.y;
    doublet_vel->z = vel_1.z + vel_2.z + vel_3.z + vel_4.z;

}

void get_panel_parameters(double p1x, double p1y, double p1z,
                          double p2x, double p2y, double p2z,
                          double p3x, double p3y, double p3z,
                          double p4x, double p4y, double p4z,
                          Vec3D *e1, Vec3D *e2, Vec3D *e3,
                          Vec2D *p1, Vec2D *p2, Vec2D *p3, Vec2D *p4)
{

    Vec3D vec1, vec2, vec3, vec4, vec5, vec6, vec7, vec8;
    double scalar;
    Vec3D p_avg;

    // Panel center
    p_avg.x = 0.25 * (p1x + p2x + p3x + p4x);
    p_avg.y = 0.25 * (p1y + p2y + p3y + p4y);
    p_avg.z = 0.25 * (p1z + p2z + p3z + p4z);

    // Orthogonal base (e3)
    vec1.x = p2x - p4x; vec1.y = p2y - p4y; vec1.z = p2z - p4z;
    vec2.x = p3x - p1x; vec2.y = p3y - p1y; vec2.z = p3z - p1z;

    vec3 = cross_func(vec1, vec2);
    scalar = norm_func(vec3);

    e3->x = vec3.x / scalar;
    e3->y = vec3.y / scalar;
    e3->z = vec3.z / scalar;

    // Vertices in global system
    vec1.x = p1x - p_avg.x; vec1.y = p1y - p_avg.y; vec1.z = p1z - p_avg.z;
    vec2.x = p2x - p_avg.x; vec2.y = p2y - p_avg.y; vec2.z = p2z - p_avg.z;
    vec3.x = p3x - p_avg.x; vec3.y = p3y - p_avg.y; vec3.z = p3z - p_avg.z;
    vec4.x = p4x - p_avg.x; vec4.y = p4y - p_avg.y; vec4.z = p4z - p_avg.z;

    vec5.x = vec1.x - e3->x * (vec1.x * e3->x + vec1.y * e3->y + vec1.z * e3->z); vec5.y = vec1.y - e3->y * (vec1.x * e3->x + vec1.y * e3->y + vec1.z * e3->z); vec5.z = vec1.z - e3->z * (vec1.x * e3->x + vec1.y * e3->y + vec1.z * e3->z);
    vec6.x = vec2.x - e3->x * (vec2.x * e3->x + vec2.y * e3->y + vec2.z * e3->z); vec6.y = vec2.y - e3->y * (vec2.x * e3->x + vec2.y * e3->y + vec2.z * e3->z); vec6.z = vec2.z - e3->z * (vec2.x * e3->x + vec2.y * e3->y + vec2.z * e3->z);
    vec7.x = vec3.x - e3->x * (vec3.x * e3->x + vec3.y * e3->y + vec3.z * e3->z); vec7.y = vec3.y - e3->y * (vec3.x * e3->x + vec3.y * e3->y + vec3.z * e3->z); vec7.z = vec3.z - e3->z * (vec3.x * e3->x + vec3.y * e3->y + vec3.z * e3->z);
    vec8.x = vec4.x - e3->x * (vec4.x * e3->x + vec4.y * e3->y + vec4.z * e3->z); vec8.y = vec4.y - e3->y * (vec4.x * e3->x + vec4.y * e3->y + vec4.z * e3->z); vec8.z = vec4.z - e3->z * (vec4.x * e3->x + vec4.y * e3->y + vec4.z * e3->z);

    // Orthogonal base (e1)
    vec1.x = vec6.x - vec5.x; vec1.y = vec6.y - vec5.y; vec1.z = vec6.z - vec5.z;
    scalar = norm_func(vec1);
    e1->x = vec1.x / scalar; e1->y = vec1.y / scalar; e1->z = vec1.z / scalar;

    // Orthogonal base (e2)
    vec1.x = e1->x; vec1.y = e1->y; vec1.z = e1->z;
    vec3.x = e3->x; vec3.y = e3->y; vec3.z = e3->z;
    vec2 = cross_func(vec3, vec1);
    scalar = norm_func(vec2);
    vec2.x = vec2.x / scalar; vec2.y = vec2.y / scalar; vec2.z = vec2.z / scalar;
    e2->x = vec2.x; e2->y = vec2.y; e2->z = vec2.z;

    // Vertices in local system
    p1->x = dot_func(vec1, vec5); p1->y = dot_func(vec2, vec5);
    p2->x = dot_func(vec1, vec6); p2->y = dot_func(vec2, vec6);
    p3->x = dot_func(vec1, vec7); p3->y = dot_func(vec2, vec7);
    p4->x = dot_func(vec1, vec8); p4->y = dot_func(vec2, vec8);

}

/*-------------------------------------------
    SURFACE COEFS
--------------------------------------------*/
void set_surface_coefs(int nf,
                       int n_sides[],
                       double p_avg[],
                       double p_ctrl[],
                       double e1[], double e2[], double e3[],
                       double p1[], double p2[], double p3[], double p4[],
                       double a_ij[],
                       double b_ij[])
{
    int i, j;

    Vec3D p_vec;
    Vec2D p1_vec, p2_vec, p3_vec, p4_vec;
    Vec3D e1_vec, e2_vec, e3_vec;

    double a, b;

    for (i = 0; i < nf; i++)
    {

        for (j = 0; j < nf; j++)
        {

            e1_vec.x = e1[3 * j]; e1_vec.y = e1[3 * j + 1]; e1_vec.z = e1[3 * j + 2];
            e2_vec.x = e2[3 * j]; e2_vec.y = e2[3 * j + 1]; e2_vec.z = e2[3 * j + 2];
            e3_vec.x = e3[3 * j]; e3_vec.y = e3[3 * j + 1]; e3_vec.z = e3[3 * j + 2];

            p_vec.x = e1_vec.x * (p_ctrl[3 * i] - p_avg[3 * j]) + e1_vec.y * (p_ctrl[3 * i + 1] - p_avg[3 * j + 1]) + e1_vec.z * (p_ctrl[3 * i + 2] - p_avg[3 * j + 2]);
            p_vec.y = e2_vec.x * (p_ctrl[3 * i] - p_avg[3 * j]) + e2_vec.y * (p_ctrl[3 * i + 1] - p_avg[3 * j + 1]) + e2_vec.z * (p_ctrl[3 * i + 2] - p_avg[3 * j + 2]);
            p_vec.z = e3_vec.x * (p_ctrl[3 * i] - p_avg[3 * j]) + e3_vec.y * (p_ctrl[3 * i + 1] - p_avg[3 * j + 1]) + e3_vec.z * (p_ctrl[3 * i + 2] - p_avg[3 * j + 2]);

            p1_vec.x = p1[2 * j]; p1_vec.y = p1[2 * j + 1];
            p2_vec.x = p2[2 * j]; p2_vec.y = p2[2 * j + 1];
            p3_vec.x = p3[2 * j]; p3_vec.y = p3[2 * j + 1];
            p4_vec.x = p4[2 * j]; p4_vec.y = p4[2 * j + 1];

            planar_potential_panel(n_sides[j], p_vec, p1_vec, p2_vec, p3_vec, p4_vec, &a, &b);

            a_ij[i * nf + j] = a / (4.0 * M_PI);
            b_ij[i * nf + j] = b / (4.0 * M_PI);
        }
    }

}

/*-------------------------------------------
    WAKE VELOCITY
--------------------------------------------*/
void wake_point_velocity(int nte,
                         int section,
                         double vertices[],
                         int faces[],
                         double areas[],
                         double circulations[],
                         double p_ctrl[],
                         double vel[])
{
    
    int i;

    Vec3D p1, p2, p3, p4, p12, p13, p14;
    Vec3D v, p;

    double area;

    p.x = p_ctrl[0]; p.y = p_ctrl[1]; p.z = p_ctrl[2];

    vel[0] = 0.0; vel[1] = 0.0; vel[2] = 0.0;

    for (i = 0; i < nte * (section - 1); i++)
    {

        p1.x = vertices[3 * faces[5 * i + 1]]; p1.y = vertices[3 * faces[5 * i + 1] + 1]; p1.z = vertices[3 * faces[5 * i + 1] + 2];
        p2.x = vertices[3 * faces[5 * i + 2]]; p2.y = vertices[3 * faces[5 * i + 2] + 1]; p2.z = vertices[3 * faces[5 * i + 2] + 2];
        p3.x = vertices[3 * faces[5 * i + 3]]; p3.y = vertices[3 * faces[5 * i + 3] + 1]; p3.z = vertices[3 * faces[5 * i + 3] + 2];
        p4.x = vertices[3 * faces[5 * i + 4]]; p4.y = vertices[3 * faces[5 * i + 4] + 1]; p4.z = vertices[3 * faces[5 * i + 4] + 2];

        non_planar_quad_doublet_panel(p, p1, p2, p3, p4, &v);

        p12.x = p2.x - p1.x; p12.y = p2.y - p1.y; p12.z = p2.z - p1.z;
        p13.x = p3.x - p1.x; p13.y = p3.y - p1.y; p13.z = p3.z - p1.z;
        p14.x = p4.x - p1.x; p14.y = p4.y - p1.y; p14.z = p4.z - p1.z;

        area = 0.5 * (norm_func(cross_func(p12, p13)) + norm_func(cross_func(p13, p14)));
        
        vel[0] = vel[0] + (circulations[i] * areas[i] / area) * v.x;
        vel[1] = vel[1] + (circulations[i] * areas[i] / area) * v.y;
        vel[2] = vel[2] + (circulations[i] * areas[i] / area) * v.z;

    }

}

/*-------------------------------------------
    SURFACE VELOCITY
--------------------------------------------*/
void surface_point_velocity(int nf,
                            int n_sides[],
                            double p_avg[],
                            double e1[], double e2[], double e3[],
                            double p1[], double p2[], double p3[], double p4[],
                            double source[],
                            double doublet[],
                            double p_ctrl[],
                            double vel[])
{
    int i, j;

    Vec3D p_vec;
    Vec2D p1_vec, p2_vec, p3_vec, p4_vec;
    Vec3D e1_vec, e2_vec, e3_vec;

    Vec3D source_vel, doublet_vel;

    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;


    for (j = 0; j < nf; j++)
    {

        e1_vec.x = e1[3 * j]; e1_vec.y = e1[3 * j + 1]; e1_vec.z = e1[3 * j + 2];
        e2_vec.x = e2[3 * j]; e2_vec.y = e2[3 * j + 1]; e2_vec.z = e2[3 * j + 2];
        e3_vec.x = e3[3 * j]; e3_vec.y = e3[3 * j + 1]; e3_vec.z = e3[3 * j + 2];

        p_vec.x = e1_vec.x * (p_ctrl[0] - p_avg[3 * j]) + e1_vec.y * (p_ctrl[1] - p_avg[3 * j + 1]) + e1_vec.z * (p_ctrl[2] - p_avg[3 * j + 2]);
        p_vec.y = e2_vec.x * (p_ctrl[0] - p_avg[3 * j]) + e2_vec.y * (p_ctrl[1] - p_avg[3 * j + 1]) + e2_vec.z * (p_ctrl[2] - p_avg[3 * j + 2]);
        p_vec.z = e3_vec.x * (p_ctrl[0] - p_avg[3 * j]) + e3_vec.y * (p_ctrl[1] - p_avg[3 * j + 1]) + e3_vec.z * (p_ctrl[2] - p_avg[3 * j + 2]);

        p1_vec.x = p1[2 * j]; p1_vec.y = p1[2 * j + 1];
        p2_vec.x = p2[2 * j]; p2_vec.y = p2[2 * j + 1];
        p3_vec.x = p3[2 * j]; p3_vec.y = p3[2 * j + 1];
        p4_vec.x = p4[2 * j]; p4_vec.y = p4[2 * j + 1];

        planar_velocity_panel(n_sides[j], p_vec, e1_vec, e2_vec, e3_vec, p1_vec, p2_vec, p3_vec, p4_vec, &source_vel, &doublet_vel);

        vel[0] = source[j] * source_vel.x + doublet[j] * doublet_vel.x;
        vel[1] = source[j] * source_vel.y + doublet[j] * doublet_vel.y;
        vel[2] = source[j] * source_vel.z + doublet[j] * doublet_vel.z;
            
    }

}

/*-------------------------------------------
    WAKE POTENTIAL
--------------------------------------------*/
void calculate_wake_potential(int nf,
                              int nte,
                              int nw,
                              int faces[],
                              double vertices[],
                              double areas[],
                              double circulations[],
                              double p_ctrl[],
                              double c_ik[],
                              double d_i[])
{
    int i, j;

    int id1, id2, id3, id4;

    Vec3D e1, e2, e3;
    Vec3D p;
    Vec2D p1, p2, p3, p4;
    double area;

    double phi, phi_none;

    for (i = 0; i < nf; i++)
    {
        
        p.x = p_ctrl[3 * i];
        p.y = p_ctrl[3 * i + 1];
        p.z = p_ctrl[3 * i + 2];

        d_i[i] = 0.0;

        for (j = 0; j < nte * nw; j++)
        {

            // Panel parameters
            id1 = faces[5 * j + 1];
            id2 = faces[5 * j + 2];
            id3 = faces[5 * j + 3];
            id4 = faces[5 * j + 4];

            get_panel_parameters(vertices[3 * id1], vertices[3 * id1 + 1], vertices[3 * id1 + 2],
                                 vertices[3 * id2], vertices[3 * id2 + 1], vertices[3 * id2 + 2],
                                 vertices[3 * id3], vertices[3 * id3 + 1], vertices[3 * id3 + 2],
                                 vertices[3 * id4], vertices[3 * id4 + 1], vertices[3 * id4 + 2],
                                 &e1, &e2, &e3,
                                 &p1, &p2, &p3, &p4);
            
            // Potential
            quad_planar_panel_potential(p, p1, p2, p3, p4, &phi, &phi_none);

            // Area
            area = tri_plane_area(p1, p2, p3) + tri_plane_area(p1, p3, p4);

            if (j < nte) { // unknown wake circulation
                c_ik[i * nte + j] = phi;
            } else { // known wake circulation
                phi = phi * circulations[j] * areas[j] / area;
                d_i[i] = d_i[i] + phi;
            }
        }
    }
}

