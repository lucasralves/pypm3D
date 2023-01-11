#include <stdio.h>
#include <math.h>

/*-------------------------------------------
    STRUCTS
--------------------------------------------*/
typedef struct {
    double x;
    double y;
    double z;
} Vec3D;

void print_vec(Vec3D a)
{
    printf("[%.2e, %.2e, %.2e]\n", a.x, a.y, a.z);
}

/*-------------------------------------------
    MATH
--------------------------------------------*/
double norm_func(Vec3D a)
{
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

double dot_func(Vec3D a, Vec3D b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double abs_func(double a)
{
    if (a < 0) {
        return - a;
    } else {
        return a;
    }
}

double division_func(double a, double b)
{
    double b_abs = abs_func(b);

    if (b_abs < 1e-12) {
        return a * b / (b_abs * (b_abs + 1e-12));
    } else {
        return a / b;
    }
}

double distance_func(Vec3D a, Vec3D b)
{
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

void unary_func(Vec3D *a)
{
    double norm = sqrt(a->x * a->x + a->y * a->y + a->z * a->z);
    a->x = a->x / norm;
    a->y = a->y / norm;
    a->z = a->z / norm;
}

Vec3D cross_func(Vec3D a, Vec3D b)
{
    Vec3D out;
    out.x = a.y * b.z - a.z * b.y;
    out.y = a.z * b.x - a.x * b.z;
    out.z = a.x * b.y - a.y * b.x;
    return out;
}

double fabs_func(double a)
{
    if (a < 0) return - a;
    return a;
}

/*-------------------------------------------
    INDUCED VELOCITY
--------------------------------------------*/
const double ZERO_ERROR = 1e-12;
const double FACTOR = 0.25 / M_PI;

Vec3D tri_source_panel(Vec3D p, Vec3D p1, Vec3D p2, Vec3D p3, Vec3D e1, Vec3D e2, Vec3D e3)
{

    Vec3D vel;
    double u, v, w;

    double r1, r2, r3;
    double d12, d23, d31;
    double S12, S23, S31;
    double C12, C23, C31;
    double Q12, Q23, Q31;
    double R12, R23, R31;
    double J12, J23, J31;

    double s12_1, s12_2;
    double s23_2, s23_3;
    double s31_3, s31_1;

    double sign, delta;

    r1 = sqrt((p.x - p1.x) * (p.x - p1.x) + (p.y - p1.y) * (p.y - p1.y) + p.z * p.z);
    r2 = sqrt((p.x - p2.x) * (p.x - p2.x) + (p.y - p2.y) * (p.y - p2.y) + p.z * p.z);
    r3 = sqrt((p.x - p3.x) * (p.x - p3.x) + (p.y - p3.y) * (p.y - p3.y) + p.z * p.z);

    d12 = sqrt((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y));
    d23 = sqrt((p3.x - p2.x) * (p3.x - p2.x) + (p3.y - p2.y) * (p3.y - p2.y));
    d31 = sqrt((p1.x - p3.x) * (p1.x - p3.x) + (p1.y - p3.y) * (p1.y - p3.y));

    S12 = (p2.y - p1.y) / d12;
    S23 = (p3.y - p2.y) / d23;
    S31 = (p1.y - p3.y) / d31;

    C12 = (p2.x - p1.x) / d12;
    C23 = (p3.x - p2.x) / d23;
    C31 = (p1.x - p3.x) / d31;

    Q12 = log((r1 + r2 + d12) / (r1 + r2 - d12));
    Q23 = log((r2 + r3 + d23) / (r2 + r3 - d23));
    Q31 = log((r3 + r1 + d31) / (r3 + r1 - d31));

    R12 = (p.x - p1.x) * S12 - (p.y - p1.y) * C12;
    R23 = (p.x - p2.x) * S23 - (p.y - p2.y) * C23;
    R31 = (p.x - p3.x) * S31 - (p.y - p3.y) * C31;

    s12_1 = (p1.x - p.x) * C12 + (p1.y - p.y) * S12;
    s12_2 = (p2.x - p.x) * C12 + (p2.y - p.y) * S12;

    s23_2 = (p2.x - p.x) * C23 + (p2.y - p.y) * S23;
    s23_3 = (p3.x - p.x) * C23 + (p3.y - p.y) * S23;

    s31_3 = (p3.x - p.x) * C31 + (p3.y - p.y) * S31;
    s31_1 = (p1.x - p.x) * C31 + (p1.y - p.y) * S31;

    J12 = atan2( (R12 * fabs(p.z) * (r1 * s12_2 - r2 * s12_1)) , (r1 * r2 * R12 * R12 + p.z * p.z * s12_1 * s12_2) );
    J23 = atan2( (R23 * fabs(p.z) * (r2 * s23_3 - r3 * s23_2)) , (r2 * r3 * R23 * R23 + p.z * p.z * s23_2 * s23_3) );
    J31 = atan2( (R31 * fabs(p.z) * (r3 * s31_1 - r1 * s31_3)) , (r3 * r1 * R31 * R31 + p.z * p.z * s31_3 * s31_1) );

    if (p.z < 0) {
        sign = - 1.0;
    } else {
        sign = 1.0;
    }

    if ((R12 < 0) && (R23 < 0) && (R31 < 0)) {
        delta = 2 * M_PI;
    } else {
        delta = 0.0;
    }
        
    u = FACTOR * (S12 * Q12 + S23 * Q23 + S31 * Q31);
    v = FACTOR * (- C12 * Q12 - C23 * Q23 - C31 * Q31);
    w = FACTOR * sign * (delta + J12 + J23 + J31);
    

    vel.x = u * e1.x + v * e2.x + w * e3.x;
    vel.y = u * e1.y + v * e2.y + w * e3.y;
    vel.z = u * e1.z + v * e2.z + w * e3.z;

    return vel;
}

Vec3D quad_source_panel(Vec3D p, Vec3D p1, Vec3D p2, Vec3D p3, Vec3D p4, Vec3D e1, Vec3D e2, Vec3D e3)
{

    Vec3D vel;
    double u, v, w;

    double r1, r2, r3, r4;
    double d12, d23, d34, d41;
    double S12, S23, S34, S41;
    double C12, C23, C34, C41;
    double Q12, Q23, Q34, Q41;
    double R12, R23, R34, R41;
    double J12, J23, J34, J41;

    double s12_1, s12_2;
    double s23_2, s23_3;
    double s34_3, s34_4;
    double s41_4, s41_1;

    double sign, delta;

    r1 = sqrt((p.x - p1.x) * (p.x - p1.x) + (p.y - p1.y) * (p.y - p1.y) + p.z * p.z);
    r2 = sqrt((p.x - p2.x) * (p.x - p2.x) + (p.y - p2.y) * (p.y - p2.y) + p.z * p.z);
    r3 = sqrt((p.x - p3.x) * (p.x - p3.x) + (p.y - p3.y) * (p.y - p3.y) + p.z * p.z);
    r4 = sqrt((p.x - p4.x) * (p.x - p4.x) + (p.y - p4.y) * (p.y - p4.y) + p.z * p.z);

    d12 = sqrt((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y));
    d23 = sqrt((p3.x - p2.x) * (p3.x - p2.x) + (p3.y - p2.y) * (p3.y - p2.y));
    d34 = sqrt((p4.x - p3.x) * (p4.x - p3.x) + (p4.y - p3.y) * (p4.y - p3.y));
    d41 = sqrt((p1.x - p4.x) * (p1.x - p4.x) + (p1.y - p4.y) * (p1.y - p4.y));

    S12 = (p2.y - p1.y) / d12;
    S23 = (p3.y - p2.y) / d23;
    S34 = (p4.y - p3.y) / d34;
    S41 = (p1.y - p4.y) / d41;

    C12 = (p2.x - p1.x) / d12;
    C23 = (p3.x - p2.x) / d23;
    C34 = (p4.x - p3.x) / d34;
    C41 = (p1.x - p4.x) / d41;

    Q12 = log((r1 + r2 + d12) / (r1 + r2 - d12));
    Q23 = log((r2 + r3 + d23) / (r2 + r3 - d23));
    Q34 = log((r3 + r4 + d34) / (r3 + r4 - d34));
    Q41 = log((r4 + r1 + d41) / (r4 + r1 - d41));

    R12 = (p.x - p1.x) * S12 - (p.y - p1.y) * C12;
    R23 = (p.x - p2.x) * S23 - (p.y - p2.y) * C23;
    R34 = (p.x - p3.x) * S34 - (p.y - p3.y) * C34;
    R41 = (p.x - p4.x) * S41 - (p.y - p4.y) * C41;

    s12_1 = (p1.x - p.x) * C12 + (p1.y - p.y) * S12;
    s12_2 = (p2.x - p.x) * C12 + (p2.y - p.y) * S12;

    s23_2 = (p2.x - p.x) * C23 + (p2.y - p.y) * S23;
    s23_3 = (p3.x - p.x) * C23 + (p3.y - p.y) * S23;

    s34_3 = (p3.x - p.x) * C34 + (p3.y - p.y) * S34;
    s34_4 = (p4.x - p.x) * C34 + (p4.y - p.y) * S34;

    s41_4 = (p4.x - p.x) * C41 + (p4.y - p.y) * S41;
    s41_1 = (p1.x - p.x) * C41 + (p1.y - p.y) * S41;

    J12 = atan2( (R12 * fabs(p.z) * (r1 * s12_2 - r2 * s12_1)) , (r1 * r2 * R12 * R12 + p.z * p.z * s12_1 * s12_2) );
    J23 = atan2( (R23 * fabs(p.z) * (r2 * s23_3 - r3 * s23_2)) , (r2 * r3 * R23 * R23 + p.z * p.z * s23_2 * s23_3) );
    J34 = atan2( (R34 * fabs(p.z) * (r3 * s34_4 - r4 * s34_3)) , (r3 * r4 * R34 * R34 + p.z * p.z * s34_3 * s34_4) );
    J41 = atan2( (R41 * fabs(p.z) * (r4 * s41_1 - r1 * s41_4)) , (r4 * r1 * R41 * R41 + p.z * p.z * s41_4 * s41_1) );

    if (p.z < 0) {
        sign = - 1.0;
    } else {
        sign = 1.0;
    }

    if ((R12 < 0) && (R23 < 0) && (R34 < 0) && (R41 < 0)) {
        delta = 2 * M_PI;
    } else {
        delta = 0.0;
    }
        
    u = FACTOR * (S12 * Q12 + S23 * Q23 + S34 * Q34 + S41 * Q41);
    v = FACTOR * (- C12 * Q12 - C23 * Q23 - C34 * Q34 - C41 * Q41);
    w = FACTOR * sign * (delta + J12 + J23 + J34 + J41);

    vel.x = u * e1.x + v * e2.x + w * e3.x;
    vel.y = u * e1.y + v * e2.y + w * e3.y;
    vel.z = u * e1.z + v * e2.z + w * e3.z;

    return vel;
}

Vec3D line_vortex(Vec3D p, Vec3D p1, Vec3D p2)
{
    Vec3D vel;

    Vec3D r1 = {p.x - p1.x, p.y - p1.y, p.z - p1.z};
    Vec3D r2 = {p.x - p2.x, p.y - p2.y, p.z - p2.z};
    Vec3D r0 = {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};

    Vec3D r1xr2 = {r1.y * r2.z - r1.z * r2.y, r1.z * r2.x - r1.x * r2.z, r1.x * r2.y - r1.y * r2.x};

    double r1xr2_square = r1xr2.x * r1xr2.x + r1xr2.y * r1xr2.y + r1xr2.z * r1xr2.z;

    double r1_norm = sqrt(r1.x * r1.x + r1.y * r1.y + r1.z * r1.z);
    double r2_norm = sqrt(r2.x * r2.x + r2.y * r2.y + r2.z * r2.z);

    if ((r1_norm < ZERO_ERROR) || (r2_norm < ZERO_ERROR) || (r1xr2_square < ZERO_ERROR))
    {
        vel.x = 0.0;
        vel.y = 0.0;
        vel.z = 0.0;
    }
    else
    {
        double r0r1 = r0.x * r1.x + r0.y * r1.y + r0.z * r1.z;
        double r0r2 = r0.x * r2.x + r0.y * r2.y + r0.z * r2.z;
        
        double k = (1 / (4 * M_PI * r1xr2_square)) * (r0r1 / r1_norm - r0r2 / r2_norm);

        vel.x = k * r1xr2.x;
        vel.y = k * r1xr2.y;
        vel.z = k * r1xr2.z;
    }

    return vel;
}

Vec3D tri_doublet_panel(Vec3D p, Vec3D p1, Vec3D p2, Vec3D p3, Vec3D e1, Vec3D e2, Vec3D e3, double area, double max_distance)
{
    Vec3D vel, p_local;

    p_local.x = dot_func(p, e1);
    p_local.y = dot_func(p, e2);
    p_local.z = dot_func(p, e3);

    double distance = norm_func(p);
    double u, v, w;

    if (distance > max_distance) {

        double den = pow(p.x * p.x + p.y * p.y + p.z * p.z, 2.5);

        u = 0.75 * FACTOR * area * p.z * p.x / den;
        v = 0.75 * FACTOR * area * p.z * p.y / den;
        w = - FACTOR * area * (p.x * p.x + p.y * p.y - 2 * p.z * p.z) / den;
    
    } else {

        Vec3D vel1, vel2, vel3, vel4;

        vel1 = line_vortex(p, p1, p2);
        vel2 = line_vortex(p, p2, p3);
        vel3 = line_vortex(p, p3, p1);

        u = vel1.x + vel2.x + vel3.x;
        v = vel1.y + vel2.y + vel3.y;
        w = vel1.z + vel2.z + vel3.z;
    }

    vel.x = u * e1.x + v * e2.x + w * e3.x;
    vel.y = u * e1.y + v * e2.y + w * e3.y;
    vel.z = u * e1.z + v * e2.z + w * e3.z;

    return vel;
}

Vec3D quad_doublet_panel(Vec3D p, Vec3D p1, Vec3D p2, Vec3D p3, Vec3D p4, Vec3D e1, Vec3D e2, Vec3D e3, double area, double max_distance)
{
    Vec3D vel, p_local;

    p_local.x = dot_func(p, e1);
    p_local.y = dot_func(p, e2);
    p_local.z = dot_func(p, e3);

    double distance = norm_func(p);
    double u, v, w;

    if (distance > max_distance) {

        double den = pow(p.x * p.x + p.y * p.y + p.z * p.z, 2.5);

        u = 0.75 * FACTOR * area * p.z * p.x / den;
        v = 0.75 * FACTOR * area * p.z * p.y / den;
        w = - FACTOR * area * (p.x * p.x + p.y * p.y - 2 * p.z * p.z) / den;
    
    } else {

        Vec3D vel1, vel2, vel3, vel4;

        vel1 = line_vortex(p, p1, p2);
        vel2 = line_vortex(p, p2, p3);
        vel3 = line_vortex(p, p3, p4);
        vel4 = line_vortex(p, p4, p1);

        u = vel1.x + vel2.x + vel3.x + vel4.x;
        v = vel1.y + vel2.y + vel3.y + vel4.y;
        w = vel1.z + vel2.z + vel3.z + vel4.z;
    }

    vel.x = u * e1.x + v * e2.x + w * e3.x;
    vel.y = u * e1.y + v * e2.y + w * e3.y;
    vel.z = u * e1.z + v * e2.z + w * e3.z;

    return vel;
}

Vec3D vel_doublet_sheet(int n_sides, Vec3D p, Vec3D e1, Vec3D e2, Vec3D e3, Vec3D p1, Vec3D p2, Vec3D p3, Vec3D p4, double area, double max_distance)
{

    Vec3D vel, p_local;

    p_local.x = dot_func(p, e1);
    p_local.y = dot_func(p, e2);
    p_local.z = dot_func(p, e3);

    if (n_sides == 4) {
        vel = quad_doublet_panel(p_local, p1, p2, p3, p4, e1, e2, e3, area, max_distance);
    } else {
        vel = tri_doublet_panel(p_local, p1, p2, p3, e1, e2, e3, area, max_distance);
    }

    return vel;
}

Vec3D vel_source_sheet(int n_sides, Vec3D p, Vec3D e1, Vec3D e2, Vec3D e3, Vec3D p1, Vec3D p2, Vec3D p3, Vec3D p4)
{

    Vec3D vel, p_local;

    p_local.x = dot_func(p, e1);
    p_local.y = dot_func(p, e2);
    p_local.z = dot_func(p, e3);

    if (n_sides == 4) {
        vel = quad_source_panel(p_local, p1, p2, p3, p4, e1, e2, e3);
    } else {
        vel = tri_source_panel(p_local, p1, p2, p3, e1, e2, e3);
    }

    return vel;
}

/*-------------------------------------------
    Calculate a_ij coefs
--------------------------------------------*/
void get_a_ij_coefs(int nf,
                    int n_sides[],
                    double p_avg[],
                    double p_ctrl[],
                    double e1[], double e2[], double e3[],
                    double p1[], double p2[], double p3[], double p4[],
                    double a_ij_x[], double a_ij_y[], double a_ij_z[])
{
    int i, j;

    int n_sides_j;
    Vec3D p_ij;
    Vec3D e1_j, e2_j, e3_j;
    Vec3D p1_j, p2_j, p3_j, p4_j;

    Vec3D vel_source;

    for (i = 0; i < nf; i++)
    {

        for (j = 0; j < nf; j++)
        {

            n_sides_j = n_sides[j];

            p_ij.x = p_ctrl[3 * i] - p_avg[3 * j]; p_ij.y = p_ctrl[3 * i + 1] - p_avg[3 * j + 1]; p_ij.z = p_ctrl[3 * i + 2] - p_avg[3 * j + 2];

            e1_j.x = e1[3 * j]; e1_j.y = e1[3 * j + 1]; e1_j.z = e1[3 * j + 2];
            e2_j.x = e2[3 * j]; e2_j.y = e2[3 * j + 1]; e2_j.z = e2[3 * j + 2];
            e3_j.x = e3[3 * j]; e3_j.y = e3[3 * j + 1]; e3_j.z = e3[3 * j + 2];

            p1_j.x = p1[2 * j]; p1_j.y = p1[2 * j + 1]; p1_j.z = 0.0;
            p2_j.x = p2[2 * j]; p2_j.y = p2[2 * j + 1]; p2_j.z = 0.0;
            p3_j.x = p3[2 * j]; p3_j.y = p3[2 * j + 1]; p3_j.z = 0.0;
            p4_j.x = p4[2 * j]; p4_j.y = p4[2 * j + 1]; p4_j.z = 0.0;

            vel_source = vel_source_sheet(n_sides_j, p_ij, e1_j, e2_j, e3_j, p1_j, p2_j, p3_j, p4_j);

            a_ij_x[i * nf + j] = vel_source.x;
            a_ij_y[i * nf + j] = vel_source.y;
            a_ij_z[i * nf + j] = vel_source.z;
        
        }
    }

}

/*-------------------------------------------
    Calculate b_kj coefs
--------------------------------------------*/
void get_b_kj_coefs(int nf,
                    int nte,
                    int faces[],
                    double vertices[],
                    double scale[],
                    double p_ctrl[],
                    double b_kj_x[], double b_kj_y[], double b_kj_z[])
{
    int i, j;

    int id1, id2, id3, id4;

    Vec3D vel1, vel2, vel3, vel4;
    Vec3D p, p1, p2, p3, p4;

    for (i = 0; i < nf; i++)
    {

        p.x = p_ctrl[3 * i];
        p.y = p_ctrl[3 * i + 1];
        p.z = p_ctrl[3 * i + 2];

        for (j = 0; j < nte; j++)
        {

            id1 = faces[j * 5 + 1]; id2 = faces[j * 5 + 2]; id3 = faces[j * 5 + 3];
            p1.x = vertices[3 * id1]; p1.y = vertices[3 * id1 + 1]; p1.z = vertices[3 * id1 + 2];
            p2.x = vertices[3 * id2]; p2.y = vertices[3 * id2 + 1]; p2.z = vertices[3 * id2 + 2];
            p3.x = vertices[3 * id3]; p3.y = vertices[3 * id3 + 1]; p3.z = vertices[3 * id3 + 2];

            if (faces[j * 5] == 3) {
                vel1 = line_vortex(p, p1, p2);
                vel2 = line_vortex(p, p2, p3);
                vel3 = line_vortex(p, p3, p1);
                vel4.x = 0.0; vel4.y = 0.0; vel4.z = 0.0;
            } else {
                id4 = faces[j * 5 + 4];
                p4.x = vertices[3 * id4]; p4.y = vertices[3 * id4 + 1]; p4.z = vertices[3 * id4 + 2];
                vel1 = line_vortex(p, p1, p2);
                vel2 = line_vortex(p, p2, p3);
                vel3 = line_vortex(p, p3, p4);
                vel4 = line_vortex(p, p4, p1);
            }

            b_kj_x[i * nte + j] = scale[j] * (vel1.x + vel2.x + vel3.x + vel4.x);
            b_kj_y[i * nte + j] = scale[j] * (vel1.y + vel2.y + vel3.y + vel4.y);
            b_kj_z[i * nte + j] = scale[j] * (vel1.z + vel2.z + vel3.z + vel4.z);
        
        }
    }

}

/*-------------------------------------------
    Calculate c_kj and d_j coefs
--------------------------------------------*/
void get_c_kj_vortex_line_coefs(int nf,
                                int nte,
                                int te[],
                                double vertices[],
                                double scale[],
                                double p_ctrl[],
                                double c_kj_x[], double c_kj_y[], double c_kj_z[])
{
    int i, j;

    int id1, id2;

    Vec3D vel;
    Vec3D p, p1, p2;

    for (i = 0; i < nf; i++)
    {

        p.x = p_ctrl[3 * i];
        p.y = p_ctrl[3 * i + 1];
        p.z = p_ctrl[3 * i + 2];

        for (j = 0; j < nte; j++)
        {
            
            id1 = te[2 * j]; id2 = te[2 * j + 1];

            p1.x = vertices[3 * id1]; p1.y = vertices[3 * id1 + 1]; p1.z = vertices[3 * id1 + 2];
            p2.x = vertices[3 * id2]; p2.y = vertices[3 * id2 + 1]; p2.z = vertices[3 * id2 + 2];

            vel = line_vortex(p, p1, p2);

            c_kj_x[i * nte + j] = scale[j] * vel.x;
            c_kj_y[i * nte + j] = scale[j] * vel.y;
            c_kj_z[i * nte + j] = scale[j] * vel.z;
        
        }
    }

}

void get_c_kj_dj_coefs(int nf,
                       int nte,
                       int nw,
                       int section,
                       int faces[],
                       double vertices[],
                       double areas[],
                       double circulation[],
                       double scale[],
                       double p_ctrl[],
                       double c_kj_x[], double c_kj_y[], double c_kj_z[],
                       double d_j_x[], double d_j_y[], double d_j_z[])
{
    int i, j, k;

    int id1, id2, id3, id4;

    double area;

    Vec3D vel1, vel2, vel3, vel4;
    Vec3D p, p1, p2, p3, p4;
    Vec3D v1, v2, v3;

    for (i = 0; i < nf; i++)
    {

        p.x = p_ctrl[3 * i];
        p.y = p_ctrl[3 * i + 1];
        p.z = p_ctrl[3 * i + 2];

        d_j_x[i] = 0.0;
        d_j_y[i] = 0.0;
        d_j_z[i] = 0.0;

        for (j = 0; j < nte; j++)
        {

            for (k = 0; k < section; k++)
            {

                id1 = faces[j * (nw * 2) + k * 2];
                id2 = faces[j * (nw * 2) + k * 2 + 1];
                id3 = faces[j * (nw * 2) + (k + 1) * 2 + 1];
                id4 = faces[j * (nw * 2) + (k + 1) * 2];

                p1.x = vertices[3 * id1]; p1.y = vertices[3 * id1 + 1]; p1.z = vertices[3 * id1 + 2];
                p2.x = vertices[3 * id2]; p2.y = vertices[3 * id2 + 1]; p2.z = vertices[3 * id2 + 2];
                p3.x = vertices[3 * id3]; p3.y = vertices[3 * id3 + 1]; p3.z = vertices[3 * id3 + 2];
                p4.x = vertices[3 * id4]; p4.y = vertices[3 * id4 + 1]; p4.z = vertices[3 * id4 + 2];
            
                vel1 = line_vortex(p, p1, p2);
                vel2 = line_vortex(p, p2, p3);
                vel3 = line_vortex(p, p3, p4);
                vel4 = line_vortex(p, p4, p1);

                if (k == 0) {

                    c_kj_x[i * nte + j] = scale[j] * (vel1.x + vel2.x + vel3.x + vel4.x);
                    c_kj_y[i * nte + j] = scale[j] * (vel1.y + vel2.y + vel3.y + vel4.y);
                    c_kj_z[i * nte + j] = scale[j] * (vel1.z + vel2.z + vel3.z + vel4.z);

                } else {

                    v1.x = p2.x - p1.x; v1.y = p2.y - p1.y; v1.z = p2.z - p1.z;
                    v2.x = p3.x - p1.x; v2.y = p3.y - p1.y; v2.z = p3.z - p1.z;
                    v3.x = p4.x - p1.x; v3.y = p4.y - p1.y; v3.z = p4.z - p1.z;

                    area = 0.5 * (norm_func(cross_func(v1, v2)) + norm_func(cross_func(v2, v3)));

                    d_j_x[i] = d_j_x[i] + (circulation[j * nw + k] * areas[j * nw + k] / area) * scale[j] * (vel1.x + vel2.x + vel3.x + vel4.x);
                    d_j_y[i] = d_j_y[i] + (circulation[j * nw + k] * areas[j * nw + k] / area) * scale[j] * (vel1.y + vel2.y + vel3.y + vel4.y);
                    d_j_z[i] = d_j_z[i] + (circulation[j * nw + k] * areas[j * nw + k] / area) * scale[j] * (vel1.z + vel2.z + vel3.z + vel4.z);
                
                }

            }
        
        }
    }

}

/*-------------------------------------------
    Point velocity
--------------------------------------------*/