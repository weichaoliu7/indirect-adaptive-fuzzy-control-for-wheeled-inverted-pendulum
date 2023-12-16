#include <webots/motor.h>
#include <webots/robot.h>
#include <webots/gps.h>
#include <webots/inertial_unit.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "sine.h"
#include "cosine.h"

// global variables declaration
#define TIME_STEP 10
#define PI 3.14159
#define ARRAY_SIZE 3000                    // sampling times
#define n_membership 5                     // number of membership functions
#define dimension 25                       // dimension of fuzzy controller parameter
#define l 1.0                              // distance from the body center of gravity to the wheel axis
static double t0 = 0.0;                    // start time
static double t1 = 30.0;                   // end time
static double R = 2.0;                     // radius of circular desired trajectory
static double Ts = 0.010;                  // sampling period
static double rotational_velocity_d = 1.0; // desird rotational velocity of the vehicle
static double rotational_angle_d = 0.010;  // desird rotational angle of the vehicle

/* reference: [1]Yue M, An C, Du Y, et al. Indirect adaptive fuzzy control for a nonholonomic/underactuated
wheeled inverted pendulum vehicle based on a data-driven trajectory planner[J]. Fuzzy Sets and Systems, 2016, 290: 158-177. */

struct _archive{
    double x_coordinate_archive[ARRAY_SIZE];
    double y_coordinate_archive[ARRAY_SIZE];
    double x_coordinate_right_archive[ARRAY_SIZE];
    double y_coordinate_right_archive[ARRAY_SIZE];
    double x_coordinate_left_archive[ARRAY_SIZE];
    double y_coordinate_left_archive[ARRAY_SIZE];
    double rotational_angle_archive[ARRAY_SIZE];
    double rotational_velocity_archive[ARRAY_SIZE];
    double error_x_coordinate_archive[ARRAY_SIZE];
    double error_y_coordinate_archive[ARRAY_SIZE];
    double error_rotational_angle_archive[ARRAY_SIZE];
    double velocity_plan_archive[ARRAY_SIZE];
    double rotational_velocity_plan_archive[ARRAY_SIZE];
    double velocity_archive[ARRAY_SIZE];
    double tilt_angle_archive[ARRAY_SIZE];
    double error_velocity_archive[ARRAY_SIZE];
    double error_rotational_velocity_archive[ARRAY_SIZE];
    double torque_velocity_archive[ARRAY_SIZE];
    double torque_rotational_velocity_archive[ARRAY_SIZE];
} archive;

struct _system_state{
    double x_coordinate;        // x-coordinate of the vehicle
    double y_coordinate;        // y-coordinate of the vehicle
    double x_velocity;          // vehicle velocity in x-direction
    double y_velocity;          // vehicle velocity in y-direction
    double rotational_angle;    // rotational angle of the vehicle
    double rotational_velocity; // rotational velocity of the vehicle
    double velocity;            // longitudinal velocity of the vehicle
    double tilt_angle;          // tilt angle of the vehicle body
    double tilt_velocity;       // tilt velocity of the vehicle body
} system_state;

// gaussian function
double gaussian(double x, double mean, double sigma){
    return exp(-pow((x - mean) / sigma, 2));
}

typedef struct{
    double negative_large;
    double negative_small;
    double zero;
    double positive_small;
    double positive_large;
} membership;

// membership function of input vectors in fuzzy sets
typedef struct{
    membership fx; // membership function of unknown smooth function f(x)
    membership gx; // membership function of unknown smooth function g(x)
} membershipfunction;

// pointer to access membership grade
double membershipgrade(membership *member, int index){
    switch (index){
    case 0:
        return member->negative_large;
    case 1:
        return member->negative_small;
    case 2:
        return member->zero;
    case 3:
        return member->positive_small;
    case 4:
        return member->positive_large;
    default:
        return 0;
    }
}

// symbolic function
double sign(double a){
    if (a > 0){
        return 1.0;
    }
    else if (a < 0){
        return -1.0;
    }
    else{
        return 0.0;
    }
}

struct _planner{
    double planner_u[6];
    double planner_out[2];
    double x_coordinate_d;                  // desired x-cordinate of midpoint of the vehicle
    double y_coordinate_d;                  // desired y-coordinate of midpoint of the vehicle
    double x_coordinate_right;              // x-coordinate of right wheel
    double y_coordinate_right;              // y-coordinate of right wheel
    double x_coordinate_left;               // x-coordinate of left wheel
    double y_coordinate_left;               // y-coordinate of left wheel
    double velocity_right;                  // velocity of right wheel
    double velocity_left;                   // velocity  of left wheel
    double error_x_coordinate;              // error of x-coordinate of midpoint of the vehicle
    double error_y_coordinate;              // error of y-coordinate of midpoint of the vehicle
    double error_rotational_angle;          // error of rotational angle of the vehicle
    double tracking_error_x;                // tracking error of x-coordinate
    double tracking_error_y;                // tracking error of y-coordinate
    double tracking_error_rotational_angle; // tracking error of rotational angle of the vehicle
    double velocity_desired;                // desired longitudinal velocity of the vehicle
    double velocity_plan;                   // planned longitudinal velocity of the vehicle
    double rotational_velocity_plan;        // planned rotational velocity of the vehicle
    double lambda1;                         // positive constant
    double lambda2;                         // positive constant
    double lambda3;                         // positive constant
} planner;

void PLANNER_init(){
    wb_robot_init();
    WbDeviceTag gps1 = wb_robot_get_device("gps1");
    wb_gps_enable(gps1, TIME_STEP);
    WbDeviceTag gps3 = wb_robot_get_device("gps3");
    wb_gps_enable(gps3, TIME_STEP);
    WbDeviceTag gps4 = wb_robot_get_device("gps4");
    wb_gps_enable(gps4, TIME_STEP);
    planner.lambda1 = 0.4;
    planner.lambda2 = 0.4;
    planner.lambda3 = 0.4;
    system_state.x_coordinate = 0.0;
    system_state.y_coordinate = 0.0;
    system_state.rotational_angle = 0.0;
    planner.planner_u[0] = 0.0;                           // desired x-coordinate of midpoint of the vehicle
    planner.planner_u[1] = 0.0;                           // desired y-coordinate of midpoint of the vehicle
    planner.planner_u[2] = rotational_angle_d;            // desired rotational angle of the vehicle
    planner.planner_u[3] = system_state.x_coordinate;     // actual x-coordinate of midpoint of the vehicle
    planner.planner_u[4] = system_state.y_coordinate;     // actual y-coordinate of midpoint of the vehicle
    planner.planner_u[5] = system_state.rotational_angle; // actual rotational angle of the vehicle
}

double PLANNER_realize(int i){

    // obtain x-coordinate and x-coordinate of vehicle using GPS sensor
    WbDeviceTag gps1 = wb_robot_get_device("gps1");
    wb_gps_enable(gps1, TIME_STEP);
    const double *coordinate1 = wb_gps_get_values(gps1);
    system_state.x_coordinate = coordinate1[0]; // x-coordinate of midpoint of two wheels of wheeled inverted pendulum
    system_state.y_coordinate = coordinate1[1]; // x-coordinate of midpoint of two wheels of wheeled inverted pendulum

    printf("x_coordinate: %f\n", system_state.x_coordinate);
    printf("y_coordinate: %f\n", system_state.y_coordinate);
    // printf("z_coordinate: %f\n", coordinate1[2]);
    archive.x_coordinate_archive[i] = system_state.x_coordinate;
    archive.y_coordinate_archive[i] = system_state.y_coordinate;

    // obtain x-coordinate and x-coordinate of right and left wheel using GPS sensor
    WbDeviceTag gps3 = wb_robot_get_device("gps3");
    wb_gps_enable(gps3, TIME_STEP);
    const double *coordinate3 = wb_gps_get_values(gps3);
    WbDeviceTag gps4 = wb_robot_get_device("gps4");
    wb_gps_enable(gps4, TIME_STEP);
    const double *coordinate4 = wb_gps_get_values(gps4);

    planner.x_coordinate_right = coordinate3[0]; // x-coordinate of right wheel
    planner.y_coordinate_right = coordinate3[1]; // y-coordinate of right wheel
    planner.x_coordinate_left = coordinate4[0];  // x-coordinate of left wheel
    planner.y_coordinate_left = coordinate4[1];  // y-coordinate of left wheel

    archive.x_coordinate_right_archive[i] = planner.x_coordinate_right;
    archive.y_coordinate_right_archive[i] = planner.y_coordinate_right;
    archive.x_coordinate_left_archive[i] = planner.x_coordinate_left;
    archive.y_coordinate_left_archive[i] = planner.y_coordinate_left;

    if (i == 0){
        planner.velocity_right = 0.2; // velocity of right wheel
        planner.velocity_left = 0.4;  // velocity of left wheel
    }
    else if (i > 0){
        planner.velocity_right = sqrt(pow(planner.x_coordinate_right - archive.x_coordinate_right_archive[i - 1], 2) + pow(planner.y_coordinate_right - archive.y_coordinate_right_archive[i - 1], 2)) / Ts;
        planner.velocity_left = sqrt(pow(planner.x_coordinate_left - archive.x_coordinate_left_archive[i - 1], 2) + pow(planner.y_coordinate_left - archive.y_coordinate_left_archive[i - 1], 2)) / Ts;    
    }

    system_state.rotational_angle =  (planner.velocity_left - planner.velocity_right) * Ts / l; // rotational angle of the vehicle, positive direction to the right 
    printf("rotational_angle: %f\n", system_state.rotational_angle);
    archive.rotational_angle_archive[i] = system_state.rotational_angle;

    if (system_state.x_coordinate < R){
        planner.x_coordinate_d = R - sqrt(pow(R, 2) - pow(system_state.y_coordinate, 2));
    }
    else if (system_state.x_coordinate >= R){
        planner.x_coordinate_d = R + sqrt(pow(R, 2) - pow(system_state.y_coordinate, 2));
    }

    if (system_state.y_coordinate < 0.0){
        planner.y_coordinate_d = - sqrt(pow(R, 2) - pow(system_state.x_coordinate - R, 2));
    }
    else if (system_state.y_coordinate >= 0.0){
        planner.y_coordinate_d = sqrt(pow(R, 2) - pow(system_state.x_coordinate - R, 2));
    }

    // planner for the WIP vehicle
    planner.planner_u[0] = planner.x_coordinate_d;        // desired x-coordinate of midpoint of the vehicle
    planner.planner_u[1] = planner.y_coordinate_d;        // desired y-coordinate of midpoint of the vehicle
    planner.planner_u[2] = rotational_angle_d;            // desired rotational angle of the vehicle
    planner.planner_u[3] = system_state.x_coordinate;     // actual x-coordinate of midpoint of the vehicle
    planner.planner_u[4] = system_state.y_coordinate;     // actual y-coordinate of midpoint of the vehicle
    planner.planner_u[5] = system_state.rotational_angle; // actual rotational angle of the vehicle

    planner.error_x_coordinate = planner.x_coordinate_d - system_state.x_coordinate;     // error of x-coordinate of midpoint of the vehicle
    planner.error_y_coordinate = planner.y_coordinate_d - system_state.y_coordinate;     // error of y-coordinate of midpoint of the vehicle
    planner.error_rotational_angle = rotational_angle_d - system_state.rotational_angle; // error of rotational angle of the vehicle

    if (i == 0){
        planner.error_x_coordinate = -0.001;
        planner.error_y_coordinate = -0.001;
    }

    printf("error_x_coordinate: %f\n", planner.error_x_coordinate);
    printf("error_y_coordinate: %f\n", planner.error_y_coordinate);
    printf("error_rotational_angle: %f\n", planner.error_rotational_angle);
    archive.error_x_coordinate_archive[i] = planner.error_x_coordinate;
    archive.error_y_coordinate_archive[i] = planner.error_y_coordinate;
    archive.error_rotational_angle_archive[i] = planner.error_rotational_angle;

    // tracking error of x-coordinate
    planner.tracking_error_x = cos(system_state.rotational_angle) * planner.error_x_coordinate + sin(system_state.rotational_angle) * planner.error_y_coordinate;
    // tracking error of y-coordinate
    planner.tracking_error_y = -sin(system_state.rotational_angle) * planner.error_x_coordinate + cos(system_state.rotational_angle) * planner.error_y_coordinate;
    // tracking error of rotational angle of the vehicle
    planner.tracking_error_rotational_angle = planner.error_rotational_angle;

    planner.velocity_desired = rotational_velocity_d * R; // desired longitudinal velocity of the vehicle

    // planned longitudinal velocity of the vehicle
    planner.velocity_plan = planner.velocity_desired * cos(planner.tracking_error_rotational_angle) + planner.lambda3 * planner.tracking_error_x;
    // planned rotational velocity of the vehicle
    planner.rotational_velocity_plan = rotational_velocity_d + planner.lambda1 * planner.velocity_desired * planner.tracking_error_y + planner.lambda2 * sin(planner.tracking_error_rotational_angle);

    printf("velocity_plan: %f\n", planner.velocity_plan);
    printf("rotational_velocity_plan: %f\n", planner.rotational_velocity_plan);
    archive.velocity_plan_archive[i] = planner.velocity_plan;
    archive.rotational_velocity_plan_archive[i] = planner.rotational_velocity_plan;
    planner.planner_out[0] = planner.velocity_plan;
    planner.planner_out[1] = planner.rotational_velocity_plan;
}

struct _controller{
    double controller_u[5];
    double controller_out[6];
    double z_coordinate;                                       // z-coordinate of the vehicle
    double x1_coordinate;                                      // x-coordinate of top of inverted pendulum
    double y1_coordinate;                                      // y-coordinate of top of inverted pendulum
    double z1_coordinate;                                      // z-coordinate of top of inverted pendulum
    double sliding;                                            // sliding mode manifold of tilt angle subsystem
    membershipfunction mu_tilt_angle;                          // membership function of tilt angle
    double xi_denominator_tilt_angle;                          // denominator of fuzzy basis function xi of tilt angle membership function
    double xi_numerator_tilt_angle[dimension];                 // numerator of fuzzy basis function xi of tilt angle membership function
    double xi_tilt_angle[dimension];                           // fuzzy basis function xi of tilt angle membership function
    double eta_tilt_angle[dimension];                          // fuzzy basis function eta of tilt angle membership function
    double theta_tilt_angle_fx[dimension];                     // fuzzy controller parameter of f(x) of tilt angle
    double theta_tilt_angle_gx[dimension];                     // fuzzy controller parameter of g(x) of tilt angle
    double theta_derivative_tilt_angle_fx[dimension];          // derivative of fuzzy controller parameter of f(x) of tilt angle
    double theta_derivative_tilt_angle_gx[dimension];          // derivative of fuzzy controller parameter of g(x) of tilt angle
    double control_tilt_angle;                                 // control amount of tilt angle
    membershipfunction mu_velocity;                            // membership function of longitudinal velocity of the vehicle
    double xi_denominator_velocity;                            // denominator of fuzzy basis function xi of longitudinal velocity
    double xi_numerator_velocity[dimension];                   // numerator of fuzzy basis function xi of longitudinal velocity
    double xi_velocity[dimension];                             // fuzzy basis function xi of longitudinal velocity
    double eta_velocity[dimension];                            // fuzzy basis function eta of longitudinal velocity
    double theta_velocity_fx[dimension];                       // fuzzy controller parameter of f(x) of longitudinal velocity
    double theta_velocity_gx[dimension];                       // fuzzy controller parameter of g(x) of longitudinal velocity
    double theta_derivative_velocity_fx[dimension];            // derivative of fuzzy controller parameter of f(x) of longitudinal velocity
    double theta_derivative_velocity_gx[dimension];            // derivative of fuzzy controller parameter of g(x) of longitudinal velocity
    double control_velocity;                                   // control amount of longitudinal velocity
    membershipfunction mu_rotational_velocity;                 // membership function of rotational velocity of the vehicle
    double xi_denominator_rotational_velocity;                 // denominator of fuzzy basis function xi of rotational velocity
    double xi_numerator_rotational_velocity[dimension];        // numerator of fuzzy basis function xi of rotational velocity
    double xi_rotational_velocity[dimension];                  // fuzzy basis function xi of rotational velocity
    double eta_rotational_velocity[dimension];                 // fuzzy basis function eta of rotational velocity
    double theta_rotational_velocity_fx[dimension];            // fuzzy controller parameter of f(x) of rotational velocity
    double theta_rotational_velocity_gx[dimension];            // fuzzy controller parameter of g(x) of rotational velocity
    double theta_derivative_rotational_velocity_fx[dimension]; // derivative of fuzzy controller parameter of f(x) of rotational velocity
    double theta_derivative_rotational_velocity_gx[dimension]; // derivative of fuzzy controller parameter of g(x) of rotational velocity
    double torque_rotational_velocity;                         // control input torque of rotational velocity
    double error_velocity;                                     // error of longitudinal velocity of the vehicle
    double error_rotational_velocity;                          // error of rotational velocity of the vehicle
    double torque_velocity;                                    // control input torque of longitudinal velocity subsystem
    double torque_right;                                       // control input torque of right wheel of vehicle
    double torque_left;                                        // control input torque of left wheel of vehicle
    double c;                                                  // parameter of sliding mode manifold
    double k1;                                                 // control gain
    double k2;                                                 // control gain
    double k3;                                                 // control gain
    double k4;                                                 // control gain
    double gamma1;                                             // controller parameter update gain
    double gamma2;                                             // controller parameter update gain
    double gamma3;                                             // controller parameter update gain
    double gamma4;                                             // controller parameter update gain
    double gamma5;                                             // controller parameter update gain
    double gamma6;                                             // controller parameter update gain
} controller;

void CONTROLLER_init(){
    WbDeviceTag gps1 = wb_robot_get_device("gps1");
    wb_gps_enable(gps1, TIME_STEP);
    WbDeviceTag gps2 = wb_robot_get_device("gps2");
    wb_gps_enable(gps2, TIME_STEP);
    WbDeviceTag gps3 = wb_robot_get_device("gps3");
    wb_gps_enable(gps3, TIME_STEP);
    WbDeviceTag gps4 = wb_robot_get_device("gps4");
    wb_gps_enable(gps4, TIME_STEP);
    WbDeviceTag motor[2];
    char motor_tag[2][8] = {"motor1", "motor2"};
    for (int j = 0; j < 2; j++){
        motor[j] = wb_robot_get_device(motor_tag[j]);
    }

    system_state.rotational_velocity = 0.0;
    system_state.velocity = 0.0;
    system_state.tilt_angle = 0.0;
    controller.controller_u[0] = planner.rotational_velocity_plan;
    controller.controller_u[1] = planner.velocity_plan;
    controller.controller_u[2] = system_state.rotational_velocity;
    controller.controller_u[3] = system_state.velocity;
    controller.controller_u[4] = system_state.tilt_angle;
    controller.c = 2;       // parameter of sliding mode manifold
    controller.k1 = 4;      // control gain
    controller.k2 = 2;      // control gain
    controller.k3 = 2;      // control gain
    controller.k4 = 1.0;    // control gain
    controller.gamma1 = 5; // controller parameter update gain
    controller.gamma2 = 1;  // controller parameter update gain
    controller.gamma3 = 5; // controller parameter update gain
    controller.gamma4 = 1;  // controller parameter update gain
    controller.gamma5 = 100; // controller parameter update gain
    controller.gamma6 = 1;  // controller parameter update gain

    for (int j = 0; j < dimension; j++){
        controller.theta_tilt_angle_fx[j] = 0.1;          // fuzzy controller parameter of f(x) of tilt angle
        controller.theta_tilt_angle_gx[j] = 0.1;          // fuzzy controller parameter of g(x) of tilt angle
        controller.theta_velocity_fx[j] = 0.1;            // fuzzy controller parameter of f(x) of longitudinal velocity
        controller.theta_velocity_gx[j] = 0.1;            // fuzzy controller parameter of g(x) of longitudinal velocity
        controller.theta_rotational_velocity_fx[j] = 0.1; // fuzzy controller parameter of f(x) of rotational velocity
        controller.theta_rotational_velocity_gx[j] = 0.1; // fuzzy controller parameter of g(x) of rotational velocity
    }
}

double CONTROLLER_realize(int i){

    // obtain the horizontal and vertical coordinates of vehicle using GPS sensor
    WbDeviceTag gps1 = wb_robot_get_device("gps1");
    wb_gps_enable(gps1, TIME_STEP);
    const double *coordinate1 = wb_gps_get_values(gps1);
    system_state.x_coordinate = coordinate1[0]; // x-coordinate of midpoint of two wheels of wheeled inverted pendulum
    system_state.y_coordinate = coordinate1[1]; // y-coordinate of midpoint of two wheels of wheeled inverted pendulum
    controller.z_coordinate = coordinate1[2];   // y-coordinate of midpoint of two wheels of wheeled inverted pendulum

    // printf("x_coordinate: %f\n", system_state.x_coordinate);
    // printf("y_coordinate: %f\n", system_state.y_coordinate);
    // printf("z_coordinate: %f\n", controller.z_coordinate);

    if (i == 0){
        system_state.x_velocity = 0.1; // vehicle velocity in x-direction
        system_state.y_velocity = 0.004; // vehicle velocity in y-direction
    }
    else if (i > 0){
        system_state.x_velocity = (system_state.x_coordinate - archive.x_coordinate_archive[i - 1]) / Ts; // vehicle velocity in x-direction
        system_state.y_velocity = (system_state.y_coordinate - archive.y_coordinate_archive[i - 1]) / Ts; // vehicle velocity in y-direction
    }

    // obtain the horizontal and vertical coordinates of right and left wheel using GPS sensor
    WbDeviceTag gps3 = wb_robot_get_device("gps3");
    wb_gps_enable(gps3, TIME_STEP);
    const double *coordinate3 = wb_gps_get_values(gps3);
    WbDeviceTag gps4 = wb_robot_get_device("gps4");
    wb_gps_enable(gps4, TIME_STEP);
    const double *coordinate4 = wb_gps_get_values(gps4);

    planner.x_coordinate_right = coordinate3[0]; // x-coordinate of right wheel
    planner.y_coordinate_right = coordinate3[1]; // y-coordinate of right wheel
    planner.x_coordinate_left = coordinate4[0]; // x-coordinate of left wheel
    planner.y_coordinate_left = coordinate4[1]; // y-coordinate of left wheel

    archive.x_coordinate_right_archive[i] = planner.x_coordinate_right;
    archive.y_coordinate_right_archive[i] = planner.y_coordinate_right;
    archive.x_coordinate_left_archive[i] = planner.x_coordinate_left;
    archive.y_coordinate_left_archive[i] = planner.y_coordinate_left;

    if (i == 0){
        planner.velocity_right = 0.2; // velocity of right wheel
        planner.velocity_left = 0.4;  // velocity of left wheel
    }
    else if (i > 0){
        planner.velocity_right = sqrt(pow(planner.x_coordinate_right - archive.x_coordinate_right_archive[i - 1], 2) + pow(planner.y_coordinate_right - archive.y_coordinate_right_archive[i - 1], 2)) / Ts;
        planner.velocity_left = sqrt(pow(planner.x_coordinate_left - archive.x_coordinate_left_archive[i - 1], 2) + pow(planner.y_coordinate_left - archive.y_coordinate_left_archive[i - 1], 2)) / Ts;    
    }

    system_state.rotational_angle =  (planner.velocity_left - planner.velocity_right) * Ts / l; // rotational angle of the vehicle, positive direction to the right 
    // printf("rotational_angle: %f\n", system_state.rotational_angle);

    if (i == 0){
        system_state.rotational_velocity = 0.2; // rotational velocity of the vehicle
    }
    else if (i > 0){
        system_state.rotational_velocity = system_state.rotational_angle / Ts; // rotational velocity of the vehicle
    }

    printf("rotational_velocity: %f\n", system_state.rotational_velocity);
    archive.rotational_velocity_archive[i] = system_state.rotational_velocity;

    system_state.velocity = sqrt(pow(system_state.x_velocity, 2) + pow(system_state.y_velocity, 2)); // longitudinal velocity of the vehicle
    printf("velocity: %f\n", system_state.velocity);
    archive.velocity_archive[i] = system_state.velocity;

    // obtain the horizontal and vertical coordinates of top of inverted pendulum
    WbDeviceTag gps2 = wb_robot_get_device("gps2");
    wb_gps_enable(gps2, TIME_STEP);
    const double *coordinate2 = wb_gps_get_values(gps2);
    controller.x1_coordinate = coordinate2[0]; // x-coordinate of midpoint of top of inverted pendulum
    controller.y1_coordinate = coordinate2[1]; // y-coordinate of midpoint of top of inverted pendulum
    controller.z1_coordinate = coordinate2[2]; // y-coordinate of midpoint of top of inverted pendulum

    // tilt angle of the vehicle body
    system_state.tilt_angle = atan(sqrt(pow(controller.x1_coordinate - system_state.x_coordinate, 2) + pow(controller.y1_coordinate - system_state.y_coordinate, 2)) / (controller.z1_coordinate - controller.z_coordinate));
    printf("tilt_angle: %f\n", system_state.tilt_angle);
    archive.tilt_angle_archive[i] = system_state.tilt_angle;

    if (i == 0){
        system_state.tilt_velocity = 0.001; // tilt velocity of the vehicle body
    }
    else if (i > 0){
        system_state.tilt_velocity = (system_state.tilt_angle - archive.tilt_angle_archive[i - 1]) / Ts; // tilt velocity of the vehicle body
    }
    printf("tilt_velocity: %f\n", system_state.tilt_velocity);

    controller.controller_u[0] = planner.rotational_velocity_plan;
    controller.controller_u[1] = planner.velocity_plan;
    controller.controller_u[2] = system_state.rotational_velocity;
    controller.controller_u[3] = system_state.velocity;
    controller.controller_u[4] = system_state.tilt_angle;

    // tilt angle subsystem
    controller.sliding = system_state.tilt_velocity + controller.c * system_state.tilt_angle; // sliding mode manifold
    printf("sliding_term: %f\n", controller.sliding);

    controller.mu_tilt_angle.fx.negative_large = gaussian(system_state.tilt_angle, -PI / 6.0, PI / 24.0);  // membership function of f(x) of tilt angle is negative large
    controller.mu_tilt_angle.fx.negative_small = gaussian(system_state.tilt_angle, -PI / 12.0, PI / 24.0); // membership function of f(x) of tilt angle is negative small
    controller.mu_tilt_angle.fx.zero = gaussian(system_state.tilt_angle, 0.0, PI / 24.0);                  // membership function of f(x) of tilt angle is zero
    controller.mu_tilt_angle.fx.positive_small = gaussian(system_state.tilt_angle, PI / 12.0, PI / 24.0);  // membership function of f(x) of tilt angle is positive small
    controller.mu_tilt_angle.fx.positive_large = gaussian(system_state.tilt_angle, PI / 6.0, PI / 24.0);   // membership function of f(x) of tilt angle is positive large

    controller.mu_tilt_angle.gx.negative_large = gaussian(system_state.tilt_angle, -PI / 6.0, PI / 24.0);  // membership function of g(x) of tilt angle is negative large
    controller.mu_tilt_angle.gx.negative_small = gaussian(system_state.tilt_angle, -PI / 12.0, PI / 24.0); // membership function of g(x) of tilt angle is negative small
    controller.mu_tilt_angle.gx.zero = gaussian(system_state.tilt_angle, 0.0, PI / 24.0);                  // membership function of g(x) of tilt angle is zero
    controller.mu_tilt_angle.gx.positive_small = gaussian(system_state.tilt_angle, PI / 12.0, PI / 24.0);  // membership function of g(x) of tilt angle is positive small
    controller.mu_tilt_angle.gx.positive_large = gaussian(system_state.tilt_angle, PI / 6.0, PI / 24.0);   // membership function of g(x) of tilt angle is positive large

    double* membership1_value[] = {
        &controller.mu_tilt_angle.fx.negative_large,
        &controller.mu_tilt_angle.fx.negative_small,
        &controller.mu_tilt_angle.fx.zero,
        &controller.mu_tilt_angle.fx.positive_small,
        &controller.mu_tilt_angle.fx.positive_large,
        &controller.mu_tilt_angle.gx.negative_large,
        &controller.mu_tilt_angle.gx.negative_small,
        &controller.mu_tilt_angle.gx.zero,
        &controller.mu_tilt_angle.gx.positive_small,
        &controller.mu_tilt_angle.gx.positive_large
    };

    int num_membership1_value = sizeof(membership1_value) / sizeof(double*);

    for (int j = 0; j < num_membership1_value; j++) {
        if (*membership1_value[j] < 0.0001) {
            *membership1_value[j] = 0.0001;
        }
    }

    // printf("mu_tilt_angle_NB = %lf\n", controller.mu_tilt_angle.fx.negative_large);
    // printf("mu_tilt_angle_NS = %lf\n", controller.mu_tilt_angle.fx.negative_small);
    // printf("mu_tilt_angle_Z = %lf\n", controller.mu_tilt_angle.fx.zero);
    // printf("mu_tilt_angle_PS = %lf\n", controller.mu_tilt_angle.fx.positive_small);
    // printf("mu_tilt_angle_PB = %lf\n", controller.mu_tilt_angle.fx.positive_large);

    // calculate the numerator and denominator of fuzzy basis function xi of tilt angle
    controller.xi_denominator_tilt_angle = 0;
    int index1 = 0;
    for (int j = 0; j < n_membership; j++){
        for (int k = 0; k < n_membership; k++){
            double mu1 = membershipgrade(&controller.mu_tilt_angle.fx, j);
            double mu2 = membershipgrade(&controller.mu_tilt_angle.gx, k);
            controller.xi_numerator_tilt_angle[index1] = mu1 * mu2; // numerator of fuzzy basis function xi
            controller.xi_denominator_tilt_angle += mu1 * mu2;     // denominator of fuzzy basis function xi
            index1++;
        }
    }

    // fuzzy basis function xi of tilt angle
    for (int j = 0; j < dimension; j++){
        if (controller.xi_denominator_tilt_angle != 0){
            controller.xi_tilt_angle[j] = controller.xi_numerator_tilt_angle[j] / controller.xi_denominator_tilt_angle;
        }
        else{
            printf("error: denominator of xi_tilt_angle is zero.\n");
            break;
        }
    }

    // for (int j = 0; j < dimension; j++) {
    //     printf("xi_tilt_angle[%d]: %f\n", j, controller.xi_tilt_angle[j]);
    // }

    // fuzzy basis function eta of tilt angle
    for (int j = 0; j < dimension; j++){
        controller.eta_tilt_angle[j] = controller.xi_tilt_angle[j] + 0.0001;
    }

    // calculate parameter of f(x) multiplied by fuzzy basis function xi
    double theta_tilt_angle_f_xi = 0.0;
    for (int j = 0; j < dimension; j++){
        theta_tilt_angle_f_xi += controller.theta_tilt_angle_fx[j] * controller.xi_tilt_angle[j];
    }
    printf("theta_tilt_angle_f_xi: %f\n", theta_tilt_angle_f_xi);

    // calculate parameter of g(x) multiplied by fuzzy basis function eta
    double theta_tilt_angle_g_eta = 0.0;
    for (int j = 0; j < dimension; j++){
        theta_tilt_angle_g_eta += controller.theta_tilt_angle_gx[j] * controller.eta_tilt_angle[j];
    }

    // control amount of tilt angle
    controller.control_tilt_angle = (1.0 / theta_tilt_angle_g_eta) * (-theta_tilt_angle_f_xi + controller.c * system_state.tilt_velocity - controller.k1 * sign(controller.sliding) - controller.k2 * controller.sliding);
    printf("control_tilt_angle: %f\n", controller.control_tilt_angle);

    // derivative of fuzzy controller parameter of f(x) of tilt angle
    for (int j = 0; j < dimension; j++){
        controller.theta_derivative_tilt_angle_fx[j] = controller.gamma1 * controller.sliding * controller.xi_tilt_angle[j];
    }

    // derivative of fuzzy controller parameter of g(x) of tilt angle
    for (int j = 0; j < dimension; j++){
        controller.theta_derivative_tilt_angle_gx[j] = controller.gamma2 * controller.sliding * controller.eta_tilt_angle[j] * controller.control_tilt_angle;
    }

    // update fuzzy controller parameter of tilt angle
    for (int j = 0; j < dimension; j++){
        controller.theta_tilt_angle_fx[j] += controller.theta_derivative_tilt_angle_fx[j] * Ts;
        controller.theta_tilt_angle_gx[j] += controller.theta_derivative_tilt_angle_gx[j] * Ts;
    }

    // velocity subsystem
    // error of longitudinal velocity of the vehicle
    controller.error_velocity = system_state.velocity - planner.velocity_plan;
    printf("error_velocity: %f\n", controller.error_velocity);
    archive.error_velocity_archive[i] = controller.error_velocity;

    controller.mu_velocity.fx.negative_large = gaussian(controller.error_velocity, -PI / 2.0, PI / 8.0); // membership function of f(x) of longitudinal velocity is negative large
    controller.mu_velocity.fx.negative_small = gaussian(controller.error_velocity, -PI / 4.0, PI / 8.0); // membership function of f(x) of longitudinal velocity is negative small
    controller.mu_velocity.fx.zero = gaussian(controller.error_velocity, 0.0, PI / 8.0);                 // membership function of f(x) of longitudinal velocity is zero
    controller.mu_velocity.fx.positive_small = gaussian(controller.error_velocity, PI / 4.0, PI / 8.0);  // membership function of f(x) of longitudinal velocity is positive small
    controller.mu_velocity.fx.positive_large = gaussian(controller.error_velocity, PI / 2.0, PI / 8.0);  // membership function of f(x) of longitudinal velocity is positive large

    controller.mu_velocity.gx.negative_large = gaussian(controller.error_velocity, -PI / 2.0, PI / 8.0); // membership function of g(x) of longitudinal velocity is negative large
    controller.mu_velocity.gx.negative_small = gaussian(controller.error_velocity, -PI / 4.0, PI / 8.0); // membership function of g(x) of longitudinal velocity is negative small
    controller.mu_velocity.gx.zero = gaussian(controller.error_velocity, 0.0, PI / 8.0);                 // membership function of g(x) of longitudinal velocity is zero
    controller.mu_velocity.gx.positive_small = gaussian(controller.error_velocity, PI / 4.0, PI / 8.0);  // membership function of g(x) of longitudinal velocity is positive small
    controller.mu_velocity.gx.positive_large = gaussian(controller.error_velocity, PI / 2.0, PI / 8.0);  // membership function of g(x) of longitudinal velocity is positive large

    double* membership2_value[] = {
        &controller.mu_velocity.fx.negative_large,
        &controller.mu_velocity.fx.negative_small,
        &controller.mu_velocity.fx.zero,
        &controller.mu_velocity.fx.positive_small,
        &controller.mu_velocity.fx.positive_large,
        &controller.mu_velocity.gx.negative_large,
        &controller.mu_velocity.gx.negative_small,
        &controller.mu_velocity.gx.zero,
        &controller.mu_velocity.gx.positive_small,
        &controller.mu_velocity.gx.positive_large
    };

    int num_membership2_value = sizeof(membership2_value) / sizeof(double*);

    for (int j = 0; j < num_membership2_value; j++) {
        if (*membership2_value[j] < 0.0001) {
            *membership2_value[j] = 0.0001;
        }
    }

    // printf("mu_velocity_NB = %lf\n", controller.mu_velocity.fx.negative_large);
    // printf("mu_velocity_NS = %lf\n", controller.mu_velocity.fx.negative_small);
    // printf("mu_velocity_Z = %lf\n", controller.mu_velocity.fx.zero);
    // printf("mu_velocity_PS = %lf\n", controller.mu_velocity.fx.positive_small);
    // printf("mu_velocity_PB = %lf\n", controller.mu_velocity.fx.positive_large);

    // calculate the numerator and denominator of fuzzy basis function xi of longitudinal velocity
    controller.xi_denominator_velocity = 0;
    int index2 = 0;
    for (int j = 0; j < n_membership; j++){
        for (int k = 0; k < n_membership; k++){
            double mu1 = membershipgrade(&controller.mu_velocity.fx, j);
            double mu2 = membershipgrade(&controller.mu_velocity.gx, k);
            controller.xi_numerator_velocity[index2] = mu1 * mu2; // numerator of fuzzy basis function xi
            controller.xi_denominator_velocity += mu1 * mu2;      // denominator of fuzzy basis function xi
            index2++;
        }
    }

    // fuzzy basis function xi of longitudinal velocity
    for (int j = 0; j < dimension; j++){
        if (controller.xi_denominator_velocity != 0){
            controller.xi_velocity[j] = controller.xi_numerator_velocity[j] / controller.xi_denominator_velocity;
        }
        else{
            printf("error: denominator of xi_velocity is zero.\n");
            break;
        }
    }

    // for (int j = 0; j < dimension; j++) {
    //     printf("xi_velocity[%d]: %f\n", j, controller.xi_velocity[j]);
    // }

    // fuzzy basis function eta of longitudinal velocity
    for (int j = 0; j < dimension; j++){
        controller.eta_velocity[j] = controller.xi_velocity[j] + 0.0001;
    }

    // calculate parameter of f(x) multiplied by fuzzy basis function xi
    double theta_velocity_f_xi = 0.0;
    for (int j = 0; j < dimension; j++){
        theta_velocity_f_xi += controller.theta_velocity_fx[j] * controller.xi_velocity[j];
    }
    printf("theta_velocity_f_xi: %f\n", theta_velocity_f_xi);

    // calculate parameter of g(x) multiplied by fuzzy basis function eta
    double theta_velocity_g_eta = 0.0;
    for (int j = 0; j < dimension; j++){
        theta_velocity_g_eta += controller.theta_velocity_gx[j] * controller.eta_velocity[j];
    }

    // control amount of longitudinal velocity
    controller.control_velocity = (1.0 / theta_velocity_g_eta) * (theta_velocity_f_xi - controller.k3 * controller.error_velocity);
    printf("control_velocity: %f\n", controller.control_velocity);

    // derivative of fuzzy controller parameter of f(x) of longitudinal velocity
    for (int j = 0; j < dimension; j++){
        controller.theta_derivative_velocity_fx[j] = controller.gamma3 * controller.error_velocity * controller.xi_velocity[j];
    }

    // for (int j = 0; j < dimension; j++) {
    //     printf("theta_derivative_velocity_fx[%d]: %f\n", j, controller.theta_derivative_velocity_fx[j]);
    // }

    // derivative of fuzzy controller parameter of g(x) of longitudinal velocity
    for (int j = 0; j < dimension; j++){
        controller.theta_derivative_velocity_gx[j] = controller.gamma4 * controller.error_velocity * controller.eta_velocity[j] * controller.control_velocity;
    }

    // for (int j = 0; j < dimension; j++) {
    //     printf("theta_derivative_velocity_gx[%d]: %f\n", j, controller.theta_derivative_velocity_gx[j]);
    // }

    // update fuzzy controller parameter of longitudinal velocity
    for (int j = 0; j < dimension; j++){
        controller.theta_velocity_fx[j] += controller.theta_derivative_velocity_fx[j] * Ts;
        controller.theta_velocity_gx[j] += controller.theta_derivative_velocity_gx[j] * Ts;
    }

    // control input torque of longitudinal velocity subsystem
    controller.torque_velocity = controller.control_tilt_angle + controller.control_velocity;
    printf("torque_velocity: %f\n", controller.torque_velocity);
    archive.torque_velocity_archive[i] = controller.torque_velocity;

    // rotational velocity subsystem
    // error of rotational velocity of the vehicle
    controller.error_rotational_velocity = system_state.rotational_velocity - planner.rotational_velocity_plan;
    printf("error_rotational_velocity: %f\n", controller.error_rotational_velocity);
    archive.error_rotational_velocity_archive[i] = controller.error_rotational_velocity;

    controller.mu_rotational_velocity.fx.negative_large = gaussian(controller.error_rotational_velocity, -PI / 4.0, PI / 16.0); // membership function of f(x) of rotational velocity is negative large
    controller.mu_rotational_velocity.fx.negative_small = gaussian(controller.error_rotational_velocity, -PI / 8.0, PI / 16.0); // membership function of f(x) of rotational velocity is negative small
    controller.mu_rotational_velocity.fx.zero = gaussian(controller.error_rotational_velocity, 0.0, PI / 16.0);                 // membership function of f(x) of rotational velocity is zero
    controller.mu_rotational_velocity.fx.positive_small = gaussian(controller.error_rotational_velocity, PI / 8.0, PI / 16.0);  // membership function of f(x) of rotational velocity is positive small
    controller.mu_rotational_velocity.fx.positive_large = gaussian(controller.error_rotational_velocity, PI / 4.0, PI / 16.0);  // membership function of f(x) of rotational velocity is positive large

    controller.mu_rotational_velocity.gx.negative_large = gaussian(controller.error_rotational_velocity, -PI / 4.0, PI / 16.0); // membership function of g(x) of rotational velocity is negative large
    controller.mu_rotational_velocity.gx.negative_small = gaussian(controller.error_rotational_velocity, -PI / 8.0, PI / 16.0); // membership function of g(x) of rotational velocity is negative small
    controller.mu_rotational_velocity.gx.zero = gaussian(controller.error_rotational_velocity, 0.0, PI / 24.0);                 // membership function of g(x) of rotational velocity is zero
    controller.mu_rotational_velocity.gx.positive_small = gaussian(controller.error_rotational_velocity, PI / 8.0, PI / 16.0);  // membership function of g(x) of rotational velocity is positive small
    controller.mu_rotational_velocity.gx.positive_large = gaussian(controller.error_rotational_velocity, PI / 4.0, PI / 16.0);  // membership function of g(x) of rotational velocity is positive large

    double* membership3_value[] = {
        &controller.mu_rotational_velocity.fx.negative_large,
        &controller.mu_rotational_velocity.fx.negative_small,
        &controller.mu_rotational_velocity.fx.zero,
        &controller.mu_rotational_velocity.fx.positive_small,
        &controller.mu_rotational_velocity.fx.positive_large,
        &controller.mu_rotational_velocity.gx.negative_large,
        &controller.mu_rotational_velocity.gx.negative_small,
        &controller.mu_rotational_velocity.gx.zero,
        &controller.mu_rotational_velocity.gx.positive_small,
        &controller.mu_rotational_velocity.gx.positive_large
    };

    int num_membership3_value = sizeof(membership3_value) / sizeof(double*);

    for (int j = 0; j < num_membership3_value; j++) {
        if (*membership3_value[j] < 0.0001) {
            *membership3_value[j] = 0.0001;
        }
    }

    // printf("mu_rotational_velocity_NB = %lf\n", controller.mu_rotational_velocity.fx.negative_large);
    // printf("mu_rotational_velocity_NS = %lf\n", controller.mu_rotational_velocity.fx.negative_small);
    // printf("mu_rotational_velocity_Z = %lf\n", controller.mu_rotational_velocity.fx.zero);
    // printf("mu_rotational_velocity_PS = %lf\n", controller.mu_rotational_velocity.fx.positive_small);
    // printf("mu_rotational_velocity_PB = %lf\n", controller.mu_rotational_velocity.fx.positive_large);

    // calculate the numerator and denominator of fuzzy basis function xi of rotational velocity
    controller.xi_denominator_rotational_velocity = 0;
    int index3 = 0;
    for (int j = 0; j < n_membership; j++){
        for (int k = 0; k < n_membership; k++){
            double mu1 = membershipgrade(&controller.mu_rotational_velocity.fx, j);
            double mu2 = membershipgrade(&controller.mu_rotational_velocity.gx, k);
            controller.xi_numerator_rotational_velocity[index3] = mu1 * mu2; // numerator of fuzzy basis function xi
            controller.xi_denominator_rotational_velocity += mu1 * mu2;     // denominator of fuzzy basis function xi
            index3++;
        }
    }

    // fuzzy basis function xi of rotational velocity
    for (int j = 0; j < dimension; j++){
        if (controller.xi_denominator_rotational_velocity != 0){
            controller.xi_rotational_velocity[j] = controller.xi_numerator_rotational_velocity[j] / controller.xi_denominator_rotational_velocity;
        }
        else{
            printf("error: denominator of xi_rotational_velocity is zero.\n");
            break;
        }
    }

    // for (int j = 0; j < dimension; j++) {
    //     printf("xi_rotational_velocity[%d]: %f\n", j, controller.xi_rotational_velocity[j]);
    // }

    // fuzzy basis function eta of rotational velocity
    for (int j = 0; j < dimension; j++){
        controller.eta_rotational_velocity[j] = controller.xi_rotational_velocity[j] + 0.0001;
    }

    // calculate parameter of f(x) multiplied by fuzzy basis function xi
    double theta_rotational_velocity_f_xi = 0.0;
    for (int j = 0; j < dimension; j++){
        theta_rotational_velocity_f_xi += controller.theta_rotational_velocity_fx[j] * controller.xi_rotational_velocity[j];
    }

    // calculate parameter of g(x) multiplied by fuzzy basis function eta
    double theta_rotational_velocity_g_eta = 0.0;
    for (int j = 0; j < dimension; j++){
        theta_rotational_velocity_g_eta += controller.theta_rotational_velocity_gx[j] * controller.eta_rotational_velocity[j];
    }

    // control input torque of rotational velocity
    controller.torque_rotational_velocity = (1.0 / theta_rotational_velocity_g_eta) * (theta_rotational_velocity_f_xi - controller.k4 * controller.error_rotational_velocity);
    printf("torque_rotational_velocity: %f\n", controller.torque_rotational_velocity);
    archive.torque_rotational_velocity_archive[i] = controller.torque_rotational_velocity;

    // derivative of fuzzy controller parameter of f(x) of rotational velocity
    for (int j = 0; j < dimension; j++){
        controller.theta_derivative_rotational_velocity_fx[j] = controller.gamma5 * controller.error_rotational_velocity * controller.xi_rotational_velocity[j];
    }

    // derivative of fuzzy controller parameter of g(x) of rotational velocity
    for (int j = 0; j < dimension; j++){
        controller.theta_derivative_rotational_velocity_gx[j] = controller.gamma6 * controller.error_rotational_velocity * controller.eta_rotational_velocity[j] * controller.torque_rotational_velocity;
    }

    // update fuzzy controller parameter of rotational velocity
    for (int j = 0; j < dimension; j++){
        controller.theta_rotational_velocity_fx[j] += controller.theta_derivative_rotational_velocity_fx[j] * Ts;
        controller.theta_rotational_velocity_gx[j] += controller.theta_derivative_rotational_velocity_gx[j] * Ts;
    }

    // control input torque of right wheel of vehicle
    controller.torque_right = (controller.torque_velocity + controller.torque_rotational_velocity) / 2.0;
    // control input torque of left wheel of vehicle
    controller.torque_left = (controller.torque_velocity - controller.torque_rotational_velocity) / 2.0;

    printf("torque_right: %f\n", controller.torque_right);
    printf("torque_left: %f\n", controller.torque_left);

    WbDeviceTag motor[4];
    char motor_tag[2][8] = {"motor1", "motor2"};
    for (int j = 0; j < 2; j++){
        motor[j] = wb_robot_get_device(motor_tag[j]);
    }

    double torque[2] = {controller.torque_right, controller.torque_left};
    for (int j = 0; j < 2; j++){
        wb_motor_set_torque(motor[j], torque[j]);
    }

    controller.controller_out[0] = controller.torque_right;
    controller.controller_out[1] = controller.torque_left;
}

void saveArchiveToTxt(double *archive, int size, const char *filename){

    FILE *file = fopen(filename, "w+");

    if (file == NULL){
        perror("Failed to open file");
        exit(1);
    }
    else{
        for (int i = 0; i < size; i++){
            fprintf(file, "%lf\n", archive[i]);
        }
        fclose(file);
        printf("Saved to file %s\n", filename);
    }
}

void saveArchive(){
    saveArchiveToTxt(archive.x_coordinate_archive, ARRAY_SIZE, "../../../report/x_coordinate_archive.txt");
    saveArchiveToTxt(archive.y_coordinate_archive, ARRAY_SIZE, "../../../report/y_coordinate_archive.txt");
    saveArchiveToTxt(archive.rotational_angle_archive, ARRAY_SIZE, "../../../report/rotational_angle.txt");
    saveArchiveToTxt(archive.rotational_velocity_archive, ARRAY_SIZE, "../../../report/rotational_velocity.txt");
    saveArchiveToTxt(archive.error_x_coordinate_archive, ARRAY_SIZE, "../../../report/error_x_coordinate.txt");
    saveArchiveToTxt(archive.error_y_coordinate_archive, ARRAY_SIZE, "../../../report/error_y_coordinate.txt");
    saveArchiveToTxt(archive.error_rotational_angle_archive, ARRAY_SIZE, "../../../report/error_rotational_angle.txt");
    saveArchiveToTxt(archive.velocity_plan_archive, ARRAY_SIZE, "../../../report/velocity_plan.txt");
    saveArchiveToTxt(archive.rotational_velocity_plan_archive, ARRAY_SIZE, "../../../report/rotational_velocity_plan.txt");
    saveArchiveToTxt(archive.velocity_archive, ARRAY_SIZE, "../../../report/velocity.txt");
    saveArchiveToTxt(archive.tilt_angle_archive, ARRAY_SIZE, "../../../report/tilt_angle.txt");
    saveArchiveToTxt(archive.error_velocity_archive, ARRAY_SIZE, "../../../report/error_velocity.txt");
    saveArchiveToTxt(archive.error_rotational_velocity_archive, ARRAY_SIZE, "../../../report/error_rotational_velocity.txt");
    saveArchiveToTxt(archive.torque_velocity_archive, ARRAY_SIZE, "../../../report/torque_velocity.txt");
    saveArchiveToTxt(archive.torque_rotational_velocity_archive, ARRAY_SIZE, "../../../report/torque_rotational_velocity.txt");
}

int main(int argc, char **argv){

    PLANNER_init();    // initialize planner parameter
    CONTROLLER_init(); // initialize controller parameter
    // PLANT_init();      // initialize plant parameter

    int i = 0;
    while (wb_robot_step(TIME_STEP) != -1){
        double time = i * TIME_STEP / 1000.0 + t0;
        printf("time at step %d: %f\n", i, time);

        if (time > 30.0){
            break;
        }

        PLANNER_realize(i);
        CONTROLLER_realize(i);
        // PLANT_realize(i);
        i++;
    }

    saveArchive();

    wb_robot_cleanup();

    return 0;
}
