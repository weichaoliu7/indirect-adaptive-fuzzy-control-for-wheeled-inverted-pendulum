#include <webots/motor.h>
#include <webots/robot.h>
#include <webots/gps.h>
#include <webots/inertial_unit.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// global variables declaration
#define TIME_STEP 10
#define PI 3.14159
#define ARRAY_SIZE 3000                    // sampling times
#define n_membership 5                     // number of membership functions
#define dimension 25                       // dimension of fuzzy controller parameter
#define l 1.0                              // distance from the body center of gravity to the wheel axis
static double t0 = 0.0;                    // start time
static double t1 = 30.0;                   // end time
static double Ts = 0.010;                  // sampling period

/* reference: [1]Yue M, An C, Du Y, et al. Indirect adaptive fuzzy control for a nonholonomic/underactuated
wheeled inverted pendulum vehicle based on a data-driven trajectory planner[J]. Fuzzy Sets and Systems, 2016, 290: 158-177. */

struct _archive{
    double x_coordinate_archive[ARRAY_SIZE];
    double y_coordinate_archive[ARRAY_SIZE];
    double velocity_archive[ARRAY_SIZE];
    double tilt_angle_archive[ARRAY_SIZE];
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

struct _controller{
    double controller_u[5];
    double controller_out[6];
    double z_coordinate;                              // z-coordinate of the vehicle
    double x1_coordinate;                             // x-coordinate of top of inverted pendulum
    double y1_coordinate;                             // y-coordinate of top of inverted pendulum
    double z1_coordinate;                             // z-coordinate of top of inverted pendulum
    double sliding;                                   // sliding mode manifold of tilt angle subsystem
    membershipfunction mu_tilt_angle;                 // membership function of tilt angle
    double xi_denominator_tilt_angle;                 // denominator of fuzzy basis function xi of tilt angle membership function
    double xi_numerator_tilt_angle[dimension];        // numerator of fuzzy basis function xi of tilt angle membership function
    double xi_tilt_angle[dimension];                  // fuzzy basis function xi of tilt angle membership function
    double eta_tilt_angle[dimension];                 // fuzzy basis function eta of tilt angle membership function
    double theta_tilt_angle_fx[dimension];            // fuzzy controller parameter of f(x) of tilt angle
    double theta_tilt_angle_gx[dimension];            // fuzzy controller parameter of g(x) of tilt angle
    double theta_derivative_tilt_angle_fx[dimension]; // derivative of fuzzy controller parameter of f(x) of tilt angle
    double theta_derivative_tilt_angle_gx[dimension]; // derivative of fuzzy controller parameter of g(x) of tilt angle
    double control_tilt_angle;                        // control amount of tilt angle
    double torque_rotational_velocity;                // control input torque of rotational velocity
    double error_velocity;                            // error of longitudinal velocity of the vehicle
    double error_rotational_velocity;                 // error of rotational velocity of the vehicle
    double torque_velocity;                           // control input torque of longitudinal velocity subsystem
    double torque_right;                              // control input torque of right wheel of vehicle
    double torque_left;                               // control input torque of left wheel of vehicle
    double c;                                         // parameter of sliding mode manifold
    double k1;                                        // control gain
    double k2;                                        // control gain
    double gamma1;                                    // controller parameter update gain
    double gamma2;                                    // controller parameter update gain
} controller;

void CONTROLLER_init(){
    wb_robot_init();
    WbDeviceTag gps1 = wb_robot_get_device("gps1");
    wb_gps_enable(gps1, TIME_STEP);
    WbDeviceTag gps2 = wb_robot_get_device("gps2");
    wb_gps_enable(gps2, TIME_STEP);
    WbDeviceTag motor[2];
    char motor_tag[2][8] = {"motor1", "motor2"};
    for (int j = 0; j < 2; j++){
        motor[j] = wb_robot_get_device(motor_tag[j]);
    }

    system_state.rotational_velocity = 0.0;
    system_state.velocity = 0.0;
    system_state.tilt_angle = 0.0;
    controller.c = 2.14;       // parameter of sliding mode manifold
    controller.k1 = 2.43;      // control gain
    controller.k2 = 3.14;      // control gain

    for (int j = 0; j < dimension; j++){
        controller.theta_tilt_angle_fx[j] = 0.1;          // fuzzy controller parameter of f(x) of tilt angle
        controller.theta_tilt_angle_gx[j] = 0.1;          // fuzzy controller parameter of g(x) of tilt angle
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

    // obtain the horizontal and vertical coordinates of top of inverted pendulum
    WbDeviceTag gps2 = wb_robot_get_device("gps2");
    wb_gps_enable(gps2, TIME_STEP);
    const double *coordinate2 = wb_gps_get_values(gps2);
    controller.x1_coordinate = coordinate2[0]; // x-coordinate of midpoint of top of inverted pendulum
    controller.y1_coordinate = coordinate2[1]; // y-coordinate of midpoint of top of inverted pendulum
    controller.z1_coordinate = coordinate2[2]; // y-coordinate of midpoint of top of inverted pendulum

    // tilt angle of the vehicle body
    double tilt_angle = atan(sqrt(pow(controller.x1_coordinate - system_state.x_coordinate, 2) + pow(controller.y1_coordinate - system_state.y_coordinate, 2)) / (controller.z1_coordinate - controller.z_coordinate));

    WbDeviceTag gps5 = wb_robot_get_device("gps5");
    wb_gps_enable(gps5, TIME_STEP);
    const double *coordinate5 = wb_gps_get_values(gps5);
    WbDeviceTag gps6 = wb_robot_get_device("gps6");
    wb_gps_enable(gps6, TIME_STEP);
    const double *coordinate6 = wb_gps_get_values(gps6);

    if (coordinate5[2] >= coordinate6[2]){
        system_state.tilt_angle = tilt_angle;
    } else{
        system_state.tilt_angle = - tilt_angle;
    }

    printf("tilt_angle: %f\n", system_state.tilt_angle);
    archive.tilt_angle_archive[i] = system_state.tilt_angle;

    if (i == 0){
        system_state.tilt_velocity = 0.001; // tilt velocity of the vehicle body
    }
    else if (i > 0){
        system_state.tilt_velocity = (system_state.tilt_angle - archive.tilt_angle_archive[i - 1]) / Ts; // tilt velocity of the vehicle body
    }
    printf("tilt_velocity: %f\n", system_state.tilt_velocity);

    // controller.controller_u[0] = planner.rotational_velocity_plan;
    // controller.controller_u[1] = planner.velocity_plan;
    // controller.controller_u[2] = system_state.rotational_velocity;
    // controller.controller_u[3] = system_state.velocity;
    // controller.controller_u[4] = system_state.tilt_angle;

    // tilt angle subsystem
    controller.sliding = system_state.tilt_velocity + controller.c * system_state.tilt_angle; // sliding mode manifold
    printf("sliding_term: %f\n", controller.sliding);

    controller.mu_tilt_angle.fx.negative_large = gaussian(system_state.tilt_angle, -PI / 30.0, PI / 120.0);  // membership function of f(x) of tilt angle is negative large
    controller.mu_tilt_angle.fx.negative_small = gaussian(system_state.tilt_angle, -PI / 60.0, PI / 120.0); // membership function of f(x) of tilt angle is negative small
    controller.mu_tilt_angle.fx.zero = gaussian(system_state.tilt_angle, 0.0, PI / 120.0);                  // membership function of f(x) of tilt angle is zero
    controller.mu_tilt_angle.fx.positive_small = gaussian(system_state.tilt_angle, PI / 60.0, PI / 120.0);  // membership function of f(x) of tilt angle is positive small
    controller.mu_tilt_angle.fx.positive_large = gaussian(system_state.tilt_angle, PI / 30.0, PI / 120.0);   // membership function of f(x) of tilt angle is positive large

    controller.mu_tilt_angle.gx.negative_large = gaussian(system_state.tilt_angle, -PI / 30.0, PI / 120.0);  // membership function of g(x) of tilt angle is negative large
    controller.mu_tilt_angle.gx.negative_small = gaussian(system_state.tilt_angle, -PI / 60.0, PI / 120.0); // membership function of g(x) of tilt angle is negative small
    controller.mu_tilt_angle.gx.zero = gaussian(system_state.tilt_angle, 0.0, PI / 120.0);                  // membership function of g(x) of tilt angle is zero
    controller.mu_tilt_angle.gx.positive_small = gaussian(system_state.tilt_angle, PI / 60.0, PI / 120.0);  // membership function of g(x) of tilt angle is positive small
    controller.mu_tilt_angle.gx.positive_large = gaussian(system_state.tilt_angle, PI / 30.0, PI / 120.0);   // membership function of g(x) of tilt angle is positive large

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

    printf("mu_tilt_angle_NB = %lf\n", controller.mu_tilt_angle.fx.negative_large);
    printf("mu_tilt_angle_NS = %lf\n", controller.mu_tilt_angle.fx.negative_small);
    printf("mu_tilt_angle_Z = %lf\n", controller.mu_tilt_angle.fx.zero);
    printf("mu_tilt_angle_PS = %lf\n", controller.mu_tilt_angle.fx.positive_small);
    printf("mu_tilt_angle_PB = %lf\n", controller.mu_tilt_angle.fx.positive_large);

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

    // control input torque of right wheel of vehicle
    if (i < 25){
        controller.torque_right = 17;
        controller.torque_left = 17;
    } else{
        controller.torque_right = (controller.control_tilt_angle) / 2.0;
        controller.torque_left = (controller.control_tilt_angle) / 2.0;
    }

    if ((i >= 201 && i <= 225) || (i >= 401 && i <= 425) || (i >= 601 && i <= 625) || (i >= 801 && i <= 825)) {
        controller.torque_right = 17;
        controller.torque_left = 17;
    } else {
        controller.torque_right = (controller.control_tilt_angle) / 2.0;
        controller.torque_left = (controller.control_tilt_angle) / 2.0;
    }

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
    saveArchiveToTxt(archive.velocity_archive, ARRAY_SIZE, "../../../report/velocity.txt");
    saveArchiveToTxt(archive.tilt_angle_archive, ARRAY_SIZE, "../../../report/tilt_angle.txt");
    saveArchiveToTxt(archive.torque_velocity_archive, ARRAY_SIZE, "../../../report/torque_velocity.txt");
    saveArchiveToTxt(archive.torque_rotational_velocity_archive, ARRAY_SIZE, "../../../report/torque_rotational_velocity.txt");
}

int main(int argc, char **argv){

    CONTROLLER_init(); // initialize controller parameter
    // PLANT_init();      // initialize plant parameter

    int i = 0;
    while (wb_robot_step(TIME_STEP) != -1){
        double time = i * TIME_STEP / 1000.0 + t0;
        printf("time at step %d: %f\n", i, time);

        if (time > 30.0){
            break;
        }

        CONTROLLER_realize(i);
        // PLANT_realize(i);
        i++;
    }

    saveArchive();

    wb_robot_cleanup();

    return 0;
}
