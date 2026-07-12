void fwd_dyn_crba( double *q, double *qd, double *tau, double *out_qdd ) {
    // The outward pass of the computation of C.  This is the first for loop of Table 5.1 in RBDA. 
    
    // Body 0
    double q0 = q[ 0 ];
    double qd0 = qd[ 0 ];
    double sin_q0 = sin( q0 );
    double cos_q0 = cos( q0 );
    
    double S0_0 = 0;
    double S0_1 = 0;
    double S0_2 = 1;
    double S0_3 = 0;
    double S0_4 = 0;
    double S0_5 = 0;
    
    double vJ0_0 = 0;
    double vJ0_1 = 0;
    double vJ0_2 = qd0;
    double vJ0_3 = 0;
    double vJ0_4 = 0;
    double vJ0_5 = 0;
    
    double X0_E00 = cos_q0;
    double X0_E01 = -sin_q0;
    double X0_E02 = 0;
    double X0_E10 = sin_q0;
    double X0_E11 = cos_q0;
    double X0_E12 = 0;
    double X0_E20 = 0;
    double X0_E21 = 0;
    double X0_E22 = 1;
    double X0_r0 = 1 + cos_q0;
    double X0_r1 = -sin_q0;
    double X0_r2 = 0;
    
    double Base_E00 = 1;
    double Base_E01 = 0;
    double Base_E02 = 0;
    double Base_E10 = 0;
    double Base_E11 = 1;
    double Base_E12 = 0;
    double Base_E20 = 0;
    double Base_E21 = 0;
    double Base_E22 = 1;
    double Base_r0 = 0;
    double Base_r1 = 0;
    double Base_r2 = 0;
    
    // Body 1
    double q1 = q[ 1 ];
    double qd1 = qd[ 1 ];
    double sin_q1 = sin( q1 );
    double cos_q1 = cos( q1 );
    
    double S1_0 = 0;
    double S1_1 = 0;
    double S1_2 = 1;
    double S1_3 = 0;
    double S1_4 = 0;
    double S1_5 = 0;
    
    double vJ1_0 = 0;
    double vJ1_1 = 0;
    double vJ1_2 = qd1;
    double vJ1_3 = 0;
    double vJ1_4 = 0;
    double vJ1_5 = 0;
    
    double X1_E00 = cos_q1;
    double X1_E01 = -sin_q1;
    double X1_E02 = 0;
    double X1_E10 = sin_q1;
    double X1_E11 = cos_q1;
    double X1_E12 = 0;
    double X1_E20 = 0;
    double X1_E21 = 0;
    double X1_E22 = 1;
    double X1_r0 = 1 + cos_q1;
    double X1_r1 = -sin_q1;
    double X1_r2 = 0;
    
    double X0_1_E00 = cos_q1;
    double X0_1_E01 = -sin_q1;
    double X0_1_E02 = 0;
    double X0_1_E10 = sin_q1;
    double X0_1_E11 = cos_q1;
    double X0_1_E12 = 0;
    double X0_1_E20 = 0;
    double X0_1_E21 = 0;
    double X0_1_E22 = 1;
    double X0_1_r0 = 0;
    double X0_1_r1 = 0;
    double X0_1_r2 = 0;
    
    // Body 2
    double q2 = q[ 2 ];
    double qd2 = qd[ 2 ];
    double sin_q2 = sin( q2 );
    double cos_q2 = cos( q2 );
    
    double S2_0 = 0;
    double S2_1 = 0;
    double S2_2 = 1;
    double S2_3 = 0;
    double S2_4 = 0;
    double S2_5 = 0;
    
    double vJ2_0 = 0;
    double vJ2_1 = 0;
    double vJ2_2 = qd2;
    double vJ2_3 = 0;
    double vJ2_4 = 0;
    double vJ2_5 = 0;
    
    double X2_E00 = cos_q2;
    double X2_E01 = -sin_q2;
    double X2_E02 = 0;
    double X2_E10 = sin_q2;
    double X2_E11 = cos_q2;
    double X2_E12 = 0;
    double X2_E20 = 0;
    double X2_E21 = 0;
    double X2_E22 = 1;
    double X2_r0 = 1 + cos_q2;
    double X2_r1 = -sin_q2;
    double X2_r2 = 0;
    
    double X0_2_E00 = cos_q2*cos_q1 + -sin_q2*sin_q1;
    double X0_2_E01 = cos_q2*-sin_q1 + -sin_q2*cos_q1;
    double X0_2_E02 = 0;
    double X0_2_E10 = sin_q2*cos_q1 + cos_q2*sin_q1;
    double X0_2_E11 = sin_q2*-sin_q1 + cos_q2*cos_q1;
    double X0_2_E12 = 0;
    double X0_2_E20 = 0;
    double X0_2_E21 = 0;
    double X0_2_E22 = 1;
    double X0_2_r0 = 0;
    double X0_2_r1 = 0;
    double X0_2_r2 = 0;
    
}
