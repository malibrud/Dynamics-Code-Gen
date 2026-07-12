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
    
    double Xlink0_E00 = cos_q0;
    double Xlink0_E01 = -sin_q0;
    double Xlink0_E02 = 0;
    double Xlink0_E10 = sin_q0;
    double Xlink0_E11 = cos_q0;
    double Xlink0_E12 = 0;
    double Xlink0_E20 = 0;
    double Xlink0_E21 = 0;
    double Xlink0_E22 = 1;
    double Xlink0_r0 = 1;
    double Xlink0_r1 = 0;
    double Xlink0_r2 = 0;
    
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
    
    double Xlink1_E00 = cos_q1;
    double Xlink1_E01 = -sin_q1;
    double Xlink1_E02 = 0;
    double Xlink1_E10 = sin_q1;
    double Xlink1_E11 = cos_q1;
    double Xlink1_E12 = 0;
    double Xlink1_E20 = 0;
    double Xlink1_E21 = 0;
    double Xlink1_E22 = 1;
    double Xlink1_r0 = 1;
    double Xlink1_r1 = 0;
    double Xlink1_r2 = 0;
    
    double Xglobal1_E00 = Xlink1_E00*Xlink0_E00 + Xlink1_E01*Xlink0_E10;
    double Xglobal1_E01 = Xlink1_E00*Xlink0_E01 + Xlink1_E01*Xlink0_E11;
    double Xglobal1_E02 = 0;
    double Xglobal1_E10 = Xlink1_E10*Xlink0_E00 + Xlink1_E11*Xlink0_E10;
    double Xglobal1_E11 = Xlink1_E10*Xlink0_E01 + Xlink1_E11*Xlink0_E11;
    double Xglobal1_E12 = 0;
    double Xglobal1_E20 = 0;
    double Xglobal1_E21 = 0;
    double Xglobal1_E22 = 1;
    double Xglobal1_r0 = 1 + Xlink0_E00;
    double Xglobal1_r1 = Xlink0_E01;
    double Xglobal1_r2 = 0;
    
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
    
    double Xlink2_E00 = cos_q2;
    double Xlink2_E01 = -sin_q2;
    double Xlink2_E02 = 0;
    double Xlink2_E10 = sin_q2;
    double Xlink2_E11 = cos_q2;
    double Xlink2_E12 = 0;
    double Xlink2_E20 = 0;
    double Xlink2_E21 = 0;
    double Xlink2_E22 = 1;
    double Xlink2_r0 = 1;
    double Xlink2_r1 = 0;
    double Xlink2_r2 = 0;
    
    double Xglobal2_E00 = Xlink2_E00*Xglobal1_E00 + Xlink2_E01*Xglobal1_E10;
    double Xglobal2_E01 = Xlink2_E00*Xglobal1_E01 + Xlink2_E01*Xglobal1_E11;
    double Xglobal2_E02 = 0;
    double Xglobal2_E10 = Xlink2_E10*Xglobal1_E00 + Xlink2_E11*Xglobal1_E10;
    double Xglobal2_E11 = Xlink2_E10*Xglobal1_E01 + Xlink2_E11*Xglobal1_E11;
    double Xglobal2_E12 = 0;
    double Xglobal2_E20 = 0;
    double Xglobal2_E21 = 0;
    double Xglobal2_E22 = 1;
    double Xglobal2_r0 = Xglobal1_r0 + Xglobal1_E00;
    double Xglobal2_r1 = Xglobal1_r1 + Xglobal1_E01;
    double Xglobal2_r2 = 0;
    
}
