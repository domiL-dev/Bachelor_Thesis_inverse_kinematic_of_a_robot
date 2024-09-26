function [T6,PSI, T5,T4,T3,T2,T1, r6,r5,r4,r3,r2,r1] = TCP_PSI_zxxzxz(...
    S1,S2,S3,S4,S5,S6, Fx,Fy,Fz, phi1,phi2,phi3,phi4,phi5,phi6)

% Rot z
D1 = [cos(phi1) -sin(phi1)   0
      sin(phi1)  cos(phi1)   0
      0          0           1];

% Rot x
D2 = [1          0     0
      0     cos(phi2) -sin(phi2)
      0     sin(phi2)  cos(phi2)];


% Rot x
D3 = [1          0     0
      0     cos(phi3) -sin(phi3)
      0     sin(phi3)  cos(phi3)];
  
% Rot z 
D4 = [cos(phi4) -sin(phi4)   0
      sin(phi4)  cos(phi4)   0
      0          0           1];

% Rot x
D5 = [1          0     0
      0     cos(phi5) -sin(phi5)
      0     sin(phi5)  cos(phi5)];  
  
% Rot z 
D6 = [cos(phi6) -sin(phi6)   0
      sin(phi6)  cos(phi6)   0
      0          0           1];
 
T6  = D1*(S1 + D2*(S2 + D3*(S3 + D4*(S4 + D5*(S5 + D6*(S6)))))); 
T5  = D1*(S1 + D2*(S2 + D3*(S3 + D4*(S4 + D5*(S5)))));
T4  = D1*(S1 + D2*(S2 + D3*(S3 + D4*(S4))));
T3  = D1*(S1 + D2*(S2 + D3*(S3)));
T2  = D1*(S1 + D2*(S2));
T1  = D1*(S1);

 % Berechnung T6 FRAME
T6_Fx = D1*(S1 + D2*(S2 + D3*(S3 + D4*(S4 + D5*(S5 + D6*(S6+Fx)))))); 
T6_Fy = D1*(S1 + D2*(S2 + D3*(S3 + D4*(S4 + D5*(S5 + D6*(S6+Fy)))))); 
T6_Fz = D1*(S1 + D2*(S2 + D3*(S3 + D4*(S4 + D5*(S5 + D6*(S6+Fz)))))); 
r6    = [T6_Fx-T6  T6_Fy-T6  T6_Fz-T6];
 
% Berechnung T5 FRAME
T5_Fx = D1*(S1 + D2*(S2 + D3*(S3 + D4*(S4 + D5*(S5+Fx)))));
T5_Fy = D1*(S1 + D2*(S2 + D3*(S3 + D4*(S4 + D5*(S5+Fy)))));
T5_Fz = D1*(S1 + D2*(S2 + D3*(S3 + D4*(S4 + D5*(S5+Fz)))));
r5    = [T5_Fx-T5  T5_Fy-T5  T5_Fz-T5];

% Berechnung T4 FRAME
T4_Fx = D1*(S1 + D2*(S2 + D3*(S3 + D4*(S4+Fx))));
T4_Fy = D1*(S1 + D2*(S2 + D3*(S3 + D4*(S4+Fy))));
T4_Fz = D1*(S1 + D2*(S2 + D3*(S3 + D4*(S4+Fz))));
r4    = [T4_Fx-T4  T4_Fy-T4  T4_Fz-T4];

% Berechnung T3 FRAME
T3_Fx  = D1*(S1 + D2*(S2 + D3*(S3+Fx)));
T3_Fy  = D1*(S1 + D2*(S2 + D3*(S3+Fy)));
T3_Fz  = D1*(S1 + D2*(S2 + D3*(S3+Fz)));
r3     = [T3_Fx-T3  T3_Fy-T3  T3_Fz-T3];

% Berechnung T2 FRAME
T2_Fx  = D1*(S1 + D2*(S2+Fx));
T2_Fy  = D1*(S1 + D2*(S2+Fy));
T2_Fz  = D1*(S1 + D2*(S2+Fz));
r2     = [T2_Fx-T2  T2_Fy-T2  T2_Fz-T2];

% Berechnung T1 FRAME
T1_Fx  = D1*(S1+Fx);
T1_Fy  = D1*(S1+Fy);
T1_Fz  = D1*(S1+Fz);
r1     = [T1_Fx-T1  T1_Fy-T1  T1_Fz-T1];

r = r6;
psiy = atan2(-r(3,1), sqrt(r(1,1)^2 + r(2,1)^2));
psiz = atan2( r(2,1)/cos(psiy), r(1,1)/cos(psiy));
psix = atan2( r(3,2)/cos(psiy), r(3,3)/cos(psiy));

PSI = [psix;  psiy;  psiz];  % gamma/beta/alpha
