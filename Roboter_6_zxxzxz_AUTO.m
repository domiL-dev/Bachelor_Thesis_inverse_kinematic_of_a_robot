clear all;
clc;

% Graphik__________________________________________________________________
figure(1);
clf;

% Anlegen subplot links ----------------
subplot(1,2,1);
hold on;
grid on;
view(-40,25);

xlabel('x');
ylabel('y');
zlabel('z');

axis equal;
axis([-4 4 -4 4  0 11]);

% Anlegen rechter subplot----------------
subplot(1,2,2);
hold on;
grid on;

axis equal;
axis([-1 1 -1 1 -1 1]*pi);
view(-40,25);

xlabel('psix');
ylabel('psiy');
zlabel('psiz');
%__________________________________________________________________________

i_start = 13;    %  1,...,13 Bahnvorgabe /0 Bahnvorgabe Auto-Mode

t    = 0;       % Startzeitwert Simulation
n    = 0;       % Startwert, Berechnete Werte Array zuordnen
dphi = 1e-8;    % Differenzuotient Jacobi
h    = .0005;   % Schrittweite Euler-Verfahren
tend1= 2*pi;    % Endzeitwert der Simulation Bahnvorgabe 1-13
 
%_________________________________________________________________________
% Pick_Up_Point Koordinaten
 x_pu =3;
 y_pu =3; 
 z_pu =4;  

% Pick_Up_Point Orientierung in Grad von 0° bis 360°
 gamma_pu =20;
 beta_pu  =20;
 alpha_pu =20;
 
% Drop_Off_Point Koordinaten
 x_do = -5;
 y_do = 1;
 z_do = 6;

% Drop_Off_Point Orientierung in Grad von 0° bis 360°
 gamma_do =10;
 beta_do  =20;
 alpha_do =15;
 
%_________________________________________________________________________
%Zwischenrechnung zur Berechnung des Winkels zwischen X-Achse und
%PU-Vektor in Grad
if x_pu~=0
a=y_pu/x_pu;

    phi_z_pu_deg= atand(a);

elseif x_pu==0 && y_pu==0
    phi_z_pu_deg=0;

elseif x_pu==0 && y_pu>0
    phi_z_pu_deg=90;
    
elseif x_pu==0 && y_pu<0
    phi_z_pu_deg=270;

end;

%Berechnung des Winkels zwischen X-Achse und  PU-Vektor im Bogenmaß 
if x_pu > 0 && y_pu >= 0
phi_z_pu = phi_z_pu_deg/180*pi;

elseif x_pu == 0 && y_pu >=0
phi_z_pu = phi_z_pu_deg/180*pi;

elseif x_pu < 0 
phi_z_pu = (180+phi_z_pu_deg)/180*pi;

elseif x_pu >= 0 && y_pu < 0
phi_z_pu = (360+phi_z_pu_deg)/180*pi;

end;

%__________________________________________________________________________
%Zwischenrechnung zur Berechnung des Winkels zwischen X-Achse > 0 und
%DO-Vektor in Grad
if x_do~=0
b=y_do/x_do;
    phi_z_do_deg= atand(b);

elseif x_do==0 && y_do==0
    phi_z_do_deg=0;
    
elseif x_do<0 && y_do==0
phi_z_do_deg=180;

elseif x_do==0 && y_do>0
    phi_z_do_deg=90;
    
elseif x_do==0 && y_do<0
    phi_z_do_deg=270;

end;

%Berechnung des Winkels zwischen X-Achse > 0 und  DO-Vektor im Bogenmaß
if x_do >= 0 && y_do >= 0
phi_z_do = phi_z_do_deg/180*pi;

elseif x_do == 0 && y_do >0
phi_z_do = 90/180*pi;

elseif x_do < 0 
phi_z_do = (180+phi_z_do_deg)/180*pi;

elseif x_do >= 0 && y_do < 0
phi_z_do = (360+phi_z_do_deg)/180*pi;

end;
%__________________________________________________________________________

% Winkeldifferenz zwische PU-Vektor und DO-Vektor
phi_z_dif=(phi_z_do*180/pi) - (phi_z_pu*180/pi);

%__________________________________________________________________________
%PU-Koordinaten -> DO-Koordinaten
%Entscheidung der kürzeren Drehrichtung (dr) um die Ursprungs-Z-Achse 
%Berechnung des Drehwinkels

if phi_z_dif <= 180 && phi_z_dif >= 0 
    dr=1;
    phi_z = phi_z_dif;

elseif phi_z_dif > 180
    dr=-1;
    phi_z = 360-phi_z_dif;
    
elseif phi_z_dif >= -180 && phi_z_dif < 0
    dr=-1;
    phi_z = -phi_z_dif;

elseif phi_z_dif < -180
    dr=1;
    phi_z = 360+phi_z_dif;
    
end;
%__________________________________________________________________________

%__________________________________________________________________________
%Berechnung des Endzeitpunkts der Simulation in Abhängigkeit des
%Drehwinkels zwischen PU- und Do-Koordinaten
if i_start == 0
    
if x_do~=0 || y_do~=0
tend = 0.5*pi+abs(phi_z)/180*pi;

%Spezialfall, x_do=0 und y_do=0, willkürlich eingesetzter Winkel für phi_z
%zur vermeidung, dass tend < t ist und die Simulation zum Zeitpunkt nach 
%anfahren der PU-Koordinate abbricht
%---------------------------------
elseif x_do==0 && y_do==0
    phi_z=90;
    tend = 0.5*pi+phi_z/180*pi;
%---------------------------------
end;

elseif i_start ~= 0
    tend = tend1;
end;
%__________________________________________________________________________

%__________________________________________________________________________
%Projezierte Funktion zum Anfahren der Drop-Off-z-Koordinate, 
%Vereinfacht, Funktion als Ursprungsgerade mit Steigung, sodass beim
%Zeitpunkt tend, die Funktion den Wert delta_z annimmt
delta_z= z_do-z_pu;
delta_t= tend-0.5*pi;
m=delta_z/delta_t;
%__________________________________________________________________________



%%_________________________________________________________________________
 %Berechnung der Winkeldifferenzen und Entscheidung der Drehrichtung 
 %für eine kürzere Rotation um die Achsen des TCPs um die Orientierung 
 %des TCPs nach Vorgabe zu erreichen
 %---------------------------------------------------------------------
 i_PSI=0;
 PSI=[gamma_pu,beta_pu,alpha_pu,gamma_do,beta_do,alpha_do];
 
 %Deklarieren der Variablen zur Berechnung der Winkeldifferenz zwischen PU und DO Orientierung 
 %delta_gamm_pickup,...,delta_alpha_dropoff
 delta_gamma_pu=0; delta_beta_pu=0; delta_alpha_pu=0;
 delta_gamma_do=0; delta_beta_do=0; delta_alpha_do=0;
 delta_PSI=[delta_gamma_pu,delta_beta_pu,delta_alpha_pu,...
               delta_gamma_do,delta_beta_do,delta_alpha_do];
 
 %Deklarieren der Variablen für Drehrichtung von
 %gammarotation_pickup,...,alpharotation_dropoff
 gr_pu=0;br_pu=0;ar_pu=0; gr_do=0; br_do=0; ar_do=0;
 dr_PSI=[gr_pu,br_pu,ar_pu,gr_do,br_do,ar_do];
 
%Ausgangsorienteriung -> PU-Orientierung
%Entscheidung der kürzeren Drehrichtung (dr_PSI) um die jeweiligen Drehachsen 
%Berechnung des Drehwinkels
 while i_PSI < 3
 i_PSI=i_PSI+1;
 
if PSI(i_PSI) <= 180
dr_PSI(i_PSI)=1;
delta_PSI(i_PSI)=PSI(i_PSI)/180*pi;

elseif PSI(i_PSI) > 180
dr_PSI(i_PSI)=-1;
delta_PSI(i_PSI)=(360-PSI(i_PSI))/180*pi;
    
end;
 end;
 
%PU-Orientierung -> DO-Orientierung
%Entscheidung der kürzeren Drehrichtung (dr_PSI) um die jeweiligen Drehachsen 
%Berechnung des Drehwinkels
  while i_PSI < 6
 i_PSI=i_PSI+1;
 i_PSI_dif=i_PSI-3;
 
dif_PSI = PSI(i_PSI)-PSI(i_PSI_dif);

 
if dif_PSI == 0 || abs(dif_PSI)==180
dr_PSI(i_PSI)=1;
delta_PSI(i_PSI)=dif_PSI/180*pi;

elseif dif_PSI < 180 && dif_PSI > 0
dr_PSI(i_PSI)=1;
delta_PSI(i_PSI)=dif_PSI/180*pi;

elseif dif_PSI > 180
dr_PSI(i_PSI)=-1;
delta_PSI(i_PSI)=(360-dif_PSI)/180*pi;

elseif dif_PSI < 0 && dif_PSI > -180
dr_PSI(i_PSI)=-1;
delta_PSI(i_PSI)=-dif_PSI/180*pi;

elseif dif_PSI < -180
dr_PSI(i_PSI)=1;
delta_PSI(i_PSI)=(360+dif_PSI)/180*pi;   
end;
 end;
%__________________________________________________________________________

%__________________________________________________________________________
%Wichtige Variablen und Komponenten für die Funktionen im Auto-Mode

%Radius der PU-Koordinate in der X-Y-Ebene mit Nullpunkt im Ursprung
 r_pu=sqrt((x_pu^2)+(y_pu^2));
 
%Drehmatrize um PU-Vektor nach DO-Vektor-Ausrichtung in der X-Y-Ebene
 D = [cos(phi_z/180*pi)     -dr*sin(phi_z/180*pi) 
     dr*sin(phi_z/180*pi)  cos(phi_z/180*pi)];
  
%Rotation des PU-Vektors in DO-Vektor-Ausrichtung
do1=D*[x_pu;y_pu];
 
%Berechnug des restlichen delta_y und delta_x um DO-Koordinate zu
%erreichen
 x_do1=do1(1);
 y_do1=do1(2);
 
%Funktion zum Anfahren der Drop-Off-Koordinate in X-Y-Ebene, 
%Vereinfacht, Funktion als Ursprungsgerade mit Steigung, sodass beim
%Zeitpunkt tend, die Funktionen die Werte von delta_x und delta_y annehmen
 vx=(x_do-x_do1)/delta_t;
 vy=(y_do-y_do1)/delta_t;
 
%PU-Koordinaten und -Orientierung in Arrays abgelegt
 TCP_1= [x_pu;y_pu;z_pu];
 PSI_1 = [gamma_pu/180*pi;beta_pu/180*pi;alpha_pu/180*pi];
 

% Ausgangskonfiguration in 0
S0 = [0;0;0]; 
S1 = [0;0;2]; 
S2 = [0;0;3]; 
S3 = [0;0;1]; 
S4 = [0;0;2]; 
S5 = [0;0;1]; 
S6 = [0;0;1];   

Fx = [1;0;0];
Fy = [0;1;0];
Fz = [0;0;1];


% (phi1,...,phi6) Anfangsgelenkwinkel --> (TCP,PSI)für t=0
alpha = 45/180*pi;   % Anfangsknickwinkel : 0°, 90° singulär 
if i_start==0
    phi1= -phi_z_pu;
    
elseif i_start==7 || i_start==13
    phi1= alpha;
    
elseif i_start==1 || i_start==9
    phi1= 2*alpha;
    
elseif i_start==13
    phi1= 45/180*pi;
    
elseif i_start~=0
phi1 =  0;
end;

phi2 =  -alpha;
phi3 =  2*alpha;
phi4 = 0;
phi5 =  -alpha;

if i_start==0
    phi6= phi_z_pu;

elseif i_start==9
phi6 =  -2*alpha;

elseif i_start~=0
phi6 =  0;
end;
phi  = [phi1; phi2; phi3; phi4; phi5; phi6];

%(TCP,PSI)für t=0
[TCP_0,PSI_0, T5,T4,T3,T2,T1, r6,r5,r4,r3,r2,r1] = TCP_PSI_zxxzxz(S1,S2,...
S3,S4,S5,S6, Fx,Fy,Fz, phi1,phi2,phi3,phi4,phi5,phi6);

%Ergänzender Wert zur Berechnung der Trajektorie im Auto-Mode
%Delta der Z-Ausgangskoordinate und der Z-PU-Koordinate
 rz=TCP_0(3)-z_pu;

%__________________________________________________________________________
%Simulation
while t <= tend 
    n = n + 1;
    
 %_________________________________________________________________________
     if i_start==0          % Auto-Mode -----------------------------------
  
        if t < 0.5*pi       %Statement TCP-/PSI-PU
 
        TCP_soll     = TCP_0 + [sin(t)*x_pu; sin(t)*y_pu; (cos(t)*rz-rz)];
        TCP_soll_dot =         [cos(t)*x_pu; cos(t)*y_pu; -sin(t)*rz];
        
        PSI_soll     = PSI_0 + [dr_PSI(1)*delta_PSI(1)/(0.5*pi)*t; dr_PSI(2)*delta_PSI(2)/(0.5*pi)*t; dr_PSI(3)*delta_PSI(3)/(0.5*pi)*t]; 
        PSI_soll_dot =         [dr_PSI(1)*delta_PSI(1)/(0.5*pi)  ; dr_PSI(2)*delta_PSI(2)/(0.5*pi)  ; dr_PSI(3)*delta_PSI(3)/(0.5*pi)];
      
        
        elseif t >= 0.5*pi  %Statement TCP-/PSI-DO
       
        TCP_soll     =    [cos(dr*(t-0.5*pi)+phi_z_pu)*r_pu+vx*(t-0.5*pi); sin(phi_z_pu+dr*(t-0.5*pi))*r_pu+vy*(t-0.5*pi); TCP_1(3)+ m*(t-0.5*pi)];
        TCP_soll_dot =    [-dr*r_pu*sin(dr*(t-pi/2)+phi_z_pu)+vx         ; dr*r_pu*cos(dr*(t-pi/2)+phi_z_pu)+vy          ; m];
                        
        PSI_soll     =  PSI_1 + [dr_PSI(4)*delta_PSI(4)/(tend-0.5*pi)*(t-0.5*pi); dr_PSI(5)*delta_PSI(5)/(tend-0.5*pi)*(t-0.5*pi); ...
                                dr_PSI(6)*delta_PSI(6)/(tend-0.5*pi)*(t-0.5*pi)];
        PSI_soll_dot =          [dr_PSI(4)*delta_PSI(4)/(tend-0.5*pi)          ; dr_PSI(5)*delta_PSI(5)/(tend-0.5*pi); ...
                                 dr_PSI(6)*delta_PSI(6)/(tend-0.5*pi)];
        end;
  %________________________________________________________________________
          
    %------------------------------------------------
    elseif i_start==1      % Kellner x-Achse -------
     TCP_soll     = TCP_0 + [2*sin(t); 0; 0];
     TCP_soll_dot =         [2*cos(t); 0; 0];
     PSI_soll     = PSI_0 + [0; 0; 0];
     PSI_soll_dot =         [0; 0; 0];

    elseif i_start==2   % Kellner y-Achse -----
     TCP_soll     = TCP_0 + [0; 2*sin(t); 0];
     TCP_soll_dot =         [0; 2*cos(t); 0];
     PSI_soll     = PSI_0 + [0; 0; 0];
     PSI_soll_dot =         [0; 0; 0];

    elseif i_start==3   % Kellner z-Achse -----
     TCP_soll     = TCP_0 + [0; 0; (10-TCP_0(3))*sin(t)];
     TCP_soll_dot =         [0; 0; (10-TCP_0(3))*cos(t)];
     PSI_soll     = PSI_0 + [0; 0; 0];
     PSI_soll_dot =         [0; 0; 0];

    elseif i_start==4   % Kellner Kreis (x,y)-Ebene ----
     TCP_soll     = TCP_0 + 1.8*[ cos(t)-1; sin(t); 0];
     TCP_soll_dot =         1.8*[-sin(t);   cos(t); 0];
     PSI_soll     = PSI_0 + [0; 0; 0];
     PSI_soll_dot =         [0; 0; 0];

    elseif i_start==5   % Kellner Kreis (y,z)-Ebene ----
     TCP_soll     = TCP_0 + 1.5*[0;  cos(t)-1; sin(t)];
     TCP_soll_dot =         1.5*[0; -sin(t);   cos(t)];
     PSI_soll     = PSI_0 + [0; 0; 0];
     PSI_soll_dot =         [0; 0; 0];

    elseif i_start==6   % Kellner Kreis (x,z)-Ebene ----
     TCP_soll     = TCP_0 + 1.5*[ cos(t)-1; 0; sin(t)];
     TCP_soll_dot =         1.5*[-sin(t);   0; cos(t)];
     PSI_soll     = PSI_0 + [0; 0; 0];
     PSI_soll_dot =         [0; 0; 0];

    elseif i_start==7   % Kellner(x,y,z)-Lissajouis-Figur ----
     TCP_soll     = TCP_0 + [-sin(3*t);   3*sin(t);    sin(2*t)];
     TCP_soll_dot =         [-cos(3*t)*3; 3*cos(t);   cos(2*t)*2];
     PSI_soll     = PSI_0 + [0; 0; 0];
     PSI_soll_dot =         [0; 0; 0];

    elseif i_start==8   % Schaukelbewegung um x-Achse ----
      TCP_soll     = TCP_0 + [0; 0; 0];
      TCP_soll_dot =         [0; 0; 0];
      PSI_soll     = PSI_0 + [sin(3*t); 0; 0];
      PSI_soll_dot =         [cos(3*t)*3; 0; 0];

    elseif i_start==9   % Schaukelbewegung um y-Achse ----
      TCP_soll     = TCP_0 + [0; 0; 0];
      TCP_soll_dot =         [0; 0; 0];
      PSI_soll     = PSI_0 + [0; sin(3*t); 0];
      PSI_soll_dot =         [0; cos(3*t)*3; 0];

    elseif i_start==10   % Rotation um z-Achse ---
      TCP_soll     = TCP_0 + [0; 0; 0];
      TCP_soll_dot =         [0; 0; 0];
      PSI_soll     = PSI_0 + [0; 0; 3*sin(2*t)];
      PSI_soll_dot =         [0; 0; 3*cos(2*t)*2];

    elseif i_start==11   % Rotation/Schwenkbewegung um xz-Achse ----
      TCP_soll     = TCP_0 + [0; 0; 0];
      TCP_soll_dot =         [0; 0; 0];
      PSI_soll     = PSI_0 + [sin(t);   0; 3*sin(0.25*t)];
      PSI_soll_dot =         [cos(t); 0; 3*cos(0.25*t)*0.25];

    elseif i_start==12   % Senkrecht auf Kreis (x,y)-Ebene ----
      TCP_soll     = TCP_0 + [-cos(3*t)+1; sin(3*t);0];
      TCP_soll_dot =         [sin(3*t)*3;  cos(3*t)*3;0];
      PSI_soll     = PSI_0 + [0; 0; -3*t];
      PSI_soll_dot =         [0; 0; -3];

    elseif i_start==13   % -Lissajouis-Figur, Schwenkbewegung um xz ----
      TCP_soll     = TCP_0 + [-sin(3*t);   3*sin(t);    sin(2*t)];
      TCP_soll_dot =         [-cos(3*t)*3; 3*cos(t);   cos(2*t)*2];
      PSI_soll     = PSI_0 + [0.3*sin(2*t); 0; 2*sin(t) ]; 
      PSI_soll_dot =         [0.3*cos(2*t)*2; 0; 2*cos(t) ]; 
    end;

    % Berechnung Jacobi----------------------------------------------------
    
    [TCP,PSI, T5,T4,T3,T2,T1, r6,r5,r4,r3,r2,r1] = ...
TCP_PSI_zxxzxz(S1,S2,S3,S4,S5,S6, Fx,Fy,Fz, phi1,phi2,phi3,phi4,phi5,phi6);
    
        %Jacobi mit Diff.quot. 
        [TCP1,PSI1, T51,T41,T31,T21,T11, r61,r51,r41,r31,r21,r11] = TCP_PSI_zxxzxz(S1,S2,S3,S4,S5,S6, Fx,Fy,Fz, phi1+dphi,phi2,     phi3,     phi4,     phi5,     phi6);
        [TCP2,PSI2, T52,T42,T32,T22,T12, r62,r52,r42,r32,r22,r12] = TCP_PSI_zxxzxz(S1,S2,S3,S4,S5,S6, Fx,Fy,Fz, phi1,     phi2+dphi,phi3,     phi4,     phi5,     phi6);
        [TCP3,PSI3, T53,T43,T33,T23,T13, r63,r53,r43,r33,r23,r13] = TCP_PSI_zxxzxz(S1,S2,S3,S4,S5,S6, Fx,Fy,Fz, phi1,     phi2,     phi3+dphi,phi4,     phi5,     phi6);
        [TCP4,PSI4, T54,T44,T34,T24,T14, r64,r54,r44,r34,r24,r14] = TCP_PSI_zxxzxz(S1,S2,S3,S4,S5,S6, Fx,Fy,Fz, phi1,     phi2,     phi3,     phi4+dphi,phi5,     phi6);
        [TCP5,PSI5, T55,T45,T35,T25,T15, r65,r55,r45,r35,r25,r15] = TCP_PSI_zxxzxz(S1,S2,S3,S4,S5,S6, Fx,Fy,Fz, phi1,     phi2,     phi3,     phi4,     phi5+dphi,phi6);
        [TCP6,PSI6, T56,T46,T36,T26,T16, r66,r56,r46,r36,r26,r16] = TCP_PSI_zxxzxz(S1,S2,S3,S4,S5,S6, Fx,Fy,Fz, phi1,     phi2,     phi3,     phi4,     phi5,     phi6+dphi);

        % Diff. quot Jacobi Position TCP
        f1 = (TCP1 - TCP)/dphi; %partielle Ableitung nach phi1
        f2 = (TCP2 - TCP)/dphi; %partielle Ableitung nach phi2
        f3 = (TCP3 - TCP)/dphi; %partielle Ableitung nach phi3
        f4 = (TCP4 - TCP)/dphi; %partielle Ableitung nach phi4
        f5 = (TCP5 - TCP)/dphi; %partielle Ableitung nach phi5
        f6 = (TCP6 - TCP)/dphi; %partielle Ableitung nach phi6

        % Diff. quot Jacobi Orientierung TCP
        a1 = (PSI1 - PSI)/dphi; %partielle Ableitung nach phi1
        a2 = (PSI2 - PSI)/dphi; %partielle Ableitung nach phi2
        a3 = (PSI3 - PSI)/dphi; %partielle Ableitung nach phi3
        a4 = (PSI4 - PSI)/dphi; %partielle Ableitung nach phi4
        a5 = (PSI5 - PSI)/dphi; %partielle Ableitung nach phi5
        a6 = (PSI6 - PSI)/dphi; %partielle Ableitung nach phi6

  %Jacobi Matrix
    f = [f1  f2  f3  f4  f5  f6
         a1  a2  a3  a4  a5  a6];
   
    t       = t + h;
    phi_dot = inv(f) * [TCP_soll_dot; PSI_soll_dot];
    phi     = phi + h  * phi_dot;

    phi1 = phi(1);
    phi2 = phi(2);
    phi3 = phi(3);
    phi4 = phi(4);
    phi5 = phi(5);
    phi6 = phi(6);
    
    TCP_kurve_x_soll(n) = TCP_soll(1);
    TCP_kurve_y_soll(n) = TCP_soll(2);
    TCP_kurve_z_soll(n) = TCP_soll(3);   
    
    TCP_kurve_x_ist(n)  = TCP(1);
    TCP_kurve_y_ist(n)  = TCP(2);
    TCP_kurve_z_ist(n)  = TCP(3);
    
    PSI_kurve_x_soll(n) = PSI_soll(1);
    PSI_kurve_y_soll(n) = PSI_soll(2);
    PSI_kurve_z_soll(n) = PSI_soll(3);
    
    PSI_kurve_x_ist(n)  = PSI(1);
    PSI_kurve_y_ist(n)  = PSI(2);
    PSI_kurve_z_ist(n)  = PSI(3);
    
    
    if mod(n,50)==0 || t>=tend-5*h
 
        % (x,y,z)-Raum : subplot links -----------
        subplot(1,2,1);
        cla;
        
        % TCP Sollkurve
        plot3(TCP_kurve_x_soll,TCP_kurve_y_soll,TCP_kurve_z_soll,'-m',...
            'linesmoothing','off','linewidth',1);            
  
        % TCP Istkurve
        plot3(TCP_kurve_x_ist, TCP_kurve_y_ist, TCP_kurve_z_ist,'-b');   
        plot3(TCP(1),TCP(2),TCP(3),'ok','markerfacecolor','y','markersize',8);       
                
        % TCP FRAME_SOLL nach Vorgabe PSI aus Trajektorie      
        wx = PSI_soll(1);  wy = PSI_soll(2);  wz = PSI_soll(3);
        
        Dz = [cos(wz) -sin(wz)   0
              sin(wz)  cos(wz)   0
              0            0     1];

        Dy = [cos(wy)  0   sin(wy)
              0        1   0
             -sin(wy)  0   cos(wy)];

        Dx = [1          0     0
              0     cos(wx) -sin(wx)
              0     sin(wx)  cos(wx)];
         
        rs = Dz*Dy*Dx*eye(3);
      
        %Visualisierung der Vorgabe der Orientierung des TCP-KS
        TCP_Fx_soll = TCP + rs(:,1)*2;   
        TCP_Fy_soll = TCP + rs(:,2)*2;   
        TCP_Fz_soll = TCP + rs(:,3)*2;
        plot3([TCP(1) TCP_Fx_soll(1)],[TCP(2) TCP_Fx_soll(2)],[TCP(3) TCP_Fx_soll(3)],'-r','linewidth',2);
        plot3([TCP(1) TCP_Fy_soll(1)],[TCP(2) TCP_Fy_soll(2)],[TCP(3) TCP_Fy_soll(3)],'-g','linewidth',2);
        plot3([TCP(1) TCP_Fz_soll(1)],[TCP(2) TCP_Fz_soll(2)],[TCP(3) TCP_Fz_soll(3)],'-b','linewidth',2);
        
        
        % 3D Visualisierung des
        % 6-Achsroboters--------------------------------------------------
       
        %Daten zur Visualisierung der Glieder und Gelenke
       
        % Würfel ------------------------------------------------------
    
    ecken =  [0 0 0
              1 0 0
              1 1 0
              0 1 0
              0 0 1
              1 0 1
              1 1 1
              0 1 1];

    ecken(:,1) = (ecken(:,1) - 0.5)*0.3;
    ecken(:,2) = (ecken(:,2) - 0.5)*0.3;

    flaechen = [1 2 6 5
                2 3 7 6
                3 4 8 7
                4 1 5 8];
             
        %Daten zur Orientierung der Glieder
        
        % Glied 1 (S1)
        W_  = [ecken(:,1).'; ecken(:,2).'; S1(3)*ecken(:,3).'];
        W_1 = r1*W_;
        ecken1G = [W_1(1,:).'  W_1(2,:).'  W_1(3,:).'];
 
        % Glied 2 (S2)
        W_  = [ecken(:,1).'; ecken(:,2).'; S2(3)*ecken(:,3).'];
        W_2 = r2*W_;
        ecken2G = [W_2(1,:).'+T1(1) W_2(2,:).'+T1(2) W_2(3,:).'+T1(3)];

        % Glied 3 (S3)
        W_  = [ecken(:,1)'; ecken(:,2)'; S3(3)*ecken(:,3)'];
        W_3 = r3*W_;
        ecken3G = [W_3(1,:)'+T2(1) W_3(2,:)'+T2(2) W_3(3,:)'+T2(3)];
        
        % Glied 4 (S4)
        W_  = [ecken(:,1)'; ecken(:,2)'; S4(3)*ecken(:,3)'];
        W_4 = r4*W_;
        ecken4G = [W_4(1,:)'+T3(1) W_4(2,:)'+T3(2) W_4(3,:)'+T3(3)];
        
        % Glied 5 (S5)
        W_  = [ecken(:,1)'; ecken(:,2)'; S5(3)*ecken(:,3)'];
        W_5 = r5*W_;
        ecken5G = [W_5(1,:)'+T4(1) W_5(2,:)'+T4(2) W_5(3,:)'+T4(3)];
        
        % Glied 6 (S)
        W_  = [ecken(:,1)'; ecken(:,2)'; S6(3)*ecken(:,3)'];
        W_6 = r6*W_;
        ecken6G = [W_6(1,:)'+T5(1) W_6(2,:)'+T5(2) W_6(3,:)'+T5(3)];


        %Visualisierung der Glieder
         patch('Vertices',ecken1G,'Faces',flaechen,...
                'FaceColor','k','EdgeColor','k',...
                'Marker','none','MarkerFaceColor','y',...
                'linewidth',1,'ambientstrength',1,'FaceAlpha',0.3);
            
         patch('Vertices',ecken2G,'Faces',flaechen,...
                'FaceColor','r','EdgeColor','k',...
                'Marker','none','MarkerFaceColor','y',...
                'linewidth',0.8,'ambientstrength',1,'FaceAlpha',0.9);
        
        patch('Vertices',ecken3G,'Faces',flaechen,...
                'FaceColor','r','EdgeColor','k',...
                'Marker','none','MarkerFaceColor','y',...
                'linewidth',0.8,'ambientstrength',1,'FaceAlpha',0.9);
            
         patch('Vertices',ecken4G,'Faces',flaechen,...
                'FaceColor','k','EdgeColor','k',...
                'Marker','none','MarkerFaceColor','y',...
                'linewidth',1,'ambientstrength',1,'FaceAlpha',0.3);
        
         patch('Vertices',ecken5G,'Faces',flaechen,...
                'FaceColor','r','EdgeColor','k',...
                'Marker','none','MarkerFaceColor','y',...
                'linewidth',0.8,'ambientstrength',1,'FaceAlpha',0.9);
            
         patch('Vertices',ecken6G,'Faces',flaechen,...
                'FaceColor','k','EdgeColor','k',...
                'Marker','none','MarkerFaceColor','y',...
                'linewidth',1,'ambientstrength',1,'FaceAlpha',.3);
        
            
        %Daten für Visualisierung der Gelenke als Kugeln
             % Aufbau Kugel
        NK          = 6;
       [XK,YK,ZK]  = sphere(NK);
        r           = 0.3;
        % Kugel 1 ---------------------------------------------------------
        for i=1:NK+1
            XZ        = r*[XK(i,:); YK(i,:); ZK(i,:)];
            K         = r1*XZ;
            XK_1(i,:) = K(1,:);
            YK_1(i,:) = K(2,:);
            ZK_1(i,:) = K(3,:);
        end;
 
        % Kugel 2
        for i=1:NK+1
            XZ        = r*[XK(i,:); YK(i,:); ZK(i,:)];
            K         = r2*XZ;
            XK_2(i,:) = K(1,:)+T1(1);
            YK_2(i,:) = K(2,:)+T1(2);
            ZK_2(i,:) = K(3,:)+T1(3);
        end;
     
        % Kugel 3
        for i=1:NK+1
            XZ        = r*[XK(i,:); YK(i,:); ZK(i,:)];
            K         = r3*XZ;
            XK_3(i,:) = K(1,:)+T2(1);
            YK_3(i,:) = K(2,:)+T2(2);
            ZK_3(i,:) = K(3,:)+T2(3);
        end;

        % Kugel 4
        for i=1:NK+1
            XZ        = r*[XK(i,:); YK(i,:); ZK(i,:)];
            K         = r4*XZ;
            XK_4(i,:) = K(1,:)+T3(1);
            YK_4(i,:) = K(2,:)+T3(2);
            ZK_4(i,:) = K(3,:)+T3(3);
        end;

        % Kugel 5
        for i=1:NK+1
            XZ        = r*[XK(i,:); YK(i,:); ZK(i,:)];
            K         = r5*XZ;
            XK_5(i,:) = K(1,:)+T4(1);
            YK_5(i,:) = K(2,:)+T4(2);
            ZK_5(i,:) = K(3,:)+T4(3);
        end;

        % Kugel 6
        for i=1:NK+1
            XZ        = r*[XK(i,:); YK(i,:); ZK(i,:)];
            K         = r6*XZ;
            XK_6(i,:) = K(1,:)+T5(1);
            YK_6(i,:) = K(2,:)+T5(2);
            ZK_6(i,:) = K(3,:)+T5(3);
        end;
        
         %Visualisierung Gelenke als Kugeln
         %Gelenk 1   
            surf('xdata',XK_1,'ydata',YK_1,'zdata',ZK_1,...
                'facealpha',1,'facecolor','k','visible','on',...
                'facelighting','gouraud','edgecolor','r',...
                'ambientstrength',.1,'linewidth',1);
         %Gelenk 2
            surf('xdata',XK_2,'ydata',YK_2,'zdata',ZK_2,...
                'facealpha',1,'facecolor','k','visible','on',...
                'facelighting','gouraud','edgecolor','r',...
                'ambientstrength',.1,'linewidth',1);
         %Gelenk 3
            surf('xdata',XK_3,'ydata',YK_3,'zdata',ZK_3,...
                'facealpha',1,'facecolor','k','visible','on',...
                'facelighting','gouraud','edgecolor','r',...
                'ambientstrength',.1,'linewidth',1);
         %Gelenk 4   
            surf('xdata',XK_4,'ydata',YK_4,'zdata',ZK_4,...
                'facealpha',1,'facecolor','k','visible','on',...
                'facelighting','gouraud','edgecolor','r',...
                'ambientstrength',.1,'linewidth',1);
         %Gelenk 5 
            surf('xdata',XK_5,'ydata',YK_5,'zdata',ZK_5,...
                'facealpha',1,'facecolor','k','visible','on',...
                'facelighting','gouraud','edgecolor','r',...
                'ambientstrength',.1,'linewidth',1);
         %Gelenk 6
            surf('xdata',XK_6,'ydata',YK_6,'zdata',ZK_6,...
                'facealpha',1,'facecolor','k','visible','on',...
                'facelighting','gouraud','edgecolor','r',...
                'ambientstrength',.1,'linewidth',1);
           
            %Visualisierung Gelenkachsen durch Pfeile
           
            sr= linspace(0,2*pi,8); %Sektoren Rotationsobjekt
            b_P = .13;              %breite Pfeilspitze
            h_P = 0.3;              %höhe Pfeilspitze
            r_Z = .05;              %Radius Schaft
            
            %Gelenkachse 1
            [Xall Yall Zall] =...
                arrow3(r1(:,3),r1(:,1),r1(:,2),S0,b_P,h_P,r_Z,sr);
            
            surf('xdata',Xall,'ydata',Yall,'zdata',Zall,...
                'facealpha',1,'facecolor','y','visible','on',...
                'facelighting','gouraud','edgecolor','k',...
                'ambientstrength',.1,'linewidth',1);
            
            %Gelenkachse 2
            [Xall Yall Zall] =...
                arrow3(r2(:,1),r2(:,2),r2(:,3),T1,b_P,h_P,r_Z,sr);

            surf('xdata',Xall,'ydata',Yall,'zdata',Zall,...
                'facealpha',1,'facecolor','y','visible','on',...
                'facelighting','gouraud','edgecolor','k',...
                'ambientstrength',.1,'linewidth',1);
            
            %Gelenkachse 3
            [Xall Yall Zall] =...
                arrow3(r3(:,1),r3(:,2),r3(:,3),T2,b_P,h_P,r_Z,sr); 
            
            surf('xdata',Xall,'ydata',Yall,'zdata',Zall,...
                'facealpha',1,'facecolor','y','visible','on',...
                'facelighting','gouraud','edgecolor','k',...
                'ambientstrength',.1,'linewidth',1);
            
            %Gelenkachse 4
            [Xall Yall Zall] =...
                arrow3(r4(:,3),r4(:,1),r4(:,2),T3,b_P,h_P,r_Z,sr);
            
            surf('xdata',Xall,'ydata',Yall,'zdata',Zall,...
                'facealpha',1,'facecolor','y','visible','on',...
                'facelighting','gouraud','edgecolor','k',...
                'ambientstrength',.1,'linewidth',1);
            
            %Gelenkachse 5
            [Xall Yall Zall] =...
                arrow3(r5(:,1),r5(:,2),r5(:,3),T4,b_P,h_P,r_Z,sr);

            surf('xdata',Xall,'ydata',Yall,'zdata',Zall,...
                'facealpha',1,'facecolor','y','visible','on',...
                'facelighting','gouraud','edgecolor','k',...
                'ambientstrength',.1,'linewidth',1);
            
            %Gelenkachse 6
            [Xall Yall Zall] =...
                arrow3(r6(:,3),r6(:,1),r6(:,2),T5,b_P,h_P,r_Z,sr);
            
            surf('xdata',Xall,'ydata',Yall,'zdata',Zall,...
                'facealpha',1,'facecolor','y','visible','on',...
                'facelighting','gouraud','edgecolor','k',...
                'ambientstrength',.1,'linewidth',1);
            
            %Visualisierung TCP-KS
            b_P = .17;
            h_P = .4;
            r_Z = .08;
            
            %TCP-X-Achse
            [Xall Yall Zall] =...
                arrow3(r6(:,1)*1.5,r6(:,2),r6(:,3),TCP,b_P,h_P,r_Z,sr);

            surf('xdata',Xall,'ydata',Yall,'zdata',Zall,...
                'facealpha',1,'facecolor','r','visible','on',...
                'facelighting','gouraud','edgecolor','k',...
                'ambientstrength',.1,'linewidth',0.3);
            
            %TCP-Y-Achse
            [Xall Yall Zall] =...
                arrow3(r6(:,2)*1.5,r6(:,3),r6(:,1),TCP,b_P,h_P,r_Z,sr);

            surf('xdata',Xall,'ydata',Yall,'zdata',Zall,...
                'facealpha',1,'facecolor','g','visible','on',...
                'facelighting','gouraud','edgecolor','k',...
                'ambientstrength',.1,'linewidth',0.3);
            
            %TCP-Z-Achse
           [Xall Yall Zall] =...
                arrow3(r6(:,3)*1.5,r6(:,1),r6(:,2),TCP,b_P,h_P,r_Z,sr);

            surf('xdata',Xall,'ydata',Yall,'zdata',Zall,...
                'facealpha',1,'facecolor','b','visible','on',...
                'facelighting','gouraud','edgecolor','k',...
                'ambientstrength',.1,'linewidth',0.3);
       
            
        % PSI-Raum : subplot rechts -----------
        subplot(1,2,2);
        cla;
        
        % PSI Sollkurve
        plot3(PSI_kurve_x_soll,PSI_kurve_y_soll,PSI_kurve_z_soll,'-m','linesmoothing','off');
        plot3(PSI_kurve_x_soll(n),PSI_kurve_y_soll(n),PSI_kurve_z_soll(n),'-om','linesmoothing','off');

        % PSI Istkurve
        plot3(PSI_kurve_x_ist,PSI_kurve_y_ist,PSI_kurve_z_ist,'-b','linesmoothing','off');
        plot3(PSI_kurve_x_ist(n),PSI_kurve_y_ist(n),PSI_kurve_z_ist(n),'-ob','linesmoothing','off');
 
        drawnow;
   
    end;   % Graphik
end;
