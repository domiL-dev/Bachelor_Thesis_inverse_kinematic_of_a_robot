function [Xall Yall Zall] =...
    arrow3(a,b,c,GK,width,height,rcyl,phi)

    L = length(phi);

    X(1,:) = zeros(1,L);
    Y(1,:) = zeros(1,L);
    Z(1,:) = zeros(1,L);

    X0(2,:) = cos(phi)*b(1) + sin(phi)*c(1);
    Y0(2,:) = cos(phi)*b(2) + sin(phi)*c(2);
    Z0(2,:) = cos(phi)*b(3) + sin(phi)*c(3);

    X(2,:) = width*X0(2,:);
    Y(2,:) = width*Y0(2,:);
    Z(2,:) = width*Z0(2,:);

    X(3,:) = height*ones(1,L)*a(1);
    Y(3,:) = height*ones(1,L)*a(2);
    Z(3,:) = height*ones(1,L)*a(3);
 
    % Spitze
    XP = X + GK(1) + a(1) - height*a(1);
    YP = Y + GK(2) + a(2) - height*a(2);
    ZP = Z + GK(3) + a(3) - height*a(3);
   
    p2 = [XP(1,1) YP(1,1) ZP(1,1)];    % Basis Pfeilspitze
    
    % Schaft
    XC(1,:) = rcyl*X0(2,:) + GK(1);
    YC(1,:) = rcyl*Y0(2,:) + GK(2);
    ZC(1,:) = rcyl*Z0(2,:) + GK(3);

    XC(2,:) = rcyl*X0(2,:) + p2(1);
    YC(2,:) = rcyl*Y0(2,:) + p2(2);
    ZC(2,:) = rcyl*Z0(2,:) + p2(3); 
    
%     % Schaft + Spitze als eine Fläche
     Xall = [X(1,:)+ GK(1); XC(1,:); XC(2,:); XP(2,:); XP(3,:)];
     Yall = [Y(1,:)+ GK(2); YC(1,:); YC(2,:); YP(2,:); YP(3,:)];
     Zall = [Z(1,:)+ GK(3); ZC(1,:); ZC(2,:); ZP(2,:); ZP(3,:)];
 
   
