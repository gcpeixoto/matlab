%% 3D Finite Element Cubic Interpolation (Residual Approach / TE16) 
%
% Author: Gustavo Charles P. de Oliveira
% Date: November, 2010
% Updated: October, 2013
%
% Description: Code to develop a cubic interpolation for finite
% element on 3D-space. Basically, we fit an arbitrary function over pattern
% element under a residual approach. We have been investigated a new way to
% conceive a mix between the shape functions associated to gradients on
% TE16 3-D element. (See reports nos. AM290 and AM297 by M. Wang and J. Xu from
% PennState University, Department of Mathematics)
%
%
% Consider the following numbering for each corner point:
%    Pt | x | y | z
%     1 | 1 | 0 | 0
%     2 | 0 | 1 | 0
%     3 | 0 | 0 | 1
%     4 | 0 | 0 | 0
%
%% Defaults
clear all;
format long;
%figure 

%% Splash Screen 

% Splash screen
disp(' _________________________________________ ');
disp('|                                         |');
disp('|*** FINITE ELEMENT INTERPOLATION CODE ***|');
disp('|_________________________________________|');
disp('Description: Code to carry on a FEM Z-Type'); 
disp('interpolation over the reference tetrahedron');
disp('finite element.');
disp(' ');
disp('Authors: Gustavo Charles P. de Oliveira');
disp('         Norberto Mangiavacchi         ');

%% User Input Commands

% It starts the profiler to exhibit computational data about the code
Profile_Power_On_Off = input('Profile function may be actived to show computational data after running the code.\n Do you desire active profile? [y] for Yes / [n] for No \n','s');

if Profile_Power_On_Off == 'y'
    sprintf('%s','Profile ON...')
    profile on
elseif Profile_Power_On_Off == 'n'
    sprintf('%s','Profile OFF...')
    profile off
else 
    sprintf('%s','Value is not valid! Closing...')
    break
end

% It allows choose between two ways to calculate the tetrahedron volume
% coordinates:
Input_Volume_Coords_Choice = input('Which version to calculate Tetrahedron Volume Coordinates? [1] for New /[2] for Old \n');

if (Input_Volume_Coords_Choice ~= 1 && Input_Volume_Coords_Choice ~= 2)
    sprintf('%s','Value is not valid! Closing...')
    break
end
  
% It decides on an element translation
Choose_Translation = input('Do you desire to translate the element? [y] for Yes / [n] for No \n','s');

if Choose_Translation == 'y'
    Translation_Vector = input('Choose a 3D vector to translate: [X-value Y-value Z-value] \n');
elseif Choose_Translation == 'n'
    sprintf('%s','Setting translation vector to [0 0 0]...')
    Translation_Vector = [0,0,0];
else
    sprintf('%s','Value is not valid! Closing...')
    break
end

% It decides on an element rotation
Choose_Rotation = input('Do you desire to rotate the element? [y] for Yes / [n] for No \n','s');
if Choose_Rotation == 'y'
    Rotation_Vector = input('Choose a 3D vector for Euler proper angles: [Phi-angle Theta-angle Psi-angle]. \n Valid intervals: 0 <= Theta < pi; 0 <= Phi, Psi < 2*pi \n');
else
    sprintf('%s','No rotation! Still running...')
end


% Choose which type of function to interpolate
Function_Type = input('Choose function to interpolate: \n[p] for a plane; \n[q] for a quadratic; \n[c] for a cubic; \n[s] for a senoidal-cossenoidal \n[g] for a general curve. \n ','s');

if Function_Type == 'p'
    sprintf('%s','Interpolating plane...')
elseif Function_Type == 'q'
    sprintf('%s','Interpolating quadratic...')
elseif Function_Type == 'c'
    sprintf('%s','Interpolating cubic...')
elseif Function_Type == 's'
    sprintf('%s','Interpolating senoidal-cossenoidal...')
elseif Function_Type == 'g'
    sprintf('%s','Interpolating general curve...')
else
    sprintf('%s','Value is not valid! Closing...')
    break
end

% % Choose grid discretization
% Grid_Discretization = input('Define the 3D grid discretization vector: [dx dy dz]. \n Numbers must be within the range (0,1), but mind that small spacing can be very expensive. \n');
% 
% if Grid_Discretization <= 0 || Grid_Discretization >= 1 
%     sprintf('%s','Value is not within the range! Closing...')
%     break
% else
%     sprintf('%s','Setting grid discretization...')
% end
%% Grid

dx = 0.1; dy = 0.1; dz = 0.1; 

xx = 0:dx:1;
yy = 0:dy:1;
zz = 0:dz:1;


% xx = 0:Grid_Discretization(1):1;
% yy = 0:Grid_Discretization(2):1;
% zz = 0:Grid_Discretization(3):1;

[XX YY ZZ] = meshgrid(xx,yy,zz); % interior points

XX=XX.*(1 - YY); % projection onto canonical parallel planes 
YY=YY.*(1 - ZZ);
ZZ=ZZ.*(1 - XX);


% translation test
XX = XX + Translation_Vector(1); 
YY = YY + Translation_Vector(2);
ZZ = ZZ + Translation_Vector(3);

if Choose_Rotation == 'y'

    phi = Rotation_Vector(1);
    theta = Rotation_Vector(2);
    psi = Rotation_Vector(3);

    SP = sin(psi);
    CP = cos(psi);

    ST = sin(theta);
    CT = cos(theta);

    SF = sin(phi);
    CF = cos(phi);

    M = [CT*CP (- CF*SP + SF*ST*CP) (SF*SP + CF*ST*CP); CT*SP (CF*CP + SF*ST*SP) (- SF*CP + CF*ST*SP); -ST SF*CT CT*CF]; 

    % rotating
    XT = M(1,1)*XX + M(1,2)*YY + M(1,3)*ZZ;
    YT = M(2,1)*XX + M(2,2)*YY + M(2,3)*ZZ;
    ZT = M(3,1)*XX + M(3,2)*YY + M(3,3)*ZZ;

    XX = XT;
    YY = YT;
    ZZ = ZT;
end
%% Plot Cube Faces (to see the tetrahedron)
hold on 
%  plot3(XX(:,:,1),YY(:,:,1),ZZ(:,:,1),'ro') % (x,y,z = z_1)
%  plot3(XX(:,:,end),YY(:,:,end),ZZ(:,:,end),'bo') % (x,y,z = z_end)
%  plot3(YY(:,:,1),ZZ(:,:,1),XX(:,:,1),'yo') % (x,z,y = y_1)
%  plot3(YY(:,:,end),ZZ(:,:,end),XX(:,:,end),'mo') % (x,z,y = y_end)
%  plot3(ZZ(:,:,1),XX(:,:,1),YY(:,:,1),'ko') % (y,z,x = x_1)
%  plot3(ZZ(:,:,end),XX(:,:,end),YY(:,:,end),'go') % (y,z,x = x_end)
%  plot3(XX(:,:,1),YY(:,:,1),1 - XX(:,:,1) - YY(:,:,1),'c*') % diagonal plane

%% INDICES X and Y ARE INVERTED!!!

% Proof:
plot3(XX(1,end,1),YY(1,end,1),ZZ(1,end,1),'^k','MarkerSize',10) % Wrong! (1,0,0) and not (0,1,0)
plot3(XX(end,1,1),YY(end,1,1),ZZ(end,1,1),'^k','MarkerSize',10) % Wrong! This one is (0,1,0) and not (1,0,0) strangely...
plot3(XX(1,1,end),YY(1,1,end),ZZ(1,1,end),'^k','MarkerSize',10) % Right! (0,0,1)
plot3(XX(1,1,1),YY(1,1,1),ZZ(1,1,1),'^k','MarkerSize',10) % Right! (0,0,0)

%%

Xn = [XX(1,end,1) XX(end,1,1) XX(1,1,end) XX(1,1,1)]; % 1 2 3 4  |_\
Yn = [YY(1,end,1) YY(end,1,1) YY(1,1,end) YY(1,1,1)]; 
Zn = [ZZ(1,end,1) ZZ(end,1,1) ZZ(1,1,end) ZZ(1,1,1)];

Xc = sum(Xn)/4;
Yc = sum(Yn)/4;
Zc = sum(Zn)/4;

%% Function and Slopes at corner points

switch Function_Type
    case 'p',
% ---- Plane

 %Func = 2*(XX + 1) + 2*(YY + 1) + 2*(ZZ + 1); % given function
 %dFunc_dx = 2; % slope x
 %dFunc_dy = 2; % slope y
 %dFunc_dz = 2; % slope z

 Func = 12/13*XX + 3*YY + 3/4*ZZ + 14/5; % given function
 dFunc_dx = 12/13; % slope x
 dFunc_dy = 3; % slope y
 dFunc_dz = 3/4; % slope z
 
    case 'q',
% ---- Quadratic Curve

%  Func = (XX - YY.*ZZ).^2 + 2*ZZ.*YY - 0.5*ZZ.^2 + 2; % given function
%  dFunc_dx = 2*(XX - YY.*ZZ); % slope x
%  dFunc_dy = -2*((XX - YY.*ZZ).*ZZ - ZZ); % slope y
%  dFunc_dz = -2*(XX - YY.*ZZ).*YY + 2*YY - ZZ; % slope z
        
Func = (XX - YY).^2 + 2*ZZ.*YY - ZZ.^2 + 2; % given function
dFunc_dx = 2*(XX - YY); % slope x
dFunc_dy = - 2*(XX - YY) + 2*ZZ; % slope y
dFunc_dz = + 2*(YY - ZZ); % slope z
 
    case 'c',
% ---- Cubic Curve

% Func = 2*(XX + 1).^3 + 2*(YY + 1).^2 + 2*(ZZ + 1); % given function Good!
% dFunc_dx = 6*(XX + 1).^2; % slope x
% dFunc_dy = 4*(YY + 1); % slope y
% dFunc_dz = 2; % slope z
 

% Func = XX.^3 + YY.^2 + ZZ;
% dFunc_dx = 3*XX.^2; % slope x
% dFunc_dy = 2*YY; % slope y
% dFunc_dz = 1; % slope z

% Func = YY.^3 + XX.^2 + ZZ;
% dFunc_dx = 2*XX; % slope x
% dFunc_dy = 3*YY.^2; % slope y
% dFunc_dz = 1; % slope z
  
Func = ZZ.^3 + XX.^2 + YY;
dFunc_dx = 2*XX; % slope x
dFunc_dy = 1; % slope y
dFunc_dz = 3*ZZ.^2; % slope z

% Func = (XX.^3) + 2*XX.*(ZZ.^2) + YY.*(ZZ - 2); % given function
% dFunc_dx = 3*(XX.^2) + 2*(ZZ.^2); % slope x
% dFunc_dy = ZZ - 2; % slope y
% dFunc_dz = 4*XX.*ZZ + YY; % slope z
 
% Func = XX.^3 + 2*(XX.*YY.*ZZ).^2 + ZZ.*(YY + 1); % given function
% dFunc_dx = 3*(XX.^2) + 4*XX.*(YY.*ZZ).^2; % slope x
% dFunc_dy = 4*YY.*(XX.*ZZ).^2 + ZZ; % slope y
% dFunc_dz = 4*ZZ.*(XX.*YY).^2 + YY + 1; % slope z
 
%  Func_dx = 2*XX; % slope x
%  dFunc_dy = 3*(YY.^2).*ZZ.^3; % slope y
%  dFunc = (YY.*ZZ).^3 + XX.^2 + (ZZ + 1); % given function
%  dFunc_dz = 3*(ZZ.^2).*YY.^3 + 1; % slope z
 
%  Func = YY.^3 + XX.*ZZ - XX.^2; % given function
%  dFunc_dx = -2*XX + ZZ; % slope x
%  dFunc_dy = 3*(YY.^2);% slope y
%  dFunc_dz = XX; % slope z
 
%  Func = ZZ.^3 + YY.*ZZ.*XX + YY.^2; % given function
%  dFunc_dx = YY.*ZZ; % slope x
%  dFunc_dy = ZZ.*XX + 2*YY; % slope y
%  dFunc_dz = 3*ZZ.^2 + XX.*YY; % slope z
 
 
    case 's',
        
Func = sin(XX.*YY).*cos(YY.*ZZ);
dFunc_dx = cos(XX.*YY).*YY.*cos(YY.*ZZ); % slope x
dFunc_dy = cos(XX.*YY).*XX.*cos(YY.*ZZ) - sin(XX.*YY).*ZZ.*sin(YY.*ZZ); % slope y
dFunc_dz = - sin(YY.*ZZ).*YY; % slope z

    otherwise 'g';
        
% ---- General Curve

Func = exp(YY).*cos(XX) + YY.*sinh(XX) + sin(XX.*ZZ);
dFunc_dx = -sin(XX).*exp(YY) + YY.*cosh(XX) + ZZ.*cos(XX.*ZZ); % slope x
dFunc_dy = exp(YY).*cos(XX) + sin(XX); % slope y
dFunc_dz = XX.*cos(XX.*ZZ); % slope z

end

f = [Func(1,end,1) Func(end,1,1) Func(1,1,end) Func(1,1,1)];
sx = [dFunc_dx(1,end,1) dFunc_dx(end,1,1) dFunc_dx(1,1,end) dFunc_dx(1,1,1)];
sy = [dFunc_dy(1,end,1) dFunc_dy(end,1,1) dFunc_dy(1,1,end) dFunc_dy(1,1,1)];
sz = [dFunc_dz(1,end,1) dFunc_dz(end,1,1) dFunc_dz(1,1,end) dFunc_dz(1,1,1)];


%% Volume Coordinates 
if Input_Volume_Coords_Choice == 1
    
%% Volume Coordinates (New Version)
% 4th order determinant formulae. Set 
% + 1 for even permutations; 
% - 1 for odd permutations. Then, sum all the products a1i * a2j *a3k * a4l
% i,j,k,l = 1,2,3,4 considering the sign.

% -- EVEN indices
%
% 1 as leader: 1234, 1342, 1423;
% 2 as leader: 2143, 2314, 2431;
% 3 as leader: 3124, 3241, 3412;
% 4 as leader: 4132, 4213, 4321;
%
% -- ODD indices
%
% 1 as leader: 1243, 1324, 1432;
% 2 as leader: 2134, 2341, 2413;
% 3 as leader: 3142, 3214, 3421;
% 4 as leader: 4123, 4231, 4312;
%
% Now consider the determinant: 
%
%                 | 1 x1 y1 z1 |
% V = [a_mn]= 1/6 | 1 x2 y2 z2 |
%                 | 1 x3 y3 z3 |
%                 | 1 x4 y4 z4 | 
%
% We have for any volume coordinate - permutating the vertices coordinates -
% the rule for the determinant given by:
%
% V = (1/6)*( a11*a22*a33*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 + ...
%   + a12*a21*a34*a43 + a12*a23*a31*a44 + a12*a24*a33*a41 + ...
%   + a13*a21*a32*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 + ...
%   + a14*a21*a33*a42 + a14*a22*a31*a43 + a14*a23*a32*a41 + ... 
%    ...
%   - a11*a22*a34*a43 - a11*a23*a32*a44 - a11*a24*a33*a42 - ...
%   - a12*a21*a33*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 - ...
%   - a13*a21*a34*a42 - a13*a22*a31*a44 - a13*a24*a32*a41 - ...
%   - a14*a21*a32*a43 - a14*a22*a33*a41 - a14*a23*a31*a42);


% Calculus for V:
a11 = 1;  a12 = Xn(1);  a13 = Yn(1);  a14 = Zn(1);

a21 = 1;  a22 = Xn(2);  a23 = Yn(2);  a24 = Zn(2);

a31 = 1;  a32 = Xn(3);  a33 = Yn(3);  a34 = Zn(3);

a41 = 1;  a42 = Xn(4);  a43 = Yn(4);  a44 = Zn(4);

V = (1/6)*( a11*a22*a33*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 + ...
  + a12*a21*a34*a43 + a12*a23*a31*a44 + a12*a24*a33*a41 + ...
  + a13*a21*a32*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 + ...
  + a14*a21*a33*a42 + a14*a22*a31*a43 + a14*a23*a32*a41 - ...
  ...
  - a11*a22*a34*a43 - a11*a23*a32*a44 - a11*a24*a33*a42 - ...
  - a12*a21*a33*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 - ...
  - a13*a21*a34*a42 - a13*a22*a31*a44 - a13*a24*a32*a41 - ...
  - a14*a21*a32*a43 - a14*a22*a33*a41 - a14*a23*a31*a42);


% Calculus for V1:
a11 = 1;  a12 = XX;  a13 = YY;  a14 = ZZ;

a21 = 1;  a22 = Xn(2);  a23 = Yn(2);  a24 = Zn(2);

a31 = 1;  a32 = Xn(3);  a33 = Yn(3);  a34 = Zn(3);

a41 = 1;  a42 = Xn(4);  a43 = Yn(4);  a44 = Zn(4);

V1 = (1/6)*( a11*a22*a33*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 + ...
  + a12*a21*a34*a43 + a12*a23*a31*a44 + a12*a24*a33*a41 + ...
  + a13*a21*a32*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 + ...
  + a14*a21*a33*a42 + a14*a22*a31*a43 + a14*a23*a32*a41 - ...
  ...
  - a11*a22*a34*a43 - a11*a23*a32*a44 - a11*a24*a33*a42 - ...
  - a12*a21*a33*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 - ...
  - a13*a21*a34*a42 - a13*a22*a31*a44 - a13*a24*a32*a41 - ...
  - a14*a21*a32*a43 - a14*a22*a33*a41 - a14*a23*a31*a42);

% clear a11 a12 a13 a14 a21 a22 a23 a24 a31 a32 a33 a34 a41 a42 a43 a44

% Calculus for V2:
a11 = 1;  a12 = Xn(1);  a13 = Yn(1);  a14 = Zn(1);

a21 = 1;  a22 = XX;  a23 = YY;  a24 = ZZ;

a31 = 1;  a32 = Xn(3);  a33 = Yn(3);  a34 = Zn(3);

a41 = 1;  a42 = Xn(4);  a43 = Yn(4);  a44 = Zn(4);

V2 = (1/6)*( a11*a22*a33*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 + ...
  + a12*a21*a34*a43 + a12*a23*a31*a44 + a12*a24*a33*a41 + ...
  + a13*a21*a32*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 + ...
  + a14*a21*a33*a42 + a14*a22*a31*a43 + a14*a23*a32*a41 - ...
  ...
  - a11*a22*a34*a43 - a11*a23*a32*a44 - a11*a24*a33*a42 - ...
  - a12*a21*a33*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 - ...
  - a13*a21*a34*a42 - a13*a22*a31*a44 - a13*a24*a32*a41 - ...
  - a14*a21*a32*a43 - a14*a22*a33*a41 - a14*a23*a31*a42);


% Calculus for V3:
a11 = 1;  a12 = Xn(1);  a13 = Yn(1);  a14 = Zn(1);

a21 = 1;  a22 = Xn(2);  a23 = Yn(2);  a24 = Zn(2);

a31 = 1;  a32 = XX;  a33 = YY;  a34 = ZZ;

a41 = 1;  a42 = Xn(4);  a43 = Yn(4);  a44 = Zn(4);

V3 = (1/6)*( a11*a22*a33*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 + ...
  + a12*a21*a34*a43 + a12*a23*a31*a44 + a12*a24*a33*a41 + ...
  + a13*a21*a32*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 + ...
  + a14*a21*a33*a42 + a14*a22*a31*a43 + a14*a23*a32*a41 - ...
  ...
  - a11*a22*a34*a43 - a11*a23*a32*a44 - a11*a24*a33*a42 - ...
  - a12*a21*a33*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 - ...
  - a13*a21*a34*a42 - a13*a22*a31*a44 - a13*a24*a32*a41 - ...
  - a14*a21*a32*a43 - a14*a22*a33*a41 - a14*a23*a31*a42);


L1 = V1/V; % volume coordinates
L2 = V2/V;
L3 = V3/V;
L4 = 1 - L1 - L2 - L3;

else
%% Volume Coordinates (Old version )
V = (1/6)*(Xn(2)*Yn(3)*Zn(4) + Xn(4)*Yn(2)*Zn(3) + Xn(3)*Yn(4)*Zn(2)...
         - Xn(3)*Yn(2)*Zn(4) - Xn(2)*Yn(4)*Zn(3) - Xn(4)*Yn(3)*Zn(2)...
         ...
         + Xn(1)*Yn(2)*Zn(4) + Xn(1)*Yn(4)*Zn(3) + Xn(1)*Yn(3)*Zn(2)...
         - Xn(1)*Yn(3)*Zn(4) - Xn(1)*Yn(2)*Zn(3) - Xn(1)*Yn(4)*Zn(2)...
         ...
         + Xn(3)*Yn(1)*Zn(4) + Xn(2)*Yn(1)*Zn(3) + Xn(4)*Yn(1)*Zn(2)...
         - Xn(2)*Yn(1)*Zn(4) - Xn(4)*Yn(1)*Zn(3) - Xn(3)*Yn(1)*Zn(2)...
         ...
         + Xn(2)*Yn(4)*Zn(1) + Xn(4)*Yn(3)*Zn(1) + Xn(3)*Yn(2)*Zn(1)...
         - Xn(3)*Yn(4)*Zn(1) - Xn(2)*Yn(3)*Zn(1) - Xn(4)*Yn(2)*Zn(1)); % element area 

V1 = (1/6)*(Xn(2)*Yn(3)*Zn(4) + Xn(4)*Yn(2)*Zn(3) + Xn(3)*Yn(4)*Zn(2)...
         - Xn(3)*Yn(2)*Zn(4) - Xn(2)*Yn(4)*Zn(3) - Xn(4)*Yn(3)*Zn(2)...
         ...
         + XX*Yn(2)*Zn(4) + XX*Yn(4)*Zn(3) + XX*Yn(3)*Zn(2)...
         - XX*Yn(3)*Zn(4) - XX*Yn(2)*Zn(3) - XX*Yn(4)*Zn(2)...
         ...
         + Xn(3)*YY*Zn(4) + Xn(2)*YY*Zn(3) + Xn(4)*YY*Zn(2)...
         - Xn(2)*YY*Zn(4) - Xn(4)*YY*Zn(3) - Xn(3)*YY*Zn(2)...
         ...
         + Xn(2)*Yn(4)*ZZ + Xn(4)*Yn(3)*ZZ + Xn(3)*Yn(2)*ZZ...
         - Xn(3)*Yn(4)*ZZ - Xn(2)*Yn(3)*ZZ - Xn(4)*Yn(2)*ZZ); 

V2 = (1/6)*(XX*Yn(3)*Zn(4) + Xn(4)*YY*Zn(3) + Xn(3)*Yn(4)*ZZ...
         - Xn(3)*YY*Zn(4) - XX*Yn(4)*Zn(3) - Xn(4)*Yn(3)*ZZ...
         ...
         + Xn(1)*YY*Zn(4) + Xn(1)*Yn(4)*Zn(3) + Xn(1)*Yn(3)*ZZ...
         - Xn(1)*Yn(3)*Zn(4) - Xn(1)*YY*Zn(3) - Xn(1)*Yn(4)*ZZ...
         ...
         + Xn(3)*Yn(1)*Zn(4) + XX*Yn(1)*Zn(3) + Xn(4)*Yn(1)*ZZ...
         - ZZ*Yn(1)*Zn(4) - Xn(4)*Yn(1)*Zn(3) - Xn(3)*Yn(1)*ZZ...
         ...
         + XX*Yn(4)*Zn(1) + Xn(4)*Yn(3)*Zn(1) + Xn(3)*YY*Zn(1)...
         - Xn(3)*Yn(4)*Zn(1) - XX*Yn(3)*Zn(1) - Xn(4)*YY*Zn(1)); 
 
V3 = (1/6)*(Xn(2)*YY*Zn(4) + Xn(4)*Yn(2)*ZZ + XX*Yn(4)*Zn(2)...
         - XX*Yn(2)*Zn(4) - Xn(2)*Yn(4)*ZZ - Xn(4)*YY*Zn(2)...
         ...
         + Xn(1)*Yn(2)*Zn(4) + Xn(1)*Yn(4)*ZZ + Xn(1)*YY*Zn(2)...
         - Xn(1)*YY*Zn(4) - Xn(1)*Yn(2)*ZZ - Xn(1)*Yn(4)*Zn(2)...
         ...
         + XX*Yn(1)*Zn(4) + Xn(2)*Yn(1)*ZZ + Xn(4)*Yn(1)*Zn(2)...
         - Xn(2)*Yn(1)*Zn(4) - Xn(4)*Yn(1)*ZZ - XX*Yn(1)*Zn(2)...
         ...
         + Xn(2)*Yn(4)*Zn(1) + Xn(4)*YY*Zn(1) + XX*Yn(2)*Zn(1)...
         - XX*Yn(4)*Zn(1) - Xn(2)*YY*Zn(1) - Xn(4)*Yn(2)*Zn(1)); 
 
L1 = V1/V; % volume coordinates
L2 = V2/V;
L3 = V3/V;
L4 = 1 - L1 - L2 - L3;

end
[LL1 LL2 LL3 LL4] = funcVol(Xn,Yn,Zn,XX,YY,ZZ);

L1 = LL1;
L2 = LL2;
L3 = LL3;
L4 = LL4;

%% Method

% Linear Approximation

PN1 = f(1).*L1 + f(2).*L2 + f(3).*L3 + f(4).*L4;

% PN1 and Func
% plot3(XX(:,:,1),YY(:,:,1),PN1(:,:,1),'o',XX(:,:,1),YY(:,:,1),Func(:,:,1))

% Jacobian

J = [Xn(1) - Xn(4) Yn(1) - Yn(4) Zn(1) - Zn(4); ...
     Xn(2) - Xn(4) Yn(2) - Yn(4) Zn(2) - Zn(4); ...
     Xn(3) - Xn(4) Yn(3) - Yn(4) Zn(3) - Zn(4)];

% Inverse Jacobian

INVJ = inv(J);

% Matrix to get nodes coordinates to make up directional gradients

XCoord = zeros(3,4);

for j = 1:4
        XCoord(1,j) = Xn(j);
        XCoord(2,j) = Yn(j);
        XCoord(3,j) = Zn(j);    
end

% Get inverse Jacobian and set it along with L4 gradient, which is the opposite 
% of the sum among another gradients

INVJ2 = [INVJ -sum(INVJ,2)];

% Set up volume coordinates onto line arrays

L = [reshape(L1,1,[]); reshape(L2,1,[]); reshape(L3,1,[]); reshape(L4,1,[])];

Pf = L;

PLinear = Pf;

% Below, all sums are carried out:
% m sorts nodal values;
% i sorts space directions;
% j,k sort shape function according to the nodes

% On the r.h.s. of "sum" are the shape functions of TE16 element associated
% to the gradients.

for m = 1:4
    sum = 0;
    for i = 1:3
        for j = 1:4
            for k = 1:4
                sum = sum - 0.5*INVJ2(i,m)*( XCoord(i,k) - XCoord(i,j) )*( L(j,:).*L(k,:) + L(j,:).^2.*L(k,:) - L(j,:).*L(k,:).^2 ); 
            end
        end
    end
    Pf(m,:) = Pf(m,:) + sum;
end

Pg1 = 0*L;
Pg2 = 0*L;
Pg3 = 0*L;

for j = 1:4
    for k = 1:4
        Pg1(j,:) = Pg1(j,:) + 0.5*( XCoord(1,k) - XCoord(1,j))*(L(j,:).*L(k,:) + L(j,:).^2.*L(k,:) - L(j,:).*L(k,:).^2 ); 
        Pg2(j,:) = Pg2(j,:) + 0.5*( XCoord(2,k) - XCoord(2,j))*(L(j,:).*L(k,:) + L(j,:).^2.*L(k,:) - L(j,:).*L(k,:).^2 ); 
        Pg3(j,:) = Pg3(j,:) + 0.5*( XCoord(3,k) - XCoord(3,j))*(L(j,:).*L(k,:) + L(j,:).^2.*L(k,:) - L(j,:).*L(k,:).^2 ); 
    end
end  

%% Interpolation of given function

% Interpolation related to residual approach

Choose_Interpolation_Order = input('Choose interpolation order: [1] to first order; [3] to third order \n','s');

if Choose_Interpolation_Order == '1'
    
    sprintf('%s','Linear interpolation...')
    Func_Approx2 = f(1).*PLinear(1,:) + f(2).*PLinear(2,:) + f(3).*PLinear(3,:) + f(4).*PLinear(4,:);     
    Func_Approx = reshape(Func_Approx2, size(XX));
    Lin_Rel_Err_L2_XY = ( ( ( Func(:,:,1) - Func_Approx(:,:,1) ).^2 )./( max(max(Func(:,:,1))).^2 ) ).^(0.5);

elseif Choose_Interpolation_Order == '3'
    
    sprintf('%s','Cubic interpolation...')
    Func_Approx2 = f(1).*Pf(1,:) + f(2).*Pf(2,:) + f(3).*Pf(3,:) + f(4).*Pf(4,:) ...
               + sx(1).*Pg1(1,:) + sy(1).*Pg2(1,:) + sz(1).*Pg3(1,:) ...
               + sx(2).*Pg1(2,:) + sy(2).*Pg2(2,:) + sz(2).*Pg3(2,:) ...
               + sx(3).*Pg1(3,:) + sy(3).*Pg2(3,:) + sz(3).*Pg3(3,:) ...
               + sx(4).*Pg1(4,:) + sy(4).*Pg2(4,:) + sz(4).*Pg3(4,:);
           
else
    sprintf('%s','Value is not valid! Closing...')
    break
end
               
Func_Approx = reshape(Func_Approx2, size(XX));
                               
%% Faces Plotting



% Handles 'MMMa' are restrictions of domain over tetrahedra faces, where 
% M = X,Y,Z and a = 1,2,3,4;
% Handles 'fffp' and 'iiiq' are restrictions of original function and interpolated
% function, respectively, where p,q = 1,2,3,4.


% ------ X-Y plane

XXX1(:,:) = XX(:,:,1);
YYY1(:,:) = YY(:,:,1);
ZZZ1(:,:) = ZZ(:,:,1);
fff1(:,:) = Func(:,:,1);
iii1(:,:) = Func_Approx(:,:,1);

%  surf(XXX1,YYY1,ZZZ1); % surf only this face


% ------ Y-Z plane

XXX2(:,:) = XX(:,1,:);
YYY2(:,:) = YY(:,1,:);
ZZZ2(:,:) = ZZ(:,1,:);
fff2(:,:) = Func(:,1,:);
iii2(:,:) = Func_Approx(:,1,:);

% surf(XXX2,YYY2,ZZZ2); % surf only this face

% ------ X-Z plane


XXX3(:,:) = XX(1,:,:);
YYY3(:,:) = YY(1,:,:);
ZZZ3(:,:) = ZZ(1,:,:);
fff3(:,:) = Func(1,:,:);
iii3(:,:) = Func_Approx(1,:,:);
    
% surf(XXX3,YYY3,ZZZ3); % surf only this face
   
% ------ diagonal plane


 XXX4(:,:) = XX(:,end,:); 
 YYY4(:,:) = YY(:,end,:);
 ZZZ4(:,:) = ZZ(:,end,:);
 fff4(:,:) = Func(:,end,:);
 iii4(:,:) = Func_Approx(:,end,:);

%surf(XXX4,YYY4,ZZZ4); % surf only this face

%% Comparison Plotting by 'plot3'

% Switch faces, function and approximation among (XXXq,YYYq,fffq,iiiq) q = 1,2,3,4 here!!



 subplot(2,2,1)
 hold on
% surf(XXX1,YYY1,iii1)
 plot3(XXX1,YYY1,iii1,'or')
 hold on
 plot3(XXX1,YYY1,fff1,'.k')

 subplot(2,2,2)
% surf(XXX2,YYY2,iii2);
 plot3(XXX2,YYY2,iii2,'og');
 hold on
 plot3(XXX2,YYY2,fff2,'.k');
% 
 subplot(2,2,3)
% %surf(XXX3,YYY3,iii3);
  plot3(XXX3,YYY3,iii3,'ob');
  hold on
  plot3(XXX3,YYY3,fff3,'.k');
%      
 subplot(2,2,4)
% %surf(XXX4,YYY4,iii4);
  plot3(XXX4,YYY4,iii4,'om');
  hold on
  plot3(XXX4,YYY4,fff4,'.k');

% subplot(1,2,2)
% 
% hold on
% surf(XXX1,YYY1,iii1)
% surf(XXX2,YYY2,iii2);
% surf(XXX3,YYY3,iii3);
% surf(XXX4,YYY4,iii4);
  
% Plots of function values over edges. Switch among (XXXq,YYYq,iiiq) q = 1,2,3,4.


%   plot3(XXX4(:,end),YYY4(:,end),iii4(:,end), 'ok'); % plot points over diagonal edge
%   plot3(XXX4(1,:),YYY4(1,:),iii4(1,:), 'ok'); 
%   plot3(XXX4(:,1),YYY4(:,1),iii4(:,1), 'ok'); 
% 
%   plot3(XXX1(:,end),YYY1(:,end),iii1(:,end), 'or'); % plot points over diagonal edge
%   plot3(XXX1(1,:),YYY1(1,:),iii1(1,:), 'or'); 
%   plot3(XXX1(:,1),YYY1(:,1),iii1(:,1), 'or'); 
% 
%   plot3(XXX2(:,end),YYY2(:,end),iii2(:,end), 'oy'); % plot points over diagonal edge
%   plot3(XXX2(1,:),YYY2(1,:),iii2(1,:), 'oy'); 
%   plot3(XXX2(:,1),YYY2(:,1),iii2(:,1), 'oy'); 
% 
%   plot3(XXX3(:,end),YYY3(:,end),iii3(:,end), 'og'); % plot points over diagonal edge
%   plot3(XXX3(1,:),YYY3(1,:),iii3(1,:), 'og'); 
%   plot3(XXX3(:,1),YYY3(:,1),iii3(:,1), 'og'); 

%  subplot(1,2,1)
%  colormap jet
%  shading interp
% %xlabel 'L1'; ylabel 'L2'; % axis labels
%  surf(XXX1,YYY1,fff1);
%  
%  subplot(1,2,2)
%  grid on
%  plot3(XXX1,YYY1,iii1);



% ------ end of plotting of faces and functions -------

%% Error Analysis

% Calculus of L2-norm relative error for each face:

% X-Y face

%Rel_Err_L2_XY = ( ( ( Func(:,:,1) - Func_Approx(:,:,1) ).^2 )./( max(max(Func(:,:,1))).^2 ) ).^(0.5);
%surf(XXX1,YYY1,Rel_Err_L2_XY);
%shading interp

%Rel_Err_L2_XZ = ( ( ( Func(1,:,:) - Func_Approx(1,:,:) ).^2 )./( max(max(Func(1,:,:))).^2 ) ).^(0.5);
%surf(XXX3,ZZZ3,Rel_Err_L2_XZ);
%shading interp

% for i = 1:size(Func,1)
%     for j = 1:size(Func,2)
%         Err_2 = ( Func(i,j,1) - Func_Approx(i,j,1) ).^2;
%         Func_2 = Func(i,j,1).^2;
%         Norm_Rel_Err_L2 = Err_2/Func_2;
%     end
% end

%% Comparison Plotting by "isosurface"


%  iso_Func_Approx = isosurface(XX,YY,ZZ,Func_Approx,7);
%  isosurface(XX,YY,ZZ,Func_Approx,7);
%  iso_Func = isosurface(XX,YY,ZZ,Func,7);
%  isosurface(XX,YY,ZZ,Func,7);
% % 
%  subplot(1,2,1);
%  pt = patch(iso_Func_Approx);
%  subplot(1,2,2);
%  pt2 = patch(iso_Func);

%% Profile Data

if Profile_Power_On_Off == 'y'
    sprintf('%s','Showing profile...')
    profile viewer
    sprintf('%s','Code run successfully!')
else 
    sprintf('%s','Code run successfully!')
    break
end

%______** end of code **______%
