% File: perturbation_built_3D.m
% Author: Gustavo Charles P. de Oliveira
% Updated: Oct, 2013
%
% Description: Rebuilds the perturbation given by the paper of 
% Phys. of Fluids. Perturbations are summed to the base state. It is
% possible to run the code loading chemical modes and hydrodynamic modes.

%% Defaults

clear all;
format long;

%% Grid 

x_min =  - sqrt(1.25e-3); % I use this value so that 5e-3 be the diameter
x_max =    sqrt(1.25e-3);

y_min =  - sqrt(1.25e-3);
y_max =    sqrt(1.25e-3);

z_min = 0;
z_max = 30;

npx = 41; npy = 41; npz = 2401;

% r coordinates
r_min = 0; 
r_max = 0.1; % 5e-3;

% theta coordinates
factor_reducing = 1; % good result with 20! 

t_min = 0; 
t_max = 2*pi/factor_reducing;

% vectors
xx = linspace(x_min,x_max,npx);
yy = linspace(y_min,y_max,npy);
zz = linspace(z_min,z_max,npz);

dz = (z_max - z_min)/(npz - 1);

% discretization points
npr = 31; npt = 31;

rr = linspace(r_min,r_max,npr);
tt = linspace(t_min,t_max,npt);

% dr = (r_max - r_min)/(npr - 1);
% dt = (t_max - t_min)/(npt - 1);

% [XX YY ZZ] = meshgrid(xx,yy,zz);

% [T R Z] = cart2pol(XX,YY,ZZ); % transforming cartesian coordinates to cylindrical coordinates

[TT RR ZZ] = meshgrid(tt,rr,zz); % mesh generation

[X Y Z] = pol2cart(TT,RR,ZZ); % transforming polar coordinates to cartesian coordinates

Omega = 900; % angular velocity
nu_infty = 30*10^(-3); % bulk viscosity
nu_zero = nu_infty*2.255; % viscosity near to electrode surface

Re = RR*(Omega/nu_infty)^(.5); % r adimensionalization = Reynolds

%% Perturbation Framework

% components of perturbation wave vector 

% wave vector 1 
% alpha = 0.50869; 
% beta = 0.10206;

% wave vector 2
alpha = 0.5;
beta = 0.33949;

Reynolds = 60.394;
% Reynolds = 208.77; %mean(mean(mean(Re))); % specific point at which I analyze. I chose the mean because 
                           % I did not obtain great values for Reynolds
                           % according to Re vector! Verify!!!

Perturb = 0*X;
Exp_Perturb = 0*Perturb; % exponential basis

for j = 1:length(zz)
    Perturb(:,:,j) = i*(alpha*Re(:,:,j) + beta*TT(:,:,j)*Reynolds);
    Exp_Perturb(:,:,j) = exp(Perturb(:,:,j));
end

%% Base State

% loading base state

% Both operations below load the same F,G,H and C profiles for chemical as for hydrodynamic
% region

% hydrodynamic region
fid_F = load('Perfis de simulação/dados_estado_base_e_perturbacao_regiao_hidro/F');
fid_G = load('Perfis de simulação/dados_estado_base_e_perturbacao_regiao_hidro/G');
fid_H = load('Perfis de simulação/dados_estado_base_e_perturbacao_regiao_hidro/H');
fid_C = load('Perfis de simulação/dados_estado_base_e_perturbacao_regiao_hidro/C');

F = fid_F(:,2);
G = fid_G(:,2);
H = fid_H(:,2);
C = fid_C(:,2);

F_Adim = zeros(size(ZZ));
G_Adim = zeros(size(ZZ));
H_Adim = zeros(size(ZZ));
C_Adim = zeros(size(ZZ));

% setting adimensional F,G,H,C 
for m = 1:size(ZZ,1)
    for jj = 1:size(ZZ,2)
        for kk = 1:size(ZZ,3)
                F_Adim(m,jj,kk) = (1/Reynolds).*Re(m,jj,kk).*F(kk,1);
                G_Adim(m,jj,kk) = (1/Reynolds).*Re(m,jj,kk).*G(kk,1);
                H_Adim(m,jj,kk) = (1/Reynolds).*H(kk,1);
                C_Adim(m,jj,kk) = C(kk,1);
        end
    end
end


if 0 == 1
%% Perturbation Profiles - Hydrodynamic Region


% loading perturbation profiles h, eta, c

% h
    fid_Real_h_Hydro = load('Perfis de simulação/dados_estado_base_e_perturbacao_regiao_hidro/Re_hHydro2.norm.dat');
    fid_Imag_h_Hydro = load('Perfis de simulação/dados_estado_base_e_perturbacao_regiao_hidro/Im_hHydro2.norm.dat');

% eta
    fid_Real_eta_Hydro = load('Perfis de simulação/dados_estado_base_e_perturbacao_regiao_hidro/Re_etaHydro2.norm.dat');
    fid_Imag_eta_Hydro = load('Perfis de simulação/dados_estado_base_e_perturbacao_regiao_hidro/Im_etaHydro2.norm.dat');

% c
    fid_Real_c_Hydro = load('Perfis de simulação/dados_estado_base_e_perturbacao_regiao_hidro/Re_CHydro2.norm.dat');
    fid_Imag_c_Hydro = load('Perfis de simulação/dados_estado_base_e_perturbacao_regiao_hidro/Im_CHydro2.norm.dat');


% extracting profiles: real and imaginary parts

% h
    h_Real_Hydro = fid_Real_h_Hydro(:,2);
    h_Imag_Hydro = fid_Imag_h_Hydro(:,2);

    h = h_Real_Hydro + i*h_Imag_Hydro; % real + imag of h

% eta
    eta_Real_Hydro = fid_Real_eta_Hydro(:,2);
    eta_Imag_Hydro = fid_Imag_eta_Hydro(:,2);

    eta = eta_Real_Hydro + i*eta_Imag_Hydro; % real + imag of eta

% c
    c_Real_Hydro = fid_Real_c_Hydro(:,2);
    c_Imag_Hydro = fid_Imag_h_Hydro(:,2);

    c = c_Real_Hydro + i*c_Imag_Hydro; % real + imag of c


% derivative of h (centered)
    dh_Real_dz_Hydro = 0*h_Real_Hydro;
    dh_Imag_dz_Hydro = 0*h_Imag_Hydro;

% boundary conditions are zero!
    for n = 2:length(dh_Real_dz_Hydro) - 1
        dh_Real_dz_Hydro(n) = (h_Real_Hydro(n+1) - h_Real_Hydro(n-1))/(2*dz);
        dh_Imag_dz_Hydro(n) = (h_Imag_Hydro(n+1) - h_Imag_Hydro(n-1))/(2*dz);
    end

    dhR = dh_Real_dz_Hydro;
    dhI = dh_Imag_dz_Hydro;

    dh = dhR + i*dhI; % real + imag of h'


% Recovering f and g profiles from continuity equation

% f_Hydro
    f = - dh./(i*(alpha - 1/Reynolds + beta*eta + beta^2/alpha));

% f_Hydro real and imaginary parts
    f_Real_Hydro = real(f);
    f_Imag_Hydro = imag(f);

% g_Hydro
    g = 1/alpha*(beta*f + eta);

% g_Hydro real and imaginary parts
    g_Real_Hydro = real(g);
    g_Imag_Hydro = imag(g);

 
% Setting Exponential Perturbation  
    Perturb_f_Hydro = zeros(size(ZZ));
    Perturb_g_Hydro = zeros(size(ZZ));
    Perturb_h_Hydro = zeros(size(ZZ));
    Perturb_c_Hydro = zeros(size(ZZ));

    for p = 1:size(ZZ,1)
        for q = 1:size(ZZ,2)
            for r = 1:size(ZZ,3)
                Perturb_f_Hydro(p,q,r) = Exp_Perturb(p,q,r)*f(r,1);
                Perturb_g_Hydro(p,q,r) = Exp_Perturb(p,q,r)*g(r,1);
                Perturb_h_Hydro(p,q,r) = Exp_Perturb(p,q,r)*h(r,1);
                Perturb_c_Hydro(p,q,r) = Exp_Perturb(p,q,r)*c(r,1);
            end
        end
    end


% Conjugates
    Perturb_f_Hydro_Conjug = conj(Perturb_f_Hydro); % conjugated of f
    Perturb_g_Hydro_Conjug = conj(Perturb_g_Hydro); % conjugated of g
    Perturb_h_Hydro_Conjug = conj(Perturb_h_Hydro); % conjugated of h
    Perturb_c_Hydro_Conjug = conj(Perturb_c_Hydro); % conjugated of c

% Complete perturbations

    Perturb_f_Hydro_Complete = Perturb_f_Hydro + Perturb_f_Hydro_Conjug; % f complete
    Perturb_g_Hydro_Complete = Perturb_g_Hydro + Perturb_g_Hydro_Conjug; % g complete
    Perturb_h_Hydro_Complete = Perturb_h_Hydro + Perturb_h_Hydro_Conjug; % h complete
    Perturb_c_Hydro_Complete = Perturb_c_Hydro + Perturb_c_Hydro_Conjug; % c complete


% Real parts
    Perturb_f_Complete_Hydro_RealPart = real(Perturb_f_Hydro_Complete); % real part of f complete
    Perturb_g_Complete_Hydro_RealPart = real(Perturb_g_Hydro_Complete); % real part of g complete
    Perturb_h_Complete_Hydro_RealPart = real(Perturb_h_Hydro_Complete); % real part of h complete
    Perturb_c_Complete_Hydro_RealPart = real(Perturb_c_Hydro_Complete); % real part of c complete


% Imaginary parts
    Perturb_f_Complete_Hydro_ImagPart = imag(Perturb_f_Hydro_Complete); % imaginary part of f complete
    Perturb_g_Complete_Hydro_ImagPart = imag(Perturb_g_Hydro_Complete); % imaginary part of g complete
    Perturb_h_Complete_Hydro_ImagPart = imag(Perturb_h_Hydro_Complete); % imaginary part of h complete
    Perturb_c_Complete_Hydro_ImagPart = imag(Perturb_c_Hydro_Complete); % imaginary part of c complete

%% Quivers - Hydrodynamic Region

% theta x z - plane

    x = Re(16,:,1:25:end/2);
    z = Z(16,:,1:25:end/2);

    u = Perturb_g_Complete_Hydro_RealPart(16,:,1:25:end/2);
    w = Perturb_h_Complete_Hydro_RealPart(16,:,1:25:end/2);
 
    u = reshape(u,size(x,2),48);
    w = reshape(w,size(z,2),48);

    x = reshape(x,size(x,2),48);
    z = reshape(z,size(z,2),48);

    quiver(x,z,u,w)

% r x z - plane
     y = TT(:,7,1:25:end/2);
     z = Z(:,7,1:25:end/2);
 
     v = Perturb_g_Complete_Hydro_RealPart(:,7,1:25:end/2);
     w = Perturb_h_Complete_Hydro_RealPart(:,7,1:25:end/2);
    
     v = reshape(v,size(y,1),48);
     w = reshape(w,size(z,1),48);
  
     y = reshape(y,size(y,1),48);
     z = reshape(z,size(z,1),48);
 
     quiver(y,z,v,w)

% ---------------------------------------------------------
% ------ End of hydrodynamic perturbations framework ------
% ---------------------------------------------------------

else
%% Perturbation Profiles - Chemical Region


% loading perturbation profiles h, eta, c

% h
    fid_Real_h_Chem = load('Perfis de simulação/dados_estado_base_e_perturbacao_regiao_quimica/Re_h_point2.norm.dat');
    fid_Imag_h_Chem = load('Perfis de simulação/dados_estado_base_e_perturbacao_regiao_quimica/Im_h_point2.norm.dat');

% eta
    fid_Real_eta_Chem = load('Perfis de simulação/dados_estado_base_e_perturbacao_regiao_quimica/Re_eta_point2.norm.dat');
    fid_Imag_eta_Chem = load('Perfis de simulação/dados_estado_base_e_perturbacao_regiao_quimica/Im_eta_point2.norm.dat');

% c
    fid_Real_c_Chem = load('Perfis de simulação/dados_estado_base_e_perturbacao_regiao_quimica/Re_C_point2.norm.dat');
    fid_Imag_c_Chem = load('Perfis de simulação/dados_estado_base_e_perturbacao_regiao_quimica/Im_C_point2.norm.dat');


% extracting profiles: real and imaginary parts

% h
    h_Real_Chem = fid_Real_h_Chem(:,2);
    h_Imag_Chem = fid_Imag_h_Chem(:,2);

    h_Chem = h_Real_Chem + i*h_Imag_Chem; % real + imag of h

% eta
    eta_Real_Chem = fid_Real_eta_Chem(:,2);
    eta_Imag_Chem = fid_Imag_eta_Chem(:,2);

    eta_Chem = eta_Real_Chem + i*eta_Imag_Chem; % real + imag of eta

% c
    c_Real_Chem = fid_Real_c_Chem(:,2);
    c_Imag_Chem = fid_Imag_c_Chem(:,2); % graphic of this part is not the same as that from paper one!

    c_Chem = c_Real_Chem + i*c_Imag_Chem; % real + imag of c


% derivative of h (centered)
    dh_Real_dz_Chem = 0*h_Real_Chem;
    dh_Imag_dz_Chem = 0*h_Imag_Chem;


% boundary conditions are zero!
    for n = 2:length(dh_Real_dz_Chem) - 1
        dh_Real_dz_Chem(n) = (h_Real_Chem(n+1) - h_Real_Chem(n-1))/(2*dz);
        dh_Imag_dz_Chem(n) = (h_Imag_Chem(n+1) - h_Imag_Chem(n-1))/(2*dz);
    end

    dhRC = dh_Real_dz_Chem;
    dhIC = dh_Imag_dz_Chem;

    dhC = dhRC + i*dhIC; % real + imag of h'


% Recovering f and g profiles from continuity equation - chemical region

% f_Chem
    f_Chem = - dhC./(i*(alpha - 1/Reynolds + beta*eta_Chem + beta^2/alpha));

% f_Chem real and imaginary parts
    f_Real_Chem = real(f_Chem);
    f_Imag_Chem = imag(f_Chem);

% g_Chem
    g_Chem = 1/alpha*(beta*f_Chem + eta_Chem);

% g_Chem real and imaginary parts
    g_Real_Chem = real(g_Chem);
    g_Imag_Chem = imag(g_Chem);

% Setting Exponential Perturbation  
    Perturb_f_Chem = zeros(size(ZZ));
    Perturb_g_Chem = zeros(size(ZZ));
    Perturb_h_Chem = zeros(size(ZZ));
    Perturb_c_Chem = zeros(size(ZZ));

    for p = 1:size(ZZ,1)
        for q = 1:size(ZZ,2)
            for r = 1:size(ZZ,3)
                Perturb_f_Chem(p,q,r) = Exp_Perturb(p,q,r)*f_Chem(r,1);
                Perturb_g_Chem(p,q,r) = Exp_Perturb(p,q,r)*g_Chem(r,1);
                Perturb_h_Chem(p,q,r) = Exp_Perturb(p,q,r)*h_Chem(r,1);
                Perturb_c_Chem(p,q,r) = Exp_Perturb(p,q,r)*c_Chem(r,1);
            end
        end
    end


% Conjugates
    Perturb_f_Chem_Conjug = conj(Perturb_f_Chem); % conjugated of f
    Perturb_g_Chem_Conjug = conj(Perturb_g_Chem); % conjugated of g
    Perturb_h_Chem_Conjug = conj(Perturb_h_Chem); % conjugated of h
    Perturb_c_Chem_Conjug = conj(Perturb_c_Chem); % conjugated of c

% Complete perturbations

    Perturb_f_Chem_Complete = Perturb_f_Chem + Perturb_f_Chem_Conjug; % f complete
    Perturb_g_Chem_Complete = Perturb_g_Chem + Perturb_g_Chem_Conjug; % g complete
    Perturb_h_Chem_Complete = Perturb_h_Chem + Perturb_h_Chem_Conjug; % h complete
    Perturb_c_Chem_Complete = Perturb_c_Chem + Perturb_c_Chem_Conjug; % c complete


% Real parts
    Perturb_f_Complete_Chem_RealPart = real(Perturb_f_Chem_Complete); % real part of f complete
    Perturb_g_Complete_Chem_RealPart = real(Perturb_g_Chem_Complete); % real part of g complete
    Perturb_h_Complete_Chem_RealPart = real(Perturb_h_Chem_Complete); % real part of h complete
    Perturb_c_Complete_Chem_RealPart = real(Perturb_c_Chem_Complete); % real part of c complete


% Imaginary parts
    Perturb_f_Complete_Chem_ImagPart = imag(Perturb_f_Chem_Complete); % imaginary part of f complete
    Perturb_g_Complete_Chem_ImagPart = imag(Perturb_g_Chem_Complete); % imaginary part of g complete
    Perturb_h_Complete_Chem_ImagPart = imag(Perturb_h_Chem_Complete); % imaginary part of h complete
    Perturb_c_Complete_Chem_ImagPart = imag(Perturb_c_Chem_Complete); % imaginary part of c complete

% ---------------------------------------------------------
% ------ End of chemical perturbations framework ------
% ---------------------------------------------------------

end

%% Quivers - Chemical Region
 
if 1 == 2
% r x z - plane (It works well!)
 
  x = Re(16,:,1:25:end/2);
  z = Z(16,:,1:25:end/2);
 
  u = Perturb_f_Complete_Chem_RealPart(16,:,1:25:end/2);
  w = Perturb_h_Complete_Chem_RealPart(16,:,1:25:end/2);
  
  u = reshape(u,size(x,2),48);
  w = reshape(w,size(z,2),48);
 
  x = reshape(x,size(x,2),48);
  z = reshape(z,size(z,2),48);
 
  quiver(x,z,u,w);

elseif 1==3

  x = TT(:,3,1:25:end/2);
  z = Z(:,3,1:25:end/2);
 
  u = Perturb_g_Complete_Chem_RealPart(:,3,1:25:end/2);
  w = Perturb_h_Complete_Chem_RealPart(:,3,1:25:end/2);
  
  u = reshape(u,size(x,1),48);
  w = reshape(w,size(z,1),48);
 
  x = reshape(x,size(x,1),48);
  z = reshape(z,size(z,1),48);
 
  
  quiver(x,z,u,w)
 
elseif 1 == 4

% quiver 3D

xf = X(20:26,5:31,1:25:150);
yf = Y(20:26,5:31,1:25:150);
zf = Z(20:26,5:31,1:25:150);
 

% vr, vt, vz - 3D
vr = Perturb_f_Complete_Chem_RealPart(20:26,5:31,1:25:150);
vt = Perturb_g_Complete_Chem_RealPart(20:26,5:31,1:25:150);
vz = Perturb_h_Complete_Chem_RealPart(20:26,5:31,1:25:150);

tt = (TT(20:26,5:31,1:25:150));
vx = vr.*cos(tt) - vt.*sin(tt);
vy = vt.*cos(tt) + vr.*sin(tt);

%quiver3(xf,yf,zf,vx,vy,vz)


else
% quiver 2D

pz = 3; % switch to vary on z of quiver 2D

xaf = X(3:31,3:31);
yaf = Y(3:31,3:31);
zaf = Z(3:31,3:31);

% vr, vt - 2D
vrf = Perturb_f_Complete_Chem_RealPart(3:31,3:31,pz);
vtf = Perturb_g_Complete_Chem_RealPart(3:31,3:31,pz);

ttf = 10*(TT(3:31,3:31,pz)); 
vxf = vrf.*cos(ttf) - vtf.*sin(ttf);
vyf = vtf.*cos(ttf) + vrf.*sin(ttf);

quiver(xaf,yaf,vxf,vyf/2)

end

profile off;

%%

%%% ------------------------------------ %%%
%%% ----------- End of code! ----------- %%%
%%% ------------------------------------ %%%
