% File: perturbations_seeing_2D.m
% Author: Gustavo Charles P. de Oliveira
% Description: Rebuilds perturbation of Phys. of Fluids and plots field
% over the disk. 

%% Defaults

clear all;
format long;

%% Grid 

% r coordinates
r_min = 0; 
r_max = 1; % 5e-3;

% theta coordinates
factor_reducing = 1;

t_min = 0; 
t_max = 2*pi/factor_reducing; % 2.18*pi fits well to the circle closing!

z_min = 0;
z_max = 50;

npr = 80; npt = 80; npz = 200;

rr = linspace(r_min,r_max,npr);
tt = linspace(t_min,t_max,npt);
zz = linspace(z_min,z_max,npz);

dz = (z_max - z_min)/(npz - 1);

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

% F = fid_F(:,2);
% G = fid_G(:,2);
% H = fid_H(:,2);
% C = fid_C(:,2);

F = fid_F(1:npz,2);
G = fid_G(1:npz,2);
H = fid_H(1:npz,2);
C = fid_C(1:npz,2);

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
%h_Real_Chem = fid_Real_h_Chem(:,2);
%h_Imag_Chem = fid_Imag_h_Chem(:,2);

h_Real_Chem = fid_Real_h_Chem(1:npz,2);
h_Imag_Chem = fid_Imag_h_Chem(1:npz,2);

h_Chem = h_Real_Chem + i*h_Imag_Chem; % real + imag of h

% eta
% eta_Real_Chem = fid_Real_eta_Chem(:,2);
% eta_Imag_Chem = fid_Imag_eta_Chem(:,2);

eta_Real_Chem = fid_Real_eta_Chem(1:npz,2);
eta_Imag_Chem = fid_Imag_eta_Chem(1:npz,2);


eta_Chem = eta_Real_Chem + i*eta_Imag_Chem; % real + imag of eta

% c
% c_Real_Chem = fid_Real_c_Chem(:,2);
% c_Imag_Chem = fid_Imag_c_Chem(:,2); % graphic of this part is not the same as that from paper one!

c_Real_Chem = fid_Real_c_Chem(1:npz,2);
c_Imag_Chem = fid_Imag_c_Chem(1:npz,2); % graphic of this part is not the same as that from paper one!

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

%% Quiver 2D

% quiver 2D

pz = 2; % switch to vary on z of quiver 2D
a = 1;
b = 80;
xaf = X(a:b,a:b,pz);
yaf = Y(a:b,a:b,pz);
zaf = Z(a:b,a:b,pz);

% vr, vt - 2D
vrf = Perturb_f_Complete_Chem_RealPart(a:b,a:b,pz);
vtf = Perturb_g_Complete_Chem_RealPart(a:b,a:b,pz);

ttf = (TT(a:b,a:b,pz)); 
vxf = vrf.*cos(ttf) - vtf.*sin(ttf);
vyf = vtf.*cos(ttf) + vrf.*sin(ttf);

 quiver(xaf,yaf,vxf,vyf,.8)


