% File: vtkToMatlab.m
% Author: Gustavo Charles P. de Oliveira
% 
% Modified: Oct 27th, 2013.
%
% Description: Routine to import data from Paraview to MATLAB. Code
% developed to read files of concentration (chemical species) values along
% several radius of disk. Finally, we obtain plots of perturbation
% oscillation and calculate the Sample Variance of data, also by plotting
% it.


%% Defaults

clear all;

%% Loading files

% raio 5
fid_raio5_iter0 = load('iteracao0/raio-5-iteracao-0.dat');
fid_raio5_iter12 = load('iteracao12/raio-5-iteracao-12.dat');
fid_raio5_iter26 = load('iteracao26/raio-5-iteracao-26.dat');
fid_raio5_iter38 = load('iteracao38/raio-5-iteracao-38.dat');

% raio 15
fid_raio15_iter0 = load('iteracao0/raio-15-iteracao-0.dat');
fid_raio15_iter12 = load('iteracao12/raio-15-iteracao-12.dat');
fid_raio15_iter26 = load('iteracao26/raio-15-iteracao-26.dat');
fid_raio15_iter38 = load('iteracao38/raio-15-iteracao-38.dat');

% raio 45
fid_raio45_iter0 = load('iteracao0/raio-45-iteracao-0.dat');
fid_raio45_iter12 = load('iteracao12/raio-45-iteracao-12.dat');
fid_raio45_iter26 = load('iteracao26/raio-45-iteracao-26.dat');
fid_raio45_iter38 = load('iteracao38/raio-45-iteracao-38.dat');

% raio 60
fid_raio60_iter0 = load('iteracao0/raio-60-iteracao-0.dat');
fid_raio60_iter12 = load('iteracao12/raio-60-iteracao-12.dat');
fid_raio60_iter26 = load('iteracao26/raio-60-iteracao-26.dat'); 
fid_raio60_iter38 = load('iteracao38/raio-60-iteracao-38.dat'); 

% raio 80
fid_raio80_iter0 = load('iteracao0/raio-80-iteracao-0.dat');
fid_raio80_iter12 = load('iteracao12/raio-80-iteracao-12.dat');
fid_raio80_iter26 = load('iteracao26/raio-80-iteracao-26.dat');
fid_raio80_iter38 = load('iteracao38/raio-80-iteracao-38.dat');

% raio 100
fid_raio100_iter0 = load('iteracao0/raio-100-iteracao-0.dat');
fid_raio100_iter12 = load('iteracao12/raio-100-iteracao-12.dat');
fid_raio100_iter26 = load('iteracao26/raio-100-iteracao-26.dat');
fid_raio100_iter38 = load('iteracao38/raio-100-iteracao-38.dat');

%% Loading concentration values along radius

% raio 5
concentracao_raio5_iter0 = fid_raio5_iter0(:,2);
concentracao_raio5_iter12 = fid_raio5_iter12(:,2);
concentracao_raio5_iter26 = fid_raio5_iter26(:,2);
concentracao_raio5_iter38 = fid_raio5_iter38(:,2);
 
% raio 15
concentracao_raio15_iter0 = fid_raio15_iter0(:,2);
concentracao_raio15_iter12 = fid_raio15_iter12(:,2);
concentracao_raio15_iter26 = fid_raio15_iter26(:,2);
concentracao_raio15_iter38 = fid_raio15_iter38(:,2);

% raio 45
concentracao_raio45_iter0 = fid_raio45_iter0(:,2);
concentracao_raio45_iter12 = fid_raio45_iter12(:,2);
concentracao_raio45_iter26 = fid_raio45_iter26(:,2);
concentracao_raio45_iter38 = fid_raio45_iter38(:,2);

% raio 60
concentracao_raio60_iter0 = fid_raio60_iter0(:,2);
concentracao_raio60_iter12 = fid_raio60_iter12(:,2);
concentracao_raio60_iter26 = fid_raio60_iter26(:,2); 
concentracao_raio60_iter38 = fid_raio60_iter38(:,2); 

% raio 80
concentracao_raio80_iter0 = fid_raio80_iter0(:,2);
concentracao_raio80_iter12 = fid_raio80_iter12(:,2);
concentracao_raio80_iter26 = fid_raio80_iter26(:,2);
concentracao_raio80_iter38 = fid_raio80_iter38(:,2); 

% raio 100
concentracao_raio100_iter0 = fid_raio100_iter0(:,2);
concentracao_raio100_iter12 = fid_raio100_iter12(:,2);
concentracao_raio100_iter26 = fid_raio100_iter26(:,2);
concentracao_raio100_iter38 = fid_raio100_iter38(:,2);

%% Mean concentration

% raio 5
concentracao_media_raio5_iter0 = mean(concentracao_raio5_iter0);
concentracao_media_raio5_iter12 = mean(concentracao_raio5_iter12);
concentracao_media_raio5_iter26 = mean(concentracao_raio5_iter26);
concentracao_media_raio5_iter38 = mean(concentracao_raio5_iter38);

% raio 15
concentracao_media_raio15_iter0 = mean(concentracao_raio15_iter0);
concentracao_media_raio15_iter12 = mean(concentracao_raio15_iter12);
concentracao_media_raio15_iter26 = mean(concentracao_raio15_iter26);
concentracao_media_raio15_iter38 = mean(concentracao_raio15_iter38);
 
% raio 45
concentracao_media_raio45_iter0 = mean(concentracao_raio45_iter0);
concentracao_media_raio45_iter12 = mean(concentracao_raio45_iter12);
concentracao_media_raio45_iter26 = mean(concentracao_raio45_iter26);
concentracao_media_raio45_iter38 = mean(concentracao_raio45_iter38);
 
% raio 60
concentracao_media_raio60_iter0 = mean(concentracao_raio60_iter0);
concentracao_media_raio60_iter12 = mean(concentracao_raio60_iter12);
concentracao_media_raio60_iter26 = mean(concentracao_raio60_iter26);
concentracao_media_raio60_iter38 = mean(concentracao_raio60_iter38);
 
% raio 80
concentracao_media_raio80_iter0 = mean(concentracao_raio80_iter0);
concentracao_media_raio80_iter12 = mean(concentracao_raio80_iter12);
concentracao_media_raio80_iter26 = mean(concentracao_raio80_iter26);
concentracao_media_raio80_iter38 = mean(concentracao_raio80_iter38);
 
% raio 100
concentracao_media_raio100_iter0 = mean(concentracao_raio100_iter0);
concentracao_media_raio100_iter12 = mean(concentracao_raio100_iter12);
concentracao_media_raio100_iter26 = mean(concentracao_raio100_iter26);
concentracao_media_raio100_iter38 = mean(concentracao_raio80_iter38);

%% Number of points of evaluation per radius

% radius 5
n_raio5 = size(concentracao_raio5_iter12,1);

%radius 15
n_raio15 = size(concentracao_raio15_iter12,1);

%radius 45
n_raio45 = size(concentracao_raio45_iter12,1);

%radius 60
n_raio60 = size(concentracao_raio60_iter12,1);

% radius 80
n_raio80 = size(concentracao_raio80_iter12,1);

% radius 100
n_raio100 = size(concentracao_raio100_iter12,1);

%% Radius

raios = [5 15 45 60 80 100];

arc_length_raio5 = linspace(0,2*pi*raios(1),n_raio5);
arc_length_raio15 = linspace(0,2*pi*raios(2),n_raio15);
arc_length_raio45 = linspace(0,2*pi*raios(3),n_raio45);
arc_length_raio60 = linspace(0,2*pi*raios(4),n_raio60);
arc_length_raio80 = linspace(0,2*pi*raios(5),n_raio80);
arc_length_raio100 = linspace(0,2*pi*raios(6),n_raio100);

%% partition for plot
dx = 0.001;

x_raio5 = 0:dx:max(arc_length_raio5);
x_raio15 = 0:dx:max(arc_length_raio15);
x_raio45 = 0:dx:max(arc_length_raio45);
x_raio60 = 0:dx:max(arc_length_raio60);
x_raio80 = 0:dx:max(arc_length_raio80);
x_raio100 = 0:dx:max(arc_length_raio100);


%% Mean concentration vector for plotting

% radius 5
ymedia_raio5_iter0 = concentracao_media_raio5_iter0*ones(1,length(x_raio5));
ymedia_raio5_iter12 = concentracao_media_raio5_iter12*ones(1,length(x_raio5));
ymedia_raio5_iter26 = concentracao_media_raio5_iter26*ones(1,length(x_raio5));
ymedia_raio5_iter38 = concentracao_media_raio5_iter38*ones(1,length(x_raio5));

% radius 15
ymedia_raio15_iter0 = concentracao_media_raio15_iter0*ones(1,length(x_raio15));
ymedia_raio15_iter12 = concentracao_media_raio15_iter12*ones(1,length(x_raio15));
ymedia_raio15_iter26 = concentracao_media_raio15_iter26*ones(1,length(x_raio15));
ymedia_raio15_iter38 = concentracao_media_raio15_iter38*ones(1,length(x_raio15));
 
% radius 45
ymedia_raio45_iter0 = concentracao_media_raio45_iter0*ones(1,length(x_raio45));
ymedia_raio45_iter12 = concentracao_media_raio45_iter12*ones(1,length(x_raio45));
ymedia_raio45_iter26 = concentracao_media_raio45_iter26*ones(1,length(x_raio45));
ymedia_raio45_iter38 = concentracao_media_raio45_iter38*ones(1,length(x_raio45));

 
% radius 60
ymedia_raio60_iter0 = concentracao_media_raio60_iter0*ones(1,length(x_raio60));
ymedia_raio60_iter12 = concentracao_media_raio60_iter12*ones(1,length(x_raio60));
ymedia_raio60_iter26 = concentracao_media_raio60_iter26*ones(1,length(x_raio60));
ymedia_raio60_iter38 = concentracao_media_raio60_iter38*ones(1,length(x_raio60));
 
% radius 80
ymedia_raio80_iter0 = concentracao_media_raio80_iter0*ones(1,length(x_raio80));
ymedia_raio80_iter12 = concentracao_media_raio80_iter12*ones(1,length(x_raio80));
ymedia_raio80_iter26 = concentracao_media_raio80_iter26*ones(1,length(x_raio80));
ymedia_raio80_iter38 = concentracao_media_raio80_iter38*ones(1,length(x_raio80));
 
% radius 100
ymedia_raio100_iter0 = concentracao_media_raio100_iter0*ones(1,length(x_raio100));
ymedia_raio100_iter12 = concentracao_media_raio100_iter12*ones(1,length(x_raio100));
ymedia_raio100_iter26 = concentracao_media_raio100_iter26*ones(1,length(x_raio100));
ymedia_raio100_iter38 = concentracao_media_raio100_iter38*ones(1,length(x_raio100));

%% Interpolation: theta x concentration plotting

% Iteration 0
interp_raio5_iter0 = interp1(arc_length_raio5,concentracao_raio5_iter0,x_raio5);
interp_raio15_iter0 = interp1(arc_length_raio15,concentracao_raio15_iter0,x_raio15);
interp_raio45_iter0 = interp1(arc_length_raio45,concentracao_raio45_iter0,x_raio45);
interp_raio60_iter0 = interp1(arc_length_raio60,concentracao_raio60_iter0,x_raio60);
interp_raio80_iter0 = interp1(arc_length_raio80,concentracao_raio80_iter0,x_raio80);
interp_raio100_iter0 = interp1(arc_length_raio100,concentracao_raio100_iter0,x_raio100);


% Iteration 12
interp_raio5_iter12 = interp1(arc_length_raio5,concentracao_raio5_iter12,x_raio5);
interp_raio15_iter12 = interp1(arc_length_raio15,concentracao_raio15_iter12,x_raio15);
interp_raio45_iter12 = interp1(arc_length_raio45,concentracao_raio45_iter12,x_raio45);
interp_raio60_iter12 = interp1(arc_length_raio60,concentracao_raio60_iter12,x_raio60);
interp_raio80_iter12 = interp1(arc_length_raio80,concentracao_raio80_iter12,x_raio80);
interp_raio100_iter12 = interp1(arc_length_raio100,concentracao_raio100_iter12,x_raio100);


% iteration 26
interp_raio5_iter26 = interp1(arc_length_raio5,concentracao_raio5_iter26,x_raio5);
interp_raio15_iter26 = interp1(arc_length_raio15,concentracao_raio15_iter26,x_raio15);
interp_raio45_iter26 = interp1(arc_length_raio45,concentracao_raio45_iter26,x_raio45);
interp_raio60_iter26 = interp1(arc_length_raio60,concentracao_raio60_iter26,x_raio60);
interp_raio80_iter26 = interp1(arc_length_raio80,concentracao_raio80_iter26,x_raio80);
interp_raio100_iter26 = interp1(arc_length_raio100,concentracao_raio100_iter26,x_raio100);

% iteration 38
interp_raio5_iter38 = interp1(arc_length_raio5,concentracao_raio5_iter38,x_raio5);
interp_raio15_iter38 = interp1(arc_length_raio15,concentracao_raio15_iter38,x_raio15);
interp_raio45_iter38 = interp1(arc_length_raio45,concentracao_raio45_iter38,x_raio45);
interp_raio60_iter38 = interp1(arc_length_raio60,concentracao_raio60_iter38,x_raio60);
interp_raio80_iter38 = interp1(arc_length_raio80,concentracao_raio80_iter38,x_raio80);
interp_raio100_iter38 = interp1(arc_length_raio100,concentracao_raio100_iter38,x_raio100);


%% Plottings concentration x theta + mean concentration - fit subplot!!! 

% Atenção!!! Remover os subplots, salvar uma por uma, colocar as legendas e
% formatar os gráficos em preto. Aumentar espessura de linha também.

EPS = 0.002; % adjust for axes


%% "If" to choose which plots to be generated

iter_escolha = 0;

if iter_escolha == 0
    
    %% Iteration 0

    figure
    % Radius = 5, Iteration = 0
    %subplot(3,2,1)
    pltRaio5Iter0 = plot(x_raio5,interp_raio5_iter0,x_raio5,ymedia_raio5_iter0);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 5, iteração 0');
    axis(gca,[0 max(x_raio5) min(interp_raio5_iter0)-EPS max(interp_raio5_iter0)+EPS]);

    figure
    % Radius = 15, Iteration = 0
    %subplot(3,2,2)
    pltRaio15Iter0 = plot(x_raio15,interp_raio15_iter0,x_raio15,ymedia_raio15_iter0);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 15, iteração 0');
    axis(gca,[0 max(x_raio15) min(interp_raio15_iter0)-EPS max(interp_raio15_iter0)+EPS]);

    figure
    % Radius = 45, Iteration = 0
    %subplot(3,2,3)
    pltRaio45Iter0 = plot(x_raio45,interp_raio45_iter0,x_raio45,ymedia_raio45_iter0);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 45, iteração 0');
    axis(gca,[0 max(x_raio45) min(interp_raio45_iter0)-EPS max(interp_raio45_iter0)+EPS]);

    figure
    % Radius = 60, Iteration = 0
    %subplot(3,2,4)
    pltRaio60Iter0 = plot(x_raio60,interp_raio60_iter0,x_raio60,ymedia_raio60_iter0);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 60, iteração 0');
    axis(gca,[0 max(x_raio60) min(interp_raio60_iter0)-EPS max(interp_raio60_iter0)+EPS]);

    figure
    % Radius = 80, Iteration = 12
    %subplot(3,2,5)
    pltRaio80Iter0 = plot(x_raio80,interp_raio80_iter0,x_raio80,ymedia_raio80_iter0);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 80, iteração 0');
    axis(gca,[0 max(x_raio80) min(interp_raio80_iter0)-EPS max(interp_raio80_iter0)+EPS]);

    figure
    % Radius = 100, Iteration = 0
    %subplot(3,2,6)
    pltRaio100Iter0 = plot(x_raio100,interp_raio100_iter0,x_raio100,ymedia_raio100_iter0);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 100, iteração 0');
    axis(gca,[0 max(x_raio100) min(interp_raio100_iter0)-EPS max(interp_raio100_iter0)+EPS]);

elseif iter_escolha == 12
    %% Iteration 12

    figure
    % Radius = 5, Iteration = 12
    %subplot(3,2,1)
    pltRaio5Iter12 = plot(x_raio5,interp_raio5_iter12,x_raio5,ymedia_raio5_iter12);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 5, iteração 12');
    axis(gca,[0 max(x_raio5) min(interp_raio5_iter12)-EPS max(interp_raio5_iter12)+EPS]);

    figure
    % Radius = 15, Iteration = 12
    %subplot(3,2,2)
    pltRaio15Iter12 = plot(x_raio15,interp_raio15_iter12,x_raio15,ymedia_raio15_iter12);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 15, iteração 12');
    axis(gca,[0 max(x_raio15) min(interp_raio15_iter12)-EPS max(interp_raio15_iter12)+EPS]);

    figure
    % Radius = 45, Iteration = 12
    %subplot(3,2,3)
    pltRaio45Iter12 = plot(x_raio45,interp_raio45_iter12,x_raio45,ymedia_raio45_iter12);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 45, iteração 12');
    axis(gca,[0 max(x_raio45) min(interp_raio45_iter12)-EPS max(interp_raio45_iter12)+EPS]);

    figure
    % Radius = 60, Iteration = 12
    %subplot(3,2,4)
    pltRaio60Iter12 = plot(x_raio60,interp_raio60_iter12,x_raio60,ymedia_raio60_iter12);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 60, iteração 12');
    axis(gca,[0 max(x_raio60) min(interp_raio60_iter12)-EPS max(interp_raio60_iter12)+EPS]);

    figure
    % Radius = 80, Iteration = 12
    %subplot(3,2,5)
    pltRaio80Iter12 = plot(x_raio80,interp_raio80_iter12,x_raio80,ymedia_raio80_iter12);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 80, iteração 12');
    axis(gca,[0 max(x_raio80) min(interp_raio80_iter12)-EPS max(interp_raio80_iter12)+EPS]);

    figure
    % Radius = 100, Iteration = 12
    %subplot(3,2,6)
    pltRaio100Iter12 = plot(x_raio100,interp_raio100_iter12,x_raio100,ymedia_raio100_iter12);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 100, iteração 12');
    axis(gca,[0 max(x_raio100) min(interp_raio100_iter12)-EPS max(interp_raio100_iter12)+EPS]);

elseif iter_escolha == 26
    %% Iteration 26

    figure
    % Radius = 5, Iteration = 26
    %subplot(3,2,1)
    pltRaio5Iter26 = plot(x_raio5,interp_raio5_iter26,x_raio5,ymedia_raio5_iter26);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 5, iteração 26');
    axis(gca,[0 max(x_raio5) min(interp_raio5_iter26)-EPS max(interp_raio5_iter26)+EPS]);

    figure
    % Radius = 15, Iteration = 26
    %subplot(3,2,2)
    pltRaio15Iter26 = plot(x_raio15,interp_raio15_iter26,x_raio15,ymedia_raio15_iter26);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 15, iteração 26');
    axis(gca,[0 max(x_raio15) min(interp_raio15_iter26)-EPS max(interp_raio15_iter26)+EPS]);

    figure
    % Radius = 45, Iteration = 26
    %subplot(3,2,3)
    pltRaio45Iter26 = plot(x_raio45,interp_raio45_iter26,x_raio45,ymedia_raio45_iter26);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 45, iteração 26');
    axis(gca,[0 max(x_raio45) min(interp_raio45_iter26)-EPS max(interp_raio45_iter26)+EPS]);

    figure
    % Radius = 60, Iteration = 26
    %subplot(3,2,4)
    pltRaio60Iter26 = plot(x_raio60,interp_raio60_iter26,x_raio60,ymedia_raio60_iter26);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 60, iteração 26');
    axis(gca,[0 max(x_raio60) min(interp_raio60_iter26)-EPS max(interp_raio60_iter26)+EPS]);

    figure
    % Radius = 80, Iteration = 26
    %subplot(3,2,5)
    pltRaio80Iter26 = plot(x_raio80,interp_raio80_iter26,x_raio80,ymedia_raio80_iter26);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 80, iteração 26');
    axis(gca,[0 max(x_raio80) min(interp_raio80_iter26)-EPS max(interp_raio80_iter26)+EPS]);

    figure
    % Radius = 100, Iteration = 26
    %subplot(3,2,6)
    pltRaio100Iter26 = plot(x_raio100,interp_raio100_iter26,x_raio100,ymedia_raio100_iter26);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 100, iteração 26');
    axis(gca,[0 max(x_raio100) min(interp_raio100_iter26)-EPS max(interp_raio100_iter26)+EPS]);

else 
    %% Iteration 38

    figure
    % Radius = 5, Iteration = 38
    %subplot(3,2,1)
    pltRaio5Iter38 = plot(x_raio5,interp_raio5_iter38,x_raio5,ymedia_raio5_iter38);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 5, iteração 38');
    axis(gca,[0 max(x_raio5) min(interp_raio5_iter38)-EPS max(interp_raio5_iter38)+EPS]);

    figure
    % Radius = 15, Iteration = 38
    %subplot(3,2,2)
    pltRaio15Iter38 = plot(x_raio15,interp_raio15_iter38,x_raio15,ymedia_raio15_iter38);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 15, iteração 38');
    axis(gca,[0 max(x_raio15) min(interp_raio15_iter38)-EPS max(interp_raio15_iter38)+EPS]);

    figure
    % Radius = 45, Iteration = 38
    %subplot(3,2,3)
    pltRaio45Iter38 = plot(x_raio45,interp_raio45_iter38,x_raio45,ymedia_raio45_iter38);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 45, iteração 38');
    axis(gca,[0 max(x_raio45) min(interp_raio45_iter38)-EPS max(interp_raio45_iter38)+EPS]);

    figure
    % Radius = 60, Iteration = 38
    %subplot(3,2,4)
    pltRaio60Iter38 = plot(x_raio60,interp_raio60_iter38,x_raio60,ymedia_raio60_iter38);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 60, iteração 38');
    axis(gca,[0 max(x_raio60) min(interp_raio60_iter38)-EPS max(interp_raio60_iter38)+EPS]);

    figure
    % Radius = 80, Iteration = 38
    %subplot(3,2,5)
    pltRaio80Iter38 = plot(x_raio80,interp_raio80_iter38,x_raio80,ymedia_raio80_iter38);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 80, iteração 38');
    axis(gca,[0 max(x_raio80) min(interp_raio80_iter38)-EPS max(interp_raio80_iter38)+EPS]);

    figure
    % Radius = 100, Iteration = 38
    %subplot(3,2,6)
    pltRaio100Iter38 = plot(x_raio100,interp_raio100_iter38,x_raio100,ymedia_raio100_iter38);
    xlabel(gca,'Pontos na direção \theta');
    ylabel(gca,'c');
    title(gca,'Plot r = 100, iteração 38');
    axis(gca,[0 max(x_raio100) min(interp_raio100_iter38)-EPS max(interp_raio100_iter38)+EPS]);

end

% --- End of "If" for plottings


%% Labelling Iterations 

iter0 = 0;
iter12 = 12;
iter26 = 26;
iter38 = 38;

%% Sample Variance
% Sample variance per iteration follow below:

%% Iteration 0

% --- radius 5
var_amostral_raio5_iter0 = sqrt( (1/n_raio5).*...
    sum(concentracao_raio5_iter0 - concentracao_media_raio5_iter0).^2 );

% --- radius 15
var_amostral_raio15_iter0 = sqrt( (1/n_raio15).*...
    sum(concentracao_raio15_iter0 - concentracao_media_raio15_iter0).^2 );

% --- radius 45
var_amostral_raio45_iter0 = sqrt( (1/n_raio45).*...
    sum(concentracao_raio45_iter0 - concentracao_media_raio45_iter0).^2 );

% --- radius 60
var_amostral_raio60_iter0 = sqrt( (1/n_raio60).*...
    sum(concentracao_raio60_iter0 - concentracao_media_raio60_iter0).^2 );

% --- radius 80
var_amostral_raio80_iter0 = sqrt( (1/n_raio80).*...
    sum(concentracao_raio80_iter0 - concentracao_media_raio80_iter0).^2 );

% --- radius 100
var_amostral_raio100_iter0 = sqrt( (1/n_raio100).*...
    sum(concentracao_raio100_iter0 - concentracao_media_raio100_iter0).^2 );

%% Iteration 12 

% --- radius 5
var_amostral_raio5_iter12 = sqrt( (1/n_raio5).*...
    sum(concentracao_raio5_iter12 - concentracao_media_raio5_iter12).^2 );

% --- radius 15
var_amostral_raio15_iter12 = sqrt( (1/n_raio15).*...
    sum(concentracao_raio15_iter12 - concentracao_media_raio15_iter12).^2 );

% --- radius 45
var_amostral_raio45_iter12 = sqrt( (1/n_raio45).*...
    sum(concentracao_raio45_iter12 - concentracao_media_raio45_iter12).^2 );

% --- radius 60
var_amostral_raio60_iter12 = sqrt( (1/n_raio60).*...
    sum(concentracao_raio60_iter12 - concentracao_media_raio60_iter12).^2 );

% --- radius 80
var_amostral_raio80_iter12 = sqrt( (1/n_raio80).*...
    sum(concentracao_raio80_iter12 - concentracao_media_raio80_iter12).^2 );

% --- radius 100
var_amostral_raio100_iter12 = sqrt( (1/n_raio100).*...
    sum(concentracao_raio100_iter12 - concentracao_media_raio100_iter12).^2 );


%% Iteration 26 

% --- radius 5
var_amostral_raio5_iter26 = sqrt( (1/n_raio5).*...
    sum(concentracao_raio5_iter26 - concentracao_media_raio5_iter26).^2 );

% --- radius 15
var_amostral_raio15_iter26 = sqrt( (1/n_raio15).*...
    sum(concentracao_raio15_iter26 - concentracao_media_raio15_iter26).^2 );

% --- radius 45
var_amostral_raio45_iter26 = sqrt( (1/n_raio45).*...
    sum(concentracao_raio45_iter26 - concentracao_media_raio45_iter26).^2 );

% --- radius 60
var_amostral_raio60_iter26 = sqrt( (1/n_raio60).*...
    sum(concentracao_raio60_iter26 - concentracao_media_raio60_iter26).^2 );

% --- radius 80
var_amostral_raio80_iter26 = sqrt( (1/n_raio80).*...
    sum(concentracao_raio80_iter26 - concentracao_media_raio80_iter26).^2 );

% --- radius 100
var_amostral_raio100_iter26 = sqrt( (1/n_raio100).*...
    sum(concentracao_raio100_iter26 - concentracao_media_raio100_iter26).^2 );


%% Iteration 38 

% --- radius 5
var_amostral_raio5_iter38 = sqrt( (1/n_raio5).*...
    sum(concentracao_raio5_iter38 - concentracao_media_raio5_iter38).^2 );

% --- radius 15
var_amostral_raio15_iter38 = sqrt( (1/n_raio15).*...
    sum(concentracao_raio15_iter38 - concentracao_media_raio15_iter38).^2 );

% --- radius 45
var_amostral_raio45_iter38 = sqrt( (1/n_raio45).*...
    sum(concentracao_raio45_iter38 - concentracao_media_raio45_iter38).^2 );

% --- radius 60
var_amostral_raio60_iter38 = sqrt( (1/n_raio60).*...
    sum(concentracao_raio60_iter38 - concentracao_media_raio60_iter38).^2 );

% --- radius 80
var_amostral_raio80_iter38 = sqrt( (1/n_raio80).*...
    sum(concentracao_raio80_iter38 - concentracao_media_raio80_iter38).^2 );

% --- radius 100
var_amostral_raio100_iter38 = sqrt( (1/n_raio100).*...
    sum(concentracao_raio100_iter38 - concentracao_media_raio100_iter38).^2 );

% ------ End of Sample Variance session ------

%% Plottings: iteration x variance (Logarithmic Scale)
% Series of plottings follow below:

%% "If" to choose which variance plots to be generated

iter_escolha_variancia_raio = 5;

% raio = 5
if iter_escolha_variancia_raio == 5;
    figure
    pltVarianceRaio5 = semilogy(iter0,var_amostral_raio5_iter0,iter12,var_amostral_raio5_iter12,...
        iter26,var_amostral_raio5_iter26,iter38,var_amostral_raio5_iter38);

% raio = 15
elseif iter_escolha_variancia_raio == 15    
    figure
    pltVarianceRaio15 = semilogy(iter0,var_amostral_raio15_iter0,iter12,var_amostral_raio15_iter12,...
        iter26,var_amostral_raio15_iter26,iter38,var_amostral_raio15_iter38);

% raio = 45
elseif iter_escolha_variancia_raio == 45    
    figure
    pltVarianceRaio45 = semilogy(iter0,var_amostral_raio45_iter0,iter12,var_amostral_raio45_iter12,...
        iter26,var_amostral_raio45_iter12,iter38,var_amostral_raio45_iter38);
 
% raio = 60
elseif iter_escolha_variancia_raio == 60
    figure
    pltVarianceRaio60 = semilogy(iter0,var_amostral_raio60_iter0,iter12,var_amostral_raio60_iter12,...
        iter26,var_amostral_raio60_iter12,iter38,var_amostral_raio60_iter38);

% raio = 80
elseif iter_escolha_variancia_raio == 80
    figure
    pltVarianceRaio80 = semilogy(iter0,var_amostral_raio80_iter0,iter12,var_amostral_raio80_iter12,...
        iter26,var_amostral_raio80_iter12,iter38,var_amostral_raio80_iter38);

% raio = 100
else 
    figure
    pltVarianceRaio100 = semilogy(iter0,var_amostral_raio100_iter0,iter12,var_amostral_raio100_iter12,...
        iter26,var_amostral_raio100_iter12,iter38,var_amostral_raio100_iter38);
    
end
