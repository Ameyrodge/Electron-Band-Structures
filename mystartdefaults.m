%% START DEFAULT COMMANDS OF *.m SCRIPTS
close all;       % Close all figures if any
clear all;       % Clear all variables/functions in memory
clc;             % Clear screen in the command window
#addpath('~/Nextcloud/Documents/bin');% Add a directory to search path 
addpath('/Users/Shared/Files From d.localized/QuanTEEM/UBFC/Numerical Methods/octave_files')
%% THE SEVEN EXACT DEFINING CONSTANTS OF THE SI UNIT SYSTEM (2019 UPDATE)
%  Link to complete list of updated SI fundamental physical constants:
%  https://physics.nist.gov/cuu/Constants/Table/allascii.txt
hyperfine = 9192631770;       % Hyperfine transition frequency of Cs_133 [Hz]
celeritas = 299792458;        % Speed of light in vacuum [m/s]
Planck    = 6.62607015E-34;   % Planck's constant [Js] 
qel       = 1.602176634E-19;  % Elementary charge (Absolute value of electron charge) [C]
kB        = 1.380640E-23;     % Boltzmann's constant [J/K]
Avogadro  = 6.02214076E23;    % Avogadro's constant [1/mole]
kcd       = 683;              % Luminous efficacy of 540 THz radiation [candela=lumen/Watt]
                              % Green light @ 555.016 nm = maximum possible luminous efficacy.
                              % Originally = peak sensitivity of "average" human eye.

%% PHYSICAL CONSTANTS OF ELECTROMAGNETISM
epsilon0  = 8.8541878128E-12; % Electrical constant (vacuum dielectric permittivity) [F/m]
mu0       = 1.25663706212E-6; % Magnetic constant   (vacuum magnetic permeability) [N/A^2]
                              % close to 4*pi*1.E-7 =1.25663706143E-6
Klitzing  = 25812.80745;      % von Klitzing's constant = Planck/qel^2 [Ohm]
                              % respecting significant digits 

%% UNITS USED IN QUANTUM PHYSICS
hbar      = 1.054571817E-34;                    % Reduced Planck's constant = Planck/(2*pi)
                                                % respecting significant digits [Js]
Angstroem = 1.0E-10;                            % Angström [m]
Dalton    = 1.66053906660E-27;                  % Atomic mass unit [kg] = 1 Dalton
                                                % = mass of Carbon_12 atom / 12
elm       = 9.1093837015E-31;                   % Electron mass [kg]
nem       = 1.67492749804E-27;                  % Neutron mass [kg]
prm       = 1.67262192369E-27;                  % Proton mass [kg]
elecint   = qel^2/(4*pi*epsilon0);              % Scale of electron-electron interaction [N*m^2]
Bohr      = hbar^2/(elm*elecint)/Angstroem;     % Bohr radius [Angström]
Rydberg   = (elm/(2*hbar^2)) * (elecint)^2/qel; % Rydberg [eV]
Hartree   = 2*Rydberg;                          % Hartree [eV]

%% WARNING
% The rest of this file can only be executed when sourced in a *.m file 
% that includes plotting commands. 

%% DEFAULT PLOT CONFIGURATION
% get(groot,'factory')           % List factory-defined plot configurations
% get(groot,'factoryObjectType') % List factory-defined properties for a specific object.
                                 % Examples of 'ObjectType' are: 
                                 % 'Axes', 'Figure', 'Image', 'Line','Surface', 'Text',
                                 % 'uicontrol' (user interface), etc.

% groot is the "handle index" (identifier) of the object 
% that is "parent" of all plot objects of the session.

% Function "get" and "set" are not case sensitive. Upper cases for readability.

%% DEFAULT FIGURE WINDOW POSITION
% default_figure_position=get(groot,'DefaultFigurePosition'); % Factory default position
golden_ratio=(1+sqrt(5))/2;
default_screen_size=get(groot,'ScreenSize'); % Get screen dimensions [pixels]
default_figure_position(4)=default_screen_size(4)/3; % Height [pixels]
default_figure_position(3)=default_figure_position(4)*golden_ratio; % Width [pixels]
set(groot,'DefaultFigurePosition',[default_screen_size(3)/2-default_figure_position(3)/2,...
                                   default_screen_size(4)/2-default_figure_position(4)/2,...
                                   default_figure_position(3), ...
                                   default_figure_position(4)]);

%% ADJUST LINE & CHARACTER PROPERTIES FOR VIDEOPROJECTION IN CLASSROOM
set(groot,'defaultLineLineWidth',1.5);  
set(groot,'defaultAxesFontSize',12); 
set(groot,'defaultAxesFontWeight','normal');
set(groot,'defaultAxesLineWidth',1.5); 
set(groot,'defaultAxesXaxisLocation','bottom'); % Location of abcissae ticks & labels 
                                                % Possible values:'top','bottom',
                                                % left', 'right', 'origin'.
                                                % 'origin' puts abcissae labels
                                                % along y=0 line

%% DEFAULT COLOR ORDER FOR COLORING CURVES
% default_colors = get(gca,'colororder'); % gca="get current axes" identifier 

% Color order of successive curves must be defined by a 7x3 matrix
% of values between 0 and 1
 
default_colors=[ % Default values in last versions of Matlab/Octave
                 0.0000   0.4470   0.7410 % (1) Tropical Blue
                 0.8500   0.3250   0.0980 % (2) Deep orange
                 0.9290   0.6940   0.1250 % (3) Deep yellow
                 0.4940   0.1840   0.5560 % (4) Violet
                 0.4660   0.6740   0.1880 % (5) Grass green
                 0.3010   0.7450   0.9330 % (6) Azure blue
                 0.6350   0.0780   0.1840 % (7) Sienna
               ];

%% OTHER COLOR ORDERS

brg_colors=[ % Reorder default values to start by blue, red, green
             0.0000   0.4470   0.7410 % (1) Tropical Blue
             0.6350   0.0780   0.1840 % (2) Sienna
             0.4660   0.6740   0.1880 % (3) Grass green
             0.3010   0.7450   0.9330 % (4) Azure blue
             0.8500   0.3250   0.0980 % (5) Deep orange
             0.9290   0.6940   0.1250 % (6) Deep yellow
             0.4940   0.1840   0.5560 % (7) Violet
           ];

crude_colors=[ % Unsoftened colors
               0    0    0      % (1) Black		    
               1    0    0      % (2) Red
               0    0    1      % (3) Blue
               0    0.5  0      % (4) Dark Green
               0.9  0.5  0.1    % (5) Orange
               0    0.75 0.75   % (6) Turquoise
               0.5  0.5  0.5    % (7) Grey
             ];

set(groot,'defaultAxesColorOrder',brg_colors); 

%% USER INTERFACE (UI) DEFAULTS
set(groot,'defaultUicontrolFontSize',13);

%% TRACING HORIZONTAL AND VERTICAL LINES ACROSS THE PLOT WINDOW

% yline(val) traces horizontal line at y(x)=yval using current xlim
yline = @(yval, varargin) line(xlim, [yval yval], varargin{:}); 

% xline(val) traces vertical   line at    x=xval using current ylim
xline = @(xval, varargin) line([xval xval], ylim, varargin{:}); 
 
%% END DEFAULT COMMANDS OF *.m SCRIPTS
 