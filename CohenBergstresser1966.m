%% INPUT PARAMETERS FOR COMPUTING
%% THE ELECTRON BAND STRUCTURE OF FCC SEMICONDUCTORS
%% USING EMPIRICAL PSEUDOPOTENTIALS 

semiconductor='Si'
tol=1e-12


%% PARAMETERS IN RECIPROCAL SPACE

BZstep=0.02 % Step along path in BZ

cutoff=21   % Deal only with |G|^2 < cutoff [2*pi/spacing]^2
            % in Hamiltonian

Gs_max=11   % |G|^2 of highest non zero Fourier coefficients in
            % expanding potential [2*pi/spacing]^2

%% FOR COMPUTING THE DISPERSION RELATION
%% PATH IN BRILLOUIN ZONE OF FCC LATTICE

if dispersion_relation
    nband=16    % Number of bands to be stored in output file
    
    qs(1:3,1)= [0.5 0.5 0.5]';    qs_str{1} = 'L';      % start point 1
    qe(1:3,1)= [0 0 0]';          qe_str{1} = '\Gamma'; % end   point 1
    
    qs(1:3,2)= [0 0 0]';          qs_str{2} = '\Gamma'; % start point 2
    qe(1:3,2)= [1 0 0]';          qe_str{2} = 'X';      % end   point 2
    
    qs(1:3,3)= [1 1 0]';          qs_str{3} = 'X';      % start point 3
    qe(1:3,3)= [0.75 0.75 0]';    qe_str{3} = 'K';      % end   point 3
    
    qs(1:3,4)= [0.75 0.75 0]';    qs_str{4} = 'K';      % start point 4
    qe(1:3,4)= [0 0 0]';          qe_str{4} = '\Gamma'; % end   point 4
end

%% FOR COMPUTING DOS
%% MONKHORST-PACK SAMPLING PARAMETERS OF BZ VOLUME 

if compute_dos
    movie=true     % Show movie of DOS construction
    mp=21          % mp^3 q vectors will sample the BZ volume if no
                   % symmetry is exploited
    foldsym=1      % th-fold symmetry to be exploited => 
                   % mp^3/foldsym q vectors will sample the BZ volume
    Energy_min=-14 % [eV]
    Energy_max=6   % [eV]
    E_step=0.005   % [eV]
   
    FWHM=E_step*5  % Full Width at Half-Maximum of the
                   % representation of the Dirac delta function
    Dirac='Gaussian' % Representation of Dirac delta function 
                     % Possible choices: 'Gaussian' or 'Lorentzian'
    
    Energy=[Energy_min:E_step:Energy_max];
    DOS = zeros(length(Energy),1);  % Initialize DOS
end

%% DATA OF EMPIRICAL PSEUDOPOTENTIALS
%% FOR COMPUTING ELECTRON ENERGY BANDS
%% OF 14 FCC SEMICONDUCTORS
%
%  from M.L. Cohen & T.K. Bergstresser, 
%       Phys. Rev. vol.141, p.789 (1966)

material_list=['Si';'Ge';'Sn';'GaP';'GaAs';'AlSb';'InP';'GaSb'; ...
               'InAs';'InSb';'ZnS';'ZnSe';'ZnTe';'CdTe','Empty lattice'];

% NB: In Matlab/Octave, a matrix of strings is unexploitable
%     => Conversion of matrix of strings into an exploitable
%        so-called "cell array" whose index is then between braces

materials=cellstr(material_list); 

%% MATERIAL IDENTIFIER

m = find(strcmp(materials,semiconductor));

if isempty(m)
    error('Semiconductor material not recognized.')
end

%% DIRECT LATTICE UNIT VECTORS AND ATOMIC POSITIONS
%% OF FCC SEMICONDUCTORS

a(1:3,1)= [0.5 0.5 0.0]' ;    % Direct lattice unit vector 1
a(1:3,2)= [0.0 0.5 0.5]' ;    % Direct lattice unit vector 2
a(1:3,3)= [0.5 0.0 0.5]' ;    % direct lattice unit vector 3
cell_volume = a(:,1)' * cross(a(:,2),a(:,3));

tau(1:3,1) = [ 0.125  0.125  0.125]' ; % Position of atom 1 in primitive cell
tau(1:3,2) = [-0.125 -0.125 -0.125]' ; % Position of atom 2 in primitive cell

%% LATTICE SPACINGS [Angström]

ls(1:14)= [5.43 5.66 6.49 5.44 5.64 6.13 5.86 6.12 6.04 6.48 5.41 ...
           5.65 6.07 6.41];

%% PSEUDOPOTENTIAL FORM FACTORS [Rydberg]

% Remark: ff(i,1) = V0 = V_{G=0} for material i
%                 = constant adjusting zero of potential

%%            V0  VS3  VS8  VS11 VA3  VA4  VA11
ff(1,:) = [ 0.00 -0.21 0.04 0.08 0.00 0.00 0.00]; % Si
ff(2,:) = [ 0.00 -0.23 0.01 0.06 0.00 0.00 0.00]; % Ge
ff(3,:) = [ 0.00 -0.20 0.00 0.04 0.00 0.00 0.00]; % Sn
ff(4,:) = [ 0.00 -0.22 0.03 0.07 0.12 0.07 0.02]; % GaP
ff(5,:) = [ 0.00 -0.23 0.01 0.06 0.07 0.05 0.01]; % GaAs
ff(6,:) = [ 0.00 -0.21 0.02 0.06 0.06 0.04 0.02]; % AlSb
ff(7,:) = [ 0.00 -0.23 0.01 0.06 0.07 0.05 0.01]; % InP
ff(8,:) = [ 0.00 -0.22 0.00 0.05 0.06 0.05 0.01]; % GaSb
ff(9,:) = [ 0.00 -0.22 0.00 0.05 0.08 0.05 0.03]; % InAs
ff(10,:)= [ 0.00 -0.20 0.00 0.04 0.06 0.05 0.01]; % InSb
ff(11,:)= [ 0.00 -0.22 0.03 0.07 0.24 0.14 0.04]; % ZnS
ff(12,:)= [ 0.00 -0.23 0.01 0.06 0.18 0.12 0.03]; % ZnSe
ff(13,:)= [ 0.00 -0.22 0.00 0.05 0.13 0.10 0.01]; % ZnTe
ff(14,:)= [ 0.00 -0.20 0.00 0.04 0.15 0.09 0.04]; % CdTe

%% CONSTANT ADJUSTING THE ZERO OF THE ENERGY SCALE
%% TO THE TOP OF THE VALENCE BAND [Rydberg]
%  For each material m, the adjusting constant
%  is deduced from the results of a prelimninary run 
%  with ff(1:m)=0

ff(1,1) = -7.704369581925116E-01 ; % Si

% ff(1:2:14) to be completed

%% WORK FUNCTIONS [eV]
% Data from various references
% Uncertainty may be large (0.1 to 1 eV)
% Sources of uncertainty: doping;
%                         quality of crystal;
%                         quality of sample surface.

wf(1)  = 4.85; % Si
wf(2)  = 4.75; % Ge
wf(3)  = 4.42; % Sn
wf(4)  = 4.34; % GaP
wf(5)  = 4.69; % GaAS
wf(6)  = 4.46; % AlSb
wf(7)  = 4.65; % InP
wf(8)  = 4.45; % GaSb  (4.1 to 4.8 eV)
wf(9)  = 4.95; % InAs
wf(10) = 4.57; % InSb
wf(11) = 7.00; % ZnS  (7 to 8 eV)
wf(12) = 5.11; % ZnSe
wf(13) = 5.80; % ZnTe (5.3 to 5.8 eV)
wf(14) = 5.70; % CdTe
