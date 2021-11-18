function [SIM_OUTPUTS]=tbcosim_uncond_func(nx,ny,model_number,scale_factor,b,covariance,seed)
% Conditional co-simulation of cross-correlated Gaussian random fields via the turning bands method
%--------------------------------------------------------------------------------------------------
% USE: [SIM_OUTPUTS]=tbcosim(nx,ny,model_number,scale_factor,b,covariance,seed)
%
% Sample call (see input definitions below), e.g.: 
% for two Gaussian zero-mean variable fields of dimension 50 (nx) by 
% 75 (ny) pixels with fields following an exponential covariance function 
% with a range (scale parameter) of 5 pixels, both with unit variance and 
% a correlation coefficient of 0.5, and a random seed input of 1 (note 
% using a fixed seed with return the same fields, to randomize, pass 
% different seeds):
%
% [SIM_OUTPUTS]=tbcosim_uncond_func(50,75,2,5,1,[1 0.5; 0.5 1],1);
% 
% INPUTS:
%
%   LOCATIONS FOR CO-SIMULATION
%   ---------------------------
%   nx,ny        :          number of grid nodes along x, y directions
%
%   COREGIONALIZATION MODEL
%   -----------------------
%   model        : covariance model for the Gaussian random fields (nst * 7 matrix, where nst is the number of nested structures)
%                  Each row refers to a nested structure and is codified as: [type, scale factors, angles]
%                  There are three scale factors (along the rotated y, x and z axes) and three angles
%                  to define the coordinate rotation (azimuth, dip and plunge), see Deutsch and Journel, 1992, p. 25
%                  More info. on the covariance functions are in the original TBSIM paper (Emery and Lantuejoul, 2006).
% Steve Margulis:  In this version of the code there are not multiple (nested) models so that only one of those below is used.           
%                  Available types (if not indicated parameter b = 1), where each function is pre-multiplied by the 
%                  variance at lag-0:
%                    1: spherical: C(r) ~ 1-3/2*(r/a) + 1/2*(r/a)^3 for r<=a; 0 for r>a
%                    2: exponential: C(r) ~ exp(-r/a)
%                    3: gamma (parameter b > 0): C(r) ~ (1+r/a)^(-b)
%                    4: stable (parameter b between 0 and 2): C(r) ~ exp((-r/a)^b)
%                    5: cubic: C(r) ~ (1-7*r^2/a^2 +35/4*r^3/a^3 - 7/2*r^5/a^5 + 3/4*r^7/a^7)
%                    6: Gaussian: C(r) ~ exp(-(r/a)^2)  
%                    7: cardinal sine: C(r) ~ a/r*sin(r/a)
%                    8: Bessel J (parameter b > 0): C(r) ~ 2^b*GAMMA(b+1)*(r/a)^(-b) * J_b(r/a)
%                    9: Bessel K (parameter b > 0): C(r) ~ 1/GAMMA(b)*2^(1-b)*(r/a)^(b) * K_b(r/a)
%                   10: generalized Cauchy (parameter b): C(r) ~ (1 + (r/a)^2)^(-b)
%                   11: exponential sine: C(r) ~ sin(pi/2*exp(-r/a))
%                   12: linear generalized covariance: C(r) ~ -r/a
%                   13: power generalized covariance (parameter b > 0): C(r) ~ (-1)^(k+1)(r/a)^b
%                   14: mixed power generalized covariance (parameter b between 0 and 2): C(r) ~ 1/log(r/a)*(1 - (r/a)^b)
%                   15: spline generalized covariance (parameter b = even integer): (-1)^(k+1)(r/a)^b*log(r/a)
%
%   Computed parameters for coregionalization model (specified as inputs in old version of code):
%   cc           : sills or slopes of the nested structures (nst * nvar^2 matrix)
%   b            : third parameters for covariance definition (nst * 1 vector) (used for covariance types 3,4,8,9,10,13,14 and 15)
%   nugget       : nugget effect variance-covariance matrix (1 * nvar^2 vector)
%
%   SIMULATION PARAMETERS
%   ---------------------
%   seed         : seed number for generating random values
%
%   OUTPUT OPTIONS
%   --------------
%   nbdecimal    : number of decimals for the values in output file
%
% OUTPUTS:
%
% SIM_OUTPUTS    : array of dimension ny, nx, nvar with the realization of the random field
%
%-----------------------------------------------------------------------------------------------------------------------------------
% This program is based on the computer code TBSIM.M formerly published in Computers & Geosciences:
%    X. Emery, C. Lantuéjoul, 2006. TBSIM: A computer program for conditional simulation of three-dimensional Gaussian random fields
%    via the turning bands method. Computers & Geosciences 32 (10), 1615-1628.
%
% It uses the following subroutines:
%   setrot.m           : set up matrix for rotation and reduction of coordinates
%   tbmain.m           : main routine for simulation along the lines
%   vdc.m              : create equidistributed directions on the sphere (Van der Corput sequence)
%
% It also requires the availability of the Matlab Statistics Toolbox, for random number generation:
%   betarnd.m          : simulate beta values
%   gamrnd.m           : simulate gamma values
%   rand.m             : simulate uniform values
%   randn.m            : simulate normal values
%-----------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------
% Author: Xavier Emery
% Paper: "A turning bands program for conditional co-simulation of cross-correlated Gaussian random fields"
%----------------------------------------------------------------------------------------------------------
% Revised by: Steve Margulis
% The original TBCOSIM code (which was more general for conditioning on 
% data), has been stripped down to use solely for co-simulation of 
% cross-correlated Gaussian random fields. Many of the original parameters 
% are no longer used and the code has been simplified and reformulated as 
% a function (rather than writing data to output file). The simplifying
% cases/assumptions used by this version of the TBCOSIM code include:
%
% -- The original function was capable of simulating three-dimensional 
% anisotropic fields. Instead parameters are hard-wired to simulate 
% two-dimensional isotropic fields.
% -- The original code include the capability of nested set of variogram
% models that are capable of representing multiple correlation length
% scales in a given field. This version assumes only a single model.
% -- The original code allowed for specification of a nugget effect. Here
% the nugget is assumed zero.
% -- The original code could generate more than one realization. This 
% version generates a single realization. To generate an ensemble of
% realizations the code can be called multiple times.

% Setup model: (this selects model type and assume isotropic with same
% scale factor in each dimension.
model = [model_number scale_factor scale_factor scale_factor 0 0 0];
% Determine number of variable fields to generate
nvar = size(covariance,1);
% Assume no nugget effect
nugget = zeros([nvar nvar]);
% Covariance (in single row format)
cc=covariance(:)'; 

%% SM: Hard-wired these parameters to inside function
% grid parameters
% x0,y0,z0: minimum grid coordinates along x, y and z directions
x0 = 1; y0 = 1; z0 = 0;
% nz: number of grid nodes along z direction (limits to 2D fields)
nz = 1;
% dx,dy,dz: grid meshes along x, y and z directions
dx = 1; dy = 1; dz = 1;
% nd: block discretization along x, y and z directions (1 * 3 vector) ([1 1 1] for point-support simulation)
nd = [1 1 1];

% nrealiz: number of realizations to draw
nrealiz=1;

% nbdecimal: number of decimals for the values in output file
nbdecimal = 4;
ntok = 1000; % ntok: maximum number of points to keep in memory and simulate simultaneously (optional)
% This number defines how many locations are projected onto the lines at each step of the simulation
% nlines: number of lines to use for each nested structure (nst * 1 vector)
% Using criteria from Emery and Lantuejoul (2006)
nlines =max(100,20*cc(1,1)*max(model(1,2:4)')/dx);

% Pre-allocate output matrix
SIM_OUTPUTS=NaN(ny,nx,nvar);

%%
% Define default values
%----------------------
warning('off','all');
max_intervalnumber = 1e4;
block = [dx dy dz];      % block dimensions
ng = prod(nd);           % number of points discretizing the block
nst = size(model,1);     % number of nested structures
nrealiz = nvar*nrealiz;  % number of realizations
b = b(:);
nlines = nlines(:);
order = -ones(nst,1);
if (ng == 1), block = [0 0 0]; end

% Define the grid that discretizes the block
%-------------------------------------------
t2 = [];
for i = 1:3
    nl = prod(nd(1:i-1));
    nr = prod(nd(i+1:3));
    t1  = [0.5*(1/nd(i)-1):1/nd(i):0.5*(1-1/nd(i))]';
    t2 = [t2,kron(ones(nl,1),kron(t1,ones(nr,1)))];
end
grid = t2.*(ones(ng,1)*block);

% Check input
%------------
if isempty(ntok), ntok = 1000; end

if length(nlines) ~= nst, nlines = nlines(1)*ones(nst,1); end

if nvar > floor(nvar), error('The number of columns in the sill matrix (cc) is inconsistent'); end

sill = zeros(nvar,nvar,nst);
A1 = zeros(nvar,nvar,nst);
for i = 1:nst
    sill(:,:,i) = reshape(cc(i,:),nvar,nvar);
    if max(abs(sill(:,:,i)-sill(:,:,i)'))>100*eps, error(['The sill matrix for structure nº',num2str(i),' is not symmetric']); end
    [eigenvectors,eigenvalues] = eig(sill(:,:,i));
    A1(:,:,i) = sqrt(eigenvalues)*eigenvectors';
    if min(diag(eigenvalues))<0, error(['The sill matrix for structure nº',num2str(i),' is not positive semi-definite']); end
end

sillnugget = reshape(nugget,nvar,nvar);
if max(abs(sillnugget-sillnugget'))>100*eps, error(['The sill matrix for the nugget effect is not symmetric']); end
[eigenvectors,eigenvalues] = eig(sillnugget);
A0 = sqrt(eigenvalues)*eigenvectors';
if min(diag(eigenvalues))<0, error(['The sill matrix for the nugget effect is not positive semi-definite']); end

for i = 1:nst
    
    if (model(i,1) == 3) % Gamma model
        if (b(i) <= 0)
            error('The parameter b of the gamma model must be positive');
        end
        
    elseif (model(i,1) == 4) % Stable model
        if (b(i) > 2)
            error('The parameter b of the stable model must be less than 2');
        elseif (b(i) <= 0)
            error('The parameter b of the stable model must be positive');
        elseif (b(i) == 2) % Gaussian model
            model(i,1) = 6;
        elseif (b(i) == 1) % Exponential model
            model(i,1) = 2;
        elseif (b(i) > 1) % Simulation via spectral method
            model(i,1) = 4.5;
        end
        
    elseif (model(i,1) == 8) % Bessel-J model
        if (b(i) < 0.5)
            error('The parameter b of the Bessel-J model must be greater than 0.5');
        elseif (b(i) == 0.5) % Cardinal sine model
            model(i,1) = 7;
        end
        
    elseif (model(i,1) == 9) % Bessel-K model
        if (b(i) <= 0)
            error('The parameter b of the Bessel-K model must be positive');
        elseif (b(i) == 0.5) % Exponential model
            model(i,1) = 2;
        elseif (b(i) > 0.5) % Simulation via spectral method
            model(i,1) = 9.5;
        end
        
    elseif (model(i,1) == 10) % Generalized Cauchy model
        if (b(i) <= 0)
            error('The parameter b of the generalized Cauchy model must be positive');
        end
        
    elseif (model(i,1) == 12) % Linear model
        order(i) = 0;
        
    elseif (model(i,1) == 13) % Power model
        if (b(i) <= 0)
            error('The exponent b of the power model must be positive');
        end
        if (b(i)/2 == floor(b(i)/2))
            model(i,1) = 13.1;
        elseif (b(i) > 1)
            model(i,1) = 13.5;
        end
        order(i) = ceil(b(i)/2)-1;
        
    elseif (model(i,1) == 14) % Mixed power model
        if (b(i) <= 0)
            error('The exponent b of the mixed power model must be positive');
        end
        if (b(i) > 2)
            error('The exponent b of the mixed power model must be less than 2');
        end
        order(i) = 0;
        
    elseif (model(i,1) == 15) % Spline model
        if (b(i) < 0)
            error('The exponent b of the spline model must be positive');
        end
        if (b(i)/2 ~= floor(b(i)/2))
            error('The exponent b of the spline model must an even integer');
        end
        order(i) = b(i)/2;
    end
    
end

%SM: Added to replace original code
m0=0;

% Create seed numbers
%--------------------

rand('state',seed); %#ok<RAND>
randn('state',seed);
seed_vdc = ceil(1e7*rand(1,nst));
seed_line = ceil(1e7*rand(nst*nrealiz,max(nlines)));


% Extremal coordinates to simulate in the original referential
%----------------------------------------------------------------------------------------
    radius = [0 0 0];
    minx = x0-block(1)/2;
    miny = y0-block(2)/2;
    minz = z0-block(3)/2;
    maxx = x0+(nx-1)*dx+block(1)/2;
    maxy = y0+(ny-1)*dy+block(2)/2;
    maxz = z0+(nz-1)*dz+block(3)/2;
    
extreme_coord = [minx miny minz; minx miny maxz; minx maxy minz; minx maxy maxz; ...
    maxx miny minz; maxx miny maxz; maxx maxy minz; maxx maxy maxz];


%--------------------------------------------------------------------------------------------

% PREPARE THE LINES
%------------------
% Initialization
all_lines = [];
all_offset = [];
all_r = [];
all_phi = [];
all_theta = [];
valid_lines = [];


for i = 1:nst
    % Line generation (Van der Corput sequence)
    %------------------------------------------
    V = vdc(nlines(i),nrealiz,seed_vdc(i));
    R = setrot(model,i);
    model_rotationmatrix(:,:,i) = R;
    lines = V*R';
    all_lines(1:nrealiz*nlines(i),:,i) = lines;
    
    % Spherical covariance model?
    %----------------------------
    if (model(i,1) == 1)
        x = extreme_coord*lines';
        sigma(i) = 2*sqrt(3./nlines(i));
        all_offset(1:nrealiz*nlines(i),i) = min(x)' - rand(nrealiz*nlines(i),1);
        
        % Identify the lines for which the scale factor is too small
        % The simulation along these lines will be replaced by a nugget effect
        interval_number = max(x)'-min(x)';
        nottoosmall = (interval_number < max_intervalnumber);
        valid_lines(1:nrealiz*nlines(i),i) = nottoosmall;     
        
        % Exponential covariance model?
        %------------------------------
    elseif (model(i,1) == 2)
        v = randn(6,nlines(i)*nrealiz);
        w = 3*rand(1,nlines(i)*nrealiz);
        v(5,:) = v(5,:).*(w>1);
        v(6,:) = v(6,:).*(w>1);
        G = 0.5*sum(v.^2)' * ones(1,3);
        lines = lines./G;
        x = extreme_coord*lines';
        sigma(i) = 2*sqrt(3./nlines(i));
        all_lines(1:nrealiz*nlines(i),:,i) = lines;
        all_offset(1:nrealiz*nlines(i),i) = min(x)' - rand(nrealiz*nlines(i),1);
        
        % Identify the lines for which the scale factor is too small
        % The simulation along these lines will be replaced by a nugget effect
        interval_number = max(x)'-min(x)';
        nottoosmall = (interval_number < max_intervalnumber);
        valid_lines(1:nrealiz*nlines(i),i) = nottoosmall;
        
        % Gamma covariance model?
        %------------------------
    elseif (model(i,1) == 3)
        v = randn(6,nlines(i)*nrealiz);
        w = 3*rand(1,nlines(i)*nrealiz);
        v(5,:) = v(5,:).*(w>1);
        v(6,:) = v(6,:).*(w>1);
        G = 0.5*(sum(v.^2)'./gamrnd(b(i),1,nlines(i)*nrealiz,1)) * ones(1,3);
        lines = lines./G;
        x = extreme_coord*lines';
        sigma(i) = 2*sqrt(3./nlines(i));
        all_lines(1:nrealiz*nlines(i),:,i) = lines;
        all_offset(1:nrealiz*nlines(i),i) = min(x)' - rand(nrealiz*nlines(i),1);
        
        % Identify the lines for which the scale factor is too small
        % The simulation along these lines will be replaced by a nugget effect
        interval_number = max(x)'-min(x)';
        nottoosmall = (interval_number < max_intervalnumber);
        valid_lines(1:nrealiz*nlines(i),i) = nottoosmall;
        
        % Stable covariance model with parameter <1?
        %-------------------------------------------
    elseif (model(i,1) == 4)
        v = randn(6,nlines(i)*nrealiz);
        w = 3*rand(1,nlines(i)*nrealiz);
        v(5,:) = v(5,:).*(w>1);
        v(6,:) = v(6,:).*(w>1);
        e = -log(rand(nlines(i)*nrealiz,1));
        f = pi*rand(nlines(i)*nrealiz,1)-pi/2;
        stable = abs( sin(b(i)*(f-pi/2))./(cos(f).^(1./b(i))).*(cos(f-b(i)*(f-pi/2))./e).^((1-b(i))./b(i)) );
        G = 0.5*(sum(v.^2)'./stable) * ones(1,3);
        lines = lines./G;
        x = extreme_coord*lines';
        sigma(i) = 2*sqrt(3./nlines(i));
        all_lines(1:nrealiz*nlines(i),:,i) = lines;
        all_offset(1:nrealiz*nlines(i),i) = min(x)' - rand(nrealiz*nlines(i),1);
        
        % Identify the lines for which the scale factor is too small
        % The simulation along these lines will be replaced by a nugget effect
        interval_number = max(x)'-min(x)';
        nottoosmall = (interval_number < max_intervalnumber);
        valid_lines(1:nrealiz*nlines(i),i) = nottoosmall;
        
        % Stable model with parameter >1?
        %--------------------------------
    elseif (model(i,1) == 4.5)
        sigma(i) = sqrt(2./nlines(i));
        t = randn(nlines(i),nrealiz).^2 + randn(nlines(i),nrealiz).^2 + randn(nlines(i),nrealiz).^2;
        all_r(1:nlines(i),1:nrealiz,i) = sqrt(6*t);
        e = -log(rand(nlines(i)*nrealiz,1));
        f = pi*rand(nlines(i)*nrealiz,1)-pi/2;
        G = abs( sin(b(i)/2*(f-pi/2))./(cos(f).^(2./b(i))).*(cos(f-b(i)/2*(f-pi/2))./e).^((1-b(i)/2)./b(i)*2) )*ones(1,3);
        lines = lines.*sqrt(G/3);
        all_lines(1:nrealiz*nlines(i),:,i) = lines;
        all_phi(1:nlines(i),1:nrealiz,i) = 2*pi*rand(nlines(i),nrealiz);
        all_theta(1:nlines(i),1:nrealiz,i) = sqrt(-log(rand(nlines(i),nrealiz)));
        
        % Cubic covariance model?
        %------------------------
    elseif (model(i,1) == 5)
        x = extreme_coord*lines';
        sigma(i) = 2*sqrt(210./nlines(i));
        all_lines(1:nrealiz*nlines(i),:,i) = lines;
        all_offset(1:nrealiz*nlines(i),i) = min(x)' - rand(nrealiz*nlines(i),1);
        
        % Identify the lines for which the scale factor is too small
        % The simulation along these lines will be replaced by a nugget effect
        interval_number = max(x)'-min(x)';
        nottoosmall = (interval_number < max_intervalnumber);
        valid_lines(1:nrealiz*nlines(i),i) = nottoosmall;
        
        % Gaussian covariance model?
        %---------------------------
    elseif (model(i,1) == 6)
        sigma(i) = sqrt(2./nlines(i));
        t = randn(nlines(i),nrealiz).^2 + randn(nlines(i),nrealiz).^2 + randn(nlines(i),nrealiz).^2;
        all_r(1:nlines(i),1:nrealiz,i) = sqrt(2*t);
        all_phi(1:nlines(i),1:nrealiz,i) = 2*pi*rand(nlines(i),nrealiz);
        all_theta(1:nlines(i),1:nrealiz,i) = sqrt(-log(rand(nlines(i),nrealiz)));
        
        % Cardinal sine covariance model?
        %--------------------------------
    elseif (model(i,1) == 7)
        sigma(i) = sqrt(2./nlines(i));
        t = 2*rand(nlines(i),nrealiz)-1;
        all_r(1:nlines(i),1:nrealiz,i) = sign(t);
        all_phi(1:nlines(i),1:nrealiz,i) = 2*pi*rand(nlines(i),nrealiz);
        all_theta(1:nlines(i),1:nrealiz,i) = sqrt(-log(rand(nlines(i),nrealiz)));
        
        % Bessel-J covariance model?
        %---------------------------
    elseif (model(i,1) == 8)
        sigma(i) = sqrt(2./nlines(i));
        t = betarnd(1.5,b(i)-0.5,nlines(i),nrealiz);
        all_r(1:nlines(i),1:nrealiz,i) = sqrt(t);
        all_phi(1:nlines(i),1:nrealiz,i) = 2*pi*rand(nlines(i),nrealiz);
        all_theta(1:nlines(i),1:nrealiz,i) = sqrt(-log(rand(nlines(i),nrealiz)));
        
        % Bessel-K covariance model (b<0.5)?
        %-----------------------------------
    elseif (model(i,1) == 9)
        v = randn(6,nlines(i)*nrealiz);
        w = 3*rand(1,nlines(i)*nrealiz);
        v(5,:) = v(5,:).*(w>1);
        v(6,:) = v(6,:).*(w>1);
        G = 0.5*(sum(v.^2)'.*sqrt(betarnd(b(i),0.5-b(i),nlines(i)*nrealiz,1))) * ones(1,3);
        lines = lines./G;
        x = extreme_coord*lines';
        sigma(i) = 2*sqrt(3./nlines(i));
        all_lines(1:nrealiz*nlines(i),:,i) = lines;
        all_offset(1:nrealiz*nlines(i),i) = min(x)' - rand(nrealiz*nlines(i),1);
        
        % Identify the lines for which the scale factor is too small
        % The simulation along these lines will be replaced by a nugget effect
        interval_number = max(x)'-min(x)';
        nottoosmall = (interval_number < max_intervalnumber);
        valid_lines(1:nrealiz*nlines(i),i) = nottoosmall;
        
        % Bessel-K covariance model (b>0.5)?
        %-----------------------------------
    elseif (model(i,1) == 9.5)
        sigma(i) = sqrt(2./nlines(i));
        t = randn(nlines(i),nrealiz).^2 + randn(nlines(i),nrealiz).^2 + randn(nlines(i),nrealiz).^2;
        all_r(1:nlines(i),1:nrealiz,i) = sqrt(6*t);
        G = gamrnd(b(i),1,nlines(i)*nrealiz,1)*ones(1,3);
        lines = lines./sqrt(G*12);
        all_lines(1:nrealiz*nlines(i),:,i) = lines;
        all_phi(1:nlines(i),1:nrealiz,i) = 2*pi*rand(nlines(i),nrealiz);
        all_theta(1:nlines(i),1:nrealiz,i) = sqrt(-log(rand(nlines(i),nrealiz)));
        
        % Generalized Cauchy model?
        %--------------------------
    elseif (model(i,1) == 10)
        sigma(i) = sqrt(2./nlines(i));
        t = randn(nlines(i),nrealiz).^2 + randn(nlines(i),nrealiz).^2 + randn(nlines(i),nrealiz).^2;
        all_r(1:nlines(i),1:nrealiz,i) = sqrt(6*t);
        G = gamrnd(b(i),1,nlines(i)*nrealiz,1)*ones(1,3);
        lines = lines.*sqrt(G/3);
        all_lines(1:nrealiz*nlines(i),:,i) = lines;
        all_phi(1:nlines(i),1:nrealiz,i) = 2*pi*rand(nlines(i),nrealiz);
        all_theta(1:nlines(i),1:nrealiz,i) = sqrt(-log(rand(nlines(i),nrealiz)));
        
        % Exponential sine model?
        %------------------------
    elseif (model(i,1) == 11)
        sigma(i) = sqrt(2./nlines(i));
        for k = 1:nrealiz
            for j = 1:nlines(i)
                iflag = -1;
                while (iflag < 0)
                    t2(j,k) = gamrnd(0.5,1);
                    u = rand;
                    c = 1;
                    puis = 1;
                    n = 0;
                    iflag = 0;
                    while (iflag == 0)
                        n = n+1;
                        puis = puis * pi*pi/4 / (2*n*(2*n-1));
                        if (n/2 == round(n/2))
                            c = c+puis*exp(-4*t2(j,k)*n*(n+1));
                            if (u > c), iflag = -1; end
                        else
                            c = c-puis*exp(-4*t2(j,k)*n*(n+1));
                            if (u < c), iflag = 1; end
                        end
                    end
                end
            end
        end
        t1 = randn(nlines(i),nrealiz).^2 + randn(nlines(i),nrealiz).^2 + randn(nlines(i),nrealiz).^2;
        all_r(1:nlines(i),1:nrealiz,i) = sqrt(t1/2./t2);
        all_phi(1:nlines(i),1:nrealiz,i) = 2*pi*rand(nlines(i),nrealiz);
        all_theta(1:nlines(i),1:nrealiz,i) = sqrt(-log(rand(nlines(i),nrealiz)));
               
        % Linear model?
        %--------------
    elseif (model(i,1) == 12)
        x = extreme_coord*lines';
        sigma(i) = 2./sqrt(nlines(i));
        
        % Modify the slope and scale factor if the latter is too small
        delta = max(x) - min(x);
        factor = max(delta)/max_intervalnumber;
        if (factor > 1)
            lines = lines./factor;
            sigma(i) = sigma(i)*sqrt(factor);
            x = extreme_coord*lines';
        end
        
        all_lines(1:nrealiz*nlines(i),:,i) = lines;
        all_offset(1:nrealiz*nlines(i),i) = min(x)' - rand(nrealiz*nlines(i),1);
             
        % Power model with exponent <=1?
        %-------------------------------
    elseif (model(i,1) == 13)
        % Set the scale factor equal to the diameter of the domain to simulate
        x = extreme_coord*lines';
        delta = max(x)-min(x);
        factor = max(delta);
        model(i,2:4) = factor*model(i,2:4);
        cc(i,:) = (factor.^b(i))*cc(i,:);
        sill(:,:,i) = reshape(cc(i,:),nvar,nvar);
        [eigenvectors,eigenvalues] = eig(sill(:,:,i));
        A1(:,:,i) = sqrt(eigenvalues)*eigenvectors';
        R = setrot(model,i);
        model_rotationmatrix(:,:,i) = R;
        lines = V*R';
        
        % The power covariance is simulated as a mixture of triangular models
        sigma(i) = sqrt(2*gamma(b(i)/2+1.5)/gamma(b(i)/2+0.5)./nlines(i));
        u = rand(nlines(i)*nrealiz,1);
        G = max(0.001,min(1,(u/(1-b(i))).^(1/b(i)))) * ones(1,3);
        lines = lines./G;
        x = extreme_coord*lines';
        
        % Modify the slope and scale factor if the latter is too small
        delta = max(x) - min(x);
        factor = max(delta)/max_intervalnumber;
        if (factor > 1)
            model(i,2:4) = factor*model(i,2:4);
            cc(i,:) = (factor.^b(i))*cc(i,:);
            sill(:,:,i) = reshape(cc(i,:),nvar,nvar);
            [eigenvectors,eigenvalues] = eig(sill(:,:,i));
            A1(:,:,i) = sqrt(eigenvalues)*eigenvectors';
            R = setrot(model,i);
            model_rotationmatrix(:,:,i) = R;
            lines = V*R';
            lines = lines./G;
            sigma(i) = sqrt(2*gamma(b(i)/2+1.5)/gamma(b(i)/2+0.5)./nlines(i));
            x = extreme_coord*lines';
        end
        
        all_lines(1:nrealiz*nlines(i),:,i) = lines;
        all_offset(1:nrealiz*nlines(i),i) = min(x)' - rand(nrealiz*nlines(i),1);     
        
        % Pure drift model?
        %------------------
    elseif (model(i,1) == 13.1)
        sigma(i) = sqrt(1./nlines(i));
        all_r(1:nlines(i),1:nrealiz,i) = ...
            exp(0.5*gammaln(2*order(i)+4)-gammaln(order(i)+2))*randn(nlines(i),nrealiz);
        all_lines(1:nrealiz*nlines(i),:,i) = lines;
        all_phi(1,1,i) = order(i)+1;      
        
        % Power model with exponent >1?
        %------------------------------
    elseif (model(i,1) == 13.5)
        sigma(i) = sqrt(gamma(b(i)+2)./nlines(i)./gamma(b(i)+2-2*order(i)));
        t = (randn(nlines(i),nrealiz)./randn(nlines(i),nrealiz)).^2;
        all_r(1:nlines(i),1:nrealiz,i) = 2*pi*t;
        all_lines(1:nrealiz*nlines(i),:,i) = lines;
        all_phi(1:nlines(i),1:nrealiz,i) = 2*pi*rand(nlines(i),nrealiz);
        all_theta(1:nlines(i),1:nrealiz,i) = sqrt( -pi.^(0.5-b(i)).*gamma(b(i)/2+1.5-order(i))./ ...
            2.^(2*order(i)-3)./gamma(order(i)-b(i)/2).*(1+t)./t.^(b(i)+0.5) );      
        
        % Mixed power model?
        %-------------------
    elseif (model(i,1) == 14)
        sigma(i) = sqrt(b(i)./nlines(i));
        t = (randn(nlines(i),nrealiz)./randn(nlines(i),nrealiz)).^2;
        all_r(1:nlines(i),1:nrealiz,i) = 2*pi*t;
        all_lines(1:nrealiz*nlines(i),:,i) = lines;
        all_phi(1:nlines(i),1:nrealiz,i) = 2*pi*rand(nlines(i),nrealiz);
        pow = rand(nlines(i),nrealiz)*b(i);
        all_theta(1:nlines(i),1:nrealiz,i) = sqrt( -8*pi.^(0.5-pow).*gamma(pow/2+1.5)./ ...
            gamma(-pow/2).*(1+t)./exp((pow+0.5).*log(t)) );       
        
        % Spline model?
        %--------------
    elseif (model(i,1) == 15)
        sigma(i) = sqrt(1./nlines(i));
        t = (randn(nlines(i),nrealiz)./randn(nlines(i),nrealiz)).^2;
        all_r(1:nlines(i),1:nrealiz,i) = 2*pi*t;
        all_lines(1:nrealiz*nlines(i),:,i) = lines;
        all_phi(1:nlines(i),1:nrealiz,i) = 2*pi*rand(nlines(i),nrealiz);
        all_theta(1:nlines(i),1:nrealiz,i) = sqrt( pi.^(1-b(i)).*gamma(b(i)+2)./ ...
            2.^(b(i)-1).*(1+t)./t.^(b(i)+0.5) );
        
    end % end model if statement
    
end % end loop over nested models structures

max_nugget = max(abs(nugget));

%--------------------------------------------------------------------------------------------
% NON CONDITIONAL SIMULATION AT THE DATA LOCATIONS
%-------------------------------------------------
simudata = zeros(m0,nrealiz);
weights = zeros(1,nrealiz);

%--------------------------------------------------------------------------------------------
nbdigit = nbdecimal(1)+3;    
    
% How many rows can be simulated simultaneously?
%-----------------------------------------------
ntok = ceil(ntok/ng);
m1 = max(1,min(ny,floor(ntok/nx)));
rows = [[0:m1:ny-0.5] ny];
m2 = length(rows)-1;
m3 = m2*nz;
seed_nugget = ceil(1e7*rand(1,m3));

% SM: Pre-allocate temporary matrix of simulated outputs
SIM_OUTPUTS_TMP=NaN(ny*nx,nrealiz); 

% Loop over the grid rows
%------------------------
for n = 1:m3

    nnz = floor((n-1)/m2);
    nny = n-nnz*m2;
    index = [rows(nny):rows(nny+1)-1]';
    m4 = length(index); % number of rows to simulate
    m5 = m4*nx; % number of blocks to simulate
    m6 = m5*ng; % number of points discretizing the blocks

    % Coordinates of the points to simulate
    %--------------------------------------
    coord0 = [ones(m4,1)*x0 y0+index*dy ones(m4,1)*(z0+nnz*dz)];
    blk_coord = kron(coord0,ones(nx,1)) + kron(ones(m4,1),[0:nx-1]'*[dx 0 0]); % coordinates of the block centres
    coord = kron(blk_coord,ones(ng,1))-kron(ones(m5,1),grid);

    % Non-conditional simulation
    %---------------------------
    simu = tbmain(coord,model,sigma,A1,nvar,nlines,nrealiz,seed_line,all_lines,all_offset, ...
        all_r,all_phi,all_theta,valid_lines);

    % Add nugget effect
    %------------------
    if max_nugget > eps
        randn('state',seed_nugget(n));
        simunug = randn(m6*nrealiz/nvar,nvar)*A0;
        simunug = reshape(simunug',nrealiz,m6);
        simu = simu + simunug';
    end

    % Average the values over the blocks
    %-----------------------------------
    if ng > 1
        simu = reshape(simu,ng,m5*nrealiz);
        simu = mean(simu);
        simu = reshape(simu,m5,nrealiz);
    end
    
    % store simulated values
    SIM_OUTPUTS_TMP(rows(n)*nx+1:rows(n+1)*nx,:)=simu;

end

%--------------------------------------------------------------------------------------------
% Store final outputs in array of dimension ny, nx, nvar
for i=1:nvar
   SIM_OUTPUTS(:,:,i)=reshape(SIM_OUTPUTS_TMP(:,i),nx,ny)'; 
end

return

%% Nested functions (for portability convenience)
%----------------------------------------------------------------------------------------------------------
% Author: Xavier Emery
% Paper: "A turning bands program for conditional co-simulation of cross-correlated Gaussian random fields"
%----------------------------------------------------------------------------------------------------------
function simu = tbmain(coord,model,sigma,A1,nvar,nlines,nrealiz,seed,all_lines,all_offset,all_r,all_phi,all_theta,valid_lines)
%--------------------------------------------------------------
% Non conditional simulation by the turning bands method
% Main routine for program tbcosim (simulation along the lines)
%--------------------------------------------------------------
% Define parameters
%------------------
nst = size(model,1);
m = size(coord,1);

% Initialize output
%------------------
simu = zeros(m,nrealiz);

% Loop over each nested structure
%--------------------------------
for i = 1:nst

  sim = zeros(m,nrealiz);

  if (model(i,1) < 4.5) | (model(i,1) == 9)
      % Spherical, exponential, gamma, stable (b<1) and K-Bessel (b<0.5) models

    % Loop over the realizations
    %---------------------------
    for k = 1:nrealiz
      
      index = [(k-1)*nlines(i)+1:k*nlines(i)];
      lines = all_lines(index,:,i);
      valid = find(valid_lines(index,i) > 0);
      nbvalid = size(valid,1);
      nbinvalid = nlines(i)-nbvalid;
      
      % Project the points to simulate over the lines of the i-th nested structure
      %---------------------------------------------------------------------------
      x = coord*lines'; 
      offset = ones(m,1)*all_offset(index,i)';
      interval_number = floor(x-offset) + 1;
      position = x-offset - interval_number + 0.5;
      state = seed((k-1)*nst+i,1:nlines(i));
      
      % Loop over the lines
      %--------------------
      for j = 1:nbvalid
       
        % Simulate the values by a partition method
        %------------------------------------------
        randn('state',state(valid(j)));
        maxnum = max(interval_number(:,valid(j)));            
        slope = sigma(i)*(2*(randn(maxnum,1)>0)-1);    
        sim(:,k) = sim(:,k) + slope(interval_number(:,valid(j))).*position(:,valid(j));
       
      end

      sim(:,k) = sim(:,k) + sqrt(nbinvalid/nlines(i))*randn(m,1);

    end 
    
  elseif (model(i,1) == 5) % Cubic model
      
    % Loop over the realizations
    %---------------------------
    for k = 1:nrealiz
      
      index = [(k-1)*nlines(i)+1:k*nlines(i)];
      lines = all_lines(index,:,i);
      valid = find(valid_lines(index,i) > 0);
      nbvalid = size(valid,1);
      nbinvalid = nlines(i)-nbvalid;

      % Project the points to simulate over the lines of the i-th nested structure
      %---------------------------------------------------------------------------        
      x = coord*lines'; 
      offset = ones(m,1)*all_offset(index,i)';
      interval_number = floor(x-offset) + 1;
      position = x-offset - interval_number + 0.5;
      state = seed((k-1)*nst+i,1:nlines(i));
      
      % Loop over the lines
      %--------------------
      for j = 1:nbvalid

        % Simulate the values by a partition method
        %------------------------------------------
        randn('state',state(valid(j)));
        maxnum = max(interval_number(:,valid(j)));
        slope = sigma(i)*(2*(randn(maxnum,1)>0)-1);                      
        sim(:,k) = sim(:,k) + slope(interval_number(:,valid(j))).*(0.25*position(:,valid(j)) - position(:,valid(j)).^3);

      end

      sim(:,k) = sim(:,k) + sqrt(nbinvalid/nlines(i))*randn(m,1);

    end

   
  elseif (model(i,1) == 12) % Linear model

    % Loop over the realizations
    %---------------------------
    for k = 1:nrealiz
      
      index = [(k-1)*nlines(i)+1:k*nlines(i)];
      lines = all_lines(index,:,i);
      
      % Project the points to simulate over the lines of the i-th nested structure
      %---------------------------------------------------------------------------       
      x = coord*lines'; 
      offset = ones(m,1)*all_offset(index,i)';
      interval_number = floor(x-offset) + 1;
      state = seed((k-1)*nst+i,1:nlines(i));
        
      % Loop over the lines
      %--------------------
      for j = 1:nlines(i)

        % Simulate the values by a partition method
        %------------------------------------------
        randn('state',state(j));
        maxnum = max(interval_number(:,j));
        stairs = sigma(i)*cumsum(2*(randn(maxnum,1)>0)-1);                      
        sim(:,k) = sim(:,k) + stairs(interval_number(:,j));

      end
          
    end


  elseif (model(i,1) == 13) % Power model (b<=1)
    % Loop over the realizations
    %---------------------------
    for k = 1:nrealiz
      
      index = [(k-1)*nlines(i)+1:k*nlines(i)];
      lines = all_lines(index,:,i);

      % Project the points to simulate over the lines of the i-th nested structure
      %---------------------------------------------------------------------------
      x = coord*lines'; 
      offset = ones(m,1)*all_offset(index,i)';
      interval_number = floor(x-offset) + 1;
      state = seed((k-1)*nst+i,1:nlines(i));

      % Loop over the lines
      %--------------------
      for j = 1:nlines(i)

        % Simulate the values by a partition method
        %------------------------------------------
        randn('state',state(j));
        maxnum = max(interval_number(:,j));            
        value = sigma(i)*randn(maxnum,1);
        sim(:,k) = sim(:,k) + value(interval_number(:,j));

      end

    end


  elseif (model(i,1) == 13.1) % Pure drift model

    % Loop over the realizations
    %---------------------------
    for k = 1:nrealiz
      
      % Project the points to simulate over the lines of the i-th nested structure
      %---------------------------------------------------------------------------
      index = [(k-1)*nlines(i)+1:k*nlines(i)];
      lines = all_lines(index,:,i);
      x = coord*lines';    

      % Simulate the values by a monomial
      %----------------------------------
      x = x.^all_phi(1,1,i);
      r = all_r(:,k,i) * ones(1,m);
      sim(:,k) = sim(:,k) + sigma(i)*sum(r.*x')';
  
    end  

  else % Stable (b>1), Gaussian, cardinal sine, Bessel J, Bessel K (b>0.5), generalized Cauchy, 
       % exponential sine, power (b>1), mixed power and spline models

    % Loop over the realizations
    %---------------------------
    for k = 1:nrealiz
      
      % Project the points to simulate over the lines of the i-th nested structure
      %---------------------------------------------------------------------------
      index = [(k-1)*nlines(i)+1:k*nlines(i)];
      lines = all_lines(index,:,i);
      x = coord*lines';    

      % Simulate the values by a continuous spectral method
      %----------------------------------------------------
      r = all_r(:,k,i) * ones(1,m);
      phi = all_phi(:,k,i) * ones(1,m);
      theta = all_theta(:,k,i) * ones(1,m);
      sim(:,k) = sim(:,k) + sigma(i)* sum( theta.*cos(r.*x'+phi)  )';
  
    end
      
  end
  
  sim = reshape(sim',nvar,nrealiz/nvar*m);
  sim = reshape(A1(:,:,i)'*sim,nrealiz,m);
  simu = simu + sim';
  
end
%%
%----------------------------------------------------------------------------------------------------------
% Author: Xavier Emery
% Paper: "A turning bands program for conditional co-simulation of cross-correlated Gaussian random fields"
%----------------------------------------------------------------------------------------------------------
function lines = vdc(nlines,nrealiz,seed)
%------------------------------------------------------------
% Generation of equidistributed lines over the unit 3D sphere
% according to a van der Corput sequence
%------------------------------------------------------------
%
% USE: lines = vdc(nlines,nrealiz,seed)
%
% INPUT:
%   nlines: number of lines to generate for each realization
%   nrealiz: number of realizations
%   seed: seed for generation of random numbers
%
% OUTPUT:
%   lines: n*3 matrix representing the lines

% This program uses the following subroutine:
%     setrot.m: set up matrix for rotation and reduction of coordinates

lines = zeros(nlines*nrealiz,3);
rand('state',seed);
seed2 = ceil(1e7*rand);
randn('state',seed2);

i = [1:nlines]';

% binary decomposition of i
j = i;
u = 0;
p = 0;
while (max(j)>0)
  p = p+1;
  t = fix(j/2);
  u = u+2*(j/2-t)./(2.^p);
  j = t; 
end

% ternary decomposition of i
j = i; 
v = 0; 
p = 0;
while (max(j)>0)
  p = p+1;
  t = fix(j/3);
  v = v+3*(j/3-t)./(3.^p);
  j = t; 
end

% directing vector of the i-th line
x  = [cos(2*pi*u).*sqrt(1-v.*v) sin(2*pi*u).*sqrt(1-v.*v) v]; 

% random rotation
for k = 1:nrealiz
  angles = 360*rand(1,3);
  R = setrot([1 1 1 1 angles],1);
  lines((k-1)*nlines+1:k*nlines,:) = x*R;
end
%%
%----------------------------------------------------------------------------------------------------------
% Author: Xavier Emery
% Paper: "A turning bands program for conditional co-simulation of cross-correlated Gaussian random fields"
%----------------------------------------------------------------------------------------------------------
function rotred_matrix = setrot(model,it);
%------------------------------------------------------------
% Set up the matrix to transform the Cartesian coordinates to 
% coordinates that account for angles and anisotropy 
%------------------------------------------------------------
%
% INPUT:
%   model: nested variogram model
%   it: index of the nested structure
%   The rotation is performed according to the GSLIB conventions

deg2rad = pi/180;
ranges = model(it,2:4);
angles = model(it,5:7);

% matrix of coordinate reduction
redmat = diag(1./(eps+ranges));

a = (90-angles(1))*deg2rad;
b = -angles(2)*deg2rad;
c = angles(3)*deg2rad;

cosa = cos(a);
sina = sin(a);
cosb = cos(b);
sinb = sin(b);
cosc = cos(c);
sinc = sin(c);

rotmat = zeros(3,3);
rotmat(1,1) = cosb * cosa;
rotmat(1,2) = cosb * sina;
rotmat(1,3) = -sinb;
rotmat(2,1) = -cosc*sina + sinc*sinb*cosa;
rotmat(2,2) = cosc*cosa + sinc*sinb*sina;
rotmat(2,3) =  sinc * cosb;
rotmat(3,1) = sinc*sina + cosc*sinb*cosa;
rotmat(3,2) = -sinc*cosa + cosc*sinb*sina;
rotmat(3,3) = cosc * cosb;

rotred_matrix = inv(rotmat)*redmat;
%%

