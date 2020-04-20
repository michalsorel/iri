function [IRI] = iri_v2(Y,METHOD,VERBOSE,segmPar,avgPar)
%
% calculate IRI 
%
% IRI = iri_v2(Y, METHOD, avg)
%
% Input:
% Y ... road profile, Nx2 vector [x,y], 
% Both x and y must be in meters.
% METHOD ... 0 - Sayers' implementation, 1 - numerical solver ode45,
%               2(default) - our semi-analytical solution
% VERBOSE ... 0(default) - no plot; 1 - plot IRI graph 
% segmPar ... parameters of segments; vector containing 3 parameters: 
%               [segment_length(m), overlapping_segments(0/1), starting_position(m)] 
%             default: [20,0,0]
% avgPar ... averaging parameter for smoothing profile; convolution with
%           Gaussian function of variance proportional to 1/avgPar
%            default: None 
% Output:
% IRI ... IRI vector calculated in segments; matrix with rows
%           [star_pos, end_pos, IRI, std_of_IRI]
% IRI is in mm/m [equivalent to m/km]
%
% Note: 
% similar as iri.m with:
% - calculation of S and p for arbitrary dx in Sayers' method
% - setting parameters for segments (segmPar)
% - profile averaging (avgPar)



if nargin<2 || isempty(METHOD)
    METHOD = 2;
end
if nargin<3 || isempty(VERBOSE)
    VERBOSE = 0;
end
if nargin<4 || isempty(segmPar)
    segmPar = [20 0 0];
end
if nargin<5 || isempty(avgPar)
    avgPar = 10;
end

% Standard car parameters
K1 = 653; 
K2 = 63.3;
C = 6;
u = 0.15;

A = [ 0 1 0 0; -K2 -C K2 C; 0 0 0 1; K2/u C/u -(K1+K2)/u -C/u ];
b = [ 0 0 0 K1/u ].';

% standard car speed 80km/h = 22.2m/s
v = 80/3.6;

% remove duplicated points (e.g. Kunraticka data has duplicated points)
DX = diff(Y(:,1));
mask = [true; DX~=0];
YC = Y(mask,:);
DX = DX(DX~=0);
if sum(~mask) ~= 0
    disp(['Removing ',num2str(sum(~mask)),' duplicated points.']);
end
% check if data are equidistant
dx = median(DX);
if sum(abs(DX - dx) < dx/1000) ~= length(DX)
    disp('Data are not equidistant within 1 permille tolerance.');
end 

% Calculate numerically S and p for Sayers' method
if METHOD == 0
    % Default for 0.25cm
    %S = [0.9966071, .01091514, -0.002083274, 0.0003190145;...
    %    -0.5563044, 0.9438768 -0.8324718, 0.05064701; ...
    %    0.02153176, 0.002126763, 0.7508714, 0.008221888; ...
    %    3.335013 , 0.3376467 ,-39.12762, 0.4347564];
    %p = [0.005476107; 1.388776; 0.2275968; 35.79262];
    N = 100;
    t = dx/v;
    S = eye(4);
    tA = A;
    for i = 1:N
        S = S + tA*(t^i/factorial(i));
        tA = tA*A;
    end
    p = A\((S-eye(4))*b);
end


% average profile with a rectangular mask of length 0.25m as done in
% Sayers' paper
%for i=1:size(YC,1)
%    ind = find(YC(:,1)>=YC(i,1) & YC(:,1)<=YC(i,1)+0.25);
%    YC(i,2) = mean(YC(ind,2));
%end

% average profile with a Gaussian mask 
if avgPar<1
    disp(['Averaging the profile.']);
    h = gaussian(avgPar);
    y = padarray(YC(:,2),(length(h)-1)/2,'symmetric','both');
    YC(:,2) = conv(y,h,'valid');
end



% calculate IRI in 20m non-overlapping windows
SX = segmPar(1);
floating = segmPar(2);
start_pos = segmPar(3);

si = find(YC(:,1)<=start_pos,1,'last');
if isempty(si)
    si = 1;
end
if floating
    end_X = SX+YC(si,1):dx:YC(end,1);
else
    end_X = SX+YC(si,1):SX:YC(end,1);
end
maxs = length(end_X);
IRI = zeros(maxs,4);

disp('========================================');
disp('IRI computation');
disp(['Method: ' num2str(METHOD)]);
disp(['Starting position: ',num2str(YC(si,1)),'m']);
disp(['Segment length: ',num2str(SX),'m']); 
disp(['Number of segments: ',num2str(maxs)]); 


% initialization of Z according to the IRI book
t0 = 0.5;
Y = [YC(si:end,1) - YC(si,1), YC(si:end,2)];
Y = Y.*[1/v,1000];
y0 = Y(1,2);
y0p = (interp1(Y(:,1),Y(:,2),t0,'linear','extrap') - Y(1,2))/(t0);
if METHOD == 0
    z0 = [y0p/v; 0; y0p/v; 0];
else
    z0 = [y0; y0p; y0; y0p];
end
    
% main loop;take profile in SX chunks
for s = 1:maxs
ei = find(YC(:,1)<=end_X(s),1,'last');
IRI(s,1:2) = [YC(si,1),YC(ei,1)];

Y = [YC(si:ei,1)-YC(si,1), YC(si:ei,2)]; % set the first point to start at distance 0
   
% convert distance along the profile to time for car speed 80km/h
ft = Y(:,1)/v; % taking distance as projected to a flat surface
%ft = cumsum([0; sqrt(sum(diff(Y).^2,2))])/v; % taking distance along the surface
% convert y to [mm]
Y(:,2) = Y(:,2)*1000;

%if s == 1
    %t0 = 0.5;
    %y0 = Y(1,2);
    %y0p = (interp1(ft,Y(:,2),t0,'linear','extrap') - Y(1,2))/(t0);
    %z0 = [y0; y0p; y0; y0p];
%else
%    z0 = Z(end,:).';
%end

switch METHOD
case 0 % M.Sayers 
    %y0 = (interp1(ft,Y(:,2),t0,'linear','extrap') - Y(1,2))/(v*t0); % Initialization according to IRI book
    %z0 = [y0; 0; y0; 0];
    %dx = 0.25;
    % resample profile with step dx
    X = 0:dx:Y(end,1);
    Y = [X.', interp1(Y(:,1),Y(:,2),X.')];
    DY = diff(Y(:,2))/dx;
    Z = zeros(length(X),4);
    Z(1,:) = z0;
    for i = 1:length(X)-1
        Z(i+1,:) = S*Z(i,:).' + p*DY(i);
    end
    T = X/v;
    D = abs(Z(2:end,1)-Z(2:end,3));
case 1 % Numerical solver ode45
    %options = odeset('Stats','on', 'AbsTol', 1e-9, 'RelTol', 1e-5);
    options = odeset('Stats','on');
    %ft2 = linspace(ft(1),ft(end),length(ft)*10);
    [T,Z] = ode45(@myode, ft, z0, options);
    X = T*v;
    D = abs(Z(2:end,2)-Z(2:end,4)).*diff(T)/X(end)*(size(Z,1)-1);
case 2 % Semi-analytical solution, Sroubek
    fintM = @iri_iMy;
    fM = @iri_M;
    Z = zeros(size(Y,1),4);
    Z(1,:) = z0;
    iZ = inv(fM(0))*z0;
    cS = 0;
    for i = 2:size(Z,1)
        cS = cS + fintM(ft(i-1),ft(i),Y(i-1,2),Y(i,2));
        if sum(isnan(cS))>0
            warning('NaN value!');
            keyboard;
        end
        Z(i,:) = fM(ft(i))*(iZ + cS);
    end
    X = Y(:,1);
    D = abs(Z(2:end,2)-Z(2:end,4)).*diff(ft)/X(end)*(size(Z,1)-1);
    %D = abs(Z(2:end,1)-Z(2:end,3)).^2;
end
IRI(s,3) = mean(D);
IRI(s,4) = std(D);
if floating
    si = si+1;
    z0 = Z(2,:).';
else
    si = ei;
    z0 = Z(end,:).';
end

end % end of the main loop

if VERBOSE
    figure;enlarge_figure(1,2); 
    %yyaxis left; 
    %plot(DX*(1:maxs),IRI);
    %yyaxis right;
    %plot(YC(:,1),YC(:,2));
    h = plotyy(IRI(:,1),IRI(:,3),YC(:,1),YC(:,2));
    legend({'IRI','Road profile'});
    ylabel(h(1),'IRI');
    ylabel(h(2),'Elevation [m]');
    xlabel({'Stationing [m]',' '});
    title(sprintf('IRI: %0.2f \\pm %0.2f',mean(IRI(:,3)),std(IRI(:,3))));
end

%figure; plot(X(2:end), D,X, Y(:,2));
%figure; plot(X,Z(:,1)-Y(:,2),X,Z(:,3)-Y(:,2));
%figure; plot(T,Z(:,1),T,Z(:,3));

    function dzdt = myode(t,z)
        y = interp1(ft,Y(:,2),t, 'linear'); % Interpolate the input profile (ft,Y(:,2)) at time t
        dzdt = A*z + b*y;
    end

end

function iMy = iri_iMy(t1,t2,y1,y2)
%IRI_IMY
%    IMY = IRI_IMY(T1,T2,Y1,Y2)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    24-May-2018 13:37:03

t4 = t1.*2.040992621332219e1;
t5 = exp(t4);
t6 = t1.*6.442262610541198e1;
t7 = t2.*2.040992621332219e1;
t8 = exp(t7);
t9 = t2.*6.442262610541198e1;
t10 = cos(t6);
t11 = t1-t2;
t12 = 1.0./t11;
t13 = cos(t9);
t14 = sin(t6);
t15 = sin(t9);
t16 = t1.*2.590073786677789;
t17 = exp(t16);
t18 = t1.*7.323397385729718;
t19 = t2.*2.590073786677789;
t20 = exp(t19);
t21 = t2.*7.323397385729718;
t22 = sin(t18);
t23 = sin(t21);
t24 = cos(t18);
t25 = cos(t21);
iMy = [t5.*t10.*y1.*(-8.023278807758917e-1)-t5.*t14.*y1.*6.760740228704583e1+t8.*t13.*y1.*8.023278807758917e-1+t8.*t15.*y1.*6.760740228704583e1-t5.*t10.*t12.*y1.*9.501254908070717e-1+t5.*t10.*t12.*y2.*9.501254908070717e-1+t5.*t12.*t14.*y1.*3.134662503279774e-1-t5.*t12.*t14.*y2.*3.134662503279774e-1+t8.*t12.*t13.*y1.*9.501254908070717e-1-t8.*t12.*t13.*y2.*9.501254908070717e-1-t8.*t12.*t15.*y1.*3.134662503279774e-1+t8.*t12.*t15.*y2.*3.134662503279774e-1-t1.*t8.*t12.*t13.*y1.*8.023278807758917e-1+t1.*t8.*t12.*t13.*y2.*8.023278807758917e-1+t2.*t8.*t12.*t13.*y1.*8.023278807758917e-1-t1.*t8.*t12.*t15.*y1.*6.760740228704583e1-t2.*t8.*t12.*t13.*y2.*8.023278807758917e-1+t1.*t8.*t12.*t15.*y2.*6.760740228704583e1+t2.*t8.*t12.*t15.*y1.*6.760740228704583e1-t2.*t8.*t12.*t15.*y2.*6.760740228704583e1;...
    t5.*t10.*y1.*6.760740228704583e1-t5.*t14.*y1.*8.023278807758917e-1-t8.*t13.*y1.*6.760740228704583e1+t8.*t15.*y1.*8.023278807758917e-1-t5.*t10.*t12.*y1.*3.134662503279774e-1+t5.*t10.*t12.*y2.*3.134662503279774e-1-t5.*t12.*t14.*y1.*9.501254908070717e-1+t5.*t12.*t14.*y2.*9.501254908070717e-1+t8.*t12.*t13.*y1.*3.134662503279774e-1-t8.*t12.*t13.*y2.*3.134662503279774e-1+t8.*t12.*t15.*y1.*9.501254908070717e-1-t8.*t12.*t15.*y2.*9.501254908070717e-1+t1.*t8.*t12.*t13.*y1.*6.760740228704583e1-t1.*t8.*t12.*t13.*y2.*6.760740228704583e1-t2.*t8.*t12.*t13.*y1.*6.760740228704583e1-t1.*t8.*t12.*t15.*y1.*8.023278807758917e-1+t2.*t8.*t12.*t13.*y2.*6.760740228704583e1+t1.*t8.*t12.*t15.*y2.*8.023278807758917e-1+t2.*t8.*t12.*t15.*y1.*8.023278807758917e-1-t2.*t8.*t12.*t15.*y2.*8.023278807758917e-1;...
    t17.*t22.*y1.*(-6.843274120037415)+t17.*t24.*y1.*5.58745119585928+t20.*t23.*y1.*6.843274120037415-t20.*t25.*y1.*5.58745119585928-t12.*t17.*t22.*y1.*3.843934011732957e-1+t12.*t17.*t22.*y2.*3.843934011732957e-1-t12.*t17.*t24.*y1.*1.070388643317791+t12.*t17.*t24.*y2.*1.070388643317791+t12.*t20.*t23.*y1.*3.843934011732957e-1-t12.*t20.*t23.*y2.*3.843934011732957e-1+t12.*t20.*t25.*y1.*1.070388643317791-t12.*t20.*t25.*y2.*1.070388643317791-t1.*t12.*t20.*t23.*y1.*6.843274120037415+t1.*t12.*t20.*t23.*y2.*6.843274120037415+t2.*t12.*t20.*t23.*y1.*6.843274120037415+t1.*t12.*t20.*t25.*y1.*5.58745119585928-t2.*t12.*t20.*t23.*y2.*6.843274120037415-t1.*t12.*t20.*t25.*y2.*5.58745119585928-t2.*t12.*t20.*t25.*y1.*5.58745119585928+t2.*t12.*t20.*t25.*y2.*5.58745119585928;...
    t17.*t22.*y1.*5.58745119585928+t17.*t24.*y1.*6.843274120037415-t20.*t23.*y1.*5.58745119585928-t20.*t25.*y1.*6.843274120037415-t12.*t17.*t22.*y1.*1.070388643317791+t12.*t17.*t22.*y2.*1.070388643317791+t12.*t17.*t24.*y1.*3.843934011732957e-1-t12.*t17.*t24.*y2.*3.843934011732957e-1+t12.*t20.*t23.*y1.*1.070388643317791-t12.*t20.*t23.*y2.*1.070388643317791-t12.*t20.*t25.*y1.*3.843934011732957e-1+t12.*t20.*t25.*y2.*3.843934011732957e-1+t1.*t12.*t20.*t23.*y1.*5.58745119585928-t1.*t12.*t20.*t23.*y2.*5.58745119585928-t2.*t12.*t20.*t23.*y1.*5.58745119585928+t1.*t12.*t20.*t25.*y1.*6.843274120037415+t2.*t12.*t20.*t23.*y2.*5.58745119585928-t1.*t12.*t20.*t25.*y2.*6.843274120037415-t2.*t12.*t20.*t25.*y1.*6.843274120037415+t2.*t12.*t20.*t25.*y2.*6.843274120037415];
end

function M = iri_M(t)
%IRI_M
%    M = IRI_M(T)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    24-May-2018 13:37:03

t2 = t.*6.442262610541198e1;
t10 = t.*2.040992621332219e1;
t3 = exp(-t10);
t4 = cos(t2);
t5 = sin(t2);
t6 = t.*7.323397385729718;
t11 = t.*2.590073786677789;
t7 = exp(-t11);
t8 = cos(t6);
t9 = sin(t6);
M = reshape([-t3.*(t4.*1.009488570317789e-3+t5.*8.315147588488827e-4),-t3.*(t4.*3.296477717707572e-2-t5.*8.200505959666402e-2),-t3.*(t4.*4.45117663500498e-3-t5.*1.404985422724839e-2),t3.*t4.*9.959766924004811e-1,t3.*(t4.*8.315147588488827e-4-t5.*1.009488570317789e-3),-t3.*(t4.*8.200505959666402e-2+t5.*3.296477717707572e-2),-t3.*(t4.*1.404985422724839e-2+t5.*4.45117663500498e-3),t3.*t5.*9.959766924004811e-1,-t7.*(t8.*4.238836503712554e-2-t9.*1.19852508949723e-1),t7.*t8.*9.875165438583369e-1,t7.*(t8.*4.057897959535129e-3+t9.*1.115993043207963e-2),t7.*(t8.*7.121835021721261e-2-t9.*5.862264258169474e-2),-t7.*(t8.*1.19852508949723e-1+t9.*4.238836503712554e-2),t7.*t9.*9.875165438583369e-1,-t7.*(t8.*1.115993043207963e-2-t9.*4.057897959535129e-3),t7.*(t8.*5.862264258169474e-2+t9.*7.121835021721261e-2)],[4,4]);
end

function h = gaussian(vx)
% h = gaussian(vx)
% return a vector that contains 1D Gaussian with such variation,
% that vx ratio of Fourier coef. along x are above exp(-0.5)
% the vector size is such that on boundarie the values are 1/20 of the 
% middle value
v = 1./(2*pi*vx);
j = ceil(sqrt(-2*log(0.05))./(2*pi*vx));
x = -j:j;
h = 1/(sqrt(2*pi)*v)*exp(-0.5*((x.^2)/v^2));
h = h/sum(h(:));
end