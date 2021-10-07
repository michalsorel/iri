function IRI = iri(Y, segment_length, start_pos, overlap, box_filter, method)
%IRI calculate IRI 
%
%   function IRI = iri(Y, segment_length, start_pos, overlap, box_filter, method)
% 
% Input:
% Y ... road profile, Nx2 vector [x,y], both x and y must be in meters.
% segment_length ... length of the IRI segment (typically 20 or 100 meters)
% start_pos ... starting position of IRI segments, if empty, set to the first sample
% overlap ... if empty, false or 0, IRI is computed in segments without overlap
%               otherwise, the overlap is overlap of IRI segments (can be
%               set for example to sampling step for regular profile sampling);
%               overlap cannot be logical true.
% box_filter ... true/false (impl. value true), if true,
%               profile is first averaged with a box filter of length 0.25m 
%               as recommended in Sayers' paper for sampling intervals
%               shorter than 0.25m to better represent the way in which 
%               the tire of a vehicle envelops the ground.
% method ...    0 - Sayers' implementation, for irregular sampling first resampled                
%               1 - numerical solver ode45,
%               2(default) - our semi-analytical solution (Sroubek & Sorel)
% 
% Output:
% IRI ... IRI vector calculated in segments; matrix with rows
%           [star_pos, end_pos, IRI, std_of_IRI]
%          std_of_IRI is standard deviation of IRI within each segment
% IRI is in mm/m [equivalent to m/km]
% 

if nargin < 1
    error('Function needs at least two parameters.');
end
if ~exist('segment_length','var') || isempty(segment_length)
    segment_length = 20; 
    disp('Parameter segment_length does not exist, set to 20m.');
end
if ~exist('start_pos','var') || isempty(start_pos)    
    disp('Parameter start_pos does not exist, set to the first sample.');
    start_pos = Y(1,1);  
end
if ~exist('overlap','var') || isempty(overlap)
    overlap = false;
end
if islogical(overlap) && overlap == true
    error('overlap parametr cannot be logical true, must be either false or a number containing overlap size.');
end
if ~exist('box_filter','var') || isempty(box_filter)
    box_filter = true;
end
if ~exist('method','var') || isempty(method)
    method = 2;
end

% basic parameter checks
if any(diff(Y(:,1))<=0)
    error('First column of the first parameter (stationing) should be sorted and non-duplicated.');
end
if start_pos < Y(1,1)
    error('Starting position should be within the sampled profile.');
end
if start_pos > Y(end,1)
    error('Starting position should be within the sampled profile.');
end
if overlap > 0
    if overlap > segment_length
        error('Segment overlap must be shorter than segment length.');
    end
end

disp('========================================');
disp('IRI computation');
disp('----------------------------------------');
% Standard car parameters
K1 = 653; 
K2 = 63.3;
C = 6;
u = 0.15;

A = [ 0 1 0 0; -K2 -C K2 C; 0 0 0 1; K2/u C/u -(K1+K2)/u -C/u ];
b = [ 0 0 0 K1/u ].';

% standard car speed 80km/h = 22.2m/s
v = 80/3.6;

% remove duplicated points - should never happen but to be sure
DX = diff(Y(:,1));
mask = [true; DX~=0];
YC = Y(mask,:);
DX = DX(DX~=0);
if sum(~mask) ~= 0
    disp(['Removing ',num2str(sum(~mask)),' duplicated points.']);
end

% check if data are equidistant
dx = median(DX);
if max(abs(DX - dx)) > dx/1e5
    disp('Data are not equidistant (tolerance 1e-5).');
    equidistant = false;
else
    disp(['Regular sampling: ' num2str(dx) ' m']);
    equidistant = true;
end 

% Calculate numerically S and p for Sayers' method
if method == 0
    % Default for dx == 0.25cm
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
    if ~equidistant
        error('Sayers method for irregular sampling does not work (would take too long).');
    end    
    if isempty(find(start_pos==YC(:,1),1))
        error('Starting positions in Sayers method must coincide with input sampling.');
    end
end

% average profile with a box filter of length 0.25m as recommended in Sayers' paper,
% for regular samplings could be implemented much faster by convolution
if box_filter
    if any(diff(YC(:,1))<=0.25/2) % check if necessary
        for i = 1:size(YC,1)
           YC(i,2) = mean(YC(YC(:,1) >= YC(i,1)-0.25/2 & YC(:,1) <= YC(i,1)+0.25/2,2));
        end
    else
        disp('Samples are more than 0.125m from each other, box filter not needed, skipped...');
    end
end

% interpolate all segment lower/upper ends and add them to profile
end_X = (start_pos+segment_length:segment_length:YC(end,1))';
YY = interp1(YC(:,1),YC(:,2),[start_pos;end_X],'linear');
YC = union(YC,[[start_pos;end_X] YY],'rows'); % now including non-overlapping segment boundaries
[~,si] = intersect(YC(:,1),start_pos); % si = index of the beginning of the first segment
if isempty(si) 
    error('Assert: unique starting position should always exist at this moment.');
end
if ~overlap    
    [~,end_si] = intersect(YC(:,1),end_X);
    start_si = [si;end_si(1:end-1)];    
    end_si_o = end_si;
    start_si_o = start_si;
    nsegments = length(end_X);
    nsegments_o = nsegments;
else % overlapping intervals with step "overlap"
    overlap = double(overlap); % convert to double
    if method == 0 && mod(overlap,dx) > 0 %overlap must be a multiple of dx 
        error('Sayers method with overlapped segments requires the overlap to be a multiple of the sampling step.');
    end    
    start_X_o = (start_pos:overlap:YC(end,1)-segment_length)'; % lower ends of the segments
    end_X_o = (start_pos+segment_length:overlap:YC(end,1))'; % upper ends
    XX = union(start_X_o,end_X_o); 
    YY = interp1(YC(:,1),YC(:,2),XX,'linear');
    YC = union(YC,[XX YY],'rows');  % now including both overlapping and non-overlapping segment boundaries   
    [~,end_si] = intersect(YC(:,1),end_X);
    start_si = [si;end_si(1:end-1)];       
    [~,start_si_o] = intersect(YC(:,1),start_X_o);
    [~,end_si_o] = intersect(YC(:,1),end_X_o);    
    nsegments = length(end_X);
    nsegments_o = length(end_X_o);
    nsegments_aux = ceil((YC(end,1)-start_pos)/segment_length);
    if mod(YC(end,1)-start_pos,segment_length) > 0
        YC = [YC; start_pos+nsegments_aux*segment_length YC(end,2)]; % add one dummy sample point
        %dummy = true;
        nsegments = nsegments_aux; % increase the number of segments by one
        start_si = [start_si;end_si(end)]; % from the end to the dummy point
        end_si = [end_si;size(YC,1)]; 
    %else 
        %dummy = false;
    end    
end

% MELO BY JIT TEN ZACATEK ODSTRANIT?
%YC = YC(si:end,:); % overit vsechny dalsi vyskyty, ze to fakt nevadi

%IRI = zeros(maxs,4);
%IRI1 = zeros(maxs,4);

disp(['Segment length: ',num2str(segment_length),'m']); 
disp(['Starting position: ',num2str(YC(si,1)),'m']);
disp(['Method: ' num2str(method)]);
disp(['Number of segments: ',num2str(nsegments_o)]); 

% initialization of Z according to the IRI book
t0 = 0.5;
Y = [YC(si:end,1) - YC(si,1), YC(si:end,2)];
Y = Y.*[1/v,1000]; 
y0 = Y(1,2); % [mm]
Y = unique(Y,'rows');  % just in case /v create duplicates (not solved in YC for now)
y0p = (interp1(Y(:,1),Y(:,2),t0,'linear','extrap') - Y(1,2))/(t0); % [mm/s]
if method == 0
    z0 = [y0p/v; 0; y0p/v; 0];
else
    z0 = [y0; y0p; y0; y0p];
end

if method == 2
    fM0 = fM(0);
end

DD = zeros(end_si(end)-start_si(1)+1,1);
% main loop - take profile in segments of length segment_length
for s = 1:nsegments    
    si = start_si(s);
    ei = end_si(s);
    Y = [YC(si:ei,1)-YC(si,1), YC(si:ei,2)]; % set the first point to start at distance 0    
    ft = Y(:,1)/v; % convert top view distance along the profile to time for car speed 80km/h 
    %ft = cumsum([0; sqrt(sum(diff(Y).^2,2))])/v; % distance along the surface    
    Y(:,2) = Y(:,2)*1000; % convert y to [mm]

    switch method
    case 0 % M.Sayers 
        % resample profile with step dx 
        %X = 0:dx:Y(end,1);
        %Y = [X.', interp1(Y(:,1),Y(:,2),X.')];
        DY = diff(Y(:,2))/dx;
        Z = zeros(size(Y,1),4);
        Z(1,:) = z0;
        for i = 1:size(Y,1)-1
            Z(i+1,:) = S*Z(i,:).' + p*DY(i);
        end
        D = abs(Z(2:end,1)-Z(2:end,3))/(segment_length/dx);
        DD(start_si(s)-start_si(1)+1:end_si(s)-start_si(1)) = D;
    case 1 % Numerical solver ode45
        options = odeset('Stats','off', 'AbsTol', 1e-9, 'RelTol', 1e-5);
        [T,Z] = ode45(@myode, ft, z0, options);
        X = T*v;
        D = abs(Z(2:end,2)-Z(2:end,4)).*diff(T)/X(end);
        DD(start_si(s)-start_si(1)+1:end_si(s)-start_si(1)) = D;
    case 2 % Semi-analytical solution, Sroubek & Sorel
        Z = zeros(size(Y,1),4);
        Z(1,:) = z0;
        iZ = fM0\z0; %iZ = inv(fM0)*z0;
        cS = 0;
        for i = 2:size(Z,1)
            if ft(i)-ft(i-1) > 1e-10                
                cS = cS + fintM(ft(i-1),ft(i),Y(i-1,2),Y(i,2));
                if any(isnan(cS))
                    warning('NaN value!');
                end
                Z(i,:) = fM(ft(i))*(iZ + cS);
            else 
                Z(i,:) = Z(i-1,:); % is it more precise or not?
                %Z(i,:) = fM(ft(i))*(iZ + cS); 
            end
        end
        D = abs(Z(2:end,2)-Z(2:end,4)).*diff(ft)/Y(end,1);
        DD(start_si(s)-start_si(1)+1:end_si(s)-start_si(1)) = D;
    end    
    z0 = Z(end,:).'; 
end % end of the main loop
IRI = zeros(nsegments_o,4);
for s = 1:nsegments_o
    IRI(s,1:2) = [YC(start_si_o(s),1),YC(end_si_o(s),1)];
    IRI(s,3) = sum(DD(start_si_o(s)-start_si_o(1)+1:end_si_o(s)-start_si_o(1)));
    IRI(s,4) = std(DD(start_si_o(s)-start_si_o(1)+1:end_si_o(s)-start_si_o(1)))...
                    *(end_si_o(s)-start_si_o(s));
end

function dzdt = myode(t,z)
    y = interp1(ft,Y(:,2),t, 'linear'); % Interpolate the input profile (ft,Y(:,2)) at time t
    dzdt = A*z + b*y;
end

end

function iMy = fintM(t1,t2,y1,y2)
%fintM
%
%    function iMy = fintM(t1,t2,y1,y2)
%
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

function M = fM(t)
%fM
%
%    function M = fM(t)
%
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
