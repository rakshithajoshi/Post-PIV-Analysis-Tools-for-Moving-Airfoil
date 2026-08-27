%check the contour of v to see if the positive and negative are accurate
%rms values are calculated for normalised and instantaneuos profile
%extended to all x/c <= -1 and across all frames (not just frame segments)

%All fluxes are positive outward
%Flux calculations based on moving CV
%all dimensional variables are in SI units
%velocity profiles being written are in cm/s

%all the output as .mat file. No .dat files are written here.

clc;
input = inputdlg({'Pitching Point',...
    'Amplitude Sequence',...
    'Frequency(Hz)',...
    'FPS',...
    'Trial No',...
    'Pixel Scale(mm)'...
    'Write vel file (1 or 0)'...
    'Write RMS file(1 or 0)'...
    'Write flux file (1 or 0)'...
    'Determine wake-width (1 or 0)'...
    'e = 1 if v is +ve upwards, -1 otherwise'},'Input Parameters');
tic
pitch = str2double(input{1});
amp_seq = str2double(input{2});
f_pitch = str2double(input{3});
fps = str2double(input{4});
%use the string data of trial number directly - input{5}
pix_scale = str2double(input{6});
file = str2double(input{7});
rms = str2double(input{8});
flux = str2double(input{9});
wake = str2double(input{10});
e = str2double(input{11});
S = size(x);
N_max = S(1,1); 
N = (floor(N_max*f_pitch/fps))*(fps/f_pitch);%total number of frames
%total number of frames to be processed should be a multiple of FPS to
%ensure complete oscillation cycles are available for determining rms
%velocity

if pitch == 1
    pp = 0.125;
elseif pitch == 2
    pp = 0.3;
elseif pitch == 3
    pp = 0.48;
end

f1 = floor(f_pitch);
f2 = f_pitch - f1;
freq_deci = f2*100;
% A = amp/pix_scale; %amplitude in terms of pixels = amplitude(mm)/pixelscale(mm)
% dt = 1/fps;
% ds = U*dt/pix_scale;
C = 38.9; %chord length
rho = 1000; %density of water
dt = 1/fps;


%% Meanpath curvefit
X_pixels = 1024;
Y_pixels = 1024;
MeanPathName = strcat('P',num2str(pitch),'A',num2str(amp_seq),...
                  'freq',num2str(f1),'p',num2str(freq_deci),'MeanPos.dat');
mean_path = importdata(MeanPathName);
%edit this to write the file name automatically

%----------------------------------------------------
%curvefit 
x_curfit = mean_path(:,1)*pix_scale/C;
y_curfit = (Y_pixels-mean_path(:,2))*pix_scale/C; %shifting Y such that y increases towards top
%origin is the bottom left corner

[xData, yData] = prepareCurveData( x_curfit, y_curfit );

% Set up fittype and options.
ft = fittype( 'poly2' );
opts = fitoptions('Method','LinearLeastSquares');
opts.Robust = 'Bisquare';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
coef = coeffvalues(fitresult);
p1 = coef(1);
p2 = coef(2);
p3 = coef(3);

%% X_ref from data file
TipName = strcat('TipCoordP',...
                  num2str(pitch),'A',num2str(amp_seq),...
                  'freq',num2str(f1),'p',num2str(freq_deci),'Trial',input{5},'.dat');
tip_coord = importdata(TipName);
%autoread the file
x_tip = tip_coord(:,1)*pix_scale/C;
y_tip = (Y_pixels-tip_coord(:,2))*pix_scale/C;

%x_tip variation with time will give the forward velocity directly
time = (0:dt:dt*(length(x_tip)-1)); %in seconds

TipCoord(:,1) = time/f_pitch;
TipCoord(:,2) = x_tip;
TipCoord(:,3) = y_tip;

FileName1 = strcat('DataProcessed\TipCoordRP',...
                  num2str(pitch),'A',num2str(amp_seq),...
                  'f',num2str(f1),'p',num2str(freq_deci),'fps',input{4},'Trial',input{5},'.dat');
fileID = fopen(FileName1,'w');
       fprintf(fileID,'%5.4f %5.4f %5.4f\n',TipCoord');
       fclose(fileID);
%curvefit
[xData2, yData2] = prepareCurveData( time', x_tip );

% Set up fittype and options.
ft2 = fittype( 'poly1' );
opts2 = fitoptions( 'Method', 'LinearLeastSquares' );
opts2.Robust = 'LAR';

% Fit model to data.
[fitresult2, gof2] = fit( xData2, yData2, ft2, opts2 );
coef2 = coeffvalues(fitresult2);
U = coef2(1); %U_sp/C in mm/s
dx = U*dt; %incremental x/C

U_s = U*C/1000;%steady state velocity in m/s
C_m = C/1000;%chord length in m

X_ref = zeros(N,1);
X_ref(1:length(x_tip))= x_tip;

for n = length(x_tip)+1:N
    X_ref(n) = X_ref(n-1)+dx;
end

%% Frame segments
IntWinSize = 8; %interrogation window size in pixel
delX = 8*pix_scale/C;
delY = 8*pix_scale/C;
%need to determine how many frames are required such that increment in
%X_ref is equal to delX
%assuming ds = dx
%the assumption is true since the curvture is very small
q_frame = round(delX/dx); %q_frame is the total no if frames in a segment
num_frame_seg = floor(N/q_frame);
frame_seg = zeros(num_frame_seg,1);%this matrix is to write the actual number of frame, the first in every segment
%which from which the required velocity values will be extracted
frame_seg(1)= 1;%first element is the starting frame number


for n =1:num_frame_seg-1
    q = q_frame;
    thresh = 0.01;
    ratio = (X_ref(n*q+1)-X_ref((n-1)*q+1)-delX)/delX;
    %if the difference in the X increment to delX is less than or greater than 1% of delx, 
    %threre is no increment in q_frame
    if abs(ratio) <= thresh
        frame_seg(n+1)= q +frame_seg(n);
    elseif ratio < thresh
        frame_seg(n+1)= q +frame_seg(n)+1;
        q_frame = q+1;
    else
        frame_seg(n+1) = q + frame_seg(n)-1;
        q_frame = q-1;
    end
end

%% Velocity values at required all x at specific frame
%roughly 127 files will be written for every frame
%the file must contain u, v, vel_mag
X_pitch = X_ref + (1-pp); %the X-coordinate of the pitching point.
%the origin is still at the left bottom corner (y-increasing upwards)
%valid for small angles. For large angles, is geometric trnasformarion
%required??

X_EX = cell(length(frame_seg),1);
Y_EX = cell(length(frame_seg),1);
VEL = cell(length(frame_seg),6); % u(m/s),v(m/s),|U|(m/s),u*,v*,|U|*
FLUX = cell(length(frame_seg),8);% dM_dx,dPx_dx,dK_dx,dKx_dx,...
%dM_norm_dx,dPx_norm_dx,dK_norm_dx,dKx_norm_dx
GRADIENTS = cell(length(frame_seg),4);%w_z(1/s),divergence(1/s),w_z*,div*


for n = 1:length(frame_seg)
    I = frame_seg(n);
    % Position with respect to origin
    % Origin is the bottom left corner
    x_array = x{I,1}; % all x values as a square matrix at nth frame
    y_array = y{I,1}; %these values are in meters
    X = x_array(1,:)*1000/C; %converting to mm and non-dimensionalising with chord length
    Y = y_array(:,1)*1000/C;
    x_size = length(X);
    y_size = length(Y);
    
    %X and Y need to be transposed such that the origin is now the pitching
    %point. X_ref is the tip position. 
    x_exact = X-X_pitch(I); %x-position with x = 0 at the pitching point
    y_shift = p1*X.^2+p2*X+p3;
    
    X_EX(n,1) = {x_exact};
    
    y_exact = zeros(length(Y),x_size);
    
    for r = 1:x_size
        y_exact(:,r) = Y-y_shift(r);
    end
    Y_EX(n,1) = {y_exact};
%% conditionally writing the files   
    u = u_smoothed{I,1};
    v = e*v_smoothed{I,1};
    vel_mag = sqrt(u.^2+v.^2);
    E = (u.^2+v.^2)/2;

%Calculate from gradients
%Issue with the default calculation of PIV Lab
    [du_dx,du_dy] = gradient(u,delX*C/1000,delY*C/1000);
    [dv_dx,dv_dy] = gradient(v,delX*C/1000,delY*C/1000);
    vor = dv_dx - du_dy; 
    div = du_dx + dv_dy;
    
    u_norm = u_smoothed{I,1}/(U*C/1000);
    v_norm = e*v_smoothed{I,1}/(U*C/1000);
    vel_mag_norm = sqrt(u_norm.^2+v_norm.^2);
    E_norm = (u_norm.^2+v_norm.^2)/2;
    
    vor_norm = vor/U;
    div_norm = div/U;
    
    VEL(n,1)={u};
    VEL(n,2)={v};
    VEL(n,3)={vel_mag};
    VEL(n,4)={u_norm};
    VEL(n,5)={v_norm};
    VEL(n,6)={vel_mag_norm};
    
    GRADIENTS(n,1) = {vor};
    GRADIENTS(n,2) = {div};
    GRADIENTS(n,3) = {vor_norm};
    GRADIENTS(n,4) = {div_norm};
    
    A = size(u);
    columns = A(2);
    
    dM_dx = zeros(columns,1); %total momentum
    Mx_flux = zeros(columns,1);%momentum flux
    dK_dx = zeros(columns,1); %total Kinetic energy
    K_flux = zeros(columns,1);%kinetic energy flux
    
    dM_norm_dx = zeros(columns,1);
    Mx_fluxNorm = zeros(columns,1);
    dK_norm_dx = zeros(columns,1);
    K_fluxNorm = zeros(columns,1);
    
    for r = 1:columns
       dK_dx(r) = rho*trapz(E(:,r))*delY*C/1000; %integral over y -(u^2+v^2)/2
       dM_dx(r) = rho*trapz(u(:,r))*delY*C/1000; %mass flux (momentum)
       Mx_flux(r) = rho*trapz(U_s*u(:,r)-u(:,r).^2)*delY*C/1000; %momentum flux integral over y u^2
       K_flux(r) = rho*trapz(U_s*E(:,r)-E(:,r).*u(:,r))*delY*C/1000; %kinetic energy flux u*(u^2+v^2)/2
       %the units are SI
       
       dK_norm_dx(r) = trapz(E_norm(:,r))*delY; %integral over y -(u^2+v^2)/2
       dM_norm_dx(r) = trapz(u_norm(:,r))*delY; %mass flux
       Mx_fluxNorm(r) = trapz(u_norm(:,r)-u_norm(:,r).^2)*delY; %momentum flux integral over y u^2
       K_fluxNorm(r) = trapz(E_norm(:,r)-E_norm(:,r).*u_norm(:,r))*delY; %kinetic energy flux u*(u^2+v^2)/2
    end
    
    FLUX(n,1) = {dM_dx};
    FLUX(n,2) = {Mx_flux};
    FLUX(n,3) = {dK_dx};
    FLUX(n,4) = {K_flux};
    FLUX(n,5) = {dM_norm_dx};
    FLUX(n,6) = {Mx_fluxNorm};
    FLUX(n,7) = {dK_norm_dx};
    FLUX(n,8) = {K_fluxNorm};
    
end

%% rms and time average

L = x_size;
X_C = zeros(L,1);
delWakeAvg = zeros(L,1);
u_mAvg = zeros(L,1);
y_mAvg = zeros(L,1);

dMdx_avg = zeros(L,1);
MxFlux_avg = zeros(L,1);
dKdx_avg = zeros(L,1);
KxFlux_avg = zeros(L,1);

dMdx_rms = zeros(L,1);
MxFlux_rms = zeros(L,1);
dKdx_rms = zeros(L,1);
KxFlux_rms = zeros(L,1);

u_avg_field = zeros(L,L);
u_rms_field = zeros(L,L);
v_avg_field = zeros(L,L);
v_rms_field = zeros(L,L);
velMag_avg_field = zeros(L,L);
velMag_rms_field = zeros(L,L);

uNorm_avg_field = zeros(L,L);
uNorm_rms_field = zeros(L,L);
vNorm_avg_field = zeros(L,L);
vNorm_rms_field = zeros(L,L);
velMagNorm_avg_field = zeros(L,L);
velMagNorm_rms_field = zeros(L,L);

vor_rms_field = zeros(L,L);
div_rms_field = zeros(L,L);
vor_avg_field = zeros(L,L);
div_avg_field = zeros(L,L);

vorNorm_rms_field = zeros(L,L);
divNorm_rms_field = zeros(L,L);
vorNorm_avg_field = zeros(L,L);
divNorm_avg_field = zeros(L,L);

y_c_field = zeros(L,L);

uv_field = zeros(L,L);
uvNorm_field = zeros(L,L);

a = 1;  
x_req = -1; %initial x/c value for downstream distance
error = 0.5*delX;
eta = 0.15; %threshold for wake width calculations
u_ref = 0.2; %Background velocity in cm/s

for k = 1:L% k = 1:(x_size-10)10 being the cut-off column
%    x_req = X_C - delX; %X-C is already negative
     x_c = X_EX{a,1};
%      x_temp = round(100*X_EX{a,1});%for first frame
%      x_c = x_temp/100; % to ensure we have only upto two decimal places in x/C
     h = a;
    for r = 1:length(x_c)
       if x_c(r)<=(x_req+error)&& x_c(r)>(x_req-error)
           colmn = r;
           break
       end
    end
    
    if colmn == 9
        colmn = colmn +1;
    end
    
    if colmn == 10
        a = a+1;
    end    
    %colmn gives the column number corresponding to required x/c
    %one delX is added with every increment in frame segment i.e., the
    %column increments by 1 for every related frame
    
    X_C(k) = x_c(colmn); %exact value of x/c at which the rms is being calculated
    
    frame_max = length(frame_seg);
    dim = frame_seg(frame_max)-frame_seg(h);
    
    u_req = zeros(y_size,dim);
    v_req = zeros(y_size,dim);
    vel_mag_req = zeros(y_size,dim);
    uv_req = zeros(y_size,dim);
    y_req = zeros(y_size,dim);
    vor_req = zeros(y_size,dim);
    div_req = zeros(y_size,dim);
    m = zeros(dim,1);
    delT = zeros(dim,1);
    
    Int1 = zeros(1,dim);
    Int2 = zeros(1,dim);
    Int3 = zeros(1,dim);
    Int4 = zeros(1,dim);
    
  for n = (h+1):length(frame_seg)
      colmn_req = colmn + (n-(h+1));
          if colmn_req == (length(x_c)-1)
             break  %breaks the loop when the given position exits the field of view
          end
         
       I_seg0 = frame_seg(n-1);
       I_seg = frame_seg(n);
       Q = I_seg - I_seg0; 
       
     for  f = 1:(Q-1)  
          I0 = I_seg0 + f -1;
          I = I_seg0 + f;
        
          u = u_filtered{I,1};
          v = e*v_filtered{I,1};
          vel_mag = sqrt(u.^2+v.^2);
          uv = u.*v;
          E = (u.^2+v.^2)/2;
          
          u2 = u_smoothed{I,1};
          v2 = v_smoothed{I,1};
          [du_dx,du_dy] = gradient(u2,delX*C/1000,delY*C/1000);
          [dv_dx,dv_dy] = gradient(v2,delX*C/1000,delY*C/1000);
          vor = dv_dx - du_dy; 
          div = du_dx + dv_dy;
          
          index = n-h+f-1;
          
          y_coord = Y_EX{n,1};
          y_req(:,index) = y_coord(:,colmn_req);
          u_req(:,index) = u(:,colmn_req);
          v_req(:,index) = v(:,colmn_req);
          vel_mag_req(:,index) = vel_mag(:,colmn_req);
          uv_req(:,index) = uv(:,colmn_req);
          vor_req(:,index) = vor(:,colmn_req);
          div_req(:,index) = div(:,colmn_req);
          delT(index) = (I-I0)/fps;
        %determine the index where Y shifts from positive to negative or
        %it's equal to zero
        for p = 1:y_size
            if y_req(p,index)>=0 && y_req((p+1),index)<0
                m(index) = p;
                break
            end
        end   
     colmn_max = index; 
     % flux calculations
        %mass flux 
        Int1(1,index) = rho*trapz(u(:,colmn_req))*delY*C/1000; %y_coord is y/C
        %x-momentum flux
        Int2(1,index) = rho*trapz(U_s*u(:,colmn_req)-u(:,colmn_req).^2)*delY*C/1000;
        %kinetic energy
        Int3(1,index) = rho*trapz(E(:,colmn_req))*delY*C/1000;
        %Kinetic energy flux
        Int4(1,index) = rho*trapz(U_s*E(:,colmn_req)-E(:,colmn_req).*u(:,colmn_req))*delY*C/1000;
     end
  end
%truncation of the required arrays
    y_req = y_req(:,1:colmn_max);
    u_req = u_req(:,1:colmn_max);
    v_req = v_req(:,1:colmn_max);
    vel_mag_req = vel_mag_req(:,1:colmn_max);
    uv_req = uv_req(:,1:colmn_max);
    vor_req = vor(:,1:colmn_max);
    div_req = div(:,1:colmn_max);
    delT = delT(1:colmn_max);
    m = m(1:colmn_max);

    m_max = max(m);
    m_min = min(m);
    num = y_size +m_max - m_min;
%M Must be even
        if rem(num,2) ==0
             M = num;
        else
             M = num+1;
        end
%concatinating the matrices such that M/2the row has all the values at y~0

    u_con = zeros(M,colmn_max);
    v_con = zeros(M,colmn_max);
    vel_mag_con = zeros(M,colmn_max);
    uv_con = zeros(M,colmn_max);
    vor_con = zeros(M,colmn_max);
    div_con = zeros(M,colmn_max);
    y_con = zeros(M,colmn_max);
    for p =1:colmn_max
        for s = 1:M/2
            if m(p)-s+1 ==0
                break;
            else
                u_con((M/2+1-s),p) = u_req((m(p)-s+1),p);
                v_con((M/2+1-s),p) = v_req((m(p)-s+1),p);
                vel_mag_con((M/2+1-s),p) = vel_mag_req((m(p)-s+1),p);
                uv_con((M/2+1-s),p) = uv_req((m(p)-s+1),p);
                vor_con((M/2+1-s),p) = vor_req((m(p)-s+1),p);
                div_con((M/2+1-s),p) = div_req((m(p)-s+1),p);
                y_con((M/2+1-s),p) = y_req((m(p)-s+1),p);
            end
        end
        for s = 1:M/2
             if (m(p)+s)>y_size
                 break;
             else
                u_con((M/2+s),p) = u_req((m(p)+s),p);
                v_con((M/2+s),p) = v_req((m(p)+s),p);
                vel_mag_con((M/2+s),p) = vel_mag_req((m(p)+s),p);
                uv_con((M/2+s),p) = uv_req((m(p)+s),p);
                vor_con((M/2+s),p) = vor_req((m(p)+s),p);
                div_con((M/2+s),p) = div_req((m(p)+s),p);
                y_con((M/2+s),p) = y_req((m(p)+s),p);
             end            
        end
    end
    sum_delT = sum(delT);

    u_sq = u_con.^2;
    u_sum = u_sq*delT;
    u_rms = sqrt(u_sum/sum_delT);
    uNorm_rms = u_rms/(U*C/1000);

    v_sq = v_con.^2;
    v_sum = v_sq*delT;
    v_rms = sqrt(v_sum/sum_delT);
    vNorm_rms = v_rms/(U*C/1000);

    vel_mag_sq = vel_mag_con.^2;
    vel_mag_sum = vel_mag_sq*delT;
    vel_mag_rms = sqrt(vel_mag_sum/sum_delT);
    velMagNorm_rms = vel_mag_rms/(U*C/1000);
    
    vor_sq = vor_con.^2;
    vor_sum = vor_sq*delT;
    vor_rms = sqrt(vor_sum/sum_delT);
    vorNorm_rms = vor_rms/U;
    
    div_sq = div_con.^2;
    div_sum = div_sq*delT;
    div_rms = sqrt(div_sum/sum_delT);
    divNorm_rms = div_rms/U;

%dealing with y increments
    for p= 1:colmn_max
        delM1 = M/2 - m(p);
        for r = 1:delM1
             y_con(delM1+1-r,p)= y_req(1,p)+delX*r;
        end
        delM2 = M-y_size-delM1;
        for r = 1:delM2
        y_con(M-delM2+r,p)= y_req(y_size,p)-delX*r;
        end
    end

    y_vec = ones(colmn_max,1);
    y_avg = (y_con*y_vec)/colmn_max;


% time average velocity
    u_avg = (u_con*delT/sum_delT);
    uNorm_avg = u_avg/(U*C/1000);

    v_avg = (v_con*delT/sum_delT);
    vNorm_avg = v_avg/(U*C/1000);

    vel_mag_avg = (vel_mag_con*delT/sum_delT);
    velMagNorm_avg = vel_mag_avg/(U*C/1000);
    
    vor_avg = (vor_con*delT/sum_delT);
    vorNorm_avg = vor_avg/U;
    
    div_avg = (div_con*delT/sum_delT);
    divNorm_avg = div_avg/U;
    
    uv_avg = (uv_con*delT/sum_delT);
    uvNorm_avg = uv_avg/(U*C/1000)^2;
   
%wake width clculation    
[delWakeAvg(k),u_mAvg(k),y_mAvg(k)] = wakeVer1(y_avg,vel_mag_avg*100,eta,u_ref);

% Avg and rms of the fluxes
Int1 = Int1(1,1:colmn_max);
Int2 = Int2(1,1:colmn_max);
Int3 = Int3(1,1:colmn_max);
Int4 = Int4(1,1:colmn_max);
delT = delT(1:colmn_max,1);

dMdx_avg(k) = Int1*delT/sum(delT);
MxFlux_avg(k) = Int2*delT/sum(delT);
dKdx_avg(k) = Int3*delT/sum(delT);
KxFlux_avg(k) = Int4*delT/sum(delT);

dMdx_rms(k) = sqrt(Int1.^2*delT/sum(delT));
MxFlux_rms(k) = sqrt(Int2.^2*delT/sum(delT));
dKdx_rms(k) = sqrt(Int3.^2*delT/sum(delT));
KxFlux_rms(k) = sqrt(Int4.^2*delT/sum(delT));


%Field calculations
m1 = round((M-L)/2);
m2 = round((M+L)/2);
uv_field(:,k) = uv_avg(m1:m2-1);
uvNorm_field(:,k)= uvNorm_avg(m1:m2-1);
u_avg_field(:,k) = u_avg(m1:m2-1); %in m/s
u_rms_field(:,k) = u_rms(m1:m2-1);
v_avg_field(:,k) = v_avg(m1:m2-1);
v_rms_field(:,k) = v_rms(m1:m2-1);
velMag_avg_field(:,k) = vel_mag_avg(m1:m2-1);
velMag_rms_field(:,k) = vel_mag_rms(m1:m2-1);
uNorm_avg_field(:,k) = uNorm_avg(m1:m2-1);
uNorm_rms_field(:,k) = uNorm_rms(m1:m2-1);
vNorm_avg_field(:,k) = vNorm_avg(m1:m2-1);
vNorm_rms_field(:,k) = vNorm_rms(m1:m2-1);
velMagNorm_avg_field(:,k) = velMagNorm_avg(m1:m2-1);
velMagNorm_rms_field(:,k) = velMagNorm_avg(m1:m2-1);
vor_rms_field(:,k) = vor_rms(m1:m2-1);
div_rms_field(:,k) = div_rms(m1:m2-1);
vor_avg_field(:,k) = vor_avg(m1:m2-1);
div_avg_field(:,k) = div_avg(m1:m2-1);
vorNorm_rms_field(:,k) = vorNorm_rms(m1:m2-1);
divNorm_rms_field(:,k) = divNorm_rms(m1:m2-1);
vorNorm_avg_field(:,k) = vorNorm_avg(m1:m2-1);
divNorm_avg_field(:,k) = divNorm_avg(m1:m2-1);
y_c_field(:,k) = y_avg(m1:m2-1);

x_req = X_C(k) -delX;
end
%--------------------------------------------------------------------------------------------------
%% Writing Flux and wake files
Y_C = sum(y_c_field,2)/L;
u_mAvgNorm = u_mAvg/(U_s*100);%u_mAvg is in cm/s
%---------------------------------------------------------------------------------------------
                         
dMdx_avgNorm = dMdx_avg/U_s/C_m/rho;
MxFlux_avgNorm = MxFlux_avg/U_s^2/C_m/rho;
dKdx_avgNorm = dKdx_avg/U_s^2/C_m/rho;
KxFlux_avgNorm = KxFlux_avg/U_s^3/C_m/rho;

dMdx_rmsNorm = dMdx_rms/U_s/C_m/rho;
MxFlux_rmsNorm = MxFlux_rms/U_s^2/C_m/rho;
dKdx_rmsNorm = dKdx_rms/U_s^2/C_m/rho;
KxFlux_rmsNorm = KxFlux_rms/U_s^3/C_m/rho;

%% saving field variables for further averaging and processing
NameAvg = strcat('DataProcessed\MATfiles\AvgRMSTrial',input{5},...
               'RP', num2str(pitch),'A',num2str(amp_seq),'f',num2str(f1),...
               'p',num2str(freq_deci));
save(NameAvg,'X_C','Y_C',...
    'u_avg_field','u_rms_field',...
    'v_avg_field','v_rms_field',...
    'velMag_avg_field','velMag_rms_field',...
    'uNorm_avg_field','uNorm_rms_field',...
    'vNorm_avg_field','vNorm_rms_field',...
    'velMagNorm_avg_field','velMagNorm_rms_field',... 
    'vor_rms_field','div_rms_field','vor_avg_field','div_avg_field',...
    'vorNorm_rms_field','divNorm_rms_field','vorNorm_avg_field','divNorm_avg_field',...
    'u_mAvg','u_mAvgNorm','y_mAvg','delWakeAvg',...
    'dMdx_avg','MxFlux_avg','dKdx_avg','KxFlux_avg',...
    'dMdx_avgNorm','MxFlux_avgNorm','dKdx_avgNorm','KxFlux_avgNorm',...
    'dMdx_rms','MxFlux_rms','dKdx_rms','KxFlux_rms',...
    'dMdx_rmsNorm','MxFlux_rmsNorm','dKdx_rmsNorm','KxFlux_rmsNorm')

NameInst = strcat('DataProcessed\MATfiles\InstTrial',input{5},...
               'RP', num2str(pitch),'A',num2str(amp_seq),'f',num2str(f1),...
               'p',num2str(freq_deci));
save(NameInst,'X_EX','Y_EX','VEL','FLUX','GRADIENTS','frame_seg')
toc