%close all
% figure(1)

% simulation step size in milliseconds.
% 1 usec per step. (or 0.001 msec)
dt = 1e-3;
showEvolution = 0;
isControl = 0;
    
T1 = 1700;  %ms
T2 = 165;   %ms  from Qin paper
    
% T1 = 1400;
% T2 = 110;

homogeneity = 'perfect'; % 'perfect'  % 'bad_B0_B1'
refocus_scheme = 'MLEV'; % 'MLEV' % 'DR180'
base_pulse =  'sinc_mod'   % 'sinc'  % 'hard'; % sech, hard, BIR

batch_mode = 0;
exportPulses = 0;

if batch_mode
    homogeneity ='';
    exportPulses = 0;
end
 
switch homogeneity
    case 'bad_B0'
        off_resonance =    100e-3 / 42576 ; % (KHz --> Tesla) ... gammabar is in kHz/T
        B1fudge = 1;
        
    case 'bad_B0_B1'
        off_resonance =   100e-3 / 42576 ; % (KHz --> Tesla) ... gammabar is in kHz/T
        B1fudge = 1.3;
        
    case 'perfect'
        off_resonance = 0%  100e-3 / 42576 ; % (KHz --> Tesla) ... gammabar is in kHz/T
        B1fudge = 1; % 1.1;
end




% default movement parameters
accel = 0 ; % in cm/ms^2
vel = 0  ; % in cm/ms
pos0 = 0 ;  % cm

Mzfinal = [];
Mzfinal_ns = [];

Nsegs = 9;

Gvs =   1.250*1e-4;  % T/cm
Gvs =   1.4*1e-4;  % T/cm   % 16900 and 17268
%Gvs =   3*1e-4;  % T/cm     % 12360
Gvs = 3*1e-4;

% experimental new pulse
%Gvs = 1.2 * 1e-4;
%Nsegs = 17;

gap = round(0.4/dt);
gap = round(0.2/dt); % (12360)
% gap = round(0.5/dt); %  17268
%gap = round(0.5/dt);% 16900 and 17268

pregap = 0; % 16900
%pregap = round(0.05/dt); % 16900
pregap = round(0.2/dt); % 12360
%pregap = round(0.08/dt);  % 17268 


ramp_len = 0.3/dt ;  % 16900

flat_len = 0.4/dt;   % 16900 
flat_len = 0.3/dt;  % 17268
flat_len = 0;        % 12360 
flat_len = 0.3/dt;  % VSS_mod_sinc

B1_area180 = 0.117e-4; % Tesla*ms


switch base_pulse
    
    case 'sinc'
        
        mySinc_dur = 3.6; % ms - 17268 pulse (9 x 0.4) 
        %mySinc_dur = 4; % ms
        
        mySinc_dur = round(mySinc_dur/dt/Nsegs) * Nsegs*dt; % ms
        mySinc = sinc(linspace(-1, 1,  mySinc_dur/dt))' ;
        mySinc = mySinc .* hanning(length(mySinc));  % 17268 pulse
        mySinc_area = (sum(mySinc) * dt);
        mySinc =  B1_area180 * mySinc / mySinc_area;
        B1pulse = mySinc;
        B1seg_len = round(mySinc_dur /dt/Nsegs);
        
    case 'sinc_mod'
        
        mySinc_dur = 1.44; % ms
        %mySinc_dur = 3.6; % ms - 17268 pulse (9 x 0.4) 

        mySinc_dur = round(mySinc_dur/dt/Nsegs) * Nsegs*dt; % ms
        mySinc = sinc(linspace(-0.8, 0.8,  mySinc_dur/dt))' ;
        mySinc = mySinc .* hanning(length(mySinc));  % 17268 pulse
        mySinc_area = (sum(mySinc) * dt);
        mySinc =  B1_area180 * mySinc / mySinc_area;
        %mySinc = mySinc* 0.52; % doing a 90 degree pulse only.
        B1pulse = mySinc;
        B1seg_len = round(mySinc_dur /dt/Nsegs);
        
    case 'BIR'
        myBIR_dur = 3; % ms
        myBIRpulse = genBIRpulse(0.2 , myBIR_dur)';
        myBIRpulse_area = (sum(myBIRpulse) * dt);
        myBIRpulse = 200e-7 *  myBIRpulse/max(abs(myBIRpulse)) ;
        B1pulse = myBIRpulse;
        B1seg_len = round(myBIR_dur/dt/Nsegs);
        
    case 'sech'
        
        mysech_dur = 3 ;
        mysech = genSech180(0.1, mysech_dur)';
        mysech_area = sum(mysech * dt);
        mysech =  B1_area180 * mysech / mysech_area;
        mysech = 130e-7 *  mysech/max(abs(mysech)) ;

        B1pulse = mysech;
        B1seg_len = round( mysech_dur /dt/Nsegs);
        
    case 'hard'
        
        hard180_dur = 1.44; % (ms) this makes it an even multiple of 4 us and 9 segments
        hard180_dur = 0.16 * 9 ; % (ms) in the paper
        
        hard180_dur = round(hard180_dur/Nsegs/dt) * Nsegs*dt;
        hard180 = ones(hard180_dur/dt);
        
        hard180_area = sum(hard180)*dt;
        hard180 = B1_area180* hard180 /hard180_area;
        B1pulse = hard180;        
        B1seg_len = round(hard180_dur /dt/Nsegs);

end


% refocuser pulse is 180 hard pulse for 0.8 ms
B1ref_dur = 0.8;   % ms
B1ref_dur = 1.0;   % ms   - 12360
%B1ref_dur = 0.5;   % ms  - 16900 pulse and 17280

B1ref_len = B1ref_dur/dt;  % points
B1ref_max = (1/B1ref_dur)* 0.1175e-4; % Tesla ... amplitude of 1 ms hard 180
B1ref = B1ref_max * ones(B1ref_len,1) ; % * exp(i*pi/2);
B1ref_area = sum(abs(B1ref));

% make Refocuser 180 be a composite pulse (90x-180y-90x)
%tmp90 = B1ref(1:end/2);  % Jia's version of the 17268 uses 
                         %  B1ref(1:(end/2-1))
% the whole composite pulse:                         
%B1ref = [tmp90; B1ref*exp(i*pi/2); tmp90]; %  16900 and 17268

% Now shorten the pulses and make them taller
% B1ref = 2*B1ref(1:2:end);

B1ref_len = length(B1ref); 
B1ref_dur = B1ref_len * dt;


%%


trap = [
    zeros(pregap,1);
    [0:ramp_len-1]'/(ramp_len-1);
    ones(flat_len,1);
    [ramp_len-1:-1:0]'/(ramp_len-1);
    zeros(gap,1);
    ]*Gvs;

trap_len = length(trap);

str = sprintf('Durations: B1seg= %0.2f  ms, B1ref=  %0.2f ms , Gtrap=%0.2f ms\n', ...
    B1seg_len*dt , B1ref_len*dt, trap_len*dt);

% another trapezoid for the crusher
trap2 = [
    zeros(pregap,1);
    [0:ramp_len-1]'/(ramp_len-1);
    ones(2/dt,1);
    [ramp_len-1:-1:0]'/(ramp_len-1);
    zeros(gap,1);
    ];


switch refocus_scheme
    case 'MLEV'
        % The MLEV-16 sequence
        mlev_phases =  repmat( [1 1 -1 -1  -1  1 1 -1 -1 -1  1 1  1  -1 -1  1], [1,4]);
        %mlev_phases =  repmat( [1 1 -1 -1  -1  1 1 -1 -1 -1  1 1  1  -1 -1  -1], [1,4]);  % trying something different?
        
    case 'DR180'
        % My scheme for refocusing the off-res.
        mlev_phases = repmat( exp(-i*pi/2)*[1 -1 -1 1], [1,16]);
end
Gz = [];
B1 = [];
for n=1:Nsegs-1
    B1seg = B1pulse(1 + B1seg_len*(n-1) : n*B1seg_len);
    
    if base_pulse=='sinc_mod'
        B1seg(:) = mean(B1seg);
    end
        
    B1 = [B1 ;
        B1seg;
        zeros(size(trap));
        B1ref * mlev_phases(2*n-1);
        zeros(size(trap));
        zeros(size(trap));
        B1ref * mlev_phases(2*n) ;
        zeros(size(trap));
        ];

    Gz = [Gz;
        zeros(size(B1seg)) ;
        trap;
        zeros(size(B1ref))
        -trap;
        trap;
        zeros(size(B1ref))
        -trap;
        ];
    
end

if isControl
    Gz = abs(Gz);
    %Gz(:)  = 0;
end

Gz = [
    Gz;
    zeros(size(B1seg))
    ];

B1seg = B1pulse(1 + B1seg_len*(Nsegs-1) : Nsegs*B1seg_len);

if base_pulse=='sinc_mod'
    B1seg(:) = mean(B1seg);
end

B1 = [
    B1;
    B1seg;
    ]*B1fudge;
%
%put a crusher at the end
%{
Gz = [
        Gz;
        Gvs*trap2;
    ];

% allow for the crusher:
B1 = [B1; 
    zeros(size(trap2))
    ];
%}

% control case:
if isControl
    Gz = abs(Gz);
end

if exportPulses
% write pulses to put on the scanner
    genScannerPulses(B1, Gz, dt);
end

Bx = real(B1);
By = imag(B1);

NSTEPS = length(B1);


% total duration of the simulation interval (in ms)
duration = NSTEPS*dt;  % ms.
vel_range =[-100:100]*1e-3;
vel_range = linspace(0, 20, 200) * 1e-3;

for vel = vel_range  % cm / msec
    
    t = linspace(0,duration, NSTEPS)'; % mseconds.
    zpos = pos0 + vel*t + 0.5*accel*(t.^2);
    Bz = zpos.*Gz ;
    
    Bz = Bz + off_resonance;
        
    beff = [Bx By  Bz];
    Mi = [0 0 1]';
    M = blochsim(Mi, beff, T1, T2, dt, NSTEPS);
    if showEvolution
        plot(M(:,3))
        drawnow
    end
    
    Mzfinal=[Mzfinal; M(end,3)];
end


% now do the control case
%{
for vel = vel_range  % cm / msec
    
    Bz = zpos.*abs(Gz) ;

    Bz = Bz + off_resonance;
        
    beff = [Bx By  Bz];
    Mi = [0 0 1]';
    M = blochsim(Mi, beff, T1, T2, dt, NSTEPS);
    if showEvolution
        plot(M(:,3))
        drawnow
    end
    
    Mzfinal_ns=[Mzfinal_ns; M(end,3)];
end
%}

if batch_mode==0
    figure(3)
    t=[0:length(Gz)-1]*dt;
    
    subplot(311)
    area(t, Gz *1e4);
    grid on
    xlabel('time (ms)')
    ylabel('G_z (G/cm)')
    title('VSI Pulse')
    
    subplot(312)
    area(t, abs(B1)*1e7);
    grid on
    xlabel('time (ms)')
    ylabel('B_1 amplitude (mG)')
    
    subplot(313)
    area(t, angle(B1));
    grid on
    xlabel('time (ms)')
    ylabel('B_1 phase (rad)')
    
    figure(4)
    hold on
    plot(vel_range*1e3, Mzfinal)
    plot(vel_range*1e3, Mzfinal_ns)
    axis([min(vel_range)*1e3 max(vel_range)*1e3, 0 1])
    fatlines
    grid on
    xlabel('Velocity (cm/s)')
    ylabel('M_z')
    title('Velocity Profiles')
    %legend('Selective Case', 'Non-selective Case', 'Location', 'SouthEast')
    print -dpng ~/FTVS_profile

    %text(-90 , -0.6, {refocus_scheme ;base_pulse });
    hold off
%%
end
%{
Gvenc = 0.3*1e2; % Gauss/m
d = 0.001;
D = 0.01;
gamma = 2*pi*4.257e3; % rad/s/gauss

Beta = pi*gamma*Gvenc*d*D;
%Beta = pi/2;
MzBIR = sinc(vel_range*Beta*2);
figure(4)
hold on
plot(vel_range*1e3,MzBIR ,'--');
area(vel_range*1e3,MzBIR.*Mzfinal' );  
legend('FT-VS profile', 'BIR8-VS Profile', 'Net Profile', 'Location', 'NorthEast')
fatlines
print -dpng ~/Net_profile
hold off
%}

%%  Calculating labeling efficency over a range of values
%{
figure(8)
efficiency = (Mzfinal - Mzfinal_ns)/2;
plot(vel_range, efficiency)
inds = find((vel_range < 0.08)  & (vel_range > 0.005));
hold on
plot(vel_range(inds), efficiency(inds))
mean(efficiency(inds))
xlabel('Velocity (cm/s)')
ylabel('Efficiency')
title('Velocity Profile of VSI pulse')
legend('Selective Case', 'Non-selective Case', 'Location', 'SouthEast')
hold off

%}


