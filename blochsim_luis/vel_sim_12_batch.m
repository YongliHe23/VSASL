close all
% figure(1)

% simulation step size in milliseconds.
% 1 usec per step. (or 0.001 msec)
dt = 1e-3;
showEvolution = 0;
exportPulses = 0;

off_resonance_range =   [-200:40:200] * 1e-3/ 42576;
B1fudge_range = [0.8:0.05:1.2] ;

B0_sim = [];
B1_sim = [];

for n_omega = 1 : length(off_resonance_range)
    
    off_resonance = off_resonance_range(n_omega)
    B1fudge = 1;
    
    vel_sim_12
    
 drawnow
    
    B0_sim = [B0_sim ; Mzfinal'];
end
%
for n_b1 = 1 :length(B1fudge_range)
    
    off_resonance = 0;
    B1fudge =B1fudge_range(n_b1) ;
    vel_sim_12
    drawnow
    B1_sim = [B1_sim ; Mzfinal'];
end
%           

%%
figure(4)
subplot(211)
imagesc(B0_sim);
caxis([-1 1])
title(sprintf('Off-resonance Effects. B1 scale = %0.2f',  B1fudge_range(floor(end/2)+1)));
colormap jet
colorbar
xlabel('Velocity');
set(gca, 'XTick', [1 length(vel_range)/2+1 length(vel_range) ])
set(gca, 'XTickLabel', [vel_range(1) 0 vel_range(end)]*1e3)

ylabel('off-resonance');
set(gca, 'YTick', [1 length(off_resonance_range)/2+1 length(off_resonance_range) ])
set(gca, 'YTickLabel', [-200 0 200])

subplot(212)
imagesc(squeeze(B1_sim));
caxis([-1 1])
title(sprintf('B1 Effects. Off resonance= %0.2f',  off_resonance_range(floor(end/2)+1)));
colormap jet
colorbar
xlabel('Velocity');
set(gca, 'XTick', [1 length(vel_range)/2+1 length(vel_range) ])
set(gca, 'XTickLabel', [vel_range(1) 0 vel_range(end)]*1e3)
ylabel('B1 factor');
set(gca, 'YTick', [1 length(B1fudge_range)/2+1 length(B1fudge_range) ])
set(gca, 'YTickLabel', [B1fudge_range(1) 0 B1fudge_range(end)])

print -dpng imperfect_profiles.png

%%
figure(5)
subplot(211)
imagesc(B0_sim - repmat(B0_sim(ceil(end/2),:), size(B0_sim,1),1));
caxis([-1 1])
title(sprintf('Off-resonance Error. B1 scale = %0.2f',  B1fudge_range(floor(end/2)+1)));
colormap jet
colorbar
xlabel('Velocity');
set(gca, 'XTick', [1 length(vel_range)/2+1 length(vel_range) ])
set(gca, 'XTickLabel', [vel_range(1) 0 vel_range(end)]*1e3)

ylabel('off-resonance');
set(gca, 'YTick', [1 length(off_resonance_range)/2+1 length(off_resonance_range) ])
set(gca, 'YTickLabel', [-200 0 200])

subplot(212)
imagesc(B1_sim - repmat(B1_sim(ceil(end/2),:), size(B1_sim,1),1));
caxis([-1 1])
title(sprintf('B1 Error. Off resonance= %0.2f',  off_resonance_range(floor(end/2)+1)));
colormap jet
colorbar
xlabel('Velocity');
set(gca, 'XTick', [1 length(vel_range)/2+1 length(vel_range) ])
set(gca, 'XTickLabel', [vel_range(1) 0 vel_range(end)]*1e3)
ylabel('B1 factor');
set(gca, 'YTick', [1 length(B1fudge_range)/2+1 length(B1fudge_range) ])
set(gca, 'YTickLabel', [B1fudge_range(1) 0 B1fudge_range(end)])

print -dpng imperfect_profiles_diff.png
