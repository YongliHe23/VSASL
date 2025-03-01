function spins_dynamics(M)
for i = 1:3:size(M,1)
    mcol = M(i,:);
    plot3([0,M(i,1)],[0,M(i,2)],[0,M(i,3)])%,'Color',abs(mcol));
    hold on; plot3(M(i,1),M(i,2),M(i,3),'r*')%,'Color',abs(mcol));
    %arrow3([0,0,0],M(i,:))
    xlim([-1,1])
    ylim([-1,1])
    zlim([-1,1])
    grid on
    xlabel('Mx');ylabel("My");zlabel("Mz")
    title(sprintf("v=%d cm/s",30))
    if i==1
        pause;
    end
    pause(0.01)
end
end