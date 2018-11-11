function unew = solveWaveEqn(dims,source,c,s,t,uold,u,unew)
    %% Inject source and solve wave equation for one timestep
    u(s) = u(s) + source(t,:);
    
    rx = (dims.dt^2)/(dims.dx^2);
    ry = (dims.dt^2)/(dims.dy^2);
    
    for i=4:dims.ny-3
        for j=4:dims.nx-3
            unew(i,j) = (1/180)*ry*c(i,j)*c(i,j)*...
                     (2*(u(i-3,j)+u(i+3,j)) - ...
                     27*(u(i-2,j)+u(i+2,j)) + ...
                     270*(u(i-1,j)+u(i+1,j)) - ...
                     490*(u(i,j)))...
                     +...
                     (1/180)*rx*c(i,j)*c(i,j)*...
                     (2*(u(i,j-3)+u(i,j+3)) - ...
                     27*(u(i,j-2)+u(i,j+2)) + ...
                     270*(u(i,j-1)+u(i,j+1)) - ...
                     490*(u(i,j))) + ...
                     2*u(i,j)-uold(i,j);
        end
    end
end