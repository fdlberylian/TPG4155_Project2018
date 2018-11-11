function [gradient, err] = calculateGradient(dims,source,model,trueRec) 
%% Setup
    recording = zeros(dims.nt,length(dims.recPos),'single');
    gradient  = zeros(dims.ny,dims.nx,'single');
    forwardField = zeros(dims.my,dims.mx,dims.nt,'single'); 
    adjointField = zeros(dims.my,dims.mx,dims.nt,'single');
    err = 0;
    for s = 1:dims.ds:length(dims.srcPos)
        %% Run forward simulation on background model
        uold = zeros(dims.ny,dims.nx,'single');
        u = zeros(dims.ny,dims.nx,'single');
        unew = zeros(dims.ny,dims.nx,'single');
        for t = 1:dims.nt
            % Solve wave equation
            srcPos = dims.srcPos(s);
            unew = solveWaveEqn(dims,source,model,srcPos,t,uold,u,unew);
            % Record traces
            recording(t,:) = unew(dims.recPos);
            % Save forward field for use in correlation
            forwardField(:,:,t) = unew(dims.modely,dims.modelx);
            uold = u;
            u = unew;
        end
        %% Calculate difference and error
        chi = recording(:,:)-trueRec(:,:,s);
        err = err + norm(chi);
        %% Run adjoint simulation
        chi = flipud(chi);
        uold = zeros(dims.ny,dims.nx,'single');
        u = zeros(dims.ny,dims.nx,'single');
        unew = zeros(dims.ny,dims.nx,'single');
        for t = 1:dims.nt
            % Solve wave equation using the difference (chi) as sources
            recPos = dims.recPos;
            unew = solveWaveEqn(dims,chi,model,recPos,t,uold,u,unew);
            % Save adjoint field for use in correlation
            adjointField(:,:,dims.nt-t+1) = unew(dims.modely,dims.modelx);
            uold = u;
            u = unew;
        end
        %% Correlate
        for t = 2:dims.nt-1
            % Calculate the time derivative of the displacement to
            % gradient.
            dadj = zeros(dims.my,dims.mx,'single');
            dfor = zeros(dims.my,dims.mx,'single');
            for i=1:length(dims.modely)
                for j=1:length(dims.modelx)
                    dadj(i,j) = (adjointField(i,j,t+1)-adjointField(i,j,t))/dims.dt;
                    dfor(i,j) = (forwardField(i,j,t)-forwardField(i,j,t-1))/dims.dt;
                    ii = dims.modely(i);
                    jj = dims.modelx(j);
                    gradient(ii,jj) = gradient(ii,jj) + dadj(i,j)*dfor(i,j);
                end
            end
        end
    end
end