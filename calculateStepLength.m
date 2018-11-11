function [stepLength,newErr] = calculateStepLength(dims,oldErr,gradient,source,model,trueRec)
%% Setup                     
    recording = zeros(dims.nt,length(dims.recPos),'single');
    stepLength = 256;
    newErr = inf;
    while (newErr > oldErr)
        newErr = 0;
        stepLength = stepLength/2;
        fprintf('  Current Step Length: %3.2f \n',stepLength);
        fprintf('    Original Error: %5.4f \n',oldErr);
        % Test model update
        modelt = model + stepLength*gradient;
        for s = 1:dims.ds:length(dims.srcPos)
            uold = zeros(dims.ny,dims.nx,'single');
            u = zeros(dims.ny,dims.nx,'single');
            unew = zeros(dims.ny,dims.nx,'single');
            for t = 1:dims.nt
                %  Solve wave equation using test model update
                srcPos = dims.srcPos(s);
                unew = solveWaveEqn(dims,source,modelt,srcPos,t,uold,u,unew);
                %  Record traces
                recording(t,:) = unew(dims.recPos);
                uold = u;
                u = unew;
            end
            %% Calculate new error and check against old
            chi = recording(:,:)-trueRec(:,:,s);
            newErr = newErr + norm(chi);
        end
        fprintf('    New Error: %5.4f \n',newErr);
    end
end

