function [stepLength,newErr] = calculateStepLength_v2(dims,oldErr,gradient,source,model,trueRec)
%% Setup                     
    recordinga = zeros(dims.nt,length(dims.recPos),'single');
    recordingb = zeros(dims.nt,length(dims.recPos),'single');
    recordingc = zeros(dims.nt,length(dims.recPos),'single');
    recordingd = zeros(dims.nt,length(dims.recPos),'single');
    a = 0; b = 256; eps = 10^-2; del = 10^-2;
    r1 = (sqrt(5)-1)/2; r2 = r1^2; h = b-a;
    c = a+r2*h; d = a+r1*h;
    erra = 0; errb = 0; errc = 0; errd = 0;
    modelta = model + a*gradient;
    modeltb = model + b*gradient;
    modeltc = model + c*gradient;
    modeltd = model + d*gradient;
        for s = 1:dims.ds:length(dims.srcPos)
            uolda = zeros(dims.ny,dims.nx,'single');
            ua = zeros(dims.ny,dims.nx,'single');
            unewa = zeros(dims.ny,dims.nx,'single');
            uoldb = zeros(dims.ny,dims.nx,'single');
            ub = zeros(dims.ny,dims.nx,'single');
            unewb = zeros(dims.ny,dims.nx,'single');
            uoldc = zeros(dims.ny,dims.nx,'single');
            uc = zeros(dims.ny,dims.nx,'single');
            unewc = zeros(dims.ny,dims.nx,'single');
            uoldd = zeros(dims.ny,dims.nx,'single');
            ud = zeros(dims.ny,dims.nx,'single');
            unewd = zeros(dims.ny,dims.nx,'single');
            for t = 1:dims.nt
                %  Solve wave equation using test model update
                srcPos = dims.srcPos(s);
                unewa = solveWaveEqn(dims,source,modelta,srcPos,t,uolda,ua,unewa);
                unewb = solveWaveEqn(dims,source,modeltb,srcPos,t,uoldb,ub,unewb);
                unewc = solveWaveEqn(dims,source,modeltc,srcPos,t,uoldc,uc,unewc);
                unewd = solveWaveEqn(dims,source,modeltd,srcPos,t,uoldd,ud,unewd);
                %  Record traces
                recordinga(t,:) = unewa(dims.recPos);
                recordingb(t,:) = unewb(dims.recPos);
                recordingc(t,:) = unewc(dims.recPos);
                recordingd(t,:) = unewd(dims.recPos);
                uolda = ua; ua = unewa;
                uoldb = ub; ub = unewb;
                uoldc = uc; uc = unewc;
                uoldd = ud; ud = unewd;
            end
            %% Calculate new error and check against old
            chia = recordinga(:,:)-trueRec(:,:,s);
            chib = recordingb(:,:)-trueRec(:,:,s);
            chic = recordingc(:,:)-trueRec(:,:,s);
            chid = recordingd(:,:)-trueRec(:,:,s);
            erra = erra + norm(chia); errb = errb + norm(chib);
            errc = errc + norm(chic); errd = errd + norm(chid);
        end
       stepLength = b; newErr = errb;
    while (abs(erra-errb)>eps)||(h>del)
        stepLength = b; newErr = errb;
        fprintf('  Current Step Length: %3.2f \n',stepLength);
        fprintf('    Original Error: %5.4f \n',oldErr);
        
        if(errc < errd)
            b = d; errb = errd;
            d = c; errd = errc;
            errc = 0;
            h = b-a; c = a+r2*h;
            modeltc = model + c*gradient;
            for s = 1:dims.ds:length(dims.srcPos)
                uoldc = zeros(dims.ny,dims.nx,'single');
                uc = zeros(dims.ny,dims.nx,'single');
                unewc = zeros(dims.ny,dims.nx,'single');
                for t = 1:dims.nt
                    %  Solve wave equation using test model update
                    srcPos = dims.srcPos(s);
                    unewc = solveWaveEqn(dims,source,modeltc,srcPos,t,uoldc,uc,unewc);
                    %  Record traces
                    recordingc(t,:) = unewc(dims.recPos);
                    uoldc = uc; uc = unewc;
                end
                %% Calculate new error and check against old
                chic = recordingc(:,:)-trueRec(:,:,s);
                errc = errc + norm(chic);
            end
        else
            a = c; erra = errc;
            c = d; errc = errd;
            errd = 0;
            h = b-a; d = a+r1*h;
            modeltd = model + d*gradient;
            for s = 1:dims.ds:length(dims.srcPos)
                uoldd = zeros(dims.ny,dims.nx,'single');
                ud = zeros(dims.ny,dims.nx,'single');
                unewd = zeros(dims.ny,dims.nx,'single');
                for t = 1:dims.nt
                    %  Solve wave equation using test model update
                    srcPos = dims.srcPos(s);
                    unewd = solveWaveEqn(dims,source,modeltd,srcPos,t,uoldd,ud,unewd);
                    %  Record traces
                    recordingd(t,:) = unewd(dims.recPos);
                    uoldd = ud; ud = unewd;
                end
                %% Calculate new error and check against old
                chid = recordingd(:,:)-trueRec(:,:,s);
                errd = errd + norm(chid);
            end
        end
        fprintf('    New Error: %5.4f \n',newErr);
    end
end