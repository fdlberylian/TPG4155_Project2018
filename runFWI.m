clear;
%% Setting up dimensions
dims.dy =     10; % [m]
dims.dx =     10; % [m]
dims.dt = 1.0e-3; % [s]

dims.ny = 201; % Cells in y-direction
dims.nx = 301; % Cells in x-direction
dims.nt = 801; % Amount of time steps

%% Model dimensions
dims.modely = 100:150;
dims.modelx = 100:200;
dims.my = length(dims.modely);
dims.mx = length(dims.modelx);

%% Source locations
sx = min(dims.modelx):max(dims.modelx);
sy = min(dims.modely)*ones(1,length(sx));
dims.srcPos = sy + dims.ny*sx;

%% Receiver locations
rx = min(dims.modelx):max(dims.modelx);
ry = min(dims.modely)*ones(1,length(rx));
dims.recPos = ry+dims.ny*rx;

%% Creating background model
bg = zeros(dims.ny,dims.nx,'single');
bg(:) = 2.0e3;         % [m/s] - Background
bg(115:end,:) = 2.3e3; % [m/s] - Layer

%% Begin iteration
model = bg;     % Starting model
dims.ds = 10;   % Grid point distance between sources
maxIter = 10;   % Maximum number of iterations per frequency
freqs = [4,6,8,10,12];  % Frequencies to use in inversion

it = 1; tic;
for f = freqs
    fprintf('<<Frequency: %2.0f Hz>> \n',f);
    %% Generating ricker source signature wavelet 
    source = rickerWave(f,dims);
    %% Load true recording
    load (['trueRec_',num2str(f),'Hz.mat']);   
    for i = 1:maxIter
        %% Calculate gradient
        fprintf('Iteration Number: %2.0f \n',i);
        [gradient,err] = calculateGradient(dims,source,model,trueRec);
            
        %% Taper and plot gradient
        gradient = taperGradient(gradient);
        figure(1);
        imagesc(gradient(dims.modely,dims.modelx));
        title('Gradient');
        axis('image');
        colorbar();
        drawnow();

        %% Calculate step length
        [stepLength,err] = calculateStepLength(dims,err,gradient,source,model,trueRec);

        %% Update model
        model = model + stepLength*gradient;
        figure(2);
        imagesc(model(dims.modely,dims.modelx));
        title('Model');
        axis('image');
        colorbar();
        drawnow();
        
        errVec(it) = err;
        alpha(it) = stepLength;
        figure(3);
        subplot(2,1,1);
        semilogy(errVec,'r-');
        title('Error');
        xlim([0,length(freqs)*maxIter]);
        subplot(2,1,2);
        semilogy(alpha,'b-');
        title('Step Length');
        xlim([0,length(freqs)*maxIter]);
        drawnow();
        
        it = it + 1;
        if stepLength < 1
            toc
            fprintf('\n');
            break
        end
        toc
        fprintf('\n');
    end
end
fprintf('Run completed.');
