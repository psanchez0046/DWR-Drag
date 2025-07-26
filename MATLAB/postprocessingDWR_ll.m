function [resultStruct, timeElapsedTotal] = postprocessingDWR_ll(geom, mesh, Bulk, iteParams, IO)
timerTotal = tic;
% This code calculates the rheological properties from amplitude ratio
% and phase lag data of a rotational rheometer with a DWR fixture

% INPUTS (all units are in the International System)
% geom: geometry parameters
% mesh: mesh parameters
% bulk: bulk phases physical parameters
% iteParams: iterative scheme parameters
% iteMax: maximum number of iterations
% IO: Input/Output data parameters

% OUTPUTS(optional)
% resultStruct: strut with results
% timeElapsedTotal: time elapsed in processing experiental files

resultStruct = struct();

expFilenames = GetFilenames('_exp.txt', IO.inputFilepath);

% Initializing optional output data
iterationsTimesData = cell([1 length(expFilenames)]);
iterationsData = cell([1 length(expFilenames)]);
bouData = cell([1 length(expFilenames)]);
bouOmegaData = cell([1 length(expFilenames)]);
FreqData = cell([1 length(expFilenames)]);
AmpData = cell([1 length(expFilenames)]);
TorqData = cell([1 length(expFilenames)]);
etasData = cell([1 length(expFilenames)]);
etasData_linear = cell([1 length(expFilenames)]);
GData = cell([1 length(expFilenames)]);
GData_linear = cell([1 length(expFilenames)]);
ARcalcData = cell([1 length(expFilenames)]);
deltaARcalcData = cell([1 length(expFilenames)]);

% Prompting the number of input data files to process
fprintf('Experiment files: %s \n\n', num2str(length(expFilenames)))

% for-end looping on each input data file (*_exp.txt file)
for ite = 1:length(expFilenames)
    fprintf('Analyzing file %s...\n', num2str(ite))
    expFileData = importdata(char(expFilenames(ite)));
    % Prompting the number of lines to process in the file
    [filasexpFiledata, ~] = size(expFileData);
    fprintf('Lines: %s\n\n', num2str(filasexpFiledata))    
    % for-end looping on each line under analysis
    timeElapsedIT = zeros(filasexpFiledata, 1);
    Bou_final = zeros(filasexpFiledata, 1);
    lambda_final = zeros(filasexpFiledata, 1);
    ARcalc_final = zeros(filasexpFiledata, 1);
    errorAR_final = zeros(filasexpFiledata, 1);
    delta_AR_final = zeros(filasexpFiledata, 1);
    freq = expFileData(:, IO.colIndexFreq);
    amp = expFileData(:, IO.colIndexAmp);
    Torq = expFileData(:, IO.colIndexTorq);
    G_linear = zeros(filasexpFiledata, 1);
    etas_linear = zeros(filasexpFiledata, 1);
    for lin = 1:filasexpFiledata        
        omegarad = 2*pi*freq(lin);
        % Displaying the number of the line under analysis
        fprintf('Analyzing line %s from file %s...\n', num2str(lin), num2str(ite))
        Re1 = (Bulk.rho_bulk1*omegarad*geom.R6*geom.R6)/Bulk.eta_bulk1;
        Re2 = (Bulk.rho_bulk2*omegarad*geom.R6*geom.R6)/Bulk.eta_bulk2;
        Y = Bulk.eta_bulk1/Bulk.eta_bulk2;
        ARexp = expFileData(lin, IO.colIndexAR)*(cos(expFileData(lin, IO.colIndexDelta)) + 1i*sin(expFileData(lin, IO.colIndexDelta)));
        % Intializing variables...
        lambda=1;
        Bou=[];
        ARcalc=[];
        errorAR=[];
        Bou(lambda) = 0;
        timerVal = tic;
        % Solving the Navier-Stokes equation
        NN(lambda) = Bou(lambda)*(Bulk.eta_bulk1 + Bulk.eta_bulk2)/Bulk.eta_bulk1;
        R5_adim = geom.R5/geom.R6;
        [~, T1, T2, T3, T4, Ts_in, Ts_out] = solve_NS_DWR_ll(Re1, Re2, NN(lambda), Y, geom.H, geom.G1, geom.G2, geom.R1, geom.R2, geom.R6, geom.ringW, geom.stepW, mesh.ringSubs, mesh.upperBC, R5_adim);
        % Calculating AR
        if mesh.DOrder == 1
            T1 = T1(1); T2 = T2(1); T3 = T3(1); T4 = T4(1); Ts_in = Ts_in(1); Ts_out = Ts_out(1);
        else
            T1 = T1(2); T2 = T2(2); T3 = T3(2); T4 = T4(2); Ts_in = Ts_in(2); Ts_out = Ts_out(2);
        end
        Cll = 1i*omegarad*2*pi*geom.R6*geom.R6*geom.R6;
        Tdown = Cll*Bulk.eta_bulk1*(T1 + T2);
        Tup = Cll*Bulk.eta_bulk2*(T3 + T4);
        Ts = Cll*NN(lambda)*Bulk.eta_bulk1*(Ts_in + Ts_out);
        if geom.ICorrected
            ARcalc(lambda) = -Tdown - Tup - Ts;
        else
            ARcalc(lambda) = -Tdown - Tup - Ts - geom.inertia*omegarad*omegarad;
        end
        % Calculating first value of the relative error 
        errorAR(lambda) = abs((ARcalc(lambda))/(ARexp)-1);

        % while loop performing the successive iterations
        while (lambda < iteParams.iteMax) && (errorAR(lambda) > iteParams.tolMin) 
            lambda = lambda + 1;
            if geom.ICorrected
                Bou(lambda) = (-ARexp - Tdown - Tup)/(Cll*(Bulk.eta_bulk1 + Bulk.eta_bulk2)*(Ts_in + Ts_out));
            else
                Bou(lambda) = (-ARexp - Tdown - Tup - geom.inertia*omegarad*omegarad)/(Cll*(Bulk.eta_bulk1 + Bulk.eta_bulk2)*(Ts_in + Ts_out));
            end       
            % Solving the Navier-Stokes equation
            NN(lambda) = Bou(lambda)*(Bulk.eta_bulk1 + Bulk.eta_bulk2)/Bulk.eta_bulk1;
            % R5_adim = geom.R5/geom.R6;
            [~, T1, T2, T3, T4, Ts_in, Ts_out] = solve_NS_DWR_ll(Re1, Re2, NN(lambda), Y, geom.H, geom.G1, geom.G2, geom.R1, geom.R2, geom.R6, geom.ringW, geom.stepW, mesh.ringSubs, mesh.upperBC, R5_adim);
            % Calculating AR
            if mesh.DOrder == 1
                T1 = T1(1); T2 = T2(1); T3 = T3(1); T4 = T4(1); Ts_in = Ts_in(1); Ts_out = Ts_out(1);
            else
                T1 = T1(2); T2 = T2(2); T3 = T3(2); T4 = T4(2); Ts_in = Ts_in(2); Ts_out = Ts_out(2);
            end
            Tdown = Cll*Bulk.eta_bulk1*(T1 + T2);
            Tup = Cll*Bulk.eta_bulk2*(T3 + T4);
            Ts = Cll*NN(lambda)*Bulk.eta_bulk1*(Ts_in + Ts_out);
            if geom.ICorrected
                ARcalc(lambda) = -Tdown - Tup - Ts;
            else
                ARcalc(lambda) = -Tdown - Tup - Ts - geom.inertia*omegarad*omegarad;
            end
                
            % Calculating first value of the relative error 
            errorAR(lambda) = abs((ARcalc(lambda))/(ARexp)-1);          
        end  
        % Exit of the iterative process
        if lambda == iteParams.iteMax               
            fprintf('Limit reached\n') % Maximum # of iterations reached
        else
            fprintf('OK convergence at %s iterations\n', num2str(lambda)) % Iterations have converged!!!
        end    
        timeElapsedIT(lin) = toc(timerVal);
        % Displaying the time used in the iterative process for each line
        fprintf('Iterative process time = %s s\n\n', num2str(timeElapsedIT(lin)))
        Bou_final(lin) = Bou(lambda);
        ARcalc_final(lin) = ARcalc(lambda);
        delta_AR_final(lin) = angle(ARcalc(lambda));
        errorAR_final(lin) = errorAR(lambda);
        lambda_final(lin) = lambda;
        G_linear(lin, 1) = ARexp/(4*pi*((geom.R5^2*geom.R1^2/(geom.R5^2 - geom.R1^2)) + (geom.R6^2*geom.R3^2/(geom.R3^2 - geom.R6^2))));
        etas_linear(lin, 1) = G_linear(lin, 1)/(1i*omegarad);
        G_linear(lin, 2) = (ARexp - ARcalc(1))/(4*pi*((geom.R5^2*geom.R1^2/(geom.R5^2 - geom.R1^2)) + (geom.R6^2*geom.R3^2/(geom.R3^2 - geom.R6^2))));
        etas_linear(lin, 2) = G_linear(lin, 2)/(1i*omegarad);
    end
    % Calculating variables depending on the N number
    eta_s_final = Bou_final*geom.R6.*(Bulk.eta_bulk1 + Bulk.eta_bulk2);% converged viscoelasticity
    G_complex = 1i*Bou_final*2*pi.*freq*geom.R6.*(Bulk.eta_bulk1 + Bulk.eta_bulk2);% converged dynamic surface moduli    
    Lb1 = sqrt((Bulk.eta_bulk1/Bulk.rho_bulk1)/omegarad);
    Lb2 = sqrt((Bulk.eta_bulk2/Bulk.rho_bulk2)/omegarad);
    Bou_omega = eta_s_final./(Lb1*Bulk.eta_bulk1 + Lb2*Bulk.eta_bulk2);

    % Exporting results to the output data file
    results = [freq real(G_complex) imag(G_complex) real(eta_s_final) imag(eta_s_final) real(Bou_final) imag(Bou_final) real(Bou_omega) imag(Bou_omega) abs(ARcalc_final) delta_AR_final timeElapsedIT lambda_final];
    if ispc
        writematrix(results, strcat(IO.outputFilepath,'\',strrep(char(expFilenames(ite)),'exp','out')), 'delimiter','\t')
    else
        writematrix(results, strcat(IO.outputFilepath,'/',strrep(char(expFilenames(ite)),'exp','out')), 'delimiter','\t')
    end    
    % MATLAB output data
    FreqData{ite} = freq;
    AmpData{ite} = amp;
    TorqData{ite} = Torq;
    GData{ite} = G_complex;
    GData_linear{ite} = G_linear;
    etasData{ite} = eta_s_final;
    etasData_linear{ite} = etas_linear;
    bouData{ite} = Bou_final;
    bouOmegaData{ite} = Bou_omega;
    ARcalcData{ite} = abs(ARcalc_final);
    deltaARcalcData{ite} = delta_AR_final;
    iterationsTimesData{ite} = timeElapsedIT;
    iterationsData{ite} = lambda_final;
end
resultStruct.FreqData = FreqData;
resultStruct.AmpData = AmpData;
resultStruct.TorqData = TorqData;
resultStruct.GData = GData;
resultStruct.GData_linear = GData_linear;
resultStruct.etasData = etasData;
resultStruct.etasData_linear = etasData_linear;
resultStruct.bouData = bouData;
resultStruct.bouOmegaData = bouOmegaData;
resultStruct.ARcalcData = ARcalcData;
resultStruct.deltaARcalcData = deltaARcalcData;
resultStruct.iterationsTimesData = iterationsTimesData;
resultStruct.iterationsData = iterationsData;
resultStruct.geom = geom;
resultStruct.mesh = mesh; 
resultStruct.Bulk = Bulk; 
resultStruct.iteParams = iteParams;
timeElapsedTotal = toc(timerTotal);
fprintf('Total postprocessing program time = %s s\n', num2str(timeElapsedTotal))
end
