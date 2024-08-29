function [ g, T1, T2, T3, T4, Ts_in, Ts_out ] = solve_NS_DWR_ll(Re1, Re2, NN, Y, H, G1, G2, R1, R2, R6, ringW, stepW, ringSubs, upperBC, R5_adim)
% This function solves the Navier-Stokes equations with no-slip
% and Boussinesq_Scriven boundary conditions with second order centered 
% finite differences for the DWR configuration

% INPUTS:
% Re1: Reynold number for lower bulk phase
% Re2: Reynold number for upper bulk phase
% NN: "Boussinesq" number
% Y: Bulk Viscosity ratio
% H: Vertical gap for bulk phases
% G1: inner interfacial gap
% G2: outer interfacial gap
% R1: Cup lower internal radius
% R2: Cup upper internal radius
% R6: Ring external radius
% ringW: ring width
% stepW: step width
% ringSubs: ring subintervals
% upperBC: upper countour boundary condition
% R5_adim: adimensional ring internal radius

% OUTPUT:
% g: velocity field (nondimensional)
% T1: lower internal bulk drag
% T2: upper internal bulk drag
% T3: upper external bulk drag
% T4: lower external bulk drag
% Ts_in: internal interfacial drag
% Ts_out: external interfacial drag

Sres = ringW/ringSubs;
stepSubs = round(stepW/Sres);
ringSubs_2 = ringSubs/2;

NG1 = round(G1/Sres);% internal gap subintervals
NG2 = round(G2/Sres);% external gap subintervals
N1 = NG1 + NG2 + ringSubs + 2*stepSubs;% upper bulk horizontal subintervals
N2 = NG1 + NG2 + ringSubs;% lower bulk horizontal subintervals
M = round(H/Sres);% bulk phases vertical subintervals
dim = (N1+1)*(M+1) + (N2+1)*M;% Dimension of the square coefficients matrix A 
% Adimensional parameters
R2_adim = R2/R6;
R1_adim = R1/R6;
Sres_adim = Sres/R6;
AA = 1/Sres_adim;

%% ACUMULATING INDEX AND VALUES
%% Upper boundary condition
if strcmp(upperBC, 'ns')
    % no-slip (closed contour)
    rows = (2:N1);
    cols = (2:N1);
    coefs = ones(1, N1-1);
elseif strcmp(upperBC, 'fb')
    % free interface nodes (open contour)
    rows = [(2:N1) (2:N1) (2:N1) (2:N1)];
    cols = [(2:N1)-1 (2:N1) (2:N1)+1 (2:N1)+(N1+1)];
    rind = ((2:N1)-1);
    coefs = [1 - 1./(2*(AA*R2_adim + rind)) ....
            -1i*(Re2/(AA*AA)) - 4 - 1./((AA*R2_adim + rind).*(AA*R2_adim + rind)) ....
            1 + 1./(2*(AA*R2_adim + rind)) ....
            2*ones(1, N1-1)];
else
    warning("Choose 'il' for free interface or 'cc' for no-slip condition at the upper boundary")
    return
end

%% Walls
% Upper left wall
rows = [rows (0:M)*(N1+1)+1];
cols = [cols (0:M)*(N1+1)+1];
coefs = [coefs ones(1, M+1)];
% Upper right wall
rows = [rows (1:M+1)*(N1+1)];
cols = [cols (1:M+1)*(N1+1)];
coefs = [coefs ones(1, M+1)];
% Lower left wall
rows = [rows (0:M-1)*(N2+1)+(N1+1)*(M+1)+1];
cols = [cols (0:M-1)*(N2+1)+(N1+1)*(M+1)+1];
coefs = [coefs ones(1, M)];
% Lower right wall
rows = [rows (1:M)*(N2+1)+(N1+1)*(M+1)];
cols = [cols (1:M)*(N2+1)+(N1+1)*(M+1)];
coefs = [coefs ones(1, M)];
% Ground
rows = [rows (2:N2)+(M-1)*(N2+1)+(N1+1)*(M+1)];
cols = [cols (2:N2)+(M-1)*(N2+1)+(N1+1)*(M+1)];
coefs = [coefs ones(1, N2-1)];
% Left step
rows = [rows (1:stepSubs)+(M)*(N1+1)+1];
cols = [cols (1:stepSubs)+(M)*(N1+1)+1];
coefs = [coefs ones(1, stepSubs)];
% Right step
rows = [rows (0:stepSubs-1)+(M+1)*(N1+1)-stepSubs];
cols = [cols (0:stepSubs-1)+(M+1)*(N1+1)-stepSubs];
coefs = [coefs ones(1, stepSubs)];

%% Interface
% Inner nodes
rows = [rows (2:NG1)+M*(N1+1)+stepSubs ....
             (2:NG1)+M*(N1+1)+stepSubs ....
             (2:NG1)+M*(N1+1)+stepSubs ....
             (2:NG1)+M*(N1+1)+stepSubs ....
             (2:NG1)+M*(N1+1)+stepSubs];
cols = [cols (2:NG1)+(M-1)*(N1+1)+stepSubs ....
             (2:NG1)+M*(N1+1)+stepSubs-1 ....
             (2:NG1)+M*(N1+1)+stepSubs ....
             (2:NG1)+M*(N1+1)+stepSubs+1 ....
             (2:NG1)+(M+1)*(N1+1)];
rind = (((2:NG1)+stepSubs)-1);
coefs = [coefs (1/(AA*Y))*ones(1, NG1-1) ....
               (1 - 1./(2*(AA*R2_adim + rind))).*(NN + (1/(2*AA))*(1 + 1/Y)) ....
               -NN*(2 + 1./((AA*R2_adim + rind).*(AA*R2_adim + rind))) - (1/(2*AA)).*(1i*(Re1/(AA*AA)) + 4 + 1./((AA*R2_adim + rind).*(AA*R2_adim + rind)) + (1/Y).*(1i*(Re2/(AA*AA)) + 4 + 1./((AA*R2_adim + rind).*(AA*R2_adim + rind)))) ....
               (1 + 1./(2*(AA*R2_adim + rind))).*(NN + (1/(2*AA))*(1+1/Y)) ....
               (1./(AA))*ones(1, NG1-1)];
% Outer nodes
rows = [rows (1:NG2-1)+NG1+M*(N1+1)+stepSubs+ringSubs+1 ....
             (1:NG2-1)+NG1+M*(N1+1)+stepSubs+ringSubs+1 ....
             (1:NG2-1)+NG1+M*(N1+1)+stepSubs+ringSubs+1 ....
             (1:NG2-1)+NG1+M*(N1+1)+stepSubs+ringSubs+1 ....
             (1:NG2-1)+NG1+M*(N1+1)+stepSubs+ringSubs+1];
cols = [cols (1:NG2-1)+NG1+(M-1)*(N1+1)+stepSubs+ringSubs+1 ....
             (1:NG2-1)+NG1+M*(N1+1)+stepSubs+ringSubs ....
             (1:NG2-1)+NG1+M*(N1+1)+stepSubs+ringSubs+1 ....
             (1:NG2-1)+NG1+M*(N1+1)+stepSubs+ringSubs+2 ....
             (1:NG2-1)+NG1+(M+1)*(N1+1)+stepSubs+1];
rind = (1:NG2-1)+NG1+stepSubs+ringSubs;
coefs = [coefs (1/(AA*Y))*ones(1, NG2-1) ....
               (1 - 1./(2*(AA*R2_adim + rind))).*(NN + (1/(2*AA))*(1+1/Y)) ....
               -NN*(2 + (1./((AA*R2_adim + rind).*(AA*R2_adim + rind)))) - (1/(2*AA)).*(1i*(Re1/(AA*AA)) + 4 + 1./((AA*R2_adim + rind).*(AA*R2_adim + rind)) + (1/Y).*(1i*(Re2/(AA*AA)) + 4 + 1./((AA*R2_adim + rind).*(AA*R2_adim + rind)))) ....
               (1 + 1./(2*(AA*R2_adim + rind))).*(NN + (1/(2*AA))*(1+1/Y)) ....
               (1./(AA))*ones(1, NG2-1)];

%% Upper Bulk
[P, Q] = ndgrid(2:M-ringSubs_2, 2:N1);
rows = [rows reshape(((P-1)*(N1+1)+Q)', 1, [])....
             reshape(((P-1)*(N1+1)+Q)', 1, [])....
             reshape(((P-1)*(N1+1)+Q)', 1, [])....
             reshape(((P-1)*(N1+1)+Q)', 1, [])....
             reshape(((P-1)*(N1+1)+Q)', 1, [])];
cols = [cols reshape(((P-2)*(N1+1)+Q)', 1, [])....
             reshape(((P-1)*(N1+1)+Q-1)', 1, [])....
             reshape(((P-1)*(N1+1)+Q)', 1, [])....
             reshape(((P-1)*(N1+1)+Q+1)', 1, [])....
             reshape((P*(N1+1)+Q)', 1, [])];
rind = ((2:N1)-1);
coefs = [coefs ones(1, (N1-1)*(M-ringSubs_2-1)) ....
               repmat(1 - 1./(2*(AA*R2_adim + rind)), [1 M-ringSubs_2-1]) ....
               repmat(-1i*(Re2/(AA*AA)) - 4 - 1./((AA*R2_adim + rind).*(AA*R2_adim + rind)), [1 M-ringSubs_2-1]) ....
               repmat(1 + 1./(2*(AA*R2_adim + rind)), [1 M-ringSubs_2-1]) ....
               ones(1, (N1-1)*(M-ringSubs_2-1))];
for k = 1:ringSubs_2
    % Inner
    rows = [rows (1:NG1+stepSubs+ringSubs_2-k)+(M-ringSubs_2+(k-1))*(N1+1)+1 ....
                 (1:NG1+stepSubs+ringSubs_2-k)+(M-ringSubs_2+(k-1))*(N1+1)+1 ....
                 (1:NG1+stepSubs+ringSubs_2-k)+(M-ringSubs_2+(k-1))*(N1+1)+1 ....
                 (1:NG1+stepSubs+ringSubs_2-k)+(M-ringSubs_2+(k-1))*(N1+1)+1 ....
                 (1:NG1+stepSubs+ringSubs_2-k)+(M-ringSubs_2+(k-1))*(N1+1)+1];
    cols = [cols (1:NG1+stepSubs+ringSubs_2-k)+(M-ringSubs_2+(k-1)-1)*(N1+1)+1 ....
                 (1:NG1+stepSubs+ringSubs_2-k)+(M-ringSubs_2+(k-1))*(N1+1)....
                 (1:NG1+stepSubs+ringSubs_2-k)+(M-ringSubs_2+(k-1))*(N1+1)+1....
                 (1:NG1+stepSubs+ringSubs_2-k)+(M-ringSubs_2+(k-1))*(N1+1)+2....
                 (1:NG1+stepSubs+ringSubs_2-k)+(M-ringSubs_2+1+(k-1))*(N1+1)+1];
    rind = (1:NG1+stepSubs+ringSubs_2-k);
    coefs = [coefs ones(1, NG1+stepSubs+ringSubs_2-k) ....
                   1 - 1./(2*(AA*R2_adim + rind)) ....
                   -1i*(Re2/(AA*AA)) - 4 - 1./((AA*R2_adim + rind).*(AA*R2_adim + rind)) ....
                   1 + 1./(2*(AA*R2_adim + rind)) ....
                   ones(1, NG1+stepSubs+ringSubs_2-k)];
    % Outer
    rows = [rows (1:NG2+stepSubs+ringSubs_2-k)+(M-ringSubs_2+k-1)*(N1+1)+NG1+stepSubs+ringSubs_2+k ....
                 (1:NG2+stepSubs+ringSubs_2-k)+(M-ringSubs_2+k-1)*(N1+1)+NG1+stepSubs+ringSubs_2+k ....
                 (1:NG2+stepSubs+ringSubs_2-k)+(M-ringSubs_2+k-1)*(N1+1)+NG1+stepSubs+ringSubs_2+k ....
                 (1:NG2+stepSubs+ringSubs_2-k)+(M-ringSubs_2+k-1)*(N1+1)+NG1+stepSubs+ringSubs_2+k ....
                 (1:NG2+stepSubs+ringSubs_2-k)+(M-ringSubs_2+k-1)*(N1+1)+NG1+stepSubs+ringSubs_2+k];
    cols = [cols (1:NG2+stepSubs+ringSubs_2-k)+(M-ringSubs_2+k-2)*(N1+1)+NG1+stepSubs+ringSubs_2+k ....
                 (1:NG2+stepSubs+ringSubs_2-k)+(M-ringSubs_2+k-1)*(N1+1)+NG1+stepSubs+ringSubs_2+k-1 ....
                 (1:NG2+stepSubs+ringSubs_2-k)+(M-ringSubs_2+k-1)*(N1+1)+NG1+stepSubs+ringSubs_2+k ....
                 (1:NG2+stepSubs+ringSubs_2-k)+(M-ringSubs_2+k-1)*(N1+1)+NG1+stepSubs+ringSubs_2+k+1 ....
                 (1:NG2+stepSubs+ringSubs_2-k)+(M-ringSubs_2+k)*(N1+1)+NG1+stepSubs+ringSubs_2+k];
    rind = (1:NG2+stepSubs+ringSubs_2-k)+NG1+stepSubs+ringSubs_2;
    coefs = [coefs ones(1, NG2+stepSubs+ringSubs_2-k) ....
                   1 - 1./(2*(AA*R2_adim + rind)) ....
                   -1i*(Re2/(AA*AA)) - 4 - 1./((AA*R2_adim + rind).*(AA*R2_adim + rind)) ....
                   1 + 1./(2*(AA*R2_adim + rind)) ....
                   ones(1, NG2+stepSubs+ringSubs_2-k)];
end

%% Lower Bulk
for k = 1:ringSubs_2
    % Inner
    rows = [rows (1:NG1+k-1)+(N1+1)*(M+1)+(N2+1)*(k-1)+1 ....
                 (1:NG1+k-1)+(N1+1)*(M+1)+(N2+1)*(k-1)+1 ....
                 (1:NG1+k-1)+(N1+1)*(M+1)+(N2+1)*(k-1)+1 ....
                 (1:NG1+k-1)+(N1+1)*(M+1)+(N2+1)*(k-1)+1 ....
                 (1:NG1+k-1)+(N1+1)*(M+1)+(N2+1)*(k-1)+1];
    cols = [cols (1:NG1+k-1)+(N1+1)*M+stepSubs+(N2+1)*(k-1)+(k>1)*stepSubs+1 ....
                 (1:NG1+k-1)+(N1+1)*(M+1)+(N2+1)*(k-1) ....
                 (1:NG1+k-1)+(N1+1)*(M+1)+(N2+1)*(k-1)+1 ....
                 (1:NG1+k-1)+(N1+1)*(M+1)+(N2+1)*(k-1)+2 ....
                 (1:NG1+k-1)+(N1+1)*(M+1)+(N2+1)*(k)+1];
    rind = 1:NG1+k-1;
    coefs = [coefs ones(1, NG1+k-1) ....
                   1 - 1./(2*(AA*R1_adim + rind)) ....
                   -1i*(Re1/(AA*AA)) - 4 - 1./((AA*R1_adim + rind).*(AA*R1_adim + rind)) ....
                   1 + 1./(2*(AA*R1_adim + rind)) ....
                   ones(1, NG1+k-1)];
    % Outer
    rows = [rows (1:NG2-1+k)+(N1+1)*(M+1)+NG1+ringSubs+N2*(k-1) ....
                 (1:NG2-1+k)+(N1+1)*(M+1)+NG1+ringSubs+N2*(k-1) ....
                 (1:NG2-1+k)+(N1+1)*(M+1)+NG1+ringSubs+N2*(k-1) ....
                 (1:NG2-1+k)+(N1+1)*(M+1)+NG1+ringSubs+N2*(k-1) ....
                 (1:NG2-1+k)+(N1+1)*(M+1)+NG1+ringSubs+N2*(k-1)];
    cols = [cols (1:NG2-1+k)+(N1+1)*M+stepSubs+NG1+ringSubs+N2*(k-1)+(k>1)*stepSubs ....
                 (1:NG2-1+k)+(N1+1)*(M+1)+NG1+ringSubs-1+N2*(k-1) ....
                 (1:NG2-1+k)+(N1+1)*(M+1)+NG1+ringSubs+N2*(k-1) ....
                 (1:NG2-1+k)+(N1+1)*(M+1)+NG1+ringSubs+1+N2*(k-1) ....
                 (1:NG2-1+k)+(N1+1)*(M+1)+NG1+ringSubs+1+N2*k];
    rind = (1:NG2-1+k)+ringSubs+NG1-k;
    coefs = [coefs ones(1, NG2-1+k) ....
                   1 - 1./(2*(AA*R1_adim + rind)) ....
                   -1i*(Re1/(AA*AA)) - 4 - 1./((AA*R1_adim + rind).*(AA*R1_adim + rind)) ....
                   1 + 1./(2*(AA*R1_adim + rind)) ....
                   ones(1, NG2+k-1)];
end
[P2, Q2] = ndgrid(1+ringSubs_2:M-1, 2:N2);
rows = [rows (M+1)*(N1+1)+reshape(((P2-1)*(N2+1)+Q2)', 1, [])....
             (M+1)*(N1+1)+reshape(((P2-1)*(N2+1)+Q2)', 1, [])....
             (M+1)*(N1+1)+reshape(((P2-1)*(N2+1)+Q2)', 1, [])....
             (M+1)*(N1+1)+reshape(((P2-1)*(N2+1)+Q2)', 1, [])....
             (M+1)*(N1+1)+reshape(((P2-1)*(N2+1)+Q2)', 1, [])];
cols = [cols (M+1)*(N1+1)+reshape(((P2-2)*(N2+1)+Q2)', 1, [])....
             (M+1)*(N1+1)+reshape(((P2-1)*(N2+1)+Q2-1)', 1, [])....
             (M+1)*(N1+1)+reshape(((P2-1)*(N2+1)+Q2)', 1, [])....
             (M+1)*(N1+1)+reshape(((P2-1)*(N2+1)+Q2+1)', 1, [])....
             (M+1)*(N1+1)+reshape((P2*(N2+1)+Q2)', 1, [])];
rind = 1:N2-1;
coefs = [coefs ones(1, (N2-1)*(M-ringSubs_2-1)) ....
               repmat(1 - 1./(2*(AA*R1_adim + rind)), [1 M-ringSubs_2-1]) ....
               repmat(-1i*(Re1/(AA*AA)) - 4 - 1./((AA*R1_adim + rind).*(AA*R1_adim + rind)), [1 M-ringSubs_2-1]) ....
               repmat(1 + 1./(2*(AA*R1_adim + rind)), [1 M-ringSubs_2-1]) ....
               ones(1, (N2-1)*(M-ringSubs_2-1))];

gDWR = [];
indDWR = [];
% Ring (no-slip)
for k = 1:ringSubs_2+1
    rows = [rows (1:(2*k-1))+(M-ringSubs_2+(k-1))*(N1+1)+NG1+stepSubs+ringSubs_2-(k-1)];
    cols = [cols (1:(2*k-1))+(M-ringSubs_2+(k-1))*(N1+1)+NG1+stepSubs+ringSubs_2-(k-1)];
    coefs = [coefs ones(1, length(1:(2*k-1)))];
    gDWR = [gDWR R2_adim+((0:(2*k-2))+NG1+stepSubs+ringSubs_2-(k-1))./AA];
    indDWR = [indDWR (1:(2*k-1))+(M-ringSubs_2+(k-1))*(N1+1)+NG1+stepSubs+ringSubs_2-(k-1)];
end
for k = 1:ringSubs_2
    rows = [rows (1:(2*(ringSubs_2-k)+1))+NG1+1+(N1+1)*(M+1)+(N2+2)*(k-1)];
    cols = [cols (1:(2*(ringSubs_2-k)+1))+NG1+1+(N1+1)*(M+1)+(N2+2)*(k-1)];
    coefs = [coefs ones(1, length((1:(2*(ringSubs_2-k)+1))))];
    gDWR = [gDWR R1_adim+((0:(2*(ringSubs_2-k)))+NG1+1+(k-1))./AA];
    indDWR = [indDWR (1:(2*(ringSubs_2-k)+1))+NG1+1+(N1+1)*(M+1)+(N2+2)*(k-1)];
end
%% FILL A
A = sparse(rows, cols, coefs);
%% FILL b with no-slip boundary conditions over the DWR surfaces
b = sparse(indDWR, ones(1, length(indDWR)), gDWR, dim, 1);
% b = sparse((1:Nb+1)+(M/2)*(N+1), ones(1, Nb+1), (0:Nb)*(1/(R1_adim*N)), dim, 1);
%% Solving the linear system of equations
g = A\b;
%% Torque calculation
% indices
gpM = zeros(ringSubs + 5);
for k = 1:ringSubs_2+3
    gpM(k, :) = (1:5+ringSubs)+NG1+stepSubs+(M-ringSubs_2-2+(k-1))*(N1+1)-2;
end
for k = 1:ringSubs_2+2
    gpM(k+ringSubs_2+3, :) = (1:5+ringSubs)+NG1+(N1+1)*(M+1)+(N2+1)*(k-1)-2;
end
gM = full(g(gpM));

% Quadrants
gp1 = zeros(3, ringSubs_2+1);
gp2 = gp1; gp3 = gp1; gp4 = gp1;

ind = 1:floor(((ringSubs+2)/2+1)/2);
if rem((ringSubs+2)/2+1, 2) ~= 0
    ind = [ind floor(((ringSubs+2)/2+1)/2)];
end
ind = [ind floor(((ringSubs+2)/2+1)/2)-1:-1:1];

for k = 1:ringSubs_2+1
    d = diag(gM, -ringSubs_2 + (k-1)*2);
    gp1(:, k) = d(ind(k):ind(k)+2);
    d = diag(rot90(gM), -ringSubs_2 + (k-1)*2);
    gp2(:, k) = d(ind(k):ind(k)+2);
    d = diag(rot90(gM, 2), -ringSubs_2 + (k-1)*2);
    gp3(:, k) = d(ind(k):ind(k)+2);
    d = diag(rot90(gM, 3), -ringSubs_2 + (k-1)*2);
    gp4(:, k) = d(ind(k):ind(k)+2);
end

r1 = R2_adim + (stepSubs+NG1:stepSubs+NG1+ringSubs_2)*Sres_adim;
T1_fo = (1/(sqrt(2)))*trapz(r1.*r1.*(gp1(2, :) - gp1(3, :))); % first order
T1_so = (1/(2*sqrt(2)))*trapz(r1.*r1.*(-3*gp1(3, :) + 4*gp1(2, :) - gp1(1, :))); % second order
T1 = [T1_fo T1_so];
r2 = R2_adim + (stepSubs+NG1+ringSubs_2:stepSubs+NG1+ringSubs)*Sres_adim;
T2_fo = (1/(sqrt(2)))*trapz(r2.*r2.*(gp2(2, :) - gp2(3, :)));
T2_so = (1/(2*sqrt(2)))*trapz(r2.*r2.*(-3*gp2(3, :) + 4*gp2(2, :) - gp2(1, :)));
T2 = [T2_fo T2_so];
r3 = flip(r2);
T3_fo = (1/(sqrt(2)))*trapz(r3.*r3.*(gp3(2, :) - gp3(3, :)));
T3_so = (1/(2*sqrt(2)))*trapz(r3.*r3.*(-3*gp3(3, :) + 4*gp3(2, :) - gp3(1, :)));
T3 = [T3_fo T3_so];
r4 = flip(r1);
T4_fo = (1/(sqrt(2)))*trapz(r4.*r4.*(gp4(2, :) - gp4(3, :)));
T4_so = (1/(2*sqrt(2)))*trapz(r4.*r4.*(-3*gp4(3, :) + 4*gp4(2, :) - gp4(1, :)));
T4 = [T4_fo T4_so];

gsin = gM(ringSubs_2+3, 1:3);
rin = R2_adim + (stepSubs+NG1-2:stepSubs+NG1)*Sres_adim;
g_r_in = gsin./rin;
Ts_in_fo = R5_adim*R5_adim*R5_adim*(1/(Sres_adim))*(g_r_in(2) - g_r_in(3));
Ts_in_so = R5_adim*R5_adim*R5_adim*(1/(2*Sres_adim))*(-3*g_r_in(3) + 4*g_r_in(2) - g_r_in(1));
Ts_in = [Ts_in_fo Ts_in_so];

gsout = gM(ringSubs_2+3, ringSubs+3:ringSubs+5);
rout = R2_adim + (stepSubs+NG1+ringSubs:stepSubs+NG1+ringSubs+2)*Sres_adim;
g_r_out = gsout./rout;
Ts_out_fo = (1/(Sres_adim))*(g_r_out(2) - g_r_out(1));
Ts_out_so = (1/(2*Sres_adim))*(-3*g_r_out(1) + 4*g_r_out(2) - g_r_out(3));
Ts_out = [Ts_out_fo Ts_out_so];
end