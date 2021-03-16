clear;

% Rod Properties
D = 1.5/1000;
E = 2.07e11; %n/m^2
G = 7.93e10; %n/m^2

I = pi/4*(D/2)^4; %m^4
K = [E*I 0 0; 0 E*I 0; 0 0 2*G*I];

% Integration constants
p_0 = zeros(3,1);
R_0 = eye(3);
N = 100;
sSpan = [0, 1];

% Now we need a forcing function f : [0,1] -> R2
d = 3;
k = 80;

tData = linspace(0, 1, 100);
fData = zeros(2, 100);
[yy, C] = LsqFitVectorSpline(fData', tData, d, k);
C(:,40) = 25.*[1; 1];
C(:,80) = 7.*[-1; -1];

f = @(t) EvalVectorSpline(yy, C, d, t)';

% Solve and plot the rod
[s, state] = SolveStaticRod(p_0, R_0, D, E, G, N, sSpan, f);
PlotStaticRod(s, state, f);