function [y1] = myNeuralNetworkFunction(x1)
%MYNEURALNETWORKFUNCTION neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 04-May-2020 18:52:04.
%
% [y1] = myNeuralNetworkFunction(x1) takes these arguments:
%   x = 2xQ matrix, input #1
% and returns:
%   y = 1xQ matrix, output #1
% where Q is the number of samples.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-3.25157134765256;-4.38500173436354];
x1_step1.gain = [0.307543612943305;0.228050080838826];
x1_step1.ymin = -1;

% Layer 1
b1 = [-26.796264010132933464;-2.7262124158017551068;4.0500436341483512237;0.012877830133416071426;0.00055192618069286193249;0.0049745280048556503663;-2.7917150101672971729;0.0001297780767057100243;4.2706907451820717014;-0.0056575212854522959061];
IW1_1 = [12.749956233331735334 -28.28949130098636644;21.717952646372634007 24.108227793933217242;-5.5155504931261685186 3.7067804103857819875;-0.10596198878348789263 -7.3576993953266764947;3.7447512279103349897 8.2523822249181666422;-37.640616228743667193 -68.957651139118240735;-22.360741126044896276 -24.807211894645401884;3.6128812442480988665 7.7012385877052018657;5.8733443191535279482 -3.884889614013987913;4.0885340153417430997 5.5568074266156948227];

% Layer 2
b2 = 0.18008631786138254438;
LW2_1 = [0.18764631590283842311 0.67888575522960981079 0.15599952375755360423 0.23968105905632847152 -19.995101466178510918 0.58881143231842214547 -0.66254218256109509433 22.47003965764509914 -0.14769616473724656025 -3.0043817036038151791];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.681724448436753;
y1_step1.xoffset = -1.46686832530809;

% ===== SIMULATION ========

% Dimensions
Q = size(x1,2); % samples

% Input 1
xp1 = mapminmax_apply(x1,x1_step1);

% Layer 1
a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*xp1);

% Layer 2
a2 = repmat(b2,1,Q) + LW2_1*a1;

% Output 1
y1 = mapminmax_reverse(a2,y1_step1);
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
y = bsxfun(@minus,x,settings.xoffset);
y = bsxfun(@times,y,settings.gain);
y = bsxfun(@plus,y,settings.ymin);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
a = 2 ./ (1 + exp(-2*n)) - 1;
end

% Map Minimum and Maximum Output Reverse-Processing Function
function x = mapminmax_reverse(y,settings)
x = bsxfun(@minus,y,settings.ymin);
x = bsxfun(@rdivide,x,settings.gain);
x = bsxfun(@plus,x,settings.xoffset);
end