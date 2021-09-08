function Vec = TridiagMatSolve(MainDiag, OffDiagConst, RHSVec)
	% Implement Thomas Algorithm with constant off-diagonal elements
	%
	% INPUTS:
	% MainDiag
	% Type: Vector of complex doubles
	% This is the leading diagonal of the matrix A in the system to solve Ax = b
	%
	% OffDiagConst
	% Type: Complex double
	% This is the value of all elements in the matrix A next to the leading diagonal
    %
	% RHSVec
	% Type: Vector of complex doubles
	% This is the b-vector in the system to solve
	%
	% OUTPUTS:
    % Vec
	% Type: Vector of complex doubles
    % This is the solution x in the system Ax = b
	% 
N = length(RHSVec);
Vec = zeros(N,1);
cprime = zeros(N-1,1);
dprime = zeros(N,1);
cprime(1) = OffDiagConst / MainDiag(1);
	
for k = 2 : N - 1
    cprime(k) = OffDiagConst / (MainDiag(k) - OffDiagConst * cprime(k - 1));
end

dprime(1) = RHSVec(1) / MainDiag(1);
for k = 2 : N
	dprime(k) = (RHSVec(k) - OffDiagConst * dprime(k-1)) / (MainDiag(k) - OffDiagConst * cprime(k-1));
end

Vec(N) = dprime(N);

for k = 1 : N - 1
    Vec(N-k) = dprime(N-k) - cprime(N-k) * Vec(N-k+1);
end
end