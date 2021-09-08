function Vec = TridiagMatMul(MainDiag, OffDiagConst, RHSVec)
	% Pre-multiply a vector by a tridiagonal matrix with constant off-diagonal elements
	%
	% INPUTS:
	% MainDiag
	% Type: Vector of complex doubles
	% This is the leading diagonal of the tridiagonal matrix
	%
	% OffDiagConst
	% Type: Complex double
	% This is the value of all elements in the tridiagonal matrix next to the leading diagonal
	%
	% RHSVec
	% Type: Vector of complex doubles
	% This is the vector to be pre-multiplied by the tridiagonal matrix
	%
	% OUTPUTS:
    % Vec
	% Type: Vector of complex doubles
	% This is the product of the tridiagonal matrix and the vector
	% 
Vec = zeros(length(RHSVec),1);
Vec(1) = MainDiag(1) * RHSVec(1) + OffDiagConst * RHSVec(2);

for k = 2:length(RHSVec)-1
    Vec(k) = MainDiag(k) * RHSVec(k) + OffDiagConst * (RHSVec(k-1) + RHSVec(k+1));
end

Vec(length(RHSVec)) = MainDiag(length(RHSVec)) * RHSVec(length(RHSVec)) + OffDiagConst * RHSVec(length(RHSVec)-1);

end