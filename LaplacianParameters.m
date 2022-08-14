function [PosLapEig, LambdaMin, ZeroIndex, LapNorm] = LaplacianParameters(Lint, nCk)

% nCk = [3, 4, 6, 7];
% 
% Lint1 = weighted_directed_Lap_ran(nCk(1));
% Lint2 = weighted_directed_Lap_ran(nCk(2));
% Lint3 = weighted_directed_Lap_ran(nCk(3));
% Lint4 = weighted_directed_Lap_ran(nCk(4));
% 
% Lint = blkdiag(Lint1, Lint2, Lint3, Lint4);

LambdaMin = zeros(1,length(nCk));
PosLapEig = [];
ZeroIndex = zeros(1,length(nCk));
LapNorm = zeros(1, length(nCk));
for i = 1:length(nCk)
    if i == 1      
        Lint_i = Lint(1:nCk(i), 1:nCk(i));
        D = real(eig(Lint_i));
        for j = 1:length(D)
            if D(j) <= 1e-7
                ZeroIndex(i) = j;
            end
        end
        LapNorm(i) = norm(Lint_i,2);
        D(ZeroIndex(i)) = [];
        LambdaMin(i) = max(D);
    else
        
        st = sum(nCk(1:i-1))+1;
        ed = st+nCk(i)-1;
        Lint_i = Lint(st:ed, st:ed);
        D = real(eig(Lint_i));
        for j = 1:length(D)
            if D(j) <= 1e-7
                ZeroIndex(i) = j;
            end
        end
        LapNorm(i) = norm(Lint_i,2);
        D(ZeroIndex(i)) = [];  
        LambdaMin(i) = max(D);
    end
    PosLapEig = [PosLapEig; D];
end
LambdaMin = min(LambdaMin);
end