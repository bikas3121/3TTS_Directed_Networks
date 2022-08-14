nCk = [3, 4, 6, 7];
n = sum(nCk);

Adj = 10*sprandn(n,n,0.5);
Adj = full(Adj);
Adj = abs(Adj);

Adj_0 = ones(n,n);

for i = 1:length(nCk)
   if i == 1
        st = 1;
        ed = nCk(1);
        for j = st:ed
            Adj_0(st:ed,j)=0;
        end
   else
       st = sum(nCk(1:i-1))+1;
       ed = st + nCk(i)-1;
       for j = st:ed
            Adj_0(st:ed,j)=0;
       end
   end


end

Adj = floor(Adj.*Adj_0);

Deg = diag(Adj*ones(n,1));
L = Deg- Adj;

