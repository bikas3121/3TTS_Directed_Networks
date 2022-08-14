function [x_Comp1, x_Comp2] =  StateComponents(x)

l = size(x,2);
x_Comp1 = [];  %component 1 of the x
x_Comp2 = [];   % component 2 of the x
for i=  2:2:l
   x_Comp1 = [x_Comp1 x(:,i-1)];
   x_Comp2 = [x_Comp2 x(:,i)];
end

end