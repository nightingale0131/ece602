function vertices = poly2D(A,b)
% Return list of vertices from inequalities A*x <= b

nbr = size(b,1);
list_x = [];
for i = 1:nbr
    for j = i+1:nbr
      
       mat = [A(i,:);A(j,:)];
       v = [b(i);b(j)];
       x = mat\v;
       
       if max(isinf(x)) == 0
           if min(A*x-b <= 0.01) == 1
                list_x = [list_x, x];
           end
       end
    end
    
end

vertices = list_x;
end
