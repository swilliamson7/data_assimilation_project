% this function builds various matrices we need. It needs 
%       k - spring coefficient
%       r - dissipation constant
% and will return
%       Ac - continuous model operator
%       Kc - submatrix of Ac (description?)
%       Rc - submatrix of Ac
function [Ac, Kc, Rc] = build_matrices(k, r)
    
   I3 = eye(3);
   Rc = -r .* eye(3);
   Kc = [-2*k, k, 0;
         k, -3*k, k;
         0, k, -2*k];
    
   Ac = [zeros(3) I3
         Kc       Rc];
end