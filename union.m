function C = union(A,B)
% UNION Union of two sets of positive integers (much faster than built-in union)
% C = union(A,B)

if isempty(A)
  ma = 0;
else
  ma = max(A);
end

if isempty(B)
  mb = 0;
else
  mb = max(B);
end

if ma==0 && mb==0
  C = [];
elseif ma==0 && mb>0
  C = B;
elseif ma>0 && mb==0
  C = A;
else
  %bits = sparse(1, max(ma,mb));
  bits = zeros(1, max(ma,mb));
  bits(A) = 1;
  bits(B) = 1;
  C = find(bits);
end
