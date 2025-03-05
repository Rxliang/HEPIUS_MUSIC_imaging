function sm = sumall(mat, arg)

if nargin < 2
  arg = 'includemissing';
end
    
sm = sum(mat(:), arg);