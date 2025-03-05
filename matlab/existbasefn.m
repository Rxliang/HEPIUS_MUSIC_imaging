function [existing, val] =  existbasefn(vari)

existing = evalin('base', ['exist(''' vari ''')']);

if existing
  val =  evalin('base', vari);
else
  val = [];
end

