function [per] = PercChange(data, p)

 data = reshape(data, p, [])';
 
 baseline = mean(data, 1);
 
 per = (data ./ baseline - 1) * 100;
 
end