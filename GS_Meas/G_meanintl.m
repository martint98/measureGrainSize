function G = G_meanintl(u)
% Calculate the ASTM grain size as a function of mean intercept length

% mean intercept length was previously calculated:
    %---- calculating the average intercept length
%         ints = xync(:,5)  --extracting number of intercept lines
%         z = sum(ints)  --total number of intercept lines
% 
%         intl = linints(:,1)  --extracting the intercept lengths
%         q = sum(intl)  --total of intercept lengths
% 
%         % calculating the mean intercept length:
%         u = q / z  --total of int lengths / total number of intercepts

A = 2.0 * log2(320.0);
B = 2.0 / log10(2.0);

G = A - B * log10(u);

end

