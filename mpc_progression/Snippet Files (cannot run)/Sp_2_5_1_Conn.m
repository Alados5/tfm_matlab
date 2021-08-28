conn = zeros(nr*nc);
for i=1:nr*nc
    if(mod(i,nc) ~= 1), conn(i, i-1) = 1; end
    if(mod(i,nc) ~= 0), conn(i, i+1) = 1; end
    if(i > nc),         conn(i, i-nc) = 1; end
    if(i <= (nr-1)*nc), conn(i, i+nc) = 1; end
end
conn(nctrl,:) = 0; % Remove connections on controlled nodes