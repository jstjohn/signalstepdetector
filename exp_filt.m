function filt = exp_filt(data, start_ind, end_ind, alpha)

filt = zeros(size(data(start_ind:end_ind)));
filt(1) = data(start_ind);
ind=1;
for i=start_ind+1:end_ind
    %FROM NOAH 
    %calculate the exponentially weighted moving average
     %mean = (1-alpha).*data(i) + alpha*mean;
    filt(ind+1) = data(i)*(1-alpha) + alpha*filt(ind);
    ind = ind+1;
end