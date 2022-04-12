function [new_data] = bin_data(new_time,old_time,data);


for kk = 1:length(new_time)-1 %go through each minute interval
            mi = find(old_time >= new_time(kk) & old_time < new_time(kk+1));
            if isempty(mi) %if no data in the current minute, write as nan
                new_data(kk) = nan;
                
            elseif ~isempty(mi)
                %new_data(kk) = nanmean(data(mi));
                new_data(kk) = nanmedian(data(mi));
                
            end
end
new_data(end+1)=nan;

end