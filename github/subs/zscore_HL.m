function zy = zscore_HL(y)

zy = (y - mean(y,1,'omitnan')) ./ std(y,'omitnan');

end