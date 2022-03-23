function v = nan_interp(v)
% interpolate nan
if isnan(v(1))
    v(1) = v(find(~isnan(v), 1, 'first'));
end
if isnan(v(end))
    v(end) = v(find(~isnan(v), 1, 'last'));
end
nanx = isnan(v);
t = 1:numel(v);
v(nanx) = interp1(t(~nanx), v(~nanx), t(nanx));
