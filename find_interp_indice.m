function[indice_]=find_interp_indice(x,xp)
indice_=zeros(length(x),2);
% for x that is ramping up, interpolate from first to last element of xp
if x(end)>x(1)
    start = 1;
    stride = 1;
    end_= length(x);
% for x that is ramping down, interpolate from last to first element of xp
else
    start=length(x);
    stride = -1;
    end_= 1;
end
% if the first x value is unchanged
if xp(start)<=x(start)
    indice_(start,1) = start+stride;
    indice_(start,2) = start;
    start = start+stride;
end
% do interpolation
for ii = start:stride:end_
    xt = xp(ii);
    while xt>x(start)
        start=start+stride;
    end
    indice_(ii,1) = start-stride;
    indice_(ii,2) = start;
%     x0=x(start-stride);
%     x1=x(start);
%     y0=y(start-stride);
%     y1=y(start);
%     yp(ii)=y0+(xt-x0)*(y1-y0)/(x1-x0);
end
% toc