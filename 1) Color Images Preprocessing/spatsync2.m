function [ssimg] = spatsync2(img,xss,yss)
% ========================================================================
% ------------------------Spatial Synchronisation -------------------------
% ---------------------------- 28-OCT-2013 ------------------------------ %
% Function to spatially synchronise images that have been affected by
% camera movement.
% ------------------------------------------------------------------------
% find most middle xss and yss and adjust everything to that
% +ve means image is up and left from middle
% -ve means image is down and right from middle
sizeimg = size(img);
for i = 1:sizeimg(3)
    xoffset(i) = (round(median(xss,3)) - round(xss(:,:,i)));
    yoffset(i) = (round(median(yss,3)) - round(yss(:,:,i)));
end
% calculate size of new matrix plus zeros for offsets
xmax = max(xoffset)+abs(min(xoffset));
ymax = max(yoffset)+abs(min(yoffset));
% most negative location is always at (1,1) when added to xmid and
% ymid
xmid = abs(min(xoffset))+1;
ymid = abs(min(yoffset))+1;

% create new domain that can contain every image with respect to the median
% image
ssimg = zeros(sizeimg(1,1)+ymax,...
    sizeimg(1,2)+xmax,...
    sizeimg(1,3));
% y = rows, x = cols
% if all of the offsets are equal to zero, then mid values will also be
% zero
for i = 1:sizeimg(3)
    if xmax == 0 && ymax == 0
        ssimg(1:sizeimg(1),...
            1:sizeimg(2),...
            i) = img(:,:,i);
    elseif xmax == 0 && ymax ~= 0
        ssimg(ymid+yoffset(i):ymid+yoffset(i)+(sizeimg(1,1)-1),...
            1:sizeimg(2),...
            i) = img(:,:,i);
    elseif xmax ~=0 && ymax == 0
        ssimg(1:sizeimg(1),...
            xmid+xoffset(i):xmid+xoffset(i)+(sizeimg(1,2)-1),...
            i) = img(:,:,i);
    else
        ssimg(ymid+yoffset(i):ymid+yoffset(i)+(sizeimg(1,1)-1),...
            xmid+xoffset(i):xmid+xoffset(i)+(sizeimg(1,2)-1),...
            i) = img(:,:,i);
    end
end