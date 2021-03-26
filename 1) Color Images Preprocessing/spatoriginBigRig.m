function [tl,tr,bl,br] = spatoriginBigRig(img,pixelsizem,dh)
% =========================================================================
% --------------- Current variation 10-DEC-2020 ---------------------------
% ---------------------- Spatial Origin of Images -------------------------
% ---------------------------- 11-DEC-2013 ------------------------------ %
% Function to determine the origin of images in space, to be representative
% of the origin of numerical modelling output. Is also used in determining
% the offsets required for spatial synchronisation of images.
% Definitions tl=top left, tr= top right, bl=bottom left, br=bottom right
% and they are point dimensions!!!
% --- Ammendments ---
% 19-AUG-2014
% 1. Added in region of interest variables to easily re-define boundarys if
% necessary (eg. for other experimental setups). NOTE: CHANGE WITHIN
% FUNCTION
% --- Variables to set ---
% Lower right region of the image
roi_rows = [0.75 0.89]; % region of interest to analyse, eg. [0.25 0.65] is
%                       equivalent to between 40% and 60% of the image rows
roi_cols = [0.9 1]; %This always a percentage
% ----------------------------------------------------------------------- %
sizeimg = size(img);
% -------------------------------------------------------------------------
% isolate roi of image (RHS) and avg along cols.
avg_cols_br = mean(img((round(sizeimg(1,1)*roi_rows(1,1)):...
    round(sizeimg(1,1)*roi_rows(1,2))),...
    (round(sizeimg(1,2)*roi_cols(1,1))):...
    round(sizeimg(1,2)*roi_cols(1,2))));
% find the abs difference between adjacent profiles
for i = 1:length(avg_cols_br)-1
    diff_brx(1,i) = abs(avg_cols_br(1,i) - avg_cols_br(1,i+1));
end
[maxval_br] = max(diff_brx);
[maxloc_br] = min(find(diff_brx == maxval_br));
avg_rows_br = mean(img((round(sizeimg(1,1)*roi_rows(1,1)):...
    round(sizeimg(1,1)*roi_rows(1,2))),...
    (round(sizeimg(1,2)*roi_cols(1,1))):...
    round(sizeimg(1,2))*roi_cols(1,2)),2);
% avg_rows_br = transpose(avg_rows_br_t);
for i = 1:length(avg_rows_br)-1
    diff_bry(i,1) = abs(avg_rows_br(i,1) - avg_rows_br(i+1,1));
end
[maxval_br(2,1)] = max(diff_bry);
[maxloc_br(2,1)] = min(find(diff_bry == maxval_br(2,1)));

br(1,1) = maxloc_br(1,1)+(roi_cols(1,1)*(sizeimg(1,2)));
br(2,1) = maxloc_br(2,1)+(roi_rows(1,1)*sizeimg(1,1));
% -------------------------------------------------------------------------
% isolate first 20% of the image (LHS) & middle half and avg along cols.
avg_cols_bl = std(img((round(sizeimg(1,1)*roi_rows(1,1)):...
    round(sizeimg(1,1)*roi_rows(1,2))),...
    1:round(sizeimg(1,2)*(1-roi_rows(1,2)))));
% find the difference between adjacent profiles
for i = 1:length(avg_cols_bl)-1
    diff_blx(1,i) = abs(avg_cols_bl(1,i) - avg_cols_bl(1,i+1));
end
[maxval_bl] = max(diff_blx);
[maxloc_bl] = find(diff_blx == maxval_bl);

bl(1,1) = maxloc_bl(1,1);
bl(2,1) = maxloc_br(2,1)+(roi_rows(1,1)*sizeimg(1,1));

tr(1,1) = br(1,1);
tr(2,1) = round(br(2,1)-(dh/pixelsizem));
tl(1,1) = bl(1,1);
tl(2,1) = tr(2,1);
end