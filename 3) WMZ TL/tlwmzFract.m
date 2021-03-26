function [tl,wmz,npoints,c25,c50,c75] = tlwmzFract(xmat,ymat,zmat,cdata,edge,offset,tl_predicted)
% =========================================================================
% ---------------- Toe Length and Width of Mixing Zone-7 ------------------
% ---------------------------- 26-SEP-2014 ------------------------------ %
% Function to calculate the toe length and width of mixing zone of a
% saltwater wedge.
% --- Input Variables ---
% xmat     - 2D matrix of x co-ordinates in real/image space
% ymat     - 2D matrix of y co-ordinates in real/image space
% zmat     - 2D matrix of concentration fields
% cdata    - concentration contour values (25, 50 and 75% of max conc)
% edge     - remove edge effects when calculating toe length
%           value of 1 - yes (experimental)
%           value of 0 = no (numerical)
% offset   - offset in x and y direction from sync point when removing edges
% ao       - analysis offset (usually 8mm [5mm + 3mm], only numerical)
% tl_predicted - predict TL at bottom boundary through extrapolation
% --- Output Variables ---
% tl       - toe length
% wmz      - width of mixing zone
% npoints  - no. of points used in determining wmz
% c25      - 25% concentration contour co-ordinates
% c50      - 50% concentration contour co-ordinates
% c75      - 75% concentration contour co-ordinates
%ONLY Change lines 224-5 lolim=0.2*tl(2,1) instead of lolim=0.2*tl
% =========================================================================
[C1,h1] = contour(xmat, ymat, zmat,cdata);
axis equal
axis tight
% Determine starting locations of contour lines in C1
% C1(1,1) = concentration of first contour, ie. cdata(1)
% C2(2,1) = no. of data points in contour
% C1(1,:) = x co-ordinates
% C1(2,:) = y co-ordinates
% But some have multiple contours, where the longest one needs to be taken
% as the main contour line
sizeC1 = size(C1);
if sizeC1(1,2) == 0 ||... % ie. no contours exist, or;
        sum(C1(1,:)== cdata(1,1)) == 0 ||... % no 25% contour exists, or;
        sum(C1(1,:)== cdata(1,2)) == 0 ||... % no 50% contour exists, or;
        sum(C1(1,:)== cdata(1,3)) == 0; % no 75% contour exists.
    
    tl(1:2,1) = 0;
    wmz = 0;
    npoints = 0;
    c25(1:2,1) = 0;
    c50(1:2,1) = 0;
    c75(1:2,1) = 0;
else
    sizezmat = size(zmat);
    i=1; % counter for the number of contour lines in C1
    ndatacontours(1,1) = C1(1,1); % First contour concentration
    ndatacontours(2,1) = C1(2,1); % No. of points in first contour
    ndatacontours(3,1) = 1; % Column location of contour in C1
    while sum(ndatacontours(2,:))+i < sizeC1(1,2) % TRY USING WHILE
        % while the sum of all the contour points, plus the number of contour
        % lines is < the length of C1 (ie. exists in C1):
        i = i + 1;
        col_counter = sum(ndatacontours(2,:))+i;
        ndatacontours(1,i) = C1(1,col_counter);
        ndatacontours(2,i) = C1(2,col_counter);
        ndatacontours(3,i) = col_counter;
        
    end
   
    % ndatacontours contains the concentration and number of points in each
    % contour. Find the contour with the largest number of points for each
    % concentration. Also, at least one of their y co-ordinates must be 0, ie.
    % the bottom of the domain.
    for i = 1:length(ndatacontours)
        minspanx = min(C1(1,ndatacontours(3,i)+1:...
            ndatacontours(3,i)+ndatacontours(2,i)));
        maxspanx = max(C1(1,ndatacontours(3,i)+1:...
            ndatacontours(3,i)+ndatacontours(2,i)));
        minspany = min(C1(2,ndatacontours(3,i)+1:...
            ndatacontours(3,i)+ndatacontours(2,i)));
        maxspany = max(C1(2,ndatacontours(3,i)+1:...
            ndatacontours(3,i)+ndatacontours(2,i)));
        
        
 %% Modified part of tlwamzNumerical.m
            spanx(i,1) = maxspanx - minspanx; % span of contour in x direction
            spanx(i,2) = ndatacontours(1,i);% conc. related to span

    end
    if sum(spanx(:,2)==cdata(1)) >= 1 &&...
            sum(spanx(:,2)==cdata(2)) >= 1 &&...
            sum(spanx(:,2)==cdata(3)) >= 1 % if sum of suitable 25%,50% and 75% 
        % contour lines > 1(ie. exist) the run tl and wmz analysis
    for i = 1:length(cdata)
        conc_cols = find(spanx(:,2) == cdata(i));
        [max_conc_loc(1,i)] = find(spanx == max(spanx(conc_cols,1)));
        max_conc_loc(2,i) = ndatacontours(2,max_conc_loc(1,i));
    end
    % if no y co-ordinates exist in the lowest
    % =========================================================================
    % --- Toe Length ---
    % Calculate TL based on max_conc_loc(1,2) for 50 % concentration line, ie.
    % the last point on the longest contour line
    % (the x co-ordinate when y co-ordinate is minimum)
    c50pre = C1(:,ndatacontours(3,max_conc_loc(1,2))+1:...
        ndatacontours(3,max_conc_loc(1,2))+max_conc_loc(2,2));
    % isolate maximum length contour lines for 25% and 75%
    c25pre = C1(:,ndatacontours(3,max_conc_loc(1,1))+1:...
        ndatacontours(3,max_conc_loc(1,1))+max_conc_loc(2,1));
    c75pre = C1(:,ndatacontours(3,max_conc_loc(1,3))+1:...
        ndatacontours(3,max_conc_loc(1,3))+max_conc_loc(2,3));
   
    % Edge effects from the experiment can significantly affect the
    % determination of a toe length.
    
    if edge == 0
        
        
        c25(1:4,length(c25pre)) = zeros;
        c25(1:2,:) = c25pre;
        loc25y = (max(find(c25pre(2,:) < ao)))+1; % 1st location of a value greater than 0.006m from bottom boundary
        c25(3:4,loc25y:end) = c25pre(:,loc25y:end);
        
        c50(1:4,1:length(c50pre)) = zeros;
        c50(1:2,:) = c50pre;
        loc50y = (max(find(c50pre(2,:) < ao)))+1;
        c50(3:4,loc50y:end) = c50pre(:,loc50y:end);
        
        c75(1:4,1:length(c75pre)) = zeros;
        c75(1:2,:) = c75pre;
        loc75y = (max(find(c75pre(2,:) < ao)))+1;
        c75(3:4,loc75y:end) = c75pre(:,loc75y:end);
        % isolate maximum length contour lines for 25% and 75%
        tl(1,1) = max(c50(1,:));
        tl(2,1) = max(c50(3,:));
        
    elseif edge == 1
        % carry out linear regression on a portion of c50 and extend
        % to the bottom boundary. Use this intersection as tl.
        % find where c50 is at its max and take the previous 20% of co-ordinates to
        % calculate the tl
        if tl_predicted == 1
            for i = 1:length(cdata)
                % isolate contour data, same script as above for c50,c25 & c75
                c_conc = C1(:,ndatacontours(3,max_conc_loc(1,i))+1:...
                    ndatacontours(3,max_conc_loc(1,i))+max_conc_loc(2,i));
                size(c_conc)
                                
                % determine maximum values of x and y
                [max_valx(i,1),max_locx(i,1)] = max(c_conc(1,:));
                [max_valy(i,1),max_locy(i,1)] = max(c_conc(2,:));
                % need x co-ordinate at max_locy to trim the contour line and v.v.
                min_locx(i,1) = c_conc(1,max_locy(i,1));
                min_locy(i,1) = c_conc(2,max_locx(i,1));
                
                
                % for max_locx values
                x0 = [1;0]; % Starting guess
                options = optimset('Display','off','MaxFunEvals',5000,'MaxIter',5000);
                [x,resnorm,residual,exitflag] = lsqcurvefit(@linear1,x0,...
                    c_conc(1,max_locx(i,1):(max_locx(i,1)+round(0.2*length(c_conc)))),...
                    c_conc(2,max_locx(i,1):(max_locx(i,1)+round(0.2*length(c_conc)))),...
                    [],[],options);
                coeffsx(i,:) = x;
                
                % for max_locy values
                x0 = [1;0]; % Starting guess
                options = optimset('Display','off','MaxFunEvals',5000,'MaxIter',5000);
                [x,resnorm,residual,exitflag] = lsqcurvefit(@linear1,x0,...
                    c_conc(1,(max_locy(i,1)+round(0.2*length(c_conc))):max_locy(i,1)),...
                    c_conc(2,(max_locy(i,1)+round(0.2*length(c_conc))):max_locy(i,1)),...
                    [],[],options);
                coeffsy(i,:) = x;
                
                if exitflag < 1
                    error('Regression failed.');
                end
            end
        end
        
        
        % trim contours to these values
        c25(1:4,length(c25pre)) = zeros;
        c25(1:2,:) = c25pre; % raw contour line
        % contour line trimmed to location of maximum values
        % Any reduction in contour line co-ordinate value after the maximum
        % is assumed to be due to edge effects.
        
                
        c50(1:4,1:length(c50pre)) = zeros;
        c50(1:2,:) = c50pre;
                
        c75(1:4,1:length(c75pre)) = zeros;
        c75(1:2,:) = c75pre;
        if tl_predicted == 1            
            tl(1,1) = abs(((-1*offset(2,1))-coeffsx(2,2))/coeffsx(2,1))+offset(1,1);
            
            c25(3:4,max_locx(1,1):max_locy(1,1)) = c25pre(:,max_locx(1,1):max_locy(1,1));
            c50(3:4,max_locx(2,1):max_locy(2,1)) = c50pre(:,max_locx(2,1):max_locy(2,1));
            c75(3:4,max_locx(3,1):max_locy(3,1)) = c75pre(:,max_locx(3,1):max_locy(3,1));
        else
            tl(1,1) = NaN;
        end
        % min value of y occurs
        [closesty] = find(c50(2,:) == min(c50(2,:)));
        % closesty could have multiple readings. The minimum value of closesty will
        % be were the wedge has intruded furthest.
        
        tl(2,1) = c50(1,min(closesty));
    end
    % =========================================================================
    % --- Width of Mixing Zone ---
    % Determining width of mixing zone 'wmz'
    % First find limits 0.2*tl and 0.8*tl
    lolim =0.2*tl(2,1);
    uplim =0.8*tl(2,1);
    
    % Find the contour limits for calculating wmz using uplim and lolim
    % --- 25% ---
    P25 = 1; % Set counter = 1  
    
    while c25(1,P25) > lolim; % contours points in C1 start with y co-ords at min
        P25 = P25+1; % Increase counter by 1
    end
    
    P25(2) = 1; % Set counter = 1
    
    while c25(1,P25(2)) > uplim; % contours points in C1 start with y co-ords at min
        P25(2) = P25(2)+1; % Increase counter by 1
    end
    
    
    % --- 75% ---
    P75 = 1; % Set counter = 1
    while c75(1,P75) > lolim; % contours points in C1 start with y co-ords at min
        P75 = P75+1; % Increase counter by 1
    end
    
    P75(2) = 1; % Set counter = 1
    while c75(1,P75(2)) > uplim; % contours points in C1 start with y co-ords at min
        P75(2) = P75(2)+1; % Increase counter by 1
    end
    
    % Define 25% and 75% contours within the limits of 0.2 and 0.8 TL
    c25lim = c25(:,P25(2):P25(1));
    c75lim = c75(:,P75(2):P75(1));
    
    % Find location of matching x co-ordinates in c25lim and c75lim
    icounter = size(c25lim);
    
    for i = 1:icounter(1,2);
        
        % Find the column numbers where CC75 equals CC25 for all values of CC25
        xf = find(c25lim(1,i)==c75lim(1,:));
        
        % If xf is empty (ie. does not match) replace with a zero
        xfempty = isempty(xf);
        
        if xfempty == 1;
            xf2(1,i) = 0;
        else
            xf2(1,i) = xf(1,1) ;
        end
        
    end
    
    % if any co-ordinates match, numbers in xf2 will be greater than zero
    if max(xf2) > 0
        icounter2 = size(xf2); % Same size as CC25
        
        n = 1;
        
        for i = 1:icounter2(1,2);
            
            % Where value doesn't equal 0, place value in xf3 (column values of
            % CC75 matches)
            if xf2(1,i) ~= 0;
                xf3(1,n) = xf2(1,i);
                n = n+1;
            end
            
        end
        
        % Find column values of non zero's in xf3 to determine column numbers of
        % CC25 values, xf4
        xf4 = find(xf2~=0);
        icounter3 = size(xf3);
        icounter4 = size(xf4);
        
        for i = 1:icounter3(1,2);
            xf75(1,i) = c75lim(2,xf3(1,i));
        end
        
        for i = 1:icounter4(1,2);
            xf25(1,i) = c25lim(2,xf4(1,i));
        end
        npoints = length(xf25);
        % mean of the difference between 25% and 75% isolines between 'tl' limits
        wmz = mean(xf25-xf75);
        
    else
        
        wmz = nan;
        npoints = nan;
    end
    else % contour lines unsuitable
        tl(1:2,1) = 0;
        wmz = 0;
        npoints = 0;
        c25(1:2,1) = 0;
        c50(1:2,1) = 0;
        c75(1:2,1) = 0;
end
end
end