%Script: Snowgreen
%Date: 25/2/2022, version 1.0
%Author: Andreas Westergaard-Nielsen, awn@ign.ku.dk
%Purpose: Analyze greenness and snow cover from RGB images

%% input variables
% Close windows and clear workspace
clear
% Path to input images
path = '/path to input images';
pathOut = '/path to output';
% Path to mask file
maskfile = '/mask.mat'; %path to predefined binary image mask
% Choose graph and data output
showgraphs = true;
exportdata = false;
% Name of exported Excel file
name = 'outputImages.xlsx';
% Snow threshold if dynamic threshold fails
fixedThreshold = 240; %variance threshold. higher = more variance allowed in snow
year = 2021; % define year
nameOut = ['output.mat'];
nameLength = 8;

%% processing
files = dir(fullfile(path, ['**' filesep '*.JPG']));
metafiles = dir(fullfile(path, ['**' filesep '*.meta']));
startDoy = 333; %set to first doy in timeseries
counter = 0;

for i = 1:length(files)
    exif = imfinfo([files(i).folder filesep files(i).name]); %grab exif information from image
    inputDate = files(i).name;
    if nameLength == 7
        inputDate = datenum(str2num(inputDate(15:18)),str2num(inputDate(20:21)),str2num(inputDate(23:24)));
    elseif nameLength == 8
        inputDate = datenum(str2num(inputDate(16:19)),str2num(inputDate(21:22)),str2num(inputDate(24:25)));
    end
    t1 = datevec(inputDate);
    imgDoy = date2doy(inputDate); %derive day of year (doy) from the image
    doyNoDecimals(i) = fix(imgDoy);
    dayList(i) = datenum(t1(1),t1(2),t1(3));
    
end

dayShift = unique(dayList); %create list of unique dates

for i = 1:length(dayShift)
    
    idx = dayList==dayShift(i); %list of days to average images from
    filesShortlist = files(idx);
    metaShortlist = metafiles(idx);
    
    for ii = 1:length(filesShortlist)
        x2 = imread([filesShortlist(ii).folder filesep filesShortlist(ii).name]); %load input image as "x
        disp(['Loading ' filesShortlist(ii).folder filesep filesShortlist(ii).name]);
        xRed(:,:,ii) = x2(:,:,1);
        xGre(:,:,ii) = x2(:,:,2);
        xBlu(:,:,ii) = x2(:,:,3);
        metaData = importMeta([metaShortlist(ii).folder filesep metaShortlist(ii).name]);
        metaData = convertStringsToChars(metaData);
        findMe1 = strfind(metaData,'Camera Temperature: ');
        findMe2 = strfind(metaData,'\r\n');
        metaDataList(ii) = str2num(metaData((findMe1+20):(findMe2(2)-1)));
        
    end
    metaDataNum(i) = mean(metaDataList);
    save([pathOut 'metaData.mat'],'metaData');
    
    xRed = single(mean(xRed,3));
    xGre = single(mean(xGre,3));
    xBlu = single(mean(xBlu,3));
end

newFiles = dir([pathOut filesep '*.jpg']);

for i = 1:length(newFiles)
    i
    disp(newFiles(i).name);
    try
    if nameLength == 7
        yearOfAvgImg = newFiles(i).name(15:18)
        monthOfAvgImg = newFiles(i).name(20:21);
        dayOfAvgImg = newFiles(i).name(23:24);
    elseif nameLength == 8
        yearOfAvgImg = newFiles(i).name(16:19)
        monthOfAvgImg = newFiles(i).name(21:22);
        dayOfAvgImg = newFiles(i).name(24:25);
    end
    catch
        pause
    end
    imgTimeAxis(i) = datenum(str2num(yearOfAvgImg),str2num(monthOfAvgImg),str2num(dayOfAvgImg));
    doyNoDecimals(i) = fix(imgDoy);
    x = imread([newFiles(i).folder filesep newFiles(i).name]);
    
    load_string = strcat('load', [' ',maskfile]);
    eval(load_string);
    roi_mask=polygon_mask;
    
    % Mask image (i.e., set areas outside mask to zero)
    roi_mask(:,:,2) = roi_mask;
    roi_mask(:,:,3) = roi_mask(:,:,1);
    roi_image = x;
    roi_image(roi_mask == 0) = 0;
    
    % Analyze local image range
    redRan = roi_image(:,:,1);
    redRan(redRan==0)=nan;
    greenRan = roi_image(:,:,2);
    greenRan(greenRan==0)=nan;
    blueRan = roi_image(:,:,3);
    blueRan(blueRan==0)=nan;
    img4Range(:,:,1) = redRan;
    img4Range(:,:,2) = greenRan;
    img4Range(:,:,3) = blueRan;
    roi_rangefilt = rangefilt(img4Range);
    roi_range(i) = nanmean(roi_rangefilt(:));
    
    if rangefiltAmpl > 50
        if roi_range(i) < 19 | metaDataNum < 32
            fixedThreshold2 = 600
        elseif metaDataNum(i) > 32
            fixedThreshold2 = 80
        elseif metaDataNum(i) < 32 | roi_range(i) > 19
            fixedThreshold2 = 240
        end
        
    else
        if roi_range(i) < 19 | metaDataNum < 32
            fixedThreshold2 = 270
        elseif metaDataNum(i) > 32
            fixedThreshold2 = 70
        elseif metaDataNum(i) < 32 | roi_range(i) > 19
            fixedThreshold2 = 120
        end
    end
    
    % Detect snow (ones = snow true)
    varianceImg = var(single(roi_image),0,3);
    variance(i) = mean2(varianceImg);
    stretchedImg = imadjust(roi_image,stretchlim(roi_image),[]);
    meanImg = mean(single(stretchedImg),3);
    meanImgThreshold = zeros(size(meanImg));
    meanImgThreshold(meanImg > 256/4) = 1; %decide if mean pixel dn is above the lowest 1/4 to avoid dark pixels being recognized as snow
    [counts,x] = histcounts(varianceImg(:),256);
    x(1) = []; x(end) = [];%remove first edge
    counts(1) = []; %remove zeros generated by maskin
    counts = smooth(counts); % 5 point moving average
    
    snowMask = zeros(size(roi_image));
    varImgThreshold = single(varianceImg < fixedThreshold2);
    snowMask = varImgThreshold + meanImgThreshold;
    snowMask(snowMask ~= 2) = 0;
    snowMask(snowMask == 2) = 1;
    
    %-------------------------------------------------------------------------%
    % end of snow classificaion
    snowMask=double(snowMask);snowMask(polygon_mask==0)=0;           %preprocessing to get snow percentage
    polygon_mask=double(polygon_mask);polygon_mask(polygon_mask==0)=NaN;%preprocessing to get snow percentage
    snow(i) = nansum(nansum(snowMask)) / nansum(nansum(polygon_mask)); %generate vector with percent snow
    snow(i)
    % Mask snow from ROI
    roi_image=double(roi_image);
    roi_image_snow(:,:,1) = roi_image(:,:,1).*~snowMask;
    roi_image_snow(:,:,2) = roi_image(:,:,2).*~snowMask;
    roi_image_snow(:,:,3) = roi_image(:,:,3).*~snowMask;
    
    % Compute greenness
    R = roi_image_snow(:,:,1);
    Rall = R;
    R = R(roi_image_snow(:,:,1)~=0);
    R = mean(R(:));
    G = roi_image_snow(:,:,2);
    Gall = G;
    G = G(roi_image_snow(:,:,2)~=0);
    G = mean(G(:));
    B = roi_image_snow(:,:,3);
    Ball = B;
    B = B(roi_image_snow(:,:,3)~=0);
    B = mean(B(:));
    RGB = [R G B];
    
    GCC = (Gall./(Rall+Gall+Ball));
    GCC(isnan(GCC))=[];
    GCCstd(i) = std(GCC);
    R_mean(i) = R;
    G_mean(i) = G;
    B_mean(i) = B;
    
    % Greenness index
    green_2grbi(i) = 2*(RGB(2)) - (RGB(1)+RGB(3)); %2G_RBi
    green_gcc(i) = (RGB(2)/(RGB(1)+RGB(2)+RGB(3))); %GCC
    green_ndgi(i) = (2*(RGB(2)) - (RGB(1)+RGB(3)))/(2*(RGB(2)) - (RGB(1)-RGB(3))); %NDGI
    green(i) = (RGB(2)/(RGB(1)+RGB(2)+RGB(3))); %GCC for plot output    
    roi_image_snow = uint8(roi_image_snow);
    imshow(roi_image_snow); %show masked snow pixels
    clearvars roi_image_snow snow_mask roi_image roi_image_gray
end

clc
green(green==0)=nan;
snow(snow==0)=nan;
snowInterp = snow;
snowInterp(snowInterp >1)=1;
snowInterp(snowInterp <0)=0;
smoothGreen = smooth(green,7)';
snowInterp(end-10:end)=nan;
snowSmooth = smooth(snow);
entropySmooth = smooth(roi_range);
maxEntropy = mean(maxk(entropySmooth,4));
minEntropy = mean(mink(entropySmooth,4));
rangeEntropy = maxEntropy - minEntropy;

% show figures
if showgraphs
    
    figure(1)
    subplot(1,2,1)
    plot(imgTimeAxis,100*snow(:,:),'kx'),title('Amount of snow in the ROI'),grid on,...
        xlabel('DOY'),ylabel('Snow cover [%]');
    
    subplot(1,2,2)
    plot(imgTimeAxis,conv(green,ones(1,1)/1,'same'),'-kx'),hold on;
    plot(imgTimeAxis,smoothGreen,'-c.'),hold on;
    title('Amount of green in the ROI'),...
        legend('Original','Averaged 7'),grid on,hold off,...
        xlabel('DOY'),ylabel('Green Chromatic Coordinate');
    
end

if exportdata
    save([pathOut filesep nameOut],'files','dayShift','green','imgTimeAxis','roi_range','snow','snowSmooth','metaDataNum','maxEntropy','minEntropy','rangeEntropy');
    
end


function doy = date2doy(inputDate)
% Modified by Andreas Westergaard-Nielsen from:
% Author: Anthony Kendall
% Contact: anthony [dot] kendall [at] gmail [dot] com
% Created: 2008-03-11
% Copyright 2008 Michigan State University.
doy = deal(zeros(size(inputDate)));
inputDate = inputDate(:);
[dateVector] = datevec(inputDate);
dateVector(:,2:end) = 0;
dateYearBegin = datenum(dateVector);
doyRow = inputDate - dateYearBegin;
doy(:) = doyRow;
end

function t1 = importMeta(metaInputFile)
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [57, 57];
opts.Delimiter = "=";

% Specify column names and types
opts.VariableNames = ["Var1", "VarName2"];
opts.SelectedVariableNames = "VarName2";
opts.VariableTypes = ["string", "string"];
opts = setvaropts(opts, [1, 2], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

%Import the data
t1 = readtable(metaInputFile, opts);

% Convert to output type
t1 = table2array(t1);

% Clear temporary variables
clear opts
end

function t2 = importExposure(metaInputFile)
% Specify column names and types
opts = delimitedTextImportOptions("NumVariables", 2);
opts.Delimiter = '='
opts.DataLines = [1, Inf];
opts.VariableNames = ["agc", "VarName2"];
opts.VariableTypes = ["string", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "agc", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "agc", "EmptyFieldRule", "auto");

% Import the data
t2 = readtable(metaInputFile, opts);
t2 = table2array(t2);
clear opts
end