close all
clear all
clc

%% Insert filepath and parameters
% File location
FileLocation = %Insert path 
% input parameters
channel = [1 2];
ns=[1 2];
sz=[11 5];
pkh=[12 10]; 
channel_track = [1 2];
param.mem=1;
param.good=10;
param.dim=2;
param.quiet=1;



oldFolder = cd(FileLocation);
listing = dir('**/*.tif');
cd(oldFolder)

%% Detect and Track

for i = 1 : size(listing,1)
    fileName = FileLocation + '\' + listing(i).name;
    Img = MinStack.loadStack(convertStringsToChars(fileName));
    [filepath,name,ext] = fileparts(listing(i).name);

    Img = Img.bpass(channel,ns,sz);
    Img = Img.pkfnd(channel,pkh,sz);
    Img = Img.cntrd(channel,sz+2);
    Img = Img.track(channel_track,5,param);
    Img.drawpoint_track(channel,channel_track,name+"_track")
    Img.saveData(name+"_data");
    close all
end
