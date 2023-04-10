%% this m file is using one single DLC file to detect whisker base of both Right whiskers and Left whiskser (in Mirror R format)
data_folder = '\\dk-server.dk.ucsd.edu\data\afassihizakeri\topview data\2023_02_22_ 163923';
tracker = WhiskerTracker(data_folder);
tracker.main()
wait_bar = waitbar(0,'Tracking in progress...');
for  vidname =trial_names
    trial = str2num(vidname);
    %% using commnads from Janelia whisker tracking algorithm 
    THISx = ['trace ' data_folder 'Mask' vidname 'L.avi ' data_folder 'Mask' vidname 'L.whiskers'];
    THISy = ['trace ' data_folder   'Mirror' vidname 'R.avi ' data_folder 'Mirror' vidname 'R.whiskers'];
    THISx2 = ['       measure --face right '  data_folder  'Mask' vidname 'L.whiskers '  data_folder  'Mask' vidname 'L.measurements               '];
    THISY2 = ['       measure --face right '  data_folder  'Mirror' vidname 'R.whiskers '  data_folder  'Mirror' vidname 'R.measurements               '];
    status = system(THISx);
    status2 = system(THISy);
    status3 = system(THISx2);
    status4 = system(THISY2);
    measurementsR = LoadMeasurements([data_folder 'Mirror' vidname 'R.measurements']);
    measurementsL = LoadMeasurements([data_folder 'Mask' vidname 'L.measurements']);
    TMR = struct2table(measurementsR);
    TML = struct2table(measurementsL);
    TWR = struct2table(LoadWhiskers([data_folder 'Mirror' vidname 'R.whiskers']));
    TWL = struct2table(LoadWhiskers([data_folder 'Mask' vidname 'L.whiskers']));
    %% loading dlc files
    textfilenameR = [data_folder  'Mask1LDLC_resnet152_ar30shiwkerSep13shuffle2_100000'];
    textfilenameL = [data_folder  'Mirror1RDLC_resnet152_ar30shiwkerSep13shuffle2_100000'];
    optsR = detectImportOptions(textfilenameR);
    optsR.VariableNames={'frames', 'x1','y1','L1'...
        ,'x2','y2','L2','x3','y3','L3','x4','y4','L4','x5','y5','L5'};
    optsL = detectImportOptions(textfilenameL);
    optsL.VariableNames={'frames', 'x1','y1','L1'...
        ,'x2','y2','L2','x3','y3','L3','x4','y4','L4','x5','y5','L5'};
    %open tables with DLC whisker files
    TLeft = readtable(textfilenameL,optsL);
    TRight = readtable(textfilenameR,optsR);        
    [TLeft,TRight]=Smooth_whiskerDLCTLeft(TLeft,TRight,Numberofwhiskers);
    %%
    % open videos
    if trial ==1
        vL = VideoReader([data_folder  'Mask' vidname 'L.avi']);
        vR = VideoReader([data_folder  'Mirror' vidname 'R.avi']);
        im2R = readindex(vR,double(100));  im2R = im2R(:,:,1);
        ThissizeR(1) = size(im2R,1);
        ThissizeR(2) = size(im2R,2);
        im2 = readindex(vL,double(100));  im2 = im2(:,:,1);
        ThissizeL(1) = size(im2,1);
        ThissizeL(2) = size(im2,2);
    end
    %%
    badframesR =zeros(1,size(TLeft,1));    
    badframesL =zeros(1,size(TLeft,1));
    frame_data = readtable([data_folder  vidname 'FrameData.xlsx']);
    Frames = frame_data.goodframes-1;
    pix_1 = [TLeft.x1,TLeft.y1];
    pix_2 =  [TLeft.x2,TLeft.y2];
    distance_left = sqrt( (pix_2(:,1)-pix_1(:,1)).^2 + (pix_2(:,2)-pix_1(:,2)).^2 );
    pix_1 = [TRight.x1,TRight.y1];
    pix_2 =  [TRight.x2,TRight.y2];
    distance_right = sqrt( (pix_2(:,1)-pix_1(:,1)).^2 + (pix_2(:,2)-pix_1(:,2)).^2 );
    theseframes =ismember(frame_data.goodframes-1, Frames(TRight.frames+1));
    nose_pixel = [frame_data.Nosex(theseframes),frame_data.Nosey(theseframes)];
    snout_pixel = [frame_data.Snoutx(theseframes),frame_data.Snouty(theseframes)];
    nose_to_snout_distance = sqrt( (snout_pixel(:,1)-nose_pixel(:,1)).^2 + (snout_pixel(:,2)-nose_pixel(:,2)).^2 );
    if numel(distance_left)~=numel(distance_right)
        continue
    end
    check = (smooth(abs(zscore(distance_left))+abs(zscore(distance_right))+abs(zscore(nose_to_snout_distance)),40));
    thisdistance = NaN(size(frame_data.Nosex));
    for i=1:numel(frame_data.Nosex)
        thisdistance(i) = pdist([frame_data.Nosex(i),frame_data.Snoutx(i);frame_data.Nosey(i),frame_data.Snouty(i)]', 'euclidean');
    end
    Whisker.R.A=NaN(size(frame_data.Nosex,1),2);
    Whisker.R.A2=NaN(size(frame_data.Nosex,1),2);
    Whisker.L.A=NaN(size(frame_data.Nosex,1),2);
    Whisker.L.A2=NaN(size(frame_data.Nosex,1),2);
    Whisker.R.X=cell(size(frame_data.Nosex,1),2);
    Whisker.R.Y=cell(size(frame_data.Nosex,1),2);
    Whisker.L.X=cell(size(frame_data.Nosex,1),2);
    Whisker.L.Y=cell(size(frame_data.Nosex,1),2);
    Whisker.R.JN=NaN(size(frame_data.Nosex,1),2);
    Whisker.R.A3=NaN(size(frame_data.Nosex,1),2);
    Whisker.L.JN=NaN(size(frame_data.Nosex,1),2);
    Whisker.L.A3=NaN(size(frame_data.Nosex,1),2);
    plottype = 1;
    L=TRight.L2;
    clear F
    Whisker.R.BaseDLCX=NaN(size(frame_data.Nosex,1),Numberofwhiskers);
    Whisker.R.BaseDLCY=NaN(size(frame_data.Nosex,1),Numberofwhiskers);
    Whisker.L.BaseDLCX=NaN(size(frame_data.Nosex,1),Numberofwhiskers);
    Whisker.L.BaseDLCY=NaN(size(frame_data.Nosex,1),Numberofwhiskers);
    Whisker.R.FoX=NaN(size(frame_data.Nosex,1),Numberofwhiskers);
    Whisker.L.FoX=NaN(size(frame_data.Nosex,1),Numberofwhiskers);
    Whisker.R.FoY=NaN(size(frame_data.Nosex,1),Numberofwhiskers);
    Whisker.L.FoY=NaN(size(frame_data.Nosex,1),Numberofwhiskers);
    for i=1:numel(Frames)
        if ~any(TRight.frames==i)
            continue
        end
        if ~any(TLeft.frames==i)
            continue
        end
        if check(i)>10
            continue
        end
        arclenththreshold=10;
        XYbiasDLC(1)=1;
        XYbiasDLC(2)=1.5;
        % find the mask and whisker values from deep lab cut
        [probability,~,~,X2,Y2,boundariesR ,BW ]= findthemask_WhiskerBaseonly_2 (ThissizeR,i,TRight,GuassFilterSize,Numberofwhiskers,XYbiasDLC);
        Ron_Leftoff=1;
        [Whisker,AllwTrialR] = addwhiskertrialfromDLCandJanelia_WhiskerBaseOnly_Sep2021(Whisker,i,boundariesR,Ron_Leftoff,ifplot,vR,BW,X2,Y2,TWR,TMR,probability,badframesR,plottype,Numberofwhiskers,minimumwhiskerlenght,arclenththreshold);
        if isempty(AllwTrialR.X{1})
            continue
        end
        if ifplot
            axis image 
            hold off
            colormap gray
            set(gca,'visible','off')
            cmap2  =jet(size(X2,1));
            subplot(2,4,6);
            im2 = readindex(vR,double(i));
            imagesc((uint8(BW)+0.5).*im2,[10 250])
            hold on
            scatter(X2',Y2', 10, cmap2, 'filled')
            set(gca,'visible','off')
            axis image
            hold off
            colormap gray
            subplot(2,4,8);
        end
        if i>1&&size(X2,1)<2
            X2=  Whisker.R.BaseDLCX(end,:);
            Y2= Whisker.R.BaseDLCY(end,:);
        end
        Whisker.R.BaseDLCX(i,:)=X2;
        Whisker.R.BaseDLCY(i,:)=Y2;
        Ron_Leftoff=0;
         XYbiasDLC(1)=1;
      XYbiasDLC(2)=3;
        [probability,~,~,X2,Y2,boundariesL ,BW ]= findthemask_WhiskerBaseonly_2 (ThissizeL,i,TLeft,GuassFilterSize,Numberofwhiskers,XYbiasDLC);
        [Whisker,AllwTrialL] = addwhiskertrialfromDLCandJanelia_WhiskerBaseOnly_Sep2021(Whisker,i,boundariesL,Ron_Leftoff,ifplot,vL,BW,X2,Y2,TWL,TML,probability,badframesL,plottype,Numberofwhiskers,minimumwhiskerlenght,arclenththreshold);
        if isempty(AllwTrialL.X{1})
            continue
        end
        if i>1&&size(X2,1)<2
            X2= Whisker.L.BaseDLCX(end,:);
            Y2= Whisker.L.BaseDLCY(end,:);
        end
        Whisker.L.BaseDLCX(i,:)=X2;
        Whisker.L.BaseDLCY(i,:)=Y2;
        if ifplot
            cmap2  = jet(size(X2,1));
            colormap gray
            set(gca,'visible','off')
            hold off
            axis image
            subplot(2,4,7);
            im2 = readindex(vL,double(i));
            imagesc((uint8(BW)+0.5).*im2,[10 250])
            hold on
            scatter(X2',Y2', 10, cmap2, 'filled')
            axis image
            set(gca,'visible','off')
            set(gca,'XDir','reverse')
            axis image
            hold off
            F(i) = getframe(gcf) ;
            drawnow
            subplot(2,4,1:4)
            hold off
            plot(1,1,'k.')
            for xx=5:8
                subplot(2,4,xx)
                hold off
                plot(1,1,'k.')
            end
        end
    end
    allwhisker.R = AllwTrialR;
    allwhisker.L= AllwTrialL;  
    save([vidname ,'.mat'],allwhisker)
    clear allwhisker AllwTrialL AllwTrialR  
    Whisker.DLC =frame_data;
    Whisker.Name = vidname;
    AllW{trial} = Whisker; %#ok<SAGROW>
    
    if ifplot
        figure('units','normalized','outerposition',[0 0 1 1]) 
        v1 = VideoReader([data_folder  vidname '.avi']);
        vL = VideoReader([data_folder  'Mask' vidname 'L.avi']);
        vR = VideoReader([data_folder  'Mirror' vidname 'R.avi']);
        subplot(2,4,1:4); 
        if Frames(i)>0
        Frind = find(frame_data.goodframes-1==Frames(i));
        im2 = readindex(v1,double(Frames(i)));
        imagesc(im2,[10 250])
        text(20,400,num2str(i),'color',[1 1 1])
        set(gca,'YDir','normal')
        set(gca,'visible','off')
        hold on
        plot(frame_data.Nosex(Frind),frame_data.Nosey(Frind),'.r','markersize',20)
        plot(frame_data.Snoutx(Frind),frame_data.Snouty(Frind),'.g','markersize',20)
        text(100,100,num2str( Frames(i)))
        hold off
        end
        subplot(2,4,5);
        writerObj = VideoWriter([vidname 'Tracked.avi']);
        writerObj.FrameRate = 10;
        open(writerObj);
        for iF=1:length(F)
            if any(F(iF).cdata)
                frame = F(iF) ;
                writeVideo(writerObj, frame);
            end
        end
        close(writerObj);
        close all
    end
    toc
end
save AllW AllW
Whisker = AllW{7};
hold off
plot(Whisker.L.A3)
hold on
theseframes = ~cellfun('isempty',Whisker.R.X);
plot(find(theseframes),Whisker.R.A2(theseframes),'r-.')
thetrialnumber=[];
right_movies=dir('*.mp4');
for trial =1:numel(right_movies)
    Index=[];
    vidname = right_movies(trial).name;
    k = strfind(vidname, 'video');
    vidname = vidname(1:k-1);
    thetrialnumber =[thetrialnumber str2num(vidname)];
end
save thetrialnumber thetrialnumber
for i=1:200
    close()
    F = findall(0,'type','figure','tag','TMWWaitbar');
    delete(F)
end
