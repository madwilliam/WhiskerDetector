% this function should be rewritten to check for all lines that are close to the whisker base 1 2 3 ....
function [whisker, AllW_Trial]=addJaneliaAndDlcResults(whisker,framei,whisker_boundaries,Ron_Leftoff,X2,Y2,TW,TM,probability,badframes,plottype,numberofwhiskers,minimumwhiskerlenght,arclenththreshold)
% numberofwhiskers = 3;
% get trial information if DLC data is accurate and frame was considered
% good Theta is for DLC data with complete whisker input and is not
%%
minprob=0.01;
mapxy= [X2' Y2'];
if framei>1
    if Ron_Leftoff
        if size(whisker.R.BaseDLCX,1)>=framei-1
            
            mapxyp= [whisker.R.BaseDLCX(framei-1,:); whisker.R.BaseDLCY(framei-1,:)]';
        else
            mapxyp= [mapxy];
        end
    else
        if size(whisker.L.BaseDLCX,1)>=framei-1
            
            mapxyp= [whisker.L.BaseDLCX(framei-1,:) ;whisker.L.BaseDLCY(framei-1,:)]';
        else
            mapxyp= [mapxy];
        end
    end
else
    mapxyp= [mapxy];
    
end


if sum(mapxyp)==0
    mapxyp= [mapxy];
    
end
thistrial1 = (TW.time ==framei-1&ismember(TW.id, 0:1e5)); % find  whisker (from ganelia tracking) data in that given frame
thistrial=[];
if any(probability>minprob)&&badframes(framei)==0&&any(thistrial1)&&any(any(whisker_boundaries))
    % if ifplot
    %     im2 = readindex(ThisVideo,double(framei));
    %     if plottype
    %         imagesc(im2,[10 250])
    % 
    %     else
    %         imagesc((uint8(BW)+0.5).*im2,[10 250])
    % 
    %     end
    %     hold on
    % end
    
    [in,on] = inpolygon(TM.follicle_x,TM.follicle_y,whisker_boundaries(:,2),whisker_boundaries(:,1)); %take  whiskers base (foliccle) that are inside the
    
    % box that is defined by DLC data
    thistrial = find( (in|on)&thistrial1); %combing point that are in and on the polygon and
    arclen = zeros(size((thistrial)));
    if any(thistrial) & any(probability>minprob) & badframes(framei)==0 & numel(thistrial)>=1
        basepoints2 =NaN(numel(thistrial),2);
        
        for iii=1:numel(thistrial)
            x=(TW.x{thistrial(iii)}) ;
            y=(TW.y{thistrial(iii)}) ;
            [arclen(iii),~] = arclength(x,y);
            basepoints2(iii,:) = [TM.follicle_x(thistrial(iii)),TM.follicle_y(thistrial(iii))];
        end
        %         remove whisker lenghts smaller than 40 pixle..
        thistrial(arclen<arclenththreshold)=[];
        basepoints2(arclen<arclenththreshold,:)=[];
        [~,sortid]= sort(basepoints2(:,2));
        thistrial=thistrial(sortid);
        arclen= arclen(sortid);
        Thetameanall=[];
        WhiskpointsinDLC=[];
        for iii=1:numel(thistrial)
            % plot all whiksers inside the
            if ifplot&&plottype
                plot (TW.x{thistrial(iii)},TW.y{thistrial(iii)},'color',[0.8 0.4 0.4]);
            end
            % get whisker base infor
            x=(TW.x{thistrial(iii)}) ;
            y=(TW.y{thistrial(iii)}) ;
            % get whiskerlent
            All_whiskers_Lenght (iii,:) = norm([x(1) y(1)]'-[x(end) y(end)]');
            AllW_Trial.X {iii} =  x;
            AllW_Trial.Y {iii} =  y;
            AllW_Trial.Lenght {iii} =  All_whiskers_Lenght;
            AllW_Trial.Folicle {iii} =   [TM.follicle_x(thistrial(iii)),TM.follicle_y(thistrial(iii))];
            AllW_Trial.Measurement{iii} =   TM(thistrial(iii),:);
            AllW_Trial.arclenght{iii} =  arclen(iii);
            if All_whiskers_Lenght (iii,:) < minimumwhiskerlenght
                meandistances(iii,1:numel(X2)) = 1e5;
                meandistancesp(iii,1:numel(X2)) = 1e5;
                distance_to_the_Base(iii,1:numel(X2))= 1e5;
                Thetameanall(iii,:)=1e5;
                continue
            end
            [~,Thetamean]=findSplinefitAngel( x, y);
            Thetameanall(iii,:) = mean(Thetamean);
            [in2,~] = inpolygon( x,y,whisker_boundaries(:,2),whisker_boundaries(:,1));
            WhiskpointsinDLC(iii)= sum(in2);
            POINT1 = [TM.follicle_x(thistrial(iii)),TM.follicle_y(thistrial(iii))];
            basepoints = [TM.follicle_x(thistrial(iii)),TM.follicle_y(thistrial(iii))];
            %     get the distance between the dlc base and whisker curves
            for dlc_i = 1:numel(X2)
                POINT2 = [X2(dlc_i),Y2(dlc_i)];
                distance_to_the_Base(iii,dlc_i) = pdist([POINT1;POINT2]);
                curvexy= [x y];
                if framei>1&&size(whisker.L.Y,2)>=numel(X2)&&any(whisker.L.Y{framei-1,dlc_i})
                    follicleN_1 =[whisker.L.X{framei-1,dlc_i}(end),whisker.L.Y{framei-1,dlc_i}(end)];
                    [~,meandistancesCurve_N_1(iii,dlc_i),~] = distance2curve(curvexy,follicleN_1,'linear');
                else
                    meandistancesCurve_N_1(iii,dlc_i) =1e10;
                end
                [~,meandistances(iii,dlc_i),~] = distance2curve(curvexy,mapxy(dlc_i,:),'linear');
                [~,meandistancesp(iii,dlc_i),~] = distance2curve(curvexy,mapxyp(dlc_i,:),'linear');
            end
        end
        meandistances2=meandistances;
        meandistances2(isnan(meandistances2))=1000;
        alli=[];
        counter=0;
        MinAdaptiveLenght = [30 20 10 10];
        MinAdaptiveLenght = [median(arclen) median(arclen)*0.8 20 20];
                MinAdaptiveLenght = [30 20 10 10];
        colors=jet(3);
        while any(meandistances2)
            counter=counter+1;
            Whikserdistances = meandistances2(:,1);
            Whikserdistances(arclen<MinAdaptiveLenght(counter)) = 1e10;
            [value,minindex] = min(Whikserdistances);
            if framei>1&&size(whisker.L.Y,2)>=numel(X2)
                [value2,minindex2] = min(meandistances2);
            else
                value2=100;
            end
            meandistances2(:,1)=[];
            meandistances2(1:minindex,:)=1e10;
            if value>1e5
                continue
            end
            alli=([alli minindex]);
        end
        iii=unique(alli);
        if any(iii)
            for thisiii=1:numel(iii)
                if Ron_Leftoff
                    whisker.R.Y{framei,thisiii} = TW.y{thistrial(iii(thisiii))};
                    whisker.R.X{framei,thisiii} = TW.x{thistrial(iii(thisiii))};
                    whisker.R.FoX(framei,thisiii) = TM.follicle_x(thistrial(iii(thisiii)));
                    whisker.R.FoY(framei,thisiii) = TM.follicle_y(thistrial(iii(thisiii)));
                    whisker.R.Lenght(framei,thisiii)=TM.length(thistrial(iii(thisiii)));
                    if ifplot&&plottype
                        plot( whisker.R.X{framei,thisiii},whisker.R.Y{framei,thisiii},'color',colors(thisiii,:))
                    end
                    [Theta,Thetamean]=findSplinefitAngel(  whisker.R.X{framei,thisiii},  whisker.R.Y{framei,thisiii});
                    whisker.R.A(framei,thisiii)=deg2rad(Theta);
                    whisker.R.A2(framei,thisiii)=deg2rad(mean(Thetamean));
                    whisker.R.JN(framei,thisiii)=true;
                else
                    whisker.L.Y{framei,thisiii} = TW.y{thistrial(iii(thisiii))};
                    whisker.L.X{framei,thisiii} = TW.x{thistrial(iii(thisiii))};
                    whisker.L.FoX(framei,thisiii) = TM.follicle_x(thistrial(iii(thisiii)));
                    whisker.L.FoY(framei,thisiii) = TM.follicle_y(thistrial(iii(thisiii)));
                    if ifplot&&plottype
                        plot( whisker.L.X{framei,thisiii},whisker.L.Y{framei,thisiii},'color',colors(thisiii,:))
                    end
                    [Theta,Thetamean]=findSplinefitAngel(  whisker.L.X{framei,thisiii},  whisker.L.Y{framei,thisiii});
                    whisker.L.A(framei,thisiii)=deg2rad(Theta);
                    whisker.L.A2(framei,thisiii)=deg2rad(mean(Thetamean));
                    whisker.L.JN(framei,thisiii)=true;
                    whisker.L.Lenght(framei,thisiii)=TM.length(thistrial(iii(thisiii)));
                end
            end
        end
    else
        if Ron_Leftoff
            whisker.R.A(framei,numberofwhiskers)=NaN;
            whisker.R.A2(framei,numberofwhiskers)=NaN;
        else
            whisker.L.A(framei,numberofwhiskers)=NaN;
            whisker.L.A2(framei,numberofwhiskers)=NaN;
        end
            AllW_Trial.X{numberofwhiskers}  =  NaN;
            AllW_Trial.Y{numberofwhiskers}  =  NaN;
            AllW_Trial.X {numberofwhiskers} =  NaN;
            AllW_Trial.Y {numberofwhiskers} =  NaN;
            AllW_Trial.Lenght {numberofwhiskers} =  NaN;
            AllW_Trial.Folicle {numberofwhiskers} =   NaN;
            AllW_Trial.Measurement{numberofwhiskers} =  NaN;
            AllW_Trial.arclenght{numberofwhiskers} =  NaN;
    end
end
if ~(any(probability>minprob)&&badframes(framei)==0&&any(thistrial1)&&any(any(whisker_boundaries)))
    if Ron_Leftoff
        whisker.R.A(framei,numberofwhiskers)=NaN;
        whisker.R.A2(framei,numberofwhiskers)=NaN;
    else
        whisker.L.A(framei,numberofwhiskers)=NaN;
        whisker.L.A2(framei,numberofwhiskers)=NaN;
    end
        AllW_Trial.X {numberofwhiskers} =  NaN;
        AllW_Trial.Y {numberofwhiskers} =  NaN;
        AllW_Trial.Lenght {numberofwhiskers} =  NaN;
        AllW_Trial.Folicle {numberofwhiskers} =   NaN;
        AllW_Trial.Measurement{numberofwhiskers} =  NaN;
        AllW_Trial.arclenght{numberofwhiskers} =  NaN;
end

if ~isempty(thistrial)
    if numel(thistrial)<2
        if Ron_Leftoff
            whisker.R.A(framei,numberofwhiskers)=NaN;
            whisker.R.A2(framei,numberofwhiskers)=NaN;
        else
            whisker.L.A(framei,numberofwhiskers)=NaN;
            whisker.L.A2(framei,numberofwhiskers)=NaN;
        end
            AllW_Trial.X {numberofwhiskers} =  NaN;
            AllW_Trial.Y {numberofwhiskers} =  NaN;
            AllW_Trial.Lenght {numberofwhiskers} =  NaN;
            AllW_Trial.Folicle {numberofwhiskers} =   NaN;
            AllW_Trial.Measurement{numberofwhiskers} =  NaN;
            AllW_Trial.arclenght{numberofwhiskers} =  NaN;
    end
end