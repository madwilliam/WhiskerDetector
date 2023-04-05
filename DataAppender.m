classdef DataAppender<handle
    properties
        head_measurement_right
        head_measurement_left
        whisker_detection_right
        whisker_detection_left
    end
    methods
        function load_janelia_detection_results(self,vidname)
            measurementsR = LoadMeasurements(fullfile(self.data_folder,['Mirror' vidname 'R.measurements']));
            measurementsL = LoadMeasurements(fullfile(self.data_folder, ['Mask' vidname 'L.measurements']));
            self.head_measurement_right = struct2table(measurementsR);
            self.head_measurement_left = struct2table(measurementsL);
            self.whisker_detection_right = struct2table(LoadWhiskers(fullfile(self.data_folder, ['Mirror' vidname 'R.whiskers'])));
            self.whisker_detection_left = struct2table(LoadWhiskers(fullfile(self.data_folder, ['Mask' vidname 'L.whiskers'])));
        end

        function previous_whisker_coordinates = get_previous_whisker_coordinates(self,framei,whisker_coordinates)
            if framei>1
                if self.Ron_Leftoff
                    if size(self.Whisker.R.BaseDLCX,1)>=framei-1
                        previous_whisker_coordinates= [self.Whisker.R.BaseDLCX(framei-1,:); self.Whisker.R.BaseDLCY(framei-1,:)]';
                    else
                        previous_whisker_coordinates= whisker_coordinates;
                    end
                else
                    if size(self.Whisker.L.BaseDLCX,1)>=framei-1
                        
                        previous_whisker_coordinates= [self.Whisker.L.BaseDLCX(framei-1,:) ;self.Whisker.L.BaseDLCY(framei-1,:)]';
                    else
                        previous_whisker_coordinates= whisker_coordinates;
                    end
                end
            else
                previous_whisker_coordinates= whisker_coordinates;
            end
            if sum(previous_whisker_coordinates)==0
                previous_whisker_coordinates= whisker_coordinates;
            end
        end

        function [X2,Y2,TM,TW] = get_data(self,side,framei)
            X2 = self.Whisker.R.BaseDLCX(framei,:);
            Y2 = self.Whisker.R.BaseDLCY(framei,:);
            switch side
                case 'right'
                    TW = self.whisker_detection_right;
                    TM = self.head_measurement_right;
                case 'left'
                    TW = self.whisker_detection_left;
                    TM = self.head_measurement_left;
            end
        end

        function AllW_Trial = append_data(self,framei,side,whisker_boundaries,probability)
            [X2,Y2,TM,TW] = self.get_data(side,framei);
            minprob=0.01;
            whisker_coordinates= [X2' Y2'];
            previous_whisker_coordinates = self.get_previous_whisker_coordinates(framei,whisker_coordinates);
            thistrial1 = (TW.time ==framei-1&ismember(TW.id, 0:1e5)); % find  self.Whisker (from ganelia tracking) data in that given frame
            thistrial=[];
            disp('')
            if any(probability>minprob)&&any(thistrial1)&&any(any(whisker_boundaries))
                [in,on] = inpolygon(TM.follicle_x,TM.follicle_y,whisker_boundaries(:,2),whisker_boundaries(:,1)); %take  whiskers base (foliccle) that are inside the
                % box that is defined by DLC data
                thistrial = find( (in|on)&thistrial1); %combing point that are in and on the polygon and
                arclen = zeros(size((thistrial)));
                if any(thistrial) && any(probability>minprob) && numel(thistrial)>=1
                    basepoints2 =NaN(numel(thistrial),2);
                    for iii=1:numel(thistrial)
                        x=(TW.x{thistrial(iii)}) ;
                        y=(TW.y{thistrial(iii)}) ;
                        [arclen(iii),~] = arclength(x,y);
                        basepoints2(iii,:) = [TM.follicle_x(thistrial(iii)),TM.follicle_y(thistrial(iii))];
                    end
                    %         remove self.Whisker lenghts smaller than 40 pixle..
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
                        % get self.Whisker base infor
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
                        %     get the distance between the dlc base and self.Whisker curves
                        for dlc_i = 1:numel(X2)
                            POINT2 = [X2(dlc_i),Y2(dlc_i)];
                            distance_to_the_Base(iii,dlc_i) = pdist([POINT1;POINT2]);
                            curvexy= [x y];
                            if framei>1&&size(self.Whisker.L.Y,2)>=numel(X2)&&any(self.Whisker.L.Y{framei-1,dlc_i})
                                follicleN_1 =[self.Whisker.L.X{framei-1,dlc_i}(end),self.Whisker.L.Y{framei-1,dlc_i}(end)];
                                [~,meandistancesCurve_N_1(iii,dlc_i),~] = distance2curve(curvexy,follicleN_1,'linear');
                            else
                                meandistancesCurve_N_1(iii,dlc_i) =1e10;
                            end
                            [~,meandistances(iii,dlc_i),~] = distance2curve(curvexy,whisker_coordinates(dlc_i,:),'linear');
                            [~,meandistancesp(iii,dlc_i),~] = distance2curve(curvexy,previous_whisker_coordinates(dlc_i,:),'linear');
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
                        if framei>1&&size(self.Whisker.L.Y,2)>=numel(X2)
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
                                self.Whisker.R.Y{framei,thisiii} = TW.y{thistrial(iii(thisiii))};
                                self.Whisker.R.X{framei,thisiii} = TW.x{thistrial(iii(thisiii))};
                                self.Whisker.R.FoX(framei,thisiii) = TM.follicle_x(thistrial(iii(thisiii)));
                                self.Whisker.R.FoY(framei,thisiii) = TM.follicle_y(thistrial(iii(thisiii)));
                                self.Whisker.R.Lenght(framei,thisiii)=TM.length(thistrial(iii(thisiii)));
                                if ifplot&&plottype
                                    plot( self.Whisker.R.X{framei,thisiii},self.Whisker.R.Y{framei,thisiii},'color',colors(thisiii,:))
                                end
                                [Theta,Thetamean]=findSplinefitAngel(  self.Whisker.R.X{framei,thisiii},  self.Whisker.R.Y{framei,thisiii});
                                self.Whisker.R.A(framei,thisiii)=deg2rad(Theta);
                                self.Whisker.R.A2(framei,thisiii)=deg2rad(mean(Thetamean));
                                self.Whisker.R.JN(framei,thisiii)=true;
                            else
                                self.Whisker.L.Y{framei,thisiii} = TW.y{thistrial(iii(thisiii))};
                                self.Whisker.L.X{framei,thisiii} = TW.x{thistrial(iii(thisiii))};
                                self.Whisker.L.FoX(framei,thisiii) = TM.follicle_x(thistrial(iii(thisiii)));
                                self.Whisker.L.FoY(framei,thisiii) = TM.follicle_y(thistrial(iii(thisiii)));
                                if ifplot&&plottype
                                    plot( self.Whisker.L.X{framei,thisiii},self.Whisker.L.Y{framei,thisiii},'color',colors(thisiii,:))
                                end
                                [Theta,Thetamean]=findSplinefitAngel(  self.Whisker.L.X{framei,thisiii},  self.Whisker.L.Y{framei,thisiii});
                                self.Whisker.L.A(framei,thisiii)=deg2rad(Theta);
                                self.Whisker.L.A2(framei,thisiii)=deg2rad(mean(Thetamean));
                                self.Whisker.L.JN(framei,thisiii)=true;
                                self.Whisker.L.Lenght(framei,thisiii)=TM.length(thistrial(iii(thisiii)));
                            end
                        end
                    end
                else
                    if self.Ron_Leftoff
                        self.Whisker.R.A(framei,self.Numberofwhiskers)=NaN;
                        self.Whisker.R.A2(framei,self.Numberofwhiskers)=NaN;
                    else
                        self.Whisker.L.A(framei,self.Numberofwhiskers)=NaN;
                        self.Whisker.L.A2(framei,self.Numberofwhiskers)=NaN;
                    end
                        AllW_Trial.X{self.Numberofwhiskers}  =  NaN;
                        AllW_Trial.Y{self.Numberofwhiskers}  =  NaN;
                        AllW_Trial.X {self.Numberofwhiskers} =  NaN;
                        AllW_Trial.Y {self.Numberofwhiskers} =  NaN;
                        AllW_Trial.Lenght {self.Numberofwhiskers} =  NaN;
                        AllW_Trial.Folicle {self.Numberofwhiskers} =   NaN;
                        AllW_Trial.Measurement{self.Numberofwhiskers} =  NaN;
                        AllW_Trial.arclenght{self.Numberofwhiskers} =  NaN;
                end
            end
            if ~(any(probability>minprob)&&any(thistrial1)&&any(any(whisker_boundaries)))
                if Ron_Leftoff
                    self.Whisker.R.A(framei,self.Numberofwhiskers)=NaN;
                    self.Whisker.R.A2(framei,self.Numberofwhiskers)=NaN;
                else
                    self.Whisker.L.A(framei,self.Numberofwhiskers)=NaN;
                    self.Whisker.L.A2(framei,self.Numberofwhiskers)=NaN;
                end
                AllW_Trial.X {self.Numberofwhiskers} =  NaN;
                AllW_Trial.Y {self.Numberofwhiskers} =  NaN;
                AllW_Trial.Lenght {self.Numberofwhiskers} =  NaN;
                AllW_Trial.Folicle {self.Numberofwhiskers} =   NaN;
                AllW_Trial.Measurement{self.Numberofwhiskers} =  NaN;
                AllW_Trial.arclenght{self.Numberofwhiskers} =  NaN;
            end
            
            if ~isempty(thistrial)
                if numel(thistrial)<2
                    if Ron_Leftoff
                        self.Whisker.R.A(framei,self.Numberofwhiskers)=NaN;
                        self.Whisker.R.A2(framei,self.Numberofwhiskers)=NaN;
                    else
                        self.Whisker.L.A(framei,self.Numberofwhiskers)=NaN;
                        self.Whisker.L.A2(framei,self.Numberofwhiskers)=NaN;
                    end
                        AllW_Trial.X {self.Numberofwhiskers} =  NaN;
                        AllW_Trial.Y {self.Numberofwhiskers} =  NaN;
                        AllW_Trial.Lenght {self.Numberofwhiskers} =  NaN;
                        AllW_Trial.Folicle {self.Numberofwhiskers} =   NaN;
                        AllW_Trial.Measurement{self.Numberofwhiskers} =  NaN;
                        AllW_Trial.arclenght{self.Numberofwhiskers} =  NaN;
                end
            end
        end
        function plot(self)
            im2 = readindex(ThisVideo,double(framei));
            if plottype
                imagesc(im2,[10 250])

            else
                imagesc((uint8(BW)+0.5).*im2,[10 250])

            end
            hold on
        end
    end
end