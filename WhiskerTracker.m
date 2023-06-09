classdef WhiskerTracker <DataAppender&WhiskerDetector
    properties
        data_folder
        minimumwhiskerlenght
        thrs
        textfilenameR
        textfilenameL
        Minmaxdist
        plottype
        arclenththreshold
        Ron_Leftoff
    end
    methods
        function self = WhiskerTracker(data_folder)
            self.data_folder = data_folder;
            self.Numberofwhiskers=2;
            self.minimumwhiskerlenght = 30;
            self.thrs= 2;
            self.Minmaxdist=[15 200];
            self.GuassFilterSize = 20;
            self.arclenththreshold=10;
            self.plottype = 1;
            self.Ron_Leftoff=1;
            self.XYbiasDLC = [1,1.5];
            [fileL,fileR] = self.select_left_and_right_whisker_dlc();
            self.textfilenameR = fullfile(data_folder,fileL);
            self.textfilenameL = fullfile(data_folder,fileR);
        end

        function [fileL,fileR] = select_left_and_right_whisker_dlc(self)
            trial_name = self.get_trial_names();
            trial_number = cellfun(@(str) str2double(str) , trial_name);
            [trial_number,id] = sort(trial_number);
            trial_name = trial_name(id);
            [indx,~] = listdlg('ListString',trial_name);
            wanted_trial = trial_number(indx);
            fileL = self.get_dlc_file(wanted_trial,'Mirror*filtered.csv','RDLC');
            fileR = self.get_dlc_file(wanted_trial,'Mask*filtered.csv','LDLC');
        end

        function fileL = get_dlc_file(self,wanted_trial,csv_pattern,dlc_pattern)
            fileL = dir(fullfile(self.data_folder,csv_pattern));
            trial_number_left = cellfun(@(str) extract(str,digitsPattern+dlc_pattern) ,{fileL.name});
            trial_number_left = cellfun(@(str) str2double(str(1:end-numel(dlc_pattern))) ,trial_number_left);
            left_id = find(trial_number_left==wanted_trial);
            fileL = fileL(left_id).name;
        end

        function trial_names = get_trial_names(self)
            files = {dir(self.data_folder).name};
            pat =  digitsPattern+".avi";
            trial_names = cellfun(@(file) extract(file,pat), files,'UniformOutput',false);
            is_trial = cellfun(@(name) numel(name),trial_names)>0;
            trial_names = trial_names(is_trial);
            trial_names = cellfun(@(name) name(1),trial_names);
            trial_names = cellfun(@(file) {extractBefore(file,'.avi')}, trial_names);
        end

        function status = run_janelia_whisker_tracking(self,vidname)
            THISx = ['trace "' self.data_folder '\Mask' vidname 'L.avi ' self.data_folder '\Mask' vidname 'L.whiskers"'];
            THISy = ['trace "' self.data_folder   '\Mirror' vidname 'R.avi ' self.data_folder '\Mirror' vidname 'R.whiskers"'];
            THISx2 = ['       measure --face right "'  self.data_folder  '\Mask' vidname 'L.whiskers" "'  self.data_folder  '\Mask' vidname 'L.measurements"               '];
            THISY2 = ['       measure --face right "'  self.data_folder  '\Mirror' vidname 'R.whiskers" "'  self.data_folder  '\Mirror' vidname 'R.measurements"               '];
            status = system(THISx);
            status2 = system(THISy);
            status3 = system(THISx2);
            status4 = system(THISY2);
            status = [status,status2,status3,status4];
        end

        function [TLeft,TRight] = load_dlc_files(self)
            optsR = detectImportOptions(self.textfilenameR);
            optsR.VariableNames={'frames', 'x1','y1','L1'...
                ,'x2','y2','L2','x3','y3','L3','x4','y4','L4','x5','y5','L5'};
            optsL = detectImportOptions(self.textfilenameL);
            optsL.VariableNames={'frames', 'x1','y1','L1'...
                ,'x2','y2','L2','x3','y3','L3','x4','y4','L4','x5','y5','L5'};
            TLeft = readtable(self.textfilenameL,optsL);
            TRight = readtable(self.textfilenameR,optsR);        
            [TLeft,TRight]=Smooth_whiskerDLCTLeft(TLeft,TRight,self.Numberofwhiskers);
        end

        function Whisker = get_whisker_data(self,frame_data)
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
            Whisker.R.BaseDLCX=NaN(size(frame_data.Nosex,1),self.Numberofwhiskers);
            Whisker.R.BaseDLCY=NaN(size(frame_data.Nosex,1),self.Numberofwhiskers);
            Whisker.L.BaseDLCX=NaN(size(frame_data.Nosex,1),self.Numberofwhiskers);
            Whisker.L.BaseDLCY=NaN(size(frame_data.Nosex,1),self.Numberofwhiskers);
            Whisker.R.FoX=NaN(size(frame_data.Nosex,1),self.Numberofwhiskers);
            Whisker.L.FoX=NaN(size(frame_data.Nosex,1),self.Numberofwhiskers);
            Whisker.R.FoY=NaN(size(frame_data.Nosex,1),self.Numberofwhiskers);
            Whisker.L.FoY=NaN(size(frame_data.Nosex,1),self.Numberofwhiskers);
        end
        function create_plot(self)
            figure('units','normalized','outerposition',[0 0 1 1]) 
            v1 = VideoReader([self.data_folder  vidname '.avi']);
            vL = VideoReader([self.data_folder  'Mask' vidname 'L.avi']);
            vR = VideoReader([self.data_folder  'Mirror' vidname 'R.avi']);
            subplot(2,4,1:4); 
            if good_frames(framei)>0
                Frind = find(frame_data.goodframes-1==good_frames(framei));
                im2 = readindex(v1,double(good_frames(framei)));
                imagesc(im2,[10 250])
                text(20,400,num2str(framei),'color',[1 1 1])
                set(gca,'YDir','normal')
                set(gca,'visible','off')
                hold on
                plot(frame_data.Nosex(Frind),frame_data.Nosey(Frind),'.r','markersize',20)
                plot(frame_data.Snoutx(Frind),frame_data.Snouty(Frind),'.g','markersize',20)
                text(100,100,num2str( good_frames(framei)))
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
            cmap2  = jet(size(X2,1));
            colormap gray
            set(gca,'visible','off')
            hold off
            axis image
            subplot(2,4,7);
            im2 = readindex(vL,double(framei));
            imagesc((uint8(BW)+0.5).*im2,[10 250])
            hold on
            scatter(X2',Y2', 10, cmap2, 'filled')
            axis image
            set(gca,'visible','off')
            set(gca,'XDir','reverse')
            axis image
            hold off
            clear F
            F(framei) = getframe(gcf) ;
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
<<<<<<< HEAD

        function shared_frames = get_shared_frames(self,frame_data,left_dlc,right_dlc)
=======
      
    function main(self)
        trial_names = self.get_trial_names();
        left_video = VideoReader(fullfile(self.data_folder,['Mask' trial_names{1} 'L.avi']));
        right_video = VideoReader(fullfile(self.data_folder,['Mirror' trial_names{1} 'R.avi']));
        for  vidname =trial_names
            vidname = vidname{1};
            % self.run_janelia_whisker_tracking(vidname);
            [left_dlc,right_dlc] = self.load_dlc_files();
            [frame_size_left,frame_size_right] = self.get_frame_sizes(left_video,right_video);
>>>>>>> e78c332a227d2fe640f8f36400f3d55c08d0a80d
            [distance_left,distance_right] = self.get_distance_left_and_right(left_dlc,right_dlc);
            assert(numel(distance_left)==numel(distance_right))
            nose_to_snout_distance = self.get_nose_to_snout_distance(frame_data,right_dlc);
            face_measure = (smooth(abs(zscore(distance_left))+abs(zscore(distance_right))+abs(zscore(nose_to_snout_distance)),40));
            good_frames = frame_data.goodframes-1;
            shared_frames = intersect(right_dlc.frames,left_dlc.frames);
            shared_frames = intersect(shared_frames,good_frames);
            shared_frames = shared_frames(face_measure(shared_frames+1)<=10);
        end
      
        function main(self)
            trial_names = self.get_trial_names();
            [frame_size_left,frame_size_right] = self.get_frame_sizes(trial_names);
            [left_dlc,right_dlc] = self.load_dlc_files();
            for  vidname =trial_names
                vidname = vidname{1};
                self.run_janelia_whisker_tracking(vidname);
                try
                    frame_data = readtable(fullfile(self.data_folder ,[vidname 'FrameData.xlsx']));
                catch
                    frame_data = readtable(fullfile(self.data_folder ,[vidname 'FrameData.xls']));
                end
                self.Whisker = self.get_whisker_data(frame_data);
                shared_frames = self.get_shared_frames(frame_data,left_dlc,right_dlc);
                for framei=shared_frames'+1
                    [probability,whisker_boundaries,~]= self.findWhiskers(frame_size_right,framei,right_dlc);
                    allwhisker_right = self.append_data(framei,'right',whisker_boundaries,probability);
                    if isempty(allwhisker_right.X{1})
                        continue
                    end
                    [probability,whisker_boundaries,~]= self.findWhiskers(frame_size_left,framei,right_dlc);
                    allwhisker_left = self.append_data(framei,'left',whisker_boundaries,probability);
                    if isempty(allwhisker_left.X{1})
                        continue
                    end
                    if framei>1&&size(X,1)<2
                        Whisker.L.BaseDLCX(framei,:) = Whisker.L.BaseDLCX(end,:);
                        Whisker.L.BaseDLCY(framei,:) = Whisker.L.BaseDLCY(end,:);
                    end
                    allwhisker.R = allwhisker_right;
                    allwhisker.L = allwhisker_left;
                    save(fullfile(self.data_folder,[vidname ,'.mat']),allwhisker)
                    Whisker.DLC =frame_data;
                    Whisker.Name = vidname;
                    self.create_plot()
                    toc
                end
            end
        end

        function [ThissizeL,ThissizeR] = get_frame_sizes(self,trial_names)
            left_video = VideoReader(fullfile(self.data_folder,['Mask' trial_names{1} 'L.avi']));
            right_video = VideoReader(fullfile(self.data_folder,['Mirror' trial_names{1} 'R.avi']));
            im2R = readindex(right_video,double(100));  
            im2R = im2R(:,:,1);
            ThissizeR(1) = size(im2R,1);
            ThissizeR(2) = size(im2R,2);
            im2 = readindex(left_video,double(100));  
            im2 = im2(:,:,1);
            ThissizeL(1) = size(im2,1);
            ThissizeL(2) = size(im2,2);
        end
        
        function nose_to_snout_distance = get_nose_to_snout_distance(self,frame_data,right_dlc)
            good_frames = frame_data.goodframes-1;
            theseframes =ismember(frame_data.goodframes-1, good_frames(right_dlc.frames+1));
            nose_pixel = [frame_data.Nosex(theseframes),frame_data.Nosey(theseframes)];
            snout_pixel = [frame_data.Snoutx(theseframes),frame_data.Snouty(theseframes)];
            nose_to_snout_distance = sqrt( (snout_pixel(:,1)-nose_pixel(:,1)).^2 + (snout_pixel(:,2)-nose_pixel(:,2)).^2 );
        end

        function [distance_left,distance_right] = get_distance_left_and_right(self,left_dlc,right_dlc)
            pix_1 = [left_dlc.x1,left_dlc.y1];
            pix_2 =  [left_dlc.x2,left_dlc.y2];
            distance_left = sqrt( (pix_2(:,1)-pix_1(:,1)).^2 + (pix_2(:,2)-pix_1(:,2)).^2 );
            pix_1 = [right_dlc.x1,right_dlc.y1];
            pix_2 =  [right_dlc.x2,right_dlc.y2];
            distance_right = sqrt( (pix_2(:,1)-pix_1(:,1)).^2 + (pix_2(:,2)-pix_1(:,2)).^2 );
        end
    end
end