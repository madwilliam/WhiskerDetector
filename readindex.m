function outputFrame=readindex(videoSource,frameNumber)
info=get(videoSource);
videoSource.CurrentTime=(frameNumber-1)./info.FrameRate;
outputFrame=readFrame(videoSource);