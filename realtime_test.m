clc, clear variables

%frameLength = 1024;
Fs = 44100; %samplerate

%recorder object
AudioInput = audiorecorder(Fs,16,1);

%record for 3 seconds
record(AudioInput, 3);
disp('recording 3 seconds')

y = [];
%wait to start recording
while ~isrecording(AudioInput)
    pause(0.01)
end
%loop while recording
while isrecording(AudioInput)
  try
    y = getaudiodata(AudioInput); %save to array
  catch
  end
  length(y)
  pause(0.05)
end

%write to output file
[~,~] = mkdir('RToutput/');
audiowrite('RToutput/Y.wav',y,Fs);

%player object
%AudioOutput = audioplayer(y,Fs);

%time scope
scope = dsp.TimeScope(...
   'SampleRate',AudioInput.SampleRate,...
   'TimeSpan',16,...
   'BufferLength',2.4e6,...
   'YLimits',[-1,1]);

release(scope) %free