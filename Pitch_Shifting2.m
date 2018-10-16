%Pitch Shifting
clc, clear variables, clear sound, clear sounds, close all;

f = 261;
filename = 'File.wav'; %accept file
wp = 900;    %pass band in hz
wz = 1000;   %stop band in hz
%Filter signal : -2 = pass band gain dB , -20 = stop band gain dB
[num,den] = makefilter('lowpass','butter',-2,-20,wp,ws);
%Perform filtering
runfilter(num,den,'File.wav');
x = sin(2*pi*f);
%Center clip to reduce noise
CL = 0.3 * max(x);
    if x > CL
        C(x) = x - CL;
    else
        abs(x) <= CL;
        x = 0;
    else
        x < -CL;
        C(x) = x + CL;
    else
        x = 0;
    end
%Define best overlap offset seeking window for 15 ms
W = 0.015;  
T = f^-1;   %time period of current frame
rT =        %auto-correlated signal
Xi =        %new window
Xi_prev =       %previous window
    for i=1 : W-T
     rT = Xi - Xi_prev;
     PT = rT;    %pitch period
    end
rT = x1;     %gives best matching window for best pitch mark
%Apply median filter to reduce noise
    for k = 1:len
        y(k,1) = median(x1(k:k+L-1));
    end
%Locate maximum pitch and find global maximum of waveform
M = max(y(k,1));     %first pitch mark at "t"
%Locate pitch marks left/right of maximum pitch mark M
tm + f * To <= t <= tm + (2-f) * To;
tm - f * To <= t <= tm - (2-f) * To;
function pitch_marks = px = find_pmarks(M)
%Interpolate pitch contour
function v = interp[x,px,M]
end
%Input pitch contour , SR = search region of speech segments of input wave
SR = [tm + f. To ; tm + (2-f). To];
SR = v[x,px,M]
%Get optimal pitch mark sequence
P (k,j) = max[P(k-1;i).*log(tk(i,j))].*log(vk(j))];
TSF = 1;    %Read output file/define time scale factor
sequence = 0.025;   %Define processing sequence for 11k Hz sample rate
overlap = 0.020;    %Define overlapping size
seek_window = 0.015;    %Define best overlap offset seeking window
input_prev = P(k,j);    %Use optimal peak mark sequence
%Use crosscorrelation function to find best overlapping offset
%where input_prev and input_new best match
best_overlap(sample*input_prev, sample*input_new);
%Pre-calculate overlapping slopes with input_prev
    for i=0 : i<overlap
        temp[i] = (float)(input_prev*i*(overlap-i);
    end
%Find best overlap within (0:Seek_window)
    for i=0 : i<seek_window
        crosscor += (float)input_new[i+j]*temp[i];
    end
%Overlap input_prev with input_new by sliding the amplitude during overlap
%samples and store result to output
output[i] = (input_prev[i] * (overlap-i) + input_new[i]*i/overlap);
%Perform time scaling for sample data given in input and result to output
%Return # of output samples and copy flat mid-sequence from current
%processing sequence to output
int(num_out_samples) = 0;
memcpy(output,seq_offset,flat_duration*sizeof(sample));
%Calculate a pointer to theoretical next processing sequence begin
input += sequence_skip-overlap;
%Seek actual best matching offset using crosscor
seq_offset = input + seek_best-overlap(prev_offset,input);
%Overlap between previous and new sequence, copy result to output
overlap(output +flat_duration,prev_offset,seq_offset);
%Update input and sequence pointers by overlapping amount
sequence - offset += overlap;
input += overlap;
%Update output pointer and sample counters
ouput += sequence - overlap;
num_out_samples += sequence_overlap;
num_in_samples -= sequence_skip;
return num_out_samples - samples; 












