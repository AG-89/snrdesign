clear variables, clc

scale_name = '12edo';
scalemap = containers.Map; %matlab hashmap
%For now: must repeat at the octave, and be 12 notes C -> note before repeat.
  %12edo (default)
    scalemap('12edo') = [2^(1/12) 2^(2/12) 2^(3/12) 2^(4/12) 2^(5/12) 2^(6/12) 2^(7/12) 2^(8/12) 2^(9/12) 2^(10/12) 2^(11/12) 2];
  %5-limit just intonation
    scalemap('just5') = [16/15 9/8 6/5 5/4 4/3 45/32 3/2 8/5 5/3 16/9 15/8 2];
scale = scalemap(scale_name); %load the scale array
scalename = scale_name; %alias

A4 = 440;
C0 = 2^(-(9+12*4)/12)*A4; %C0

f = 48.999;
fetchnote_fastf(f,A4,scalemap(scalename))
fetchnote_print(f,A4,scalename,scalemap)
Hz2Cents(f,A4,scalename,scalemap)

function f_OR = octave_reduce(f)
%octave reduce f (Hz) to range [1,2) Hz
    f_OR = exp(mod(log(f),log(2)));
end

function cents = Hz2Cents(f,A4freq,scale_name,scalemap)
%change Hz to cents based off C0 (from A4)
    if(f <= 0) %handle bad input
        cents = 0; %f > 0 Hz required
        return;
    end
    scale_ARRAY = scalemap(scale_name);
    Czero = getCzero(A4freq,scale_name,scalemap); %C0
    if(A4freq < 0.01) %override
        Czero = 1;
    end
    cents = 1200.*log(f./Czero)./log(2);
end

function f = Cents2Hz(cents,A4freq,scale_name,scalemap)
%change cents to Hz based off C0 (from A4)
    scale_ARRAY = scalemap(scale_name);
    Czero = getCzero(A4freq,scale_name,scalemap); %C0
    if(A4freq < 0.01) %override
        Czero = 1;
    end
    f = 2.^(cents./1200).*Czero;
end

%fix this stuff
function Czero = getCzero(A4freq,scale_name,scalemap)
%returns a calculated Czero value
    scale_ARRAY = scalemap(scale_name);
    scalelen = length(scale_ARRAY);
    scaleA = 0;
    switch(scale_name)
        case {'major' 'minor' 'hminor' 'mminor'}
            scaleA = 5/3;
        otherwise
            scaleA = scale_ARRAY(round(scalelen*9/12)); %interpolate
    end
    Czero = A4freq/scaleA/16;
end

function note_names = notes_from_scale(scale_name,scalemap)
%returns the note names from scale, with appropriate tonic note (default C)
%todo: add scales other than tonic C later
    scale_ARRAY = scalemap(scale_name);
    scalelen = length(scale_ARRAY);
    %all the characters used
    L_i = [1 2 3 4 5 6 7] - 0; %chromatic shift. C=1
    L = ["C" "D" "E" "F" "G" "A" "B"]; %in alpha order
    %allocate zeros
        L_f = strings(1,length(L)); %flats
        L_s = strings(1,length(L)); %sharps
        L_Df = strings(1,length(L)); %double flats
        L_Ds = strings(1,length(L)); %double sharps
        L_Hf = strings(1,length(L)); %half flats
        L_Hs = strings(1,length(L)); %half sharps
        L_32s = strings(1,length(L)); %3/2 sharps
        L_32f = strings(1,length(L)); %3/2 flats
    %generate accidentals
    for t = 1:length(L)
        L_ts_acc = '';
        L(t) = strcat(L(t), L_ts_acc); %modify base letters for tonic shift
        L_f(t) = strcat(L(t), 'b'); %flats
        L_s(t) = strcat(L(t), '#'); %sharps
        L_Df(t) = strcat(L(t), 'bb'); %double flats
        L_Ds(t) = strcat(L(t), 'x'); %double sharps
        L_Hf(t) = strcat(L(t), 'd'); %half flats
        L_Hs(t) = strcat(L(t), '/'); %half sharps
        L_32s(t) = strcat(L(t), '#/'); %3/2 sharps
        L_32f(t) = strcat(L(t), 'db'); %3/2 flats
    end
    %for 12 note scales use this
    L_12tone = [L(1) strcat(L_s(1),'/',L_f(2)) L(2) strcat(L_s(2),'/',L_f(3)) L(3) L(4) strcat(L_s(4),'/',L_f(5)) L(5) strcat(L_s(5),'/',L_f(6)) L(6) strcat(L_s(6),'/',L_f(7)) L(7)];
    %L_12tone = circshift(L_12tone,tonicshift);
    switch(scale_name)
        case 'major'
            note_names = L;
        case 'minor'
            note_names = L;
            note_names(3) = L_f(3);
            note_names(6) = L_f(6);
            note_names(7) = L_f(7);
        case 'hminor'
            note_names = L;
            note_names(3) = L_f(3);
            note_names(6) = L_f(6);
        case 'mminor'
            note_names = [L(1) L(2) L_f(3) L(4) L(5) L_f(6) L(6) L_f(7) L(7)];
        case 'octatonic'
            note_names = [L(1) L(2) L_f(3) L(4) L_f(5) L_f(6) L_Df(7) L_f(1)];
        otherwise
            %build list of generic note number names
            if(scalelen == 12)
                note_names = L_12tone;
            else
                note_names = strings(1,scalelen);
                for a = 1:scalelen
                    note_names(a) = ['''N' int2str(a-1) ''''];
                end
            end
    end
    if(scalelen ~= length(note_names))
        fprintf("! Scale note names length error\n");
    end
end

function [note_name, error_cents, target_f] = fetchnote(x,A4freq,scale_name,scalemap,fastcalc)
%return closest frequency (Hz) from scale (12 notes), target_f.
%also returns error e (cents) and note name
%octave equivalency assumed
%use fastcalc = true to only calculate target_f and save time
    if(x <= 0) %handle bad input
        fprintf("! fetchnote() only to be used with positive f input\n");
        note_name = '?';
        error_cents = 0;
        target_f = 0;
        return;
    end
    
    %fetchnote_fastf
    scale_ARRAY = scalemap(scale_name);
    scale_ARRAY = [1 scale_ARRAY]; %append 1 at beginning
    Czero = getCzero(A4freq,scale_name,scalemap); %C0
    SA_cents = Hz2Cents(scale_ARRAY,1);
    x_cents = Hz2Cents(x,Czero);
    [~,note_index] = min(abs(SA_cents - mod(x_cents,1200)));

    target_f = Cents2Hz(x_cents + SA_cents(note_index) - mod(x_cents,1200),Czero);

    if(fastcalc)
        note_name = '?';
        error_cents = 0;
        return;
    end
    %the rest
    notes = notes_from_scale(scale_name,scalemap);
    notes = [notes notes(1)]; %append octave

    x_octave = floor(Hz2Cents(target_f,Czero)/1200); %get octave number
    %error from closest pitch in cents
    error_cents = mod(x_cents,1200) - SA_cents(note_index);
    if(abs(error_cents) < 0.01) %round off tiny error
        error_cents = 0;
    end
    %select from 'notes' + octave number
    note_name = strcat(notes(note_index),num2str(x_octave)); 
    
    function cents = Hz2Cents(f,Czero)
        cents = 1200.*log(f./Czero)./log(2);
    end
    function f = Cents2Hz(cents,Czero)
        f = 2.^(cents./1200).*Czero;
    end
end

function fetchnote_print(x,A4freq,scale_name,scalemap)
%print out fetchnote information
    [n, e, f] = fetchnote(x,A4freq,scale_name,scalemap,false);
    fprintf("Closest note = %s (%.2f Hz) Error (cents) = %+.1f\n",n,f,e);
end