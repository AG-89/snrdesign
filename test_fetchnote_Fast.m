clc, clear variables
scale_array = [16/15 9/8 6/5 5/4 4/3 45/32 3/2 8/5 5/3 16/9 15/8 2];
scale_array_front1 = [1 scale_array];
SA_cents = 1200.*log(scale_array_front1)./log(2);
for a = 1:20000
    fetchnote_fast(250,1,SA_cents);
    fetchnote_fastold(250,1,scale_array);
end
fetchnote_fast(200,1,SA_cents)
fetchnote_fastold(200,1,SA_cents)

function target_f = fetchnote_fast(x,Czero,SA_cents)
%return closest frequency (Hz) from scale, target_f.
%Speed designed for detection algorithm. Input Hz2Cents(scale array)
%octave equivalency assumed
    if(x <= 0) %handle bad input
        %fprintf("! fetchnote() only to be used with positive f input\n");
        target_f = 0;
        return;
    end
    
    %scale_ARRAY = [1 scale_ARRAY]; %append 1 to front
    %SA_cents = Hz2Cents(scale_ARRAY,1); %convert to cents
    x_cents = 1200.*log(x/Czero)./log(2); %inline Hz2Cents
    [~,note_index] = min(abs(SA_cents - mod(x_cents,1200))); %for loop isn't much faster

    %inline Cents2Hz
    target_f = 2.^((x_cents + SA_cents(note_index) - mod(x_cents,1200))./1200).*Czero;
end

function target_f = fetchnote_fastold(x,Czero,scale_ARRAY)
%return closest frequency (Hz) from scale, target_f.
%Speed designed for detection algorithm. Input Hz2Cents(scale array)
%octave equivalency assumed
    if(x <= 0) %handle bad input
        %fprintf("! fetchnote() only to be used with positive f input\n");
        target_f = 0;
        return;
    end
    
    scale_ARRAY = [1 scale_ARRAY]; %append 1 at beginning
    SA_cents = Hz2Cents(scale_ARRAY,1);
    x_cents = Hz2Cents(x,Czero);
    [~,note_index] = min(abs(SA_cents - mod(x_cents,1200)));

    target_f = Cents2Hz(x_cents + SA_cents(note_index) - mod(x_cents,1200),Czero);
    
    function cents = Hz2Cents(f,Czero)
        cents = 1200.*log(f./Czero)./log(2);
    end
    function f = Cents2Hz(cents,Czero)
        f = 2.^(cents./1200).*Czero;
    end
end