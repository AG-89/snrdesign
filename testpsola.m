f = 261;
samplerate = 44100;
x = 1:2000;
y = square(2*pi*f.*x/samplerate);
volume = 25;

pitch_scale = 1.25;

y_pm = PitchMarking(y, y, samplerate);
y_new = PSOLA(y, samplerate, y_pm, 1, pitch_scale)

if(true) %plays the sounds
    volpeakpoint = 100/volume;
    y_toplay = [volpeakpoint y];
    fprintf("Playing y...\n");
    soundsc(y_toplay,(samplerate))
    pause
    clear sound, clear sounds;
    y_toplay = [volpeakpoint y_new];
    fprintf("Playing y pitch shifted...\n\n");
    soundsc(y_toplay,(samplerate_tuned))
    %pause
end