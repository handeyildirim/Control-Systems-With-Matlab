# Control-Systems-With-Matlab
You can find the explanations above :
This question examines the frequency components of the audio signals obtained by pressing the phone keys.
so it is about determining which key is pressed.
On phones, a sound occurs as a result of pressing each key. These sounds are two frequencies
It consists of the component. These include dual tone multi frequency (DTMF)
It called. Search this topic by searching with these keywords. Briefly (a few paragraphs) own
Put them in your report by summarizing them with your sentences (not copy-paste).
Take, for example, keys from zero to nine. Each key (number) and when this key is pressed
The frequency components of the sound are given in the table below.

Table 1: Numbers and frequency components
Rakam Frekans 1 Frekans 2
0     941       1336
1     697       1209
2     697       1336
3     697       1477
4     770       1209
5     770       1336
6     770       1477
7     852       1209
8     852       1336
9     852       1477

A sample sound file is given on the site as sound01.wav. 
Download it and listen. You will notice that this sound is generated by pressing nine keys on the phone.
Read and learn the audioread MATLAB command documentation. 
Read the audio file with this command, get the sampling frequency (fs) and draw the signal. 
You are expected to get a graphic like the following for the sample file. 
Put the resulting graphic in your report and interpret it.

Analyze the signal to identify the parts that correspond to the key tones. 
Obtain the frequency response of each part with the fft command. 
Draw the frequency components of the signal in a single figure using the subplot command. 
On each of the sub-figures, write the two frequency components and the corresponding key (number) 
that form the corresponding sub-signal as the title. You are expected to create a graphic like the 
following for the sample audio file. Put the resulting graphics in your report and interpret them.

NOTES:
The sampling frequency of the signal in the sample audio file is 8000 Hz. 
Therefore, the maximum frequency is 8000/2 = 4000 Hz and the frequency axis in the graphics should 
be in the range of 0-4000 Hz. Your code should be as general as possible. It should work not only 
for the sample sound file provided but also for other files consisting of phone key sounds. Some information about it:
You can assume that there will be a maximum of nine key tones in the file. But it may be less.
Pressing times may not be fixed. Silence times between the two keys may also not be constant. 
These can also vary from file to file.
It may or may not be a silent track at the beginning or end of the file.
Therefore, your program should be able to automatically find where the key was pressed and where there was silence.
You just need to recognize the sounds corresponding to the numbers 0-9. 
You don't need to recognize the sounds for the # and * keys. 
The frequencies given in Table 1 are universal for key sounds. 
The audio file will not contain any other frequency.
You can generate different sound files by recording from the computer by pressing the keys of your phone to try your code. 
Or you can search for websites that you can create and download DTMF audio files for certain key sequences. 
Direct MATLAB
You can also find and use codes that generate DTMF signals under. 
However, you can only use them to create audio files; 
Of course, it is forbidden to use the codes that do the opposite (detecting keys in the audio file) correctly.
In the light of the above explanations, create two more sound files besides the sample file given from the question. 
Submit these files together with assignments. Get results by running your program for them, put them in your report and interpret them.
NOTE: Please note that we can not make any changes in your code, the explanations such as “The same code will be executed, 
only the filename in this line will be changed” are invalid. 
You can create and send different .m files for different audio files.
