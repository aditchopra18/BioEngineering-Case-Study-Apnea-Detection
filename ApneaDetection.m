%% Reading Data from the file
% Opening the header file
DataFile = fopen ('a03.hea','r');
test = fscanf (DataFile, '%s',2); % To skip values

% Opening the ECG data file
fid = fopen ('a03.dat','r');
ECG = fread (fid, 'int16');

% Extracting the Sampling Frequency
SampFreq = fscanf(DataFile, '%f', 1);

%% Extracting four hours of information and Interpolating

% Extracting 4 hours of data from the ECG data
ECG_Data_Before = ECG (1 : 1440000);

N =  SampFreq; % Storing the Sampling Frequency in N
t = 1/N : 1/N : 14400; % Creating the time instances
t = t.'; % Taking the Transpose

N2 = 500; % Creating new Sampling Frequency for Interpolation
t2 = 1/N2 : 1/N2 : 14400; % Creating the time instances
t2 = t2.'; % Taking the Transpose
ECG_Data_After = interp1(t, ECG_Data_Before, t2); % Interpolating from 100 Hz to 500 Hz

% Plotting two seconds of data before and after interpolating
figure
plot (t(1 : 200), ECG_Data_Before(1 : 200));
figure
plot (t2(1 : 1000), ECG_Data_After(1 : 1000));
figure
plot (t(1 : 200),ECG_Data_Before(1:200));
hold on;
plot (t2(1:1000),ECG_Data_After(1:1000));

% Formatting the Graphs by Axis Labels and Title
title('Test Case 1: Interpolation')
xlabel('Time (in seconds)');
ylabel('Signal Amplitude (in volts)');
legend('Sampling Frequency - 100Hz', 'Sampling Frequency - 500HZ');
hold off;

%% Filtering the data using a Low-Pass Filter

% Designing a Low Pass Filter using "designfilt"
d = designfilt('lowpassfir', 'Filterorder', 20, 'CutoffFrequency', 11, 'SampleRate', N2);
ECG_Data_Filtered = filter(d, ECG_Data_After); % Applying the low pass filter

% Plotting two seconds of data before and after filtering to see the difference
figure
plot (t2 (5000 : 6000), ECG_Data_After (5000 : 6000));

hold on;
plot (t2 (5000 : 6000), ECG_Data_Filtered (5000 : 6000));
% Note that only the new Sample Frequency Data is used in the Low-Pass Filter

% Formatting the graph with appropriate Title and Axis Labels
title('Test Case 2: Filtration')
xlabel('Time (in seconds)');
ylabel('Signal Amplitude (in volts)');
legend('Raw','Filtered');
hold off;

%% Time-Domain Analysis 

% Creating empty arrays
Average_Hour = [];
Max_Hour = [];
Min_Hour = [];

HR = [];
Peaks = [];

% Setting start and end points
Start = 1;
End = 1;

for i = 1 : 30000 : 7200000

    % Reading through each 60 seconds of data for the four hours duration
    newArr = ECG_Data_Filtered(i:i+29999);
    for j = 2 : length (newArr) - 1
        % Condition for calculating peaks
        if (newArr (j) > 100 && (abs (newArr (j) - newArr (j-1)) >= 0.5) && (newArr (j) > newArr (j+1)) && (newArr (j) > newArr (j-1)))
            Peaks = [Peaks, j];
            End = End+1;
        end
    end

    for k = Start + 1 : End - 1
        % Measuring the RR-values and the subsequent Heart Rate values
        diff = Peaks(k) - Peaks(k-1);
        HR = [HR, 60/(diff * 0.002)];

    end
    
    % Appending heart rate information into corresponding arrays
    Average_Hour = [Average_Hour, mean(HR)];
    Max_Hour = [Max_Hour, max(HR)];
    Min_Hour = [Min_Hour, min(HR)];

    % Empyting the HR array
    HR = [];

    % Assiging new start point
    Start = End;
end

% Creating time axis
t3 = [1:240];
t3 = t3.'; % Taking the Transpose

% Reading data from Annotated Files
Apnea_File = fopen ('a03.apn.txt', 'r');

if (Apnea_File == -1)
    error("File couldn't be opened"); % Printing Error Message if File not Found
end

% Skipping data
data = fscanf (Apnea_File, '%s', 3);

% Creating array to store annotated values in forms of numbers
apnea = [];
for i = 1 : 240
    data = fscanf(Apnea_File, '%s', 1);
    if (data == 'N')
        apnea = [apnea, 50]; % Assigning 50 for Non-Apnea Condition
    end
    if (data == 'A')
        apnea = [apnea, 100]; % Assiging 100 for apnea
    end
    data = fscanf(Apnea_File, '%s', 1); % Skipping values
end

% Plotting the Average Heart Rate graph and the Apnea Data graph
figure
plot(t3 (2 : 240), Average_Hour(2 : 240));

hold on;
plot(t3 (2 : 240), apnea(2 : 240));

% Formatting the Graph with Appropriate Title and Axis Labels
title('Test Case 3: Time-Domain Analysis')
xlabel('Time (in mins)');
ylabel('Average Heart Rate (in BPM)');
legend('Average Heart Rate','Presence of Apneas');
hold off;

%% Downsampling and Frequency vs Power Spectrum

% Assigning new sampling frequency of 4Hz
New_Sampling_Freq = 4;
Resampled_ECG = resample(ECG_Data_Before, New_Sampling_Freq, N); % Resampling the data
Resampled_ECG_60 = Resampled_ECG (1 : 240); % Extracting first 60 seconds of data
ff60 = fft (Resampled_ECG_60); % Applying fft
Num1 = length (ff60);
freq = -New_Sampling_Freq/2:New_Sampling_Freq/Num1:New_Sampling_Freq/2-New_Sampling_Freq/Num1; %defining the x-axis
%Storing power spectral density info
PSDX = (1/New_Sampling_Freq*Num1)*abs(ff60).^2; 
PSDX = PSDX(1:Num1/2);

%plotting the graphs
figure
plot(freq,ff60);
title('Test Case 4: Frequency vs Power Spectrum')
xlabel('Frequency (in Hz)');
ylabel('Power (in volts)');

%% Frequency-Domain Analysis 
PSDarr = [];

% Performing the analysis for every 60 seconds in the 4 hour period
for i = 1 : 240 : 57600
    newFreq = Resampled_ECG(i:i+239);
    ffn  = fft (newFreq);
    NumN = length (ffn);
    freqN = -New_Sampling_Freq / 2 : New_Sampling_Freq / NumN : New_Sampling_Freq / 2 - New_Sampling_Freq / NumN;
    PSDn = (1 / New_Sampling_Freq * NumN) * abs(ffn).^2;
    PSDn = PSDn(1 : NumN / 2);
    PSDarr = [PSDarr, [PSDn]]; % Appending the PSDn array into the PSDarr

end

% Colour Plotting the Information
figure
pcolor (PSDarr);

hold on;
plot(t3, Average_Hour, 'r');
plot(t3, apnea, 'y');

% Formatting the Graph with Appropriate Title and Axis Labels
title('Test Case 5: Frequency-Domain Analysis')
xlabel('Time (in minutes)');
ylabel('Frequency (in Hz)');
legend('Power Density Graph','Averge Heart Rate','Presence of Apneas');
hold off;

%% Algorithm Efficiency Calculation

check = [];
Threshold = 63; % Threshold value for Apnea Detection based on observation
for i=1 : 240
    if Average_Hour(i) <= Threshold
        check = [check, 100];
    else
        check = [check, 50];
    end
end

% Algorithm efficiency
Efficiency = 0;
for i = 1 : 240
    if apnea(i) == check(i)
        Efficiency = Efficiency + 1;
    end
end

% Displaying final efficiency in Percentage
finaleff = Efficiency * 100/240;