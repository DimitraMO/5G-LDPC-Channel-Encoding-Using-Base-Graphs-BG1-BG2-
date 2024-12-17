# 5G-LDPC-Channel-Encoding-Using-Base-Graphs-BG1-BG2-
This is a program designed to operate the 5G channel encoding procedure using LDPC (Low Density Parity Check Codes). 
The communication channel used in this program is the AWGN channel.  
The main purspose of the given code is to generate the througput of the channel (with BER and FER metrics with respect to different signal to noise ratio (SNR) values)) and its reliability.
The program consists of two functions, nrLDPCEncoding_3.m and ldpc_find_parameters.m, where the main function that runs the whole program is the first one (the second function given is auxiliary to the first one in order to calculate the parameters of the given BG according to the 5G 3GPP TS 38.212)
To run the program the user needs to initialize the values of the variables (initial_message, BG, rate) and then pass them to the nrLDPCEncoding_3.m function. After that only the execution of the function is needed.
