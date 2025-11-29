%% QPSK Receiver with ADALM-PLUTO Radio in Simulink
% This model shows how to use the the ADALM-PLUTO Radio with Simulink(R) to
% implement a QPSK receiver. The receiver addresses practical issues in
% wireless communications, e.g. carrier frequency and phase offset, timing
% offset and frame synchronization. This model receives the signal sent by
% the
% <matlab:openExample('plutoradio/QPSKTransmitterWithADALMPLUTORadioSimulinkExample','supportingFile','plutoradioQPSKTransmitterSimulinkExample')
% QPSK Transmitter with ADALM-PLUTO Radio> model. The receiver demodulates
% the received symbols and outputs a simple message to the MATLAB(R)
% command line.

% Copyright 2017-2023 The MathWorks, Inc.
%% Overview
% This model performs all processing at complex baseband to handle a
% time-varying frequency offset, a time-varying symbol delay, and Gaussian
% noise. To cope with the above-mentioned impairments, this example
% provides a reference design of a practical digital receiver, which
% includes correlation-based coarse frequency compensation, symbol timing
% recovery with fixed-rate resampling and bit stuffing/skipping, fine
% frequency compensation, frame synchronization and phase ambiguity
% resolution. The example uses some key algorithms in MATLAB,
% emphasizing textual algorithm expression over graphical algorithm
% expression.
%% Structure of the Example
% The top-level structure of the model is shown in the following figure,
% which includes a ADALM-PLUTO receiver block, a QPSK Receiver subsystem and a BER
% Display blocks.
modelname = 'plutoradioQPSKReceiverSimulinkExample';
open_system(modelname);
set_param(modelname, 'SimulationCommand', 'update')
%%
% The detailed structures of the QPSK Receiver
% subsystem are illustrated in the following figure.
%%
%
open_system([modelname '/QPSK Receiver']);

%%
% The components are further described in the following sections.
%
% * *AGC* - Automatic gain control 
% * *Raised Cosine Receive Filterring* - Uses a rolloff factor of 0.5
% * *Coarse Frequency Compensation* - Estimates an approximate frequency
% offset of the received signal and corrects it
% * *Symbol Synchronization* - Resamples the input signal according to a
% recovered timing strobe so that symbol decisions are made at the optimum
% sampling instants
% * *Carrier Synchronization* - Compensates for the residual frequency
% offset and the phase offset
% * *Frame Synchronization* - Detects location of the frame header and
% aligns the frame boundaries at the known frame header
% * *Decode Data* - Resolves the phase ambiguity caused by the Carrier
% Synchronizer, demodulates the signal, and decodes the text message
%% Receiver
%% AGC
% The received signal amplitude affects the accuracy of the carrier and
% symbol synchronizer. Therefore the signal amplitude should be stabilized
% to ensure an optimum loop design. The AGC output power is set to a value
% ensuring that the equivalent gains of the phase and timing error
% detectors keep constant over time. The AGC is placed before the Raised
% Cosine Receive Filter so that the signal amplitude can be measured with
% an oversampling factor of two, thus improving the accuracy of the
% estimate. You can refer to Chapter 7.2.2 and Chapter 8.4.1 of [ <#16 1> ]
% for details on how to design the phase detector gain.

%% Raised Cosine Receive Filtering
% The Raised Cosine Receive Filter provides matched filtering for 
% the transmitted waveform with a rolloff factor of 0.5.
%% Coarse Frequency Compensation
% The Coarse Frequency Compensation subsystem corrects the input signal 
% with a rough estimate of the frequency offset. The following diagram
% shows the subsystem, in which the frequency offset is estimated by
% averaging the output of the correlation-based algorithm of the Coarse
% Frequency Compensator block. The compensation is performed by the
% Phase/Frequency Offset block. There is usually a residual frequency
% offset even after the coarse frequency compensation, which would cause a
% slow rotation of the constellation.  The Carrier Synchronizer block
% compensates for this residual frequency.
%
% The accuracy of the Coarse Frequency Compensator decreases with its
% maximum frequency offset value. Ideally, this value should be set just
% above the expected frequency offset range.
%
open_system([modelname '/QPSK Receiver/Coarse Frequency Compensation']);
%% Symbol Synchronization
% The timing recovery is performed by a Symbol Synchronizer library
% block, which implements a PLL, described in Chapter 8 of [ <#16 1> ], to
% correct the timing error in the received signal. The timing error
% detector is estimated using the Gardner algorithm, which is rotationally
% invariant. In other words, this algorithm can be used before or after
% frequency offset compensation. The input to the block is oversampled by
% two. On average, the block generates one output symbol for every two
% input samples. However, when the channel timing error (delay) reaches symbol 
% boundaries, there will be one extra or missing symbol in the output frame. 
% In that case, the block implements bit stuffing/skipping and generates
% one more or less samples comparing to the desired frame size. So the output
% of this block is a variable-size signal.
%
% The _Damping factor_, _Normalized loop bandwidth_, and _Detector gain_
% parameters of the block are tunable. Their default values are set to 1
% (critical damping), 0.01 and 5.4 respectively, so that the PLL quickly
% locks to the correct timing while introducing little timing jitter.
%% Carrier Synchronization
% The fine frequency compensation is performed by a Carrier Synchronizer
% library block, which implements a phase-locked loop (PLL), described in
% Chapter 7 of [ <#16 1> ], to track the residual frequency offset and the
% phase offset in the input signal. The PLL uses a Direct Digital
% Synthesizer (DDS) to generate the compensating phase that offsets the
% residual frequency and phase offsets. The phase offset estimate from DDS
% is the integral of the phase error output of a Loop Filter.
%
% The _Damping factor_ and _Normalized loop bandwidth_ parameters of the
% block are tunable. Their default values are set to 1 (critical damping)
% and 0.01 respectively, so that the PLL quickly locks to the intended
% phase while introducing little phase noise.
%% Frame Synchronization
% The Frame Synchronizer MATLAB System block uses the known frame header
% (QPSK-modulated Barker code) to correlate against the received QPSK
% symbols to locate the frame header and align the frame boundaries. It
% also transforms the variable-size output of the Symbol Synchronizer block
% into a fixed-size frame, which is necessary for the downstream
% processing. The second output of the block is a boolean scalar indicating
% if the first output is a valid frame with the desired header and if so,
% enables the Data Decoding subsystem to run.
%% Decode Data
% The Data Decoding enabled subsystem performs phase ambiguity resolution,
% demodulation and text message decoding. The Carrier Synchronizer block
% may lock to the unmodulated carrier with a phase shift of 0, 90, 180, or
% 270 degrees, which can cause a phase ambiguity. For details of phase
% ambiguity and its resolution, please refer to Chapter 7.2.2 and 7.7 in [
% <#16 1> ]. The Phase Offset Estimator subsystem determines this phase
% shift. The Phase Ambiguity Correction & Demodulation subsystem rotates
% the input signal by the estimated phase offset and demodulates the
% corrected data. The payload bits are descrambled, and then decoded. All
% of the stored bits are converted to characters are printed in the
% Simulink Diagnostic Viewer.
close_system([modelname '/QPSK Receiver/Coarse Frequency Compensation']);
open_system([modelname '/QPSK Receiver/Data Decoding']);
%% Running the Example
% Before running this model, connect two ADALM-PLUTO Radios to the computer
% and start the
% <matlab:openExample('plutoradio/QPSKTransmitterWithADALMPLUTORadioSimulinkExample','supportingFile','plutoradioQPSKTransmitterSimulinkExample')
% QPSK Transmitter with ADALM-PLUTO Radio> model first. This receiver model
% is capable of handling a frequency offset of 12.5kHz between the
% transmitter and receiver boards. However, when the frequency offset
% exceeds this range, the Coarse Frequency Compensation subsystem cannot
% accurately determine the offset of the received signal, which is critical
% for correct timing recovery and data decoding. Due to hardware variations
% among the ADALM-PLUTO Radios, a frequency offset will likely exist
% between the transmitter hardware and the receiver hardware. In that case,
% perform a manual frequency calibration using the companion frequency
% offset calibration
% <matlab:openExample('plutoradio/FrequencyOffsetCalibrationWithADALMPLUTORadioSimulinkExample','supportingFile','plutoradiofreqcalib')
% transmitter> and
% <matlab:openExample('plutoradio/FrequencyOffsetCalibrationWithADALMPLUTORadioSimulinkExample','supportingFile','plutoradiofreqcalib_rx')
% receiver> models and examine the resulting behavior. With that frequency
% offset value, you can manually adjust the _Center frequency_ of the
% ADALM-PLUTO Receiver subsystem in the receiver model to ensure a
% residual frequency offset that the model can track.
% 
% If the received signal is too weak or too strong, you might notice some
% garbled message output. In that case, you can change the gain of either
% the ADALM-PLUTO Transmitter subsystem in the
% <matlab:openExample('plutoradio/QPSKTransmitterWithADALMPLUTORadioSimulinkExample','supportingFile','plutoradioQPSKTransmitterSimulinkExample')
% QPSK Transmitter with ADALM-PLUTO Radio> model or the ADALM-PLUTO
% Receiver subsystem in the current model for better reception.
% 
% Set the _Center frequency_ parameter of the ADALM-PLUTO Radio Receiver
% block according to the center frequency setting of the
% <matlab:openExample('plutoradio/QPSKTransmitterWithADALMPLUTORadioSimulinkExample','supportingFile','plutoradioQPSKTransmitterSimulinkExample')
% QPSK Transmitter with ADALM-PLUTO Hardware> model and the frequency
% calibration result. Then run the model. To ensure real-time processing,
% the model is by default set to run in Accelerator mode, and to remove all
% signal visualization. The received messages are decoded and printed out
% in the |View diagnostics| window while the simulation is running. To
% ensure real-time operation, run the models on two separate MATLAB
% sessions; one for transmitter and the other for receiver.
close_system([modelname '/QPSK Receiver/Data Decoding']);
%% Exploring the Example
% The example allows you to experiment with multiple system capabilities to
% examine their effect on bit error rate performance.
% 
% You can tune the _Normalized loop bandwidth_ and _Damping factor_
% parameters of the Symbol Synchronizer and Carrier Synchronizer
% blocks, to assess their convergence time and estimation accuracy. In
% addition, you can assess the pull-in range of the Carrier Synchronizer
% block. With a large _Normalized loop bandwidth_ and _Damping factor_, the
% PLL can acquire over a greater frequency offset range. However a large
% _Normalized loop bandwidth_ allows more noise, which leads to a large
% mean squared error in the phase estimation. "Underdamped systems (with
% Damping Factor less than one) have a fast settling time, but exhibit
% overshoot and oscillation; overdamped systems (with Damping Factor
% greater than one) have a slow settling time but no oscillations." [ <#16
% 1> ]. For more detail on the design of these PLL parameters, you can
% refer to Appendix C in [ <#16 1> ].
%
% You can also tune the Preamble detector threshold and see its effects on
% the output message.
close_system(modelname, 0);
%% References
% 1. Michael Rice, "Digital Communications - A Discrete-Time
% Approach", Prentice Hall, April 2008.
