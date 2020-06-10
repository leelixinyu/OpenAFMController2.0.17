# OpenAFMController-BetaVersionv2.0.17
Open Version with close-Loop, Simple Scan Engine
20190611

This is the README file for the Open AFM Controller provided to you by Charlène Brillard (Host & FPGA) and Adrian Nievergelt (FPGA) while they were working at Laboratory for Nano- and Bio- Instrumentation (EPFL, Lausanne, Switzerland).

This project started with the need to find a way to implement a faster imaging mode for AFM, not provided by AFM companies.

This mode is called Off-Resonance Tapping (ORT) (Peak-Force(TM) from Bruker) and can run at frequencies as high as the hardware allows.
Eventually, we implemented an even faster mode which is called Photothermal Off-Resonance Tapping (PORT), which is suitable for corresponding AFM Hardware.

The FPGA and Host Code is made available fully, the code is written is a comprehensive way. It can run without modifications any microscope, most of the modes (ORT, PORT, SICM (particularly Hopping Mode), Contact, plus any other mode which needs an external devices to be connected to like Intermittent Contact modes, Multifrequency, in fact all of the modes).
Other features are made available like the Data-Driven Resonance Filter Controller which makes any AFM 10 time faster (https://doi.org/10.1002/smtd.201900031).

It currently can run on at least three types of NI-FPGAs: 
  - USB-RIO : 7856R-OEM (10MHz), 7856R (80MHz)
  - Flex-RIO: PXIe-7975R
  - PXI: PXI-7851R, PXI-7842R, PXI-7841R
  
On top of all the advantages to have access to an AFM Controller which we hope you will enjoy using as all the current users are, you will have the opportunity to add you own features to the code: for the first time, the FPGA code for an AFM is fully modifiable.
  
Instructions on Labview (32bits, >=2016) install to run the Code, how to setup your system, and how to use the software are made available on this github repository and by contacting me directly. I would be happy to provide you support and explainations, which is how I usually proceed.

For any question, please contact cbrillard from github or by email if you can find my email address.
