<?xml version='1.0' encoding='UTF-8'?>
<Project Type="Project" LVVersion="15008000">
	<Item Name="My Computer" Type="My Computer">
		<Property Name="server.app.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="server.control.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="server.tcp.enabled" Type="Bool">false</Property>
		<Property Name="server.tcp.port" Type="Int">0</Property>
		<Property Name="server.tcp.serviceName" Type="Str">My Computer/VI Server</Property>
		<Property Name="server.tcp.serviceName.default" Type="Str">My Computer/VI Server</Property>
		<Property Name="server.vi.callsEnabled" Type="Bool">true</Property>
		<Property Name="server.vi.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="specify.custom.address" Type="Bool">false</Property>
		<Item Name="FPGA Target" Type="FPGA Target">
			<Property Name="AutoRun" Type="Bool">false</Property>
			<Property Name="configString.guid" Type="Str">{062F3BB5-B2C3-40A5-8181-89FABE3D9AE9}ResourceName=40 MHz Onboard Clock;TopSignalConnect=Clk40;ClockSignalName=Clk40;MinFreq=40000000.000000;MaxFreq=40000000.000000;VariableFreq=0;NomFreq=40000000.000000;PeakPeriodJitter=250.000000;MinDutyCycle=50.000000;MaxDutyCycle=50.000000;Accuracy=50.000000;RunTime=0;SpreadSpectrum=0;GenericDataHash=D41D8CD98F00B204E9800998ECF8427E;{23BEBDDA-449D-41F6-8E0E-56C12A0ED26C}kIntCalDebugEn=0;kIntCalibrateData=0;kIntInitGain=0;kIntResChannel=0;kIntSimEnable=0;resource=/Connector0/AI0;0;ReadMethodType=Targets\NI\FPGA\RIO\R Series\78XXR\resource\USB-7855R\Usb7855AiFxp.ctl{35E5F89A-EA1D-428A-9CDE-506ED10FFC82}"ControlLogic=1;NumberOfElements=513;Type=0;ReadArbs=Arbitrate if Multiple Requestors Only;ElementsPerRead=1;WriteArbs=Arbitrate if Multiple Requestors Only;ElementsPerWrite=1;Implementation=2;;DataType=1000800000000001003C005F03510012000000050001001200000005FFFFFFFFFFFE00000001001200000005000000000001FFFF00000001FFFFFFF4000000000000000100010000000000000000000000000000;DisableOnOverflowUnderflow=FALSE"{6A8E710E-FF06-4DD9-8D82-A7291CBD661F}"ControlLogic=1;NumberOfElements=513;Type=0;ReadArbs=Arbitrate if Multiple Requestors Only;ElementsPerRead=1;WriteArbs=Arbitrate if Multiple Requestors Only;ElementsPerWrite=1;Implementation=2;;DataType=1000800000000001003C005F03510012000000050001001200000005FFFFFFFFFFFE00000001001200000005000000000001FFFF00000001FFFFFFF4000000000000000100010000000000000000000000000000;DisableOnOverflowUnderflow=FALSE"{DF575759-346B-457B-B228-A11C4F64E568}kIntCalDebugEn=0;kIntCalibrateData=0;kIntResChannel=0;kIntSimEnable=0;resource=/Connector0/AO0;0;WriteMethodType=Targets\NI\FPGA\RIO\R Series\78XXR\resource\USB-7855R\Usb7855AoFxp.ctlUSB-7856R OEM/Clk40/falsefalseFPGA_EXECUTION_MODEFPGA_TARGETFPGA_TARGET_CLASSUSB_7856R_OEMFPGA_TARGET_FAMILYKINTEX7TARGET_TYPEFPGA/[rSeriesConfig.Begin]rseries.aio./Connector0/AI0=0,rseries.aio./Connector0/AI1=0,rseries.aio./Connector0/AI2=0,rseries.aio./Connector0/AI3=0,rseries.aio./Connector0/AI4=0,rseries.aio./Connector0/AI5=0,rseries.aio./Connector0/AI6=0,rseries.aio./Connector0/AI7=0,rseries.analogCalibratedType=1,rseries.hsdio=0,rseries.lsdio=0,rseries.terminalConfig=1[rSeriesConfig.End]</Property>
			<Property Name="configString.name" Type="Str">40 MHz Onboard ClockResourceName=40 MHz Onboard Clock;TopSignalConnect=Clk40;ClockSignalName=Clk40;MinFreq=40000000.000000;MaxFreq=40000000.000000;VariableFreq=0;NomFreq=40000000.000000;PeakPeriodJitter=250.000000;MinDutyCycle=50.000000;MaxDutyCycle=50.000000;Accuracy=50.000000;RunTime=0;SpreadSpectrum=0;GenericDataHash=D41D8CD98F00B204E9800998ECF8427E;Connector0/AI0kIntCalDebugEn=0;kIntCalibrateData=0;kIntInitGain=0;kIntResChannel=0;kIntSimEnable=0;resource=/Connector0/AI0;0;ReadMethodType=Targets\NI\FPGA\RIO\R Series\78XXR\resource\USB-7855R\Usb7855AiFxp.ctlConnector0/AO0kIntCalDebugEn=0;kIntCalibrateData=0;kIntResChannel=0;kIntSimEnable=0;resource=/Connector0/AO0;0;WriteMethodType=Targets\NI\FPGA\RIO\R Series\78XXR\resource\USB-7855R\Usb7855AoFxp.ctlinput"ControlLogic=1;NumberOfElements=513;Type=0;ReadArbs=Arbitrate if Multiple Requestors Only;ElementsPerRead=1;WriteArbs=Arbitrate if Multiple Requestors Only;ElementsPerWrite=1;Implementation=2;;DataType=1000800000000001003C005F03510012000000050001001200000005FFFFFFFFFFFE00000001001200000005000000000001FFFF00000001FFFFFFF4000000000000000100010000000000000000000000000000;DisableOnOverflowUnderflow=FALSE"output"ControlLogic=1;NumberOfElements=513;Type=0;ReadArbs=Arbitrate if Multiple Requestors Only;ElementsPerRead=1;WriteArbs=Arbitrate if Multiple Requestors Only;ElementsPerWrite=1;Implementation=2;;DataType=1000800000000001003C005F03510012000000050001001200000005FFFFFFFFFFFE00000001001200000005000000000001FFFF00000001FFFFFFF4000000000000000100010000000000000000000000000000;DisableOnOverflowUnderflow=FALSE"USB-7856R OEM/Clk40/falsefalseFPGA_EXECUTION_MODEFPGA_TARGETFPGA_TARGET_CLASSUSB_7856R_OEMFPGA_TARGET_FAMILYKINTEX7TARGET_TYPEFPGA/[rSeriesConfig.Begin]rseries.aio./Connector0/AI0=0,rseries.aio./Connector0/AI1=0,rseries.aio./Connector0/AI2=0,rseries.aio./Connector0/AI3=0,rseries.aio./Connector0/AI4=0,rseries.aio./Connector0/AI5=0,rseries.aio./Connector0/AI6=0,rseries.aio./Connector0/AI7=0,rseries.analogCalibratedType=1,rseries.hsdio=0,rseries.lsdio=0,rseries.terminalConfig=1[rSeriesConfig.End]</Property>
			<Property Name="NI.LV.FPGA.CompileConfigString" Type="Str">USB-7856R OEM/Clk40/falsefalseFPGA_EXECUTION_MODEFPGA_TARGETFPGA_TARGET_CLASSUSB_7856R_OEMFPGA_TARGET_FAMILYKINTEX7TARGET_TYPEFPGA</Property>
			<Property Name="NI.LV.FPGA.Version" Type="Int">6</Property>
			<Property Name="NI.SortType" Type="Int">3</Property>
			<Property Name="Resource Name" Type="Str">RIO5</Property>
			<Property Name="rseries.aio./Connector0/AI0" Type="Str">0</Property>
			<Property Name="rseries.aio./Connector0/AI1" Type="Str">0</Property>
			<Property Name="rseries.aio./Connector0/AI2" Type="Str">0</Property>
			<Property Name="rseries.aio./Connector0/AI3" Type="Str">0</Property>
			<Property Name="rseries.aio./Connector0/AI4" Type="Str">0</Property>
			<Property Name="rseries.aio./Connector0/AI5" Type="Str">0</Property>
			<Property Name="rseries.aio./Connector0/AI6" Type="Str">0</Property>
			<Property Name="rseries.aio./Connector0/AI7" Type="Str">0</Property>
			<Property Name="rseries.analogCalibratedType" Type="Str">1</Property>
			<Property Name="rseries.hsdio" Type="Str">0</Property>
			<Property Name="rseries.lsdio" Type="Str">0</Property>
			<Property Name="rseries.terminalConfig" Type="Str">1</Property>
			<Property Name="Target Class" Type="Str">USB-7856R OEM</Property>
			<Property Name="Top-Level Timing Source" Type="Str">40 MHz Onboard Clock</Property>
			<Property Name="Top-Level Timing Source Is Default" Type="Bool">true</Property>
			<Item Name="Connector0" Type="Folder">
				<Item Name="Connector0/AI0" Type="Elemental IO">
					<Property Name="eioAttrBag" Type="Xml"><AttributeSet name="">
   <Attribute name="kIntCalDebugEn">
   <Value>0</Value>
   </Attribute>
   <Attribute name="kIntCalibrateData">
   <Value>0</Value>
   </Attribute>
   <Attribute name="kIntInitGain">
   <Value>0</Value>
   </Attribute>
   <Attribute name="kIntResChannel">
   <Value>0</Value>
   </Attribute>
   <Attribute name="kIntSimEnable">
   <Value>0</Value>
   </Attribute>
   <Attribute name="resource">
   <Value>/Connector0/AI0</Value>
   </Attribute>
</AttributeSet>
</Property>
					<Property Name="FPGA.PersistentID" Type="Str">{23BEBDDA-449D-41F6-8E0E-56C12A0ED26C}</Property>
				</Item>
				<Item Name="Connector0/AO0" Type="Elemental IO">
					<Property Name="eioAttrBag" Type="Xml"><AttributeSet name="">
   <Attribute name="kIntCalDebugEn">
   <Value>0</Value>
   </Attribute>
   <Attribute name="kIntCalibrateData">
   <Value>0</Value>
   </Attribute>
   <Attribute name="kIntResChannel">
   <Value>0</Value>
   </Attribute>
   <Attribute name="kIntSimEnable">
   <Value>0</Value>
   </Attribute>
   <Attribute name="resource">
   <Value>/Connector0/AO0</Value>
   </Attribute>
</AttributeSet>
</Property>
					<Property Name="FPGA.PersistentID" Type="Str">{DF575759-346B-457B-B228-A11C4F64E568}</Property>
				</Item>
			</Item>
			<Item Name="40 MHz Onboard Clock" Type="FPGA Base Clock">
				<Property Name="FPGA.PersistentID" Type="Str">{062F3BB5-B2C3-40A5-8181-89FABE3D9AE9}</Property>
				<Property Name="NI.LV.FPGA.BaseTSConfig" Type="Str">ResourceName=40 MHz Onboard Clock;TopSignalConnect=Clk40;ClockSignalName=Clk40;MinFreq=40000000.000000;MaxFreq=40000000.000000;VariableFreq=0;NomFreq=40000000.000000;PeakPeriodJitter=250.000000;MinDutyCycle=50.000000;MaxDutyCycle=50.000000;Accuracy=50.000000;RunTime=0;SpreadSpectrum=0;GenericDataHash=D41D8CD98F00B204E9800998ECF8427E;</Property>
				<Property Name="NI.LV.FPGA.BaseTSConfig.Accuracy" Type="Dbl">50</Property>
				<Property Name="NI.LV.FPGA.BaseTSConfig.ClockSignalName" Type="Str">Clk40</Property>
				<Property Name="NI.LV.FPGA.BaseTSConfig.MaxDutyCycle" Type="Dbl">50</Property>
				<Property Name="NI.LV.FPGA.BaseTSConfig.MaxFrequency" Type="Dbl">40000000</Property>
				<Property Name="NI.LV.FPGA.BaseTSConfig.MinDutyCycle" Type="Dbl">50</Property>
				<Property Name="NI.LV.FPGA.BaseTSConfig.MinFrequency" Type="Dbl">40000000</Property>
				<Property Name="NI.LV.FPGA.BaseTSConfig.NominalFrequency" Type="Dbl">40000000</Property>
				<Property Name="NI.LV.FPGA.BaseTSConfig.PeakPeriodJitter" Type="Dbl">250</Property>
				<Property Name="NI.LV.FPGA.BaseTSConfig.ResourceName" Type="Str">40 MHz Onboard Clock</Property>
				<Property Name="NI.LV.FPGA.BaseTSConfig.SupportAndRequireRuntimeEnableDisable" Type="Bool">false</Property>
				<Property Name="NI.LV.FPGA.BaseTSConfig.TopSignalConnect" Type="Str">Clk40</Property>
				<Property Name="NI.LV.FPGA.BaseTSConfig.VariableFrequency" Type="Bool">false</Property>
				<Property Name="NI.LV.FPGA.Valid" Type="Bool">true</Property>
				<Property Name="NI.LV.FPGA.Version" Type="Int">5</Property>
			</Item>
			<Item Name="output" Type="FPGA FIFO">
				<Property Name="Actual Number of Elements" Type="UInt">513</Property>
				<Property Name="Arbitration for Read" Type="UInt">1</Property>
				<Property Name="Arbitration for Write" Type="UInt">1</Property>
				<Property Name="Control Logic" Type="UInt">1</Property>
				<Property Name="Data Type" Type="UInt">9</Property>
				<Property Name="Disable on Overflow/Underflow" Type="Bool">false</Property>
				<Property Name="fifo.configuration" Type="Str">"ControlLogic=1;NumberOfElements=513;Type=0;ReadArbs=Arbitrate if Multiple Requestors Only;ElementsPerRead=1;WriteArbs=Arbitrate if Multiple Requestors Only;ElementsPerWrite=1;Implementation=2;;DataType=1000800000000001003C005F03510012000000050001001200000005FFFFFFFFFFFE00000001001200000005000000000001FFFF00000001FFFFFFF4000000000000000100010000000000000000000000000000;DisableOnOverflowUnderflow=FALSE"</Property>
				<Property Name="fifo.configured" Type="Bool">true</Property>
				<Property Name="fifo.projectItemValid" Type="Bool">true</Property>
				<Property Name="fifo.valid" Type="Bool">true</Property>
				<Property Name="fifo.version" Type="Int">12</Property>
				<Property Name="FPGA.PersistentID" Type="Str">{35E5F89A-EA1D-428A-9CDE-506ED10FFC82}</Property>
				<Property Name="Local" Type="Bool">false</Property>
				<Property Name="Memory Type" Type="UInt">2</Property>
				<Property Name="Number Of Elements Per Read" Type="UInt">1</Property>
				<Property Name="Number Of Elements Per Write" Type="UInt">1</Property>
				<Property Name="Requested Number of Elements" Type="UInt">10</Property>
				<Property Name="Type" Type="UInt">0</Property>
				<Property Name="Type Descriptor" Type="Str">1000800000000001003C005F03510012000000050001001200000005FFFFFFFFFFFE00000001001200000005000000000001FFFF00000001FFFFFFF4000000000000000100010000000000000000000000000000</Property>
			</Item>
			<Item Name="input" Type="FPGA FIFO">
				<Property Name="Actual Number of Elements" Type="UInt">513</Property>
				<Property Name="Arbitration for Read" Type="UInt">1</Property>
				<Property Name="Arbitration for Write" Type="UInt">1</Property>
				<Property Name="Control Logic" Type="UInt">1</Property>
				<Property Name="Data Type" Type="UInt">9</Property>
				<Property Name="Disable on Overflow/Underflow" Type="Bool">false</Property>
				<Property Name="fifo.configuration" Type="Str">"ControlLogic=1;NumberOfElements=513;Type=0;ReadArbs=Arbitrate if Multiple Requestors Only;ElementsPerRead=1;WriteArbs=Arbitrate if Multiple Requestors Only;ElementsPerWrite=1;Implementation=2;;DataType=1000800000000001003C005F03510012000000050001001200000005FFFFFFFFFFFE00000001001200000005000000000001FFFF00000001FFFFFFF4000000000000000100010000000000000000000000000000;DisableOnOverflowUnderflow=FALSE"</Property>
				<Property Name="fifo.configured" Type="Bool">true</Property>
				<Property Name="fifo.projectItemValid" Type="Bool">true</Property>
				<Property Name="fifo.valid" Type="Bool">true</Property>
				<Property Name="fifo.version" Type="Int">12</Property>
				<Property Name="FPGA.PersistentID" Type="Str">{6A8E710E-FF06-4DD9-8D82-A7291CBD661F}</Property>
				<Property Name="Local" Type="Bool">false</Property>
				<Property Name="Memory Type" Type="UInt">2</Property>
				<Property Name="Number Of Elements Per Read" Type="UInt">1</Property>
				<Property Name="Number Of Elements Per Write" Type="UInt">1</Property>
				<Property Name="Requested Number of Elements" Type="UInt">10</Property>
				<Property Name="Type" Type="UInt">0</Property>
				<Property Name="Type Descriptor" Type="Str">1000800000000001003C005F03510012000000050001001200000005FFFFFFFFFFFE00000001001200000005000000000001FFFF00000001FFFFFFF4000000000000000100010000000000000000000000000000</Property>
			</Item>
			<Item Name="IP Builder" Type="IP Builder Target">
				<Item Name="Dependencies" Type="Dependencies"/>
				<Item Name="Build Specifications" Type="Build"/>
			</Item>
			<Item Name="sosfilter.vi" Type="VI" URL="../sosfilter.vi">
				<Property Name="configString.guid" Type="Str">{062F3BB5-B2C3-40A5-8181-89FABE3D9AE9}ResourceName=40 MHz Onboard Clock;TopSignalConnect=Clk40;ClockSignalName=Clk40;MinFreq=40000000.000000;MaxFreq=40000000.000000;VariableFreq=0;NomFreq=40000000.000000;PeakPeriodJitter=250.000000;MinDutyCycle=50.000000;MaxDutyCycle=50.000000;Accuracy=50.000000;RunTime=0;SpreadSpectrum=0;GenericDataHash=D41D8CD98F00B204E9800998ECF8427E;{23BEBDDA-449D-41F6-8E0E-56C12A0ED26C}kIntCalDebugEn=0;kIntCalibrateData=0;kIntInitGain=0;kIntResChannel=0;kIntSimEnable=0;resource=/Connector0/AI0;0;ReadMethodType=Targets\NI\FPGA\RIO\R Series\78XXR\resource\USB-7855R\Usb7855AiFxp.ctl{35E5F89A-EA1D-428A-9CDE-506ED10FFC82}"ControlLogic=1;NumberOfElements=513;Type=0;ReadArbs=Arbitrate if Multiple Requestors Only;ElementsPerRead=1;WriteArbs=Arbitrate if Multiple Requestors Only;ElementsPerWrite=1;Implementation=2;;DataType=1000800000000001003C005F03510012000000050001001200000005FFFFFFFFFFFE00000001001200000005000000000001FFFF00000001FFFFFFF4000000000000000100010000000000000000000000000000;DisableOnOverflowUnderflow=FALSE"{6A8E710E-FF06-4DD9-8D82-A7291CBD661F}"ControlLogic=1;NumberOfElements=513;Type=0;ReadArbs=Arbitrate if Multiple Requestors Only;ElementsPerRead=1;WriteArbs=Arbitrate if Multiple Requestors Only;ElementsPerWrite=1;Implementation=2;;DataType=1000800000000001003C005F03510012000000050001001200000005FFFFFFFFFFFE00000001001200000005000000000001FFFF00000001FFFFFFF4000000000000000100010000000000000000000000000000;DisableOnOverflowUnderflow=FALSE"{DF575759-346B-457B-B228-A11C4F64E568}kIntCalDebugEn=0;kIntCalibrateData=0;kIntResChannel=0;kIntSimEnable=0;resource=/Connector0/AO0;0;WriteMethodType=Targets\NI\FPGA\RIO\R Series\78XXR\resource\USB-7855R\Usb7855AoFxp.ctlUSB-7856R OEM/Clk40/falsefalseFPGA_EXECUTION_MODEFPGA_TARGETFPGA_TARGET_CLASSUSB_7856R_OEMFPGA_TARGET_FAMILYKINTEX7TARGET_TYPEFPGA/[rSeriesConfig.Begin]rseries.aio./Connector0/AI0=0,rseries.aio./Connector0/AI1=0,rseries.aio./Connector0/AI2=0,rseries.aio./Connector0/AI3=0,rseries.aio./Connector0/AI4=0,rseries.aio./Connector0/AI5=0,rseries.aio./Connector0/AI6=0,rseries.aio./Connector0/AI7=0,rseries.analogCalibratedType=1,rseries.hsdio=0,rseries.lsdio=0,rseries.terminalConfig=1[rSeriesConfig.End]</Property>
				<Property Name="configString.name" Type="Str">40 MHz Onboard ClockResourceName=40 MHz Onboard Clock;TopSignalConnect=Clk40;ClockSignalName=Clk40;MinFreq=40000000.000000;MaxFreq=40000000.000000;VariableFreq=0;NomFreq=40000000.000000;PeakPeriodJitter=250.000000;MinDutyCycle=50.000000;MaxDutyCycle=50.000000;Accuracy=50.000000;RunTime=0;SpreadSpectrum=0;GenericDataHash=D41D8CD98F00B204E9800998ECF8427E;Connector0/AI0kIntCalDebugEn=0;kIntCalibrateData=0;kIntInitGain=0;kIntResChannel=0;kIntSimEnable=0;resource=/Connector0/AI0;0;ReadMethodType=Targets\NI\FPGA\RIO\R Series\78XXR\resource\USB-7855R\Usb7855AiFxp.ctlConnector0/AO0kIntCalDebugEn=0;kIntCalibrateData=0;kIntResChannel=0;kIntSimEnable=0;resource=/Connector0/AO0;0;WriteMethodType=Targets\NI\FPGA\RIO\R Series\78XXR\resource\USB-7855R\Usb7855AoFxp.ctlinput"ControlLogic=1;NumberOfElements=513;Type=0;ReadArbs=Arbitrate if Multiple Requestors Only;ElementsPerRead=1;WriteArbs=Arbitrate if Multiple Requestors Only;ElementsPerWrite=1;Implementation=2;;DataType=1000800000000001003C005F03510012000000050001001200000005FFFFFFFFFFFE00000001001200000005000000000001FFFF00000001FFFFFFF4000000000000000100010000000000000000000000000000;DisableOnOverflowUnderflow=FALSE"output"ControlLogic=1;NumberOfElements=513;Type=0;ReadArbs=Arbitrate if Multiple Requestors Only;ElementsPerRead=1;WriteArbs=Arbitrate if Multiple Requestors Only;ElementsPerWrite=1;Implementation=2;;DataType=1000800000000001003C005F03510012000000050001001200000005FFFFFFFFFFFE00000001001200000005000000000001FFFF00000001FFFFFFF4000000000000000100010000000000000000000000000000;DisableOnOverflowUnderflow=FALSE"USB-7856R OEM/Clk40/falsefalseFPGA_EXECUTION_MODEFPGA_TARGETFPGA_TARGET_CLASSUSB_7856R_OEMFPGA_TARGET_FAMILYKINTEX7TARGET_TYPEFPGA/[rSeriesConfig.Begin]rseries.aio./Connector0/AI0=0,rseries.aio./Connector0/AI1=0,rseries.aio./Connector0/AI2=0,rseries.aio./Connector0/AI3=0,rseries.aio./Connector0/AI4=0,rseries.aio./Connector0/AI5=0,rseries.aio./Connector0/AI6=0,rseries.aio./Connector0/AI7=0,rseries.analogCalibratedType=1,rseries.hsdio=0,rseries.lsdio=0,rseries.terminalConfig=1[rSeriesConfig.End]</Property>
				<Property Name="NI.LV.FPGA.InterfaceBitfile" Type="Str">U:\projects\qnmController\filtering\FPGA Bitfiles\fxpFilter_FPGATarget_sosfilter_zeX6SNQLYNc.lvbitx</Property>
			</Item>
			<Item Name="Dependencies" Type="Dependencies">
				<Item Name="vi.lib" Type="Folder">
					<Item Name="FxpSim.dll" Type="Document" URL="/&lt;vilib&gt;/rvi/FXPMathLib/sim/FxpSim.dll"/>
					<Item Name="Clear Errors.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Clear Errors.vi"/>
					<Item Name="lvSimController.dll" Type="Document" URL="/&lt;vilib&gt;/rvi/Simulation/lvSimController.dll"/>
				</Item>
			</Item>
			<Item Name="Build Specifications" Type="Build">
				<Item Name="sosfilter" Type="{F4C5E96F-7410-48A5-BB87-3559BC9B167F}">
					<Property Name="AllowEnableRemoval" Type="Bool">false</Property>
					<Property Name="BuildSpecDecription" Type="Str"></Property>
					<Property Name="BuildSpecName" Type="Str">sosfilter</Property>
					<Property Name="Comp.BitfileName" Type="Str">fxpFilter_FPGATarget_sosfilter_zeX6SNQLYNc.lvbitx</Property>
					<Property Name="Comp.CustomXilinxParameters" Type="Str"></Property>
					<Property Name="Comp.MaxFanout" Type="Int">-1</Property>
					<Property Name="Comp.RandomSeed" Type="Bool">false</Property>
					<Property Name="Comp.Version.Build" Type="Int">0</Property>
					<Property Name="Comp.Version.Fix" Type="Int">0</Property>
					<Property Name="Comp.Version.Major" Type="Int">1</Property>
					<Property Name="Comp.Version.Minor" Type="Int">0</Property>
					<Property Name="Comp.VersionAutoIncrement" Type="Bool">false</Property>
					<Property Name="Comp.Vivado.EnableMultiThreading" Type="Bool">true</Property>
					<Property Name="Comp.Vivado.OptDirective" Type="Str"></Property>
					<Property Name="Comp.Vivado.PhysOptDirective" Type="Str"></Property>
					<Property Name="Comp.Vivado.PlaceDirective" Type="Str"></Property>
					<Property Name="Comp.Vivado.RouteDirective" Type="Str"></Property>
					<Property Name="Comp.Vivado.RunPowerOpt" Type="Bool">false</Property>
					<Property Name="Comp.Vivado.Strategy" Type="Str">Default</Property>
					<Property Name="Comp.Xilinx.DesignStrategy" Type="Str">balanced</Property>
					<Property Name="Comp.Xilinx.MapEffort" Type="Str">default(noTiming)</Property>
					<Property Name="Comp.Xilinx.ParEffort" Type="Str">standard</Property>
					<Property Name="Comp.Xilinx.SynthEffort" Type="Str">normal</Property>
					<Property Name="Comp.Xilinx.SynthGoal" Type="Str">speed</Property>
					<Property Name="Comp.Xilinx.UseRecommended" Type="Bool">true</Property>
					<Property Name="DefaultBuildSpec" Type="Bool">true</Property>
					<Property Name="DestinationDirectory" Type="Path">FPGA Bitfiles</Property>
					<Property Name="NI.LV.FPGA.LastCompiledBitfilePath" Type="Path">/U/projects/qnmController/filtering/FPGA Bitfiles/fxpFilter_FPGATarget_sosfilter_zeX6SNQLYNc.lvbitx</Property>
					<Property Name="NI.LV.FPGA.LastCompiledBitfilePathRelativeToProject" Type="Path">FPGA Bitfiles/fxpFilter_FPGATarget_sosfilter_zeX6SNQLYNc.lvbitx</Property>
					<Property Name="ProjectPath" Type="Path">/U/projects/qnmController/filtering/fxpFilter.lvproj</Property>
					<Property Name="RelativePath" Type="Bool">true</Property>
					<Property Name="RunWhenLoaded" Type="Bool">false</Property>
					<Property Name="SupportDownload" Type="Bool">true</Property>
					<Property Name="SupportResourceEstimation" Type="Bool">false</Property>
					<Property Name="TargetName" Type="Str">FPGA Target</Property>
					<Property Name="TopLevelVI" Type="Ref">/My Computer/FPGA Target/sosfilter.vi</Property>
				</Item>
			</Item>
		</Item>
		<Item Name="Dependencies" Type="Dependencies"/>
		<Item Name="Build Specifications" Type="Build"/>
	</Item>
</Project>
