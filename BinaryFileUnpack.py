import numpy as np
class BinaryFileUnpack:
    '''
    Class to neatly unpack binary files. Provides the header information as well as the data for the sensors.

    Methods:
        __init__: Method to initialize a BinaryFileUnpack object.
        spectra (np.ndarray): Returns the frequencies (in Hertz) and the power (in deciBels) present in the spectrum of the data. Utilizes the scipy module to return the spectrum.
        plot_static: Outputs a static plot of the sensor data against a time series. Utilizes the matplotlib module.
        plot_interactive: Outputs a plot of the sensor data against a time series. Implements a HoverTool and a CrossHairTool for the user to analyze the data.
        getDtype: Method to read the bytes of a file as a given data type. Main purpose is to be used in __init__ method.
        endOfFile (bool): Checks whether the end of file has been reached. Main purpose is to be used in __init__ method.
    '''
    def __init__(self, fileName:str):
        '''
        Method to initialize a BinaryFileUnpack object.

        Parameters:
            fileName: The file path for the binary file to be analyzed

        Instance variables:
            fileName (str): File path of the input binary file.
            fs (float): The number of samples taken per second (sampling rate). Measured in Hertz.
            header_info (dict: [str, Any]): Header info from parsing the metadata of the binary file
                Keys: Description of the data.
                Values: Value of the data (can be int, float, or ndarray depending on the data).
            num_sens (int): Number of sensors recording data.
            P (ndarray): A 2-dimensional array that returns the pressure data for each sensor.
                First axis: The sensor indices (0 to self.num_sens-1, inclusive).
                Second axis: The pressure data for the respective sensor for a time series (given by self.time).
            T (ndarray): A 2-dimensional array that returns the temperature data for each sensor.
                First axis: The sensor indices (0 to self.num_sens-1, inclusive).
                Second axis: The temperature data for the respective sensor for a time series (given by self.time).
            time (ndarray): The time series for the data given in self.P and self.T.
            offset: A variable to track the number of bytes transversed in the binary file. Primarily meant for internal use.
        '''

        self.fileName = fileName.replace('\\', '/')
        self.offset = 0

        ## Parsing Metadata
        with open(fileName, "rb") as file:
            self.fileContent = file.read()   

        fileVersion = self.getDtype('i', 4)
        self.fs = self.getDtype('f', 4)
        devCount = self.getDtype('I', 4)

        DevID = np.empty(devCount, np.int32)
        SNL = np.empty(devCount, np.uint32)
        NameL = np.empty(devCount, np.uint32)
        NumEnChan = np.empty(devCount, np.uint32)
        # Lists to be reshaped into arrays later
        SN_lst = []
        Name_lst = []
        ChanNum_lst = []
        for i in range(devCount):
            DevID[i] = self.getDtype('i', 4)
            
            SNL[i] = self.getDtype('I', 4)
            SN_lst.append(self.getDtype('b', 1, SNL[i]))
            
            NameL[i] = self.getDtype('I', 4)
            Name_lst.append(self.getDtype('b', 1, NameL[i]))

            NumEnChan[i] = self.getDtype('I', 4)
            ChanNum_lst.append(self.getDtype('i', 4, NumEnChan[i]))
        
        # Converting Lists into Arrays
        SN = np.array(SN_lst)
        Name = np.array(Name_lst)
        chanNum = np.array(ChanNum_lst)

        # Names might not be in the order of the sensors, so get order of sensors
        self.num_sens:int = 2*devCount
        order_name = np.argsort(Name[:, -1])
        order = np.empty(self.num_sens, 'int')
        for i in range(len(order_name)):
            order[2*i] = 2*order_name[i]
            order[2*i+1] = 2*order_name[i]+1
        
        SN = SN[order_name, :]
        Name = Name[order_name, :]

        # Store all header information in dictionary
        self.header_info = {
            "File Version": fileVersion, "Sampling Rate": self.fs,
            "Device Count": devCount, "Device ID": DevID,
            "Serial Number Length": SNL, "Serial Number": SN,
            "Name Length": NameL, "Name": Name,
            "Number of Enabled Channels": NumEnChan, "Channel Number": chanNum            
        }

        ## Parsing data
        status = False  # EOF marker
        c = 0
        Trel = 0
        self.data = np.empty((int(1e7), NumEnChan[0], devCount))
        NS_sum = 0
        while not status:
            for i in range(devCount):
                Tsec = self.getDtype('Q', 8)
                TNsec = self.getDtype('I', 4)
                NS = self.getDtype('I', 4)
                NS_sum+=NS
                Nt = NS // NumEnChan[i]
                d = self.getDtype('f', 4, NumEnChan[i]*Nt).reshape((Nt, NumEnChan[i]))
                if c >= self.data.shape[0]:
                    tmp = np.empty((self.data.shape[0]*10, NumEnChan[0], devCount))
                    tmp[:c, :, :] = self.data
                    self.data = tmp
                self.data[c:c+Nt, :, i] = d

            c += Nt
            status = BinaryFileUnpack.endOfFile(self)
            
        self.data = self.data[:NS_sum//12]

        ## Getting Temperature and Pressure Data
        def temp_convert(temp) -> float:
            'Converts measured temp data from volts to deg Celsius'
            return (temp - 1.478) / 0.01746 + 25


        P = np.empty((self.num_sens, self.data.shape[0]))
        T = np.empty((self.num_sens, self.data.shape[0]))
        for i in range(devCount):
            T[2*i]   = temp_convert(self.data[:, 1, i])
            T[2*i+1] = temp_convert(self.data[:, 3, i])
            P[2*i]   = self.data[:, 0, i]
            P[2*i+1] = self.data[:, 2, i]

        self.P = P[order]
        self.T = T[order]
        # Creating time series
        self.time = np.linspace(0, (self.P.shape[1])/self.fs, self.P.shape[1])

        ## Apply pressure and temperature corrections by serial number
        # Numpy array with all the sensor corrections, first column is serial number, 
        # second column is P correction, third column is P correction std. dev.
        # fourth column is T correction (volts), fifth column is T correction std. dev.
        sens_corr = np.array([
            [5122769, -0.0047, 0.0029,  0.0147, 0.0091],
            [5122770,  0.0007, 0.0024, -0.0044, 0.0185],
            [5122777,  0.0034, 0.0104, -0.0013, 0.0090],
            [5122778,  0.0182, 0.0046, -0.0015, 0.0102],
            [5940428, -0.0046, 0.0031,  0.0019,	0.0089],
            [5940430, -0.0131, 0.0031, -0.0094, 0.0133],
        ])

        # Numpy array with sensor SN with index corresponding to position
        sens_used = np.array([5122778, 5122769, 5940428, 5122770, 5122777, 5940430])

        # # Apply temperature and pressure corrections
        for i in range(len(sens_used)):
            # See if correction can be applied to the sensor
            ind_arr = np.where(sens_corr[:, 0] == sens_used[i])[0]
            if len(ind_arr) > 0:
                j = ind_arr[0]
                self.P[i] += sens_corr[j, 1]
                self.T[i] += temp_convert(sens_corr[j, 3] + 1.478) - 25

        # Futher calibration to the atmospheric pressure (linear map)
        # Apply correction as Pcorr = A*P + B
        # first column is serial number
        # Second column are A coefficents
        # Third column are B coefficents
        sens_atm = np.array([
            [5122778,  1.05029088,  -0.298257335],
            [5122769,  1.04999036,  -0.294259218],
            [5940428, 0.874466484, -0.0786998490],
            [5122770,  1.17073226,  -0.452581661],
            [5122777, 0.933776111,  -0.158336162],
            [5940430,  1.00000000,   0.000000000]
        ])

        # # Apply corrections to sensors
        # for i in range(len(sens_used)):
        #     # See if correction can be applied to the sensor
        #     ind_arr = np.where(sens_corr[:, 0] == sens_used[i])[0]
        #     if len(ind_arr) > 0:
        #         j = ind_arr[0]
        #         self.P[i] = sens_atm[j, 1] * self.P[i] + sens_atm[j, 2]
        
        # Futher calibration to the atmospheric pressure (constant map)
        # Apply correction as Pcorr = A*P + B
        # first column is serial number
        # Second column are constant corrections
        sens_atm = np.array([
            [5122778, -0.23539],
            [5122769, -0.23194],
            [5940428, -0.23564],
            [5122770, -0.23862],
            [5122777, -0.24152],
            [5940430,  0.00000],
        ])

        # Apply corrections to sensors
        for i in range(len(sens_used)):
            # See if correction can be applied to the sensor
            ind_arr = np.where(sens_corr[:, 0] == sens_used[i])[0]
            if len(ind_arr) > 0:
                j = ind_arr[0]
                self.P[i] = self.P[i] + sens_atm[j, 1]


    def spectra(self, data_spec: np.ndarray) -> np.ndarray:
        '''
        Returns the frequencies (in Hertz) and the power (in deciBels) present in the spectrum of the data. Utilizes the scipy module to return the spectrum.

        Parameters:
            data_spec (ndarray): Frequencies contatining the data for each sensor. The array is 2-dimensional: the first axis is the sensor index and the second axis is the data.
        Returns:
            Pxx (ndarray): The resulting 3-dimensional array storing frequency and power data of the complete spectrum.
                First axis: Index 0 is frequency data, index 1 is power data.
                Second axis: Contains the sensor indices.
                Third axis: The data of either frequency or power for the respective sensor.
        '''
        # Getting Frequencies for Spectral Analysis
        import scipy.signal as signal
        # Length of first axis indicates number of sensors
        # Length of second axis is number of data points

        # Reshape array to 2D if one sensor is available
        if len(data_spec.shape) == 1:
            data_spec = np.reshape(data_spec, (1, len(data_spec)))

        spec_sens = data_spec.shape[0]
        N = data_spec.shape[1]

        # Periodogram approach
        Pxx = np.empty((2, spec_sens, N//2 + 1))
        for i in range(spec_sens):
            f, Pper_spec = signal.periodogram(data_spec[i], self.fs, 'cosine', scaling='density')
            power = 10 * np.log10(Pper_spec)
            Pxx[:, i, :] = np.stack([f, power])

        return Pxx

    def plot_static(self, x:np.ndarray, y:np.ndarray, x_label:str, y_label:str, plots_shape:tuple, color:str='blue', 
                    y2:np.ndarray=None, y2_label:str=None, color2:str='red', x_axis_type:str='linear', ordering:list=None,
                    x_range:tuple=None, sharex=True, plot_type='line', t1=None, t2=None, figfilename:str=None):
        '''
        Outputs a static plot of the sensor data against a time series. Utilizes the matplotlib module.

        Parameters:
            x (ndarray): A (1-D) array of the time series or (2-D) array of the frequency or measurement data. If 2-D, axis convention follows that of the param y.
            y (ndarray): A 2-dimensional array with the sensor data or power spectrum. 
                First axis: The sensor index.
                Second axis: Contains the data for that sensor.
            x_label (str): The name of the data that defines the x-axis.
            y_label (str): The type of sensor data along with its units.
            plots_shape (tuple): The shape of the plots in the output in terms of number of rows and columns. Input should be (rows, cols).
            color (Any): Identifies the line_color feature of each line glyph for the plots. Allowed inputs are those allowed by line_color. Default is 'blue'.
            y2 (ndarray): An optional second array that is of the same dimension as :param: y. Will be plotted alongside y with an additional axis.
            y2_label (str): The type of secondary sensor data along with its units.
            color2 (str): The color of the optional second data array. Default is 'red'.
            x_axis_type ('linear', 'log'): The scale of the x-axis, either can be a linear axis (='linear;) or logarithmic axis (='log'). Default is 'linear'.
            ordering (list): A custom ordering of the sensor plot data. Indicate the sensor index (starting from 0), not the number.
                            Will throw a AssertError if len(ordering) != y.shape[0] or max(ordering) != y.shape[0] - 1.
            x_range (tuple): The lower and upper bounds of the x-axis.
            sharex (bool): If true, the plots will share the same x-scale. Otherwise, they will have independent scales.
            plot_type ('line', 'scatter'): Defines the type of plot displayed. Default is 'line'.
        
        Raises:
            ValueError:
                The plot dimension does not match the number of plots (same as the number of sensors).
        '''
        import matplotlib.pyplot as plt

        if plots_shape[0]*plots_shape[1] != y.shape[0]:
            raise ValueError(f"Plot dimension does not match number of plots. Plot Dimension: {plots_shape}, Number of plots: {y.shape[0]}")

        fig, ax = plt.subplots(plots_shape[0], plots_shape[1], sharex=sharex, figsize=(10, 5/3*self.num_sens), tight_layout=True)

        fig.text(0.5, -0.015, x_label, ha='center')
        fig.text(-0.015, 0.5, y_label, va='center', rotation='vertical')

        if plots_shape[0] == 1:
            ax = np.reshape(ax, (1, plots_shape[1]))
        elif plots_shape[1] == 1:
            ax = np.reshape(ax, (plots_shape[0], 1))
        
        # Create a seperate vertical axis if secondary data source is provided
        if y2 is not None:
            assert y2.shape == y.shape, f"Additional data shape of {y2.shape} does not match primary data shape of {y.shape}."

            if y2_label is not None:
                fig.text(1.015, 0.5, y2_label, va='center', rotation=-90)

            ax2 = np.empty(ax.shape, dtype=object)
            for i in range(ax2.shape[0]):
                for j in range(ax2.shape[1]):
                    ax2[i][j] = ax[i][j].twinx()
        
        if ordering is None:
            ordering = range(y.shape[0])
        assert len(ordering) == y.shape[0] and max(ordering) == y.shape[0]-1, f"Make sure ordering contains correct sensor numbers from 0 to {y.shape[0]-1}."
        for count, ind_sens in enumerate(ordering):
            # Determining where sensors are on plot
            i = count // ax.shape[1]
            j = count % ax.shape[1]

            if x.shape[0] == y.shape[0]:
                x_select = x[ind_sens]
            else:
                x_select = x
        
            # Plotting graphs
            if plot_type == 'line':
                ax[i][j].plot(x_select, y[ind_sens], color=color)
                if y2 is not None:
                    # Plot second data source
                    ax2[i][j].plot(x_select, y2[ind_sens], color=color2)

            elif plot_type == 'scatter':
                if color == "time":
                    # Colormap by time
                    ax[i][j].scatter(x_select, y[ind_sens], c=self.time[t1:t2], cmap="summer", s=2)
                else:
                    ax[i][j].scatter(x_select, y[ind_sens], color=color, s=2)
                if y2 is not None:
                    # Plot second data source
                    ax2[i][j].scatter(x_select, y2[ind_sens], color=color2, s=2)
            
            # Additional plot features
            ax[i][j].set_xscale(x_axis_type)
            ax[i][j].set_title(f"Sensor {ind_sens + 1}")
            if x_range is not None:
                if len(x_range) != 2: raise IndexError(f"x_range should contain only 2 values, has {len(x_range)}.")
                range_x = x_range[1] - x_range[0]
                ax[i][j].set_xlim((x_range[0] - range_x*0.1, x_range[1] + range_x*0.1))

        if figfilename is not None:
            fig.savefig(figfilename,transparent=True)

    def plot_interactive(self, x:np.ndarray, y:np.ndarray, x_label:str, y_label:str, 
                         plots_shape:tuple, color=None, x_axis_type:str='linear', output_format:str='file'):
        '''
        Outputs a plot of the sensor data against a time series. Implements a HoverTool and a CrossHairTool for the user to analyze the data.
        Parameters:
            x (ndarray): A (1-D) array of the time series or (2-D) array of the frequency data. If 2-D, axis convention follows that of the param y.
            y (ndarray): A 2-dimensional array with the sensor data or power spectrum. 
                First axis: The sensor index.
                Second axis: Contains the data for that sensor.
            x_label (str): The name of the data that defines the x-axis.
            y_label (str): The type of sensor data along with its units.
            plots_shape (tuple): The shape of the plots in the output in terms of number of rows and columns. Input should be (rows, cols).
            color (Any): Identifies the line_color feature of each line glyph for the plots. Allowed inputs are those allowed by line_color (default is element in bokeh.palettes.Turbo6).
            x_axis_type ('linear', 'log'): The scale of the x-axis, either can be a linear axis (='linear;) or logarithmic axis (='log'). Default is 'linear'. 
            output_format ('file', 'notebook'): Outputting a file in either a .html file (='file') or within the notebook (='notebook'). Default is 'file'.
        
        Raises:
            ValueError: 
                output_format is not either 'file' or 'notebook'.
                The plot dimension does not match the number of plots (same as the number of sensors).
        '''
        from bokeh.plotting import figure, output_file, output_notebook, reset_output, show, ColumnDataSource
        from bokeh.layouts import gridplot, column
        from bokeh.models import CustomJS
        from bokeh.models.tools import HoverTool, CrosshairTool, BoxZoomTool, WheelZoomTool, SaveTool, ResetTool
        from bokeh.models.widgets import RangeSlider
        
        if plots_shape[0]*plots_shape[1] != y.shape[0]:
            raise ValueError(f"Plot dimension does not match number of plots. Plot Dimension: {plots_shape}, Number of plots: {self.num_sens}")

        # Initialize the plots into a numpy array
        plots = np.empty(self.num_sens, dtype='object')
        for i in range(self.num_sens):
            # Determine the x_axis range
            x_to_use = x if len(x.shape)==1 else x[i]
            x_range = (x_to_use[0], x_to_use[-1])
            
            # Logarithm Case - make sure there are no zeros involved
            if x_axis_type == 'log':
                if x_range[0] <= 0:
                    x_range = (x_to_use[1], x_range[1])
                
                x_range = (10**np.floor(np.log10(x_range[0])), 10**np.ceil(np.log10(x_range[1])))

            plots[i] = figure (
                title = f"Sensor {i+1}",
                sizing_mode = "stretch_width",
                plot_height = 675//plots_shape[0],
                x_range=x_range,
                x_axis_type=x_axis_type,
            )

        # Create ColumnDataSources with x and y data for each sensor
        sources = np.empty(self.num_sens, dtype=ColumnDataSource)
        for i in range(self.num_sens):
            sources[i] = ColumnDataSource({'x': x[i] if x.shape[0] == y.shape[0] else x , f'Sensor {i+1}':y[i]})

        # Double-ended range slider so user can select data in a given time frame
        x_to_use = x if len(x.shape)==1 else x[0]
        time_range = RangeSlider(start=x_to_use[1], end=x_to_use[-1], value=(x_to_use[1], x_to_use[-1]), step=0.01, sizing_mode="stretch_width")
        for plot in plots:
            time_range.js_link('value', plot.x_range, 'start', attr_selector=0)
            time_range.js_link('value', plot.x_range, 'end', attr_selector=1)

        # Adjust the y-axis for each plot on change of time_range
        time_range.js_on_change('value', CustomJS(
            args=dict(plots=plots, sources=sources), 
            code="""
            // Getting x-data (shared between all sources)
            var x_arr = sources[0].data['x'];
            var x_range = cb_obj.value;
            // Getting indices of range slider
            var start_ind = x_arr.indexOf(x_range[0]);
            var end_ind = x_arr.indexOf(x_range[1]);
            if (start_ind == -1) {
                start_ind = 0;
            }
            if (end_ind == -1) {
                end_ind = x_arr.length - 1;
            }
            // Iterate through all the sensors
            for (let i = 0; i < plots.length; i++) {
                // Sensor number
                var sensNum = 'Sensor '+(i+1);
                
                // Getting y-data
                var y_arr = sources[i].data[sensNum];
                var y_range_arr = y_arr.slice(start_ind, end_ind);
                
                // Getting max, min, and range of data
                var max = Math.max(...y_range_arr.filter((el) => !isNaN(el)));
                var min = Math.min(...y_range_arr.filter((el) => !isNaN(el)));
                var range = max - min;
                
                // Dynamically changing the axis range
                plots[i].y_range.start = min - range/25;
                plots[i].y_range.end = max + range/25;
                console.log(min)
            }
            """
            ))

        if color is None:
            from bokeh.palettes import Turbo6
            import random 
            color = random.choice(Turbo6)

        hover_tools = np.empty(self.num_sens, dtype=object) 
        for i in range(self.num_sens):
            hover_tools[i] = HoverTool(
                tooltips = [
                    (x_label, f"@x"),
                    (y_label, f"@{{Sensor {i+1}}}")
                ],
                renderers = [plots[i].line(x='x', y=f"Sensor {i+1}", line_color = color, source=sources[i])]
            )
            plots[i].tools = [
                hover_tools[i], 
                CrosshairTool(dimensions="height", line_alpha=0.5, line_color="purple"), 
                BoxZoomTool(), WheelZoomTool(), SaveTool(), ResetTool()
            ]

        # Adding axis labels
        grid_layout = np.reshape(plots, plots_shape)
        # x-axis
        for i in range(grid_layout.shape[1]):
            grid_layout[-1][i].xaxis.axis_label = x_label
        # y-axis
        for i in range(grid_layout.shape[0]):
            grid_layout[i][0].yaxis.axis_label = y_label

        # Formatting the Graphs
        grid_plots = gridplot(list(grid_layout), sizing_mode = "stretch_width")
        grid = column(grid_plots, time_range, sizing_mode="stretch_width")

        reset_output()
        if output_format == 'file':
            y_strip = y_label.strip(' \\><:\"|?*$')
            x_strip = x_label.strip(' \\><:\"|?*$')
            output_file(f"{self.fileName.split('.')[0]}-{y_strip}-{x_strip}.html")
        elif output_format == 'notebook':
            output_notebook()
        else:
            raise ValueError(f"{output_format} is not a valid option for output_format. Allowed options are 'file' and 'notebook'")
        show(grid)

    def getDtype(self, typeStr:str, numBytes:int, numObjects:int=1):
        '''
        Method to read the bytes of a file as a given data type. Main purpose is to be used in __init__ method.

        Parameters:
            typeStr (str): String containing data types of requested objects. See https://docs.python.org/3/library/struct.html#format-characters
            numBytes (int): The number of bytes for all data types indicated in @typeStr
            numObjects (int): The number of the inputted typeStr to read. Default is 1.

        Returns:
            Either an object of type np.ndarray or the indicated data type (param type) containing said data types read from the file.
        '''
        import struct
        tmp = self.offset
        self.offset += numBytes*numObjects
        # Total number of requested objects
        totNumObjects = int(numObjects*len(typeStr))
        if totNumObjects == 1:
            return struct.unpack(typeStr, self.fileContent[tmp:self.offset])[0]
        return np.array(struct.unpack(typeStr*numObjects, self.fileContent[tmp:self.offset]))

    def endOfFile(self) -> bool:
        '''
        Checks whether the end of file has been reached. Main purpose is to be used in __init__ method.

        Returns:
            (bool) True if end of file is reached; otherwise False.
        '''
        return self.offset == len(np.fromfile(self.fileName, dtype='byte'))
