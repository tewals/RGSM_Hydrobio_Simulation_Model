# USU Rio Grande Silvery Minnow Alternative Hydrograph Strategies Simulation Model: Getting Started Guide

Updated Sept. 2024

## Analyzing Alternative Approaches to a Single, Specific Annual Hydrograph
### Step 1:
Extract all files from the .ZIP file (Appendix A). Save these into the desired location on your hard drive or cloud storage system.

### Step 2: 
Within the extracted “USU_RGSM_AlternativeHydrograph_SimulationModel_update2024” folder, open the R file “USU_RGSM_AlternativeManagement_SimulationModel_update2024.R”

### Step 3: 
Set up the correct working directories. On lines 43-45, you will need to replace the “…” with the correct file location on your system. Insert the correct path for the location you saved the extracted files. For example, if you extracted the files to your “Documents” folder, the pathway would be similar to: 

“C:/Users/YourName/Documents/USU_RGSM_AlternativeHydrograph_SimulationModel_update2024”

Put the correct pathway for your system in place of the ellipses, taking care not to delete the directory extensions following the ellipses.

Alternatively, if you are running R through RStudio (not necessary to run the code), you can use the “Files” panel to navigate to the correct directory, click on the gear symbol labeled “More”, select “Set Working Directory”, and then copy and paste the correct directory path into the model code.

If you have correctly set up the directories, you should be able to run lines 43-45, followed by lines 67 through 70 without any error messages being initiated.

### Step 4. 
Determine and input your baseline and alternative hydrographs. Within the folder “InputData”, open the file titled “RGSM_AlternativeHydrographInput_AllowKnownDrying.csv”. The first four columns (year, day, month, doy) should not be changed. All following columns can be changed to reflect alternative hydrographs. As many alternative hydrographs as you wish can be considered, simply add them as new columns on the right. Ideally, each alternative will have a unique identifier as its column header.

***New update***: Users can now specify the drying extent in the San Acacia and Isleta reaches rather than having the model predict the quantiles of plausible summer drying from the hydrograph. This capability can be useful when comparing hypothetical future scenarios to the predicted performance of hydrographs from historical years with known drying. If you want to specify the drying extent, enter the mile-days dry for each reach in the first two rows of the table below the headers in the appropriate column (Isleta in row #2, San Acacia in row #3). If you do not wish to specify a known extent of drying, set these rows equal to -1.

 ![image](https://github.com/user-attachments/assets/d17ac6dd-8188-4ea6-ae84-45d1cee977e4)

 
  The hydrograph submitted under the “Baseline” column should be entered as mean daily discharge in CFS. The alternative strategies may be entered in one of three ways, though all alternative strategies examined should be entered in a consistent format. (Note: The example shown above is entered as the difference in discharge from the baseline hydrograph).
	
 First, alternatives may be submitted as daily discharge values in CFS per day. If this is the case, when you are running the simulations in R, you will need to set the object named “alternative.format” equal to “cfs” on line 86.

 ![image](https://github.com/user-attachments/assets/b42c34ee-7c25-4850-a005-3bb9d548389e)

  Second, alternatives may be supplied as daily differences from the baseline hydrograph (shown in figure above), with daily values being positive or negative values indicating how many more or less CFS were discharged on a given day. If alternatice strategies are formatted in this way, you will need to set the object named “alternative.format” equal to “differences” on line 86.

 ![image](https://github.com/user-attachments/assets/40383d04-1c10-4cc3-9df8-566f04aa070d)

  Finally, alternatives may be provided as proportional changes from the baseline. For example, if an alternative strategy calls for storing 10% of the flow on January 1, the alternative strategy would have a value of 0.9 for January 1. If the strategy called for increasing flows by 10% on January 2, it would have a value of 1.1 for January 2. If alternative strategies are formatted this way, you will need to set the object named “alternative.format” equal to “proportion” on line 86.

 ![image](https://github.com/user-attachments/assets/d488570d-c8f4-40c3-b364-82ddc551e32d)

  **IMPORTANT!** If you will be switching between formats between model runs, be sure to set the “alternative.format” value to the correct option.
	Once all desired alternative strategies have been input into the file, save the file without changing its name.

### Step 5: 
Run the model. Within R, run lines 86-186. If no errors are thrown (some warning messages may be triggered when loading packages – these are ok), continue on to run lines 207-323. This will run and store the model simulations within R. The code will produce a running index in the R Console of which strategy is being simulated so you can gage progress. 

Line 323 (“save.image(…)”) will save your R Workspace so you can reload the model outputs into R at a later time. This defaults to a generically named file to your specified output directory.

### Step 6: 
Make figures. The remainder of the code calculates relative performance of management strategies (by ranking their performance relative to other strategies within each stochastic simulation and plotting the distribution of ranks across stochastic simulations) and produces output figures. By running lines 335-519, the code will produce two PDF files of figures summarizing the relative performance of the alternative strategies at the scale of the MRG (“Simulation_Figures.pdf”) and at the scale of individual reaches broken out by different drying scenarios (“Simulation_Figures_ReachSpecificProbabilities.pdf”). The PDF files will be saved in the output directory specified on line 45.


## Exploring Generalized Single Year Management Strategies Across Randomized Hydrographs
### Step 1: 
Extract all files from the .ZIP file (Appendix A). Save these into the desired location on your hard drive or cloud storage system.

### Step 2: 
Within the extracted “USU_RGSM_AlternativeHydrograph_SimulationModel_update2024” folder, open the R file “USU_RGSM_SingleYearManagement_GeneralizedSimulations_update2024.R”

### Step 3: 
Set up the correct working directories. On lines 45-47, you will need to replace the “…” with the correct file location on your system. Insert the correct path for the location you saved the extracted files. For example, if you extracted the files to your “Documents” folder, the pathway would be similar to:

“C:/Users/YourName/Documents/USU_RGSM_AlternativeHydrograph_SimulationModel_update2024”

Put this pathway in place of the ellipses, taking care not to delete the directory extensions following the ellipses.

Alternatively, if you are running R, through RStudio (not necessary to run the code), you can use the “files” panel to navigate to the correct directory, click on the gear symbol labeled “More”, select “Set Working Directory”, and then copy and paste the correct directory path into the model code.

If you have correctly set up the directories, you should be able to run lines 45-46, followed by lines 54-137 without any error messages being initiated.

### Step 4: 
Specify alternative strategies in CSV file. In Excel, open the file “ManagementStrategies_SingleYear.csv”, then specify the alternative flow management strategies you would like to examine. Important: Alternative management strategies are specified differently for the generalized simulations than for the specific hydrograph simulations described above. 

Give each strategy a number in the first column, and a short name in the second column. Specify the type of strategy you are specifying (currently only “prop” is available in the model). “prop” type strategies store a specified proportion of flow during months in which flows are stored for later use. Identify any months in which you would like to store flows (i.e., reduce discharge in river for later release) with a “-1”, any months in which you would like to supplement flows (i.e., increase discharge with stored flows from previous months and/or discretionary flows) with a “1”, and months you do not wish to alter with a “0”. If you would like to use discretionary water in your management strategy, place a “-1” in the “Discretionary” column. 

Once all desired alternative strategies have been input into the file, **save the file without changing its name**.

![image](https://github.com/user-attachments/assets/356dc6eb-efe2-40fb-8901-98570ca684d2)

### Step 5: 
Set up management scaling parameters. In R, on lines 148 through 153, specify the scale at which management actions are implemented. The “scls” parameter specifies the proportions of water that are stored during storage periods. The “damt.acft” parameter specifies the volume of discretionary water available to managers. The “sl.acft" parameter specifies the storage capacity in acre feet.

### Step 6: 
Run the model. Within R, run lines 1-205 to load all packages, data, parameters, and set management strategies. If no errors are thrown (some warning messages may be triggered when loading packages – these are ok), continue on to run lines 212-335. This will run and store the model simulations within R. Models will be run in parallel by default. ***This may consume a lot of computing power unless you specify a small number of cores to use for processing on line 215***. To specify a specific number of cores (using 2 as an example here) to run the simulations across, change line 215 from:

![image](https://github.com/user-attachments/assets/4044473d-ca41-4d74-96fe-f41ed24bfaf4)

to: ![image](https://github.com/user-attachments/assets/d48582c9-a26f-4f2e-aefd-3831bde1a77b)

### Step 7: 
Make figures. Within R, run lines 340-579. This will produce the figures found in the Walsworth and Budy (2022) report. *Note: if you examine different storage, discretionary water, and proportional flow management scenarios, you will need to change axis labels on lines 574-576.*

## Exploring Generalized Multi-Year Management Strategies Across Randomized Hydrographs
### Step 1: 
Extract all files from the .ZIP file (Appendix A). Save these into the desired location on your hard drive or cloud storage system.

### Step 2: 
Within the extracted “USU_RGSM_AlternativeHydrograph_SimulationModel_update2024” folder, open the R file “USU_RGSM_MultiYearManagement_GeneralizedSimulations_update2024.R”

### Step 3: 
Set up the correct working directories. On lines 48-51, you will need to replace the “…” with the correct file location on your system. Insert the correct path for the location you saved the extracted files. For example, if you extracted the files to your “Documents” folder, the pathway would be similar to: 

“C:/Users/YourName/Documents/USU_RGSM_AlternativeHydrograph_SimulationModel_update2024”

Put this pathway in place of the ellipses, taking care not to delete the directory extensions following the ellipses.

Alternatively, if you are running R, through RStudio (not necessary to run the code), you can use the “files” panel to navigate to the correct directory, click on the gear symbol labeled “More”, select “Set Working Directory”, and then copy and paste the correct directory path into the model code.

If you have correctly set up the directories, you should be able to run lines 48-51, followed by lines 58-90 without any error messages being initiated.

### Step 4: 
Specify alternative strategies in CSV file. In Excel, open the file “ManagementStrategies_SingleYear.csv”, then specify the alternative flow management strategies you would like to examine. Important: Alternative management strategies are specified differently for the generalized simulations than for the specific hydrograph simulations described above. 

Give each strategy a number in the first column, and a short name in the second column. Each strategy needs to have an approach to all water year types (high, low, mid, highlow -  a low year following a high year, highmid – a mid year following a high year, and midlow – a low year following a mid year), and will therefore require six lines per strategy. Specify the type of strategy you are specifying (currently only “prop” is available in the model). “prop” type strategies store a specified proportion of flow during months in which flows are stored for later use. 

Each line represents a water year type (listed above), and will need to have a monthly action specified. Identify any months in which you would like to store flows with a “-1”, any months in which you would like to supplement flows with a “1”, and months you do not wish to alter with a “0”. If you would like to use discretionary water in your management strategy, place a “-1” in the “Discretionary” column. If you would like water to be stored in a given year for use in subsequent years, enter a “1” in the “storage” column. If you would like to use water stored from previous years in a given water year type, enter “-1” in the “Storage” column.

Once all desired alternative strategies have been input into the file, **save the file without changing its name**.

![image](https://github.com/user-attachments/assets/4a993d92-6701-4501-a7f5-ca652765b2ca)

### Step 5: 
Set up management scaling parameters. In R, on lines 148 through 155, specify the scale at which management actions are implemented. The “scls” parameter specifies the proportions of water that are stored during storage periods. The “damt.acft” parameter specifies the volume of discretionary water available to managers. The “sl.acft" parameter specifies the storage capacity in acre feet.

### Step 6: 
Run the model. Within R, run lines 1-203. If no errors are thrown (some warning messages may be triggered when loading packages – these are ok), continue on to run lines 209-376. This will run and store the model simulations within R. Models will be run in parallel by default. ***This may consume a lot of computing power unless you specify a small number of cores to use for processing on line 212***. To specify a specific number of cores (using 2 as an example here) to run the simulations across, change line 21 from:

![image](https://github.com/user-attachments/assets/3dbf6a4c-a22d-432d-a763-f8c127be4b49)

to: 
![image](https://github.com/user-attachments/assets/97bdd8a2-7a9c-4ed6-9f8f-d214abfe2429)

### Step 7: 
Make figures. Within R, lines 382-653 will create two PDF files of figures comparing the performance of the alternative strategies across a range of management scale scenarios. To group similar strategies by color on the figures, open the file "ManagementStrategies_MultiYear_ColorLineIndicators.csv" in Excel and provide an integer (1-4) for each strategy. Values given a 1 will be plotted in blue, those given a 2 will be plotted in orange, those given a 3 will be plotted in green, and those given a 4 will be in black. These colors can be changed on line 408. If you do not specify colors in the csv file, the four colors will be used sequentially for each strategy. 

Run lines 382-653. This will produce the figures found in the report. *Note: if you examine different storage, discretionary water, and proportional flow management scenarios, you will need to change axis labels on lines 463-465, 517-519, 560-562, and 612-614*.




