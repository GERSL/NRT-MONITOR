# NRT-MONITOR
We developed an algorithm for Near-Real-Time MOnitoring of laNd dIsturbance based on Time-series of harmOnized Reflectance (NRT-MONITOR) from Landsats 7-8 and Sentinel-2 data at a 30-m spatial resolution. 

Increasing the temporal frequency of high-resolution observations and improving the time efficiency of the NRT monitoring algorithm are the two key solutions to reduce the time lag and provide management-related land change information. 

To increase the temporal frequency of high-resolution observations, we will explore the use of four satellite sensors, including Landsat 7, Landsat 8, Sentinel-2A, and Sentinel-2B. 

To increasing the efficiency of the NRT algorithms, we incorporates an online recursive algorithm called Forgetting Factor to improve efficiency in the determination of land disturbance to get fast detection based on the harmonized data. 

The NRT-MONITOR algorithm can monitor land disturbance in near-real-time with and without historical detection results. 

The update of NRT land disturbance monitoring starts from the M file (main_NRT_MONITOR.m). It needs to read the TRA-adjusted and outlier-removed observations (Sample_Data.mat) and the previous monitoring results ended in the year 2014 (Sample_previousMonitoring_end2014.mat). The monitored results will be saved as a MAT file (Sample_results_updatefrom2015.mat).
