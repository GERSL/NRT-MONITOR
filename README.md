# NRT-MONITOR
We developed an algorithm for Near-Real-Time MOnitoring of laNd dIsturbance based on Time-series of harmOnized Reflectance (NRT-MONITOR) [(Shang et. al., 2022)](https://www.sciencedirect.com/science/article/pii/S0034425722001870) from Landsats 7-8 and Sentinel-2 data at a 30-m spatial resolution. 

Increasing the temporal frequency of high-resolution observations, reducing the number of needed consecutive observations to confirm land disturbances, and improving the time efficiency of the NRT monitoring algorithm are the three key directions to reduce the time lag and provide timely land change information. 

To increase the temporal frequency of high-resolution observations, we will explore the use of four satellite sensors, including Landsat 7, Landsat 8, Sentinel-2A, and Sentinel-2B. 

To reduce the number of required consecutive observations to confirm land disturbances, innovative methods that could better exclude outliers and can detect land disturbance based on higher thresholds without causing omission errors of land disturbances are needed.

To increase the efficiency of the NRT algorithms, we incorporates an online recursive algorithm called Forgetting Factor to improve efficiency in the determination of land disturbance to get fast detection based on the harmonized data. 

The NRT-MONITOR algorithm can monitor land disturbance in near-real-time with or without historical detection results. 

The update of NRT land disturbance monitoring starts from the M file (***main_NRT_MONITOR.m***). It needs to read the TRA-adjusted [(Shang and Zhu, 2019)](https://www.sciencedirect.com/science/article/pii/S0034425719304584) and outlier-removed observations (***Sample_Data.mat***) and the previous monitoring results ended in the year 2014 (***Sample_previousMonitoring_end2014.mat***). The monitored results will be saved as a MAT file (***Sample_results_updatefrom2015.mat***).




## Reference

Rong Shang, Zhe Zhu, 2019. [Harmonizing Landsat 8 and Sentinel-2: A time-series-based reflectance adjustment approach](https://www.sciencedirect.com/science/article/pii/S0034425719304584). *Remote Sensing of Environment*. 235, 111439.

Rong Shang, Zhe Zhu, Junxue Zhang, Shi Qiu Zhiqiang Yang, Tian Li, Xiucheng Yang, 2022. [Near-real-time monitoring of land disturbance with harmonized Landsats 7-8 and Sentinel-2 data](https://www.sciencedirect.com/science/article/pii/S0034425722001870). *Remote Sensing of Environment*. 278, 113073.
