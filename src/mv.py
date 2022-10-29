from src.utils.monitor_values_plus import MonitorValuesPlus as MV

mv=MV(["x","y"], out_file="out/mv.tsv", overwrite=False)

for x,y in zip(range(10),range(10)):
    mv()
