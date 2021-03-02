echo "Starting up jeb"
python3 jeb.py
echo "Jeb finished successfully! Hoorah! Now let's plot everything."
python3 Plots/makePlots.py
python3 Plots/plotFeatImp.py
python3 Plots/plotPerfMetric.py