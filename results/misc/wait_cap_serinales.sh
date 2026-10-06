while squeue -h -u jstajich -o "%j" | grep -q -E "polishcap|serinales-all-882aa01"; do sleep 120; done
cd /bigdata/stajichlab/jstajich/projects/MATPredict/results
/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python 2026-09-26_polish_cap/compare.py > 2026-09-26_polish_cap/compare_output.txt 2>&1
/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python serinales_scan_analysis.py 2026-09-25_serinales_all_8dec9e1 2026-09-26_serinales_all_882aa01 > 2026-09-26_serinales_all_882aa01/analysis_vs_8dec9e1.txt 2>&1
sacct -u jstajich -S 2026-09-25 --format=JobID,JobName%24,State,Elapsed -P | grep -E "polishcap|serinales-all-882aa01" | grep -v "\."
echo finished
