from gwf import Workflow
import glob
import sys, os, glob, itertools, re

gwf = Workflow(defaults={'account': 'simons'})

#################################################################################
# Utility functions
#################################################################################

def modpath(p, parent=None, base=None, suffix=None):
    par, name = os.path.split(p)
    name_no_suffix, suf = os.path.splitext(name)
#    if suffix is not None:
    if type(suffix) is str:
        suf = suffix
    if parent is not None:
        par = parent
    if base is not None:
        name_no_suffix = base

    # return os.path.join(par, name_no_suffix + suf)
    new_path = os.path.join(par, name_no_suffix + suf)
    if type(suffix) is tuple:
        assert len(suffix) == 2
        new_path, nsubs = re.subn(r'{}$'.format(suffix[0]), suffix[1], new_path)
        assert nsubs == 1, nsubs
    return new_path


def notebook(notebook_file, html_file, prev_html_files, cores=1, memory='8g', walltime='11:00:00'):

    options = {'memory': memory,
               'walltime': walltime,
               'cores': cores
              } 

    ipcluster_str = cores > 1 and 'ipcluster start -n {} &\n sleep 10'.format(cores) or ''

    spec = """
    source activate simons
    {ipcluster}
    jupyter nbconvert --to html --execute --ExecutePreprocessor.timeout=-1  --output-dir ./html {notebook}
    """.format(notebook=notebook_file, ipcluster=ipcluster_str)

    return [notebook_file]+prev_html_files, [html_file], options, spec


#prev_files = ['analysis_globals.py']
prev_files = []

def run(notebook_file, **kwargs):
	name = modpath(notebook_file, parent='', suffix='')
	html_file = modpath(notebook_file, parent='html', suffix='.html') 
	gwf.target_from_template(name, notebook(notebook_file, html_file, prev_files, **kwargs))
	prev_files.append(html_file)

# for notebook_file in sorted(glob.glob('./[0-9]*.ipynb')):
# 	run(notebook_file)

run('nb_01_pi_data.ipynb', memory='24g')
run('nb_02_auxiliary_data.ipynb', memory='24g')
run('nb_03_admixture_data.ipynb', memory='10g')
run('nb_04_argweaver_tmrca_half.ipynb', cores=3, memory='24g')
run('nb_05_call_sweeps_from_pi.ipynb', memory='50g')
run('nb_06_runs_of_windows_with_same_sweep_calling.ipynb', memory='24g', walltime='3-00:00:00')
run('nb_07_sweep_enrichments_overlaps.ipynb')
run('nb_08_two_waves_from_pi.ipynb', memory='10g')
run('nb_09_male_dist_data.ipynb', memory='100g')
run('nb_10_pruned_tmrca_stats.ipynb', memory='10g')
run('nb_11_sweeps_from_male_dist.ipynb', cores=15, memory='100g')
run('nb_12_sweep_peaks.ipynb', memory='8g')
run('nb_12_sweep_peaks_and_rel_tmrca_half.ipynb', memory='8g')
run('nb_13_sweep_regions_and_rel_tmrca_half.ipynb', memory='8g')
run('nb_14_sweep_tmrca.ipynb', memory='8g')
run('nb_15_sweep_dating.ipynb', cores=15, memory='100g')
run('nb_16_tishkoff_plots.ipynb', cores=10, memory='50g', walltime='3-00:00:00')
run('nb_17_sweeps_and_admixture.ipynb', memory='8g')
run('nb_18_sweeps_and_dist_to_africa.ipynb', memory='8g')
run('nb_19_sweeps_and_dist_to_ust_ishim.ipynb', memory='8g')
run('nb_20_sweep_enrichments.ipynb', memory='8g')
run('nb_21_alternative_explanations.ipynb', memory='8g')
run('nb_22_circRNA.ipynb', memory='8g')
run('nb_23_iHH_analysis.ipynb', memory='8g')




# run('nb_10_sweeps_from_male_dist.ipynb', cores=15, memory='100g')
# run('nb_11_sweep_peaks_and_enrichments.ipynb', memory='24g')
# run('nb_12_sweep_dating.ipynb', cores=15, memory='100g')
# run('nb_13_tishkoff_plots.ipynb', cores=10, memory='50g', walltime='3-00:00:00')
# run('nb_14_sweeps_and_admixture.ipynb', memory='24g')
# run('nb_15_sweeps_and_dist_to_africa.ipynb', memory='24g')
# run('nb_16_two_waves.ipynb')
# run('nb_17_circRNA.ipynb')
# run('nb_18_alternative_explanations.ipynb')












