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
    source conda_init.sh
    conda activate simons
    {ipcluster}
    jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=-1 {notebook}
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

if not os.path.exists('html'):
    os.makedirs('html')

run('nb_01_pi_data.ipynb', memory='24g')
run('nb_02_auxiliary_data.ipynb', memory='24g')
run('nb_03_admixture_data.ipynb', memory='10g')
run('nb_04_argweaver_tmrca_half.ipynb', cores=10, memory='24g')
run('nb_05_call_sweeps_from_pi.ipynb', memory='50g')
run('nb_06_runs_of_windows_with_same_sweep_calling.ipynb', memory='24g', walltime='3-00:00:00')
run('nb_07_sweep_enrichments_overlaps.ipynb')
run('nb_08_two_waves_from_pi.ipynb', memory='10g')
run('nb_09_male_dist_data.ipynb', memory='100g')
run('nb_10_pruned_tmrca_stats.ipynb', memory='10g')
run('nb_11_sweeps_from_male_dist_and_pi.ipynb')
run('nb_12_sweep_peaks_ILS_and_rel_tmrca_half.ipynb')
run('nb_13_sweep_regions_and_rel_tmrca_half.ipynb')
run('nb_14_sweep_tmrca.ipynb')
#run('nb_15_sweep_dating.ipynb', cores=15, memory='100g')
run('nb_16_tishkoff_plots.ipynb', cores=10, memory='50g', walltime='3-00:00:00')
run('nb_17_sweeps_and_admixture.ipynb')
run('nb_18_sweeps_and_dist_to_africa.ipynb')
run('nb_19_sweeps_and_dist_to_ust_ishim.ipynb')
run('nb_20_sweep_enrichments.ipynb')
run('nb_21_neutral_scenarios.ipynb')
run('nb_22_slim_simulations.ipynb')
run('nb_23_male_dist_to_archaic.ipynb')
run('nb_24_derived_freqs.ipynb')
run('nb_25_selection_coeficients.ipynb')
run('nb_26_sweeps_1000genomes.ipynb')
run('nb_27_fst_1000genomes.ipynb')
run('nb_28_hapdaf_1000genomes.ipynb')
run('nb_29_circRNA.ipynb')
run('nb_30_scRNA_expression.ipynb')
run('nb_31_maria_iHH.ipynb')










