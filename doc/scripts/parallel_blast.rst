.. _parallel_blast:

.. index:: parallel_blast.py

*parallel_blast.py* -- Parallel BLAST
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script for performing blast while making use of multicore/multiprocessor environments to perform analyses in parallel.


**Usage:** :file:`parallel_blast.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-infile_path
		Path of sequences to use as queries [REQUIRED]
	-r, `-`-refseqs_path
		Path to fasta sequences to search against [REQUIRED]
	-o, `-`-output_dir
		Name of output directory for blast jobs [REQUIRED]
	
	**[OPTIONAL]**
		
	-c, `-`-disable_low_complexity_filter
		Disable filtering of low-complexity sequences (i.e., -F F is passed to blast) [default: False]
	-e, `-`-e_value
		E-value threshold for blasts [default: 1e-30]
	-n, `-`-num_hits
		Number of hits per query for blast results [default: 1]
	-w, `-`-word_size
		Word size for blast searches [default: 30]
	-D, `-`-suppress_format_blastdb
		Supress format of blastdb [default: False]
	-a, `-`-blastmat_dir
		Full path to directory containing blastmat file [default: /Applications/blast-2.2.22/data/]
	-b, `-`-blastall_fp
		Path to blastall [default: blastall]
	-O, `-`-jobs_to_start
		Number of jobs to start [default: 2]
	-P, `-`-poller_fp
		Full path to qiime/parallel/`poller.py <./poller.html>`_ [default: /Users/caporaso/code/Qiime/scripts/`poller.py <./poller.html>`_]
	-R, `-`-retain_temp_files
		Retain temporary files after runs complete (useful for debugging) [default: False]
	-S, `-`-suppress_submit_jobs
		Only split input and write commands file - don't submit jobs [default: False]
	-T, `-`-poll_directly
		Poll directly for job completion rather than running poller as a separate job. If -T is specified this script will not return until all jobs have completed. [default: False]
	-U, `-`-cluster_jobs_fp
		Path to cluster jobs script (defined in qiime_config)  [default: /Users/caporaso/code/Qiime/scripts/`start_parallel_jobs.py <./start_parallel_jobs.html>`_]
	-W, `-`-suppress_polling
		Suppress polling of jobs and merging of results upon completion [default: False]
	-X, `-`-job_prefix
		Job prefix [default: descriptive prefix + random chars]
	-Y, `-`-python_exe_fp
		Full path to python executable [default: python]
	-Z, `-`-seconds_to_sleep
		Number of seconds to sleep between checks for run  completion when polling runs [default: 6]


**Output:**

 


**Example:**

BLAST /home/qiime_user/10_seq.fasta (-i) via three (-O) independent jobs against a blast database created from /home/qiime_user/1000_seq.fasta (-r). Store the results in /home/qiime_user/bla_out/ (-o).

::

	parallel_blast.py -i /home/qiime_user/10_seq.fasta -r /home/qiime_user/1000_seq.fasta -O 3 -o /home/qiime_user/bla_out/


