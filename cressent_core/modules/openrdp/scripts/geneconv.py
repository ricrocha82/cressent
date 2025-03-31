import glob
import logging
import os
import subprocess
import sys
import tempfile
import shutil


class GeneConv:
    def __init__(self, gscale=1, ignore_indels=False, min_length=1, min_poly=2, min_score=2,
                 max_overlap=1, settings=None):
        """
        Constructs a GeneConv object
        :param gscale: mismatch penalty
        :param ignore_indels: Ignore indels or treat indels as one polymorphism (default = False)
        :param min_length: Minimum length of the fragments
        :param min_poly: Minimum number of polymorphic sites
        :param min_score: Minimum pairwise score
        :param max_overlap: Maximum number of overlapping pairs
        """
        if settings is not None:
            self.set_options_from_config(settings)
            self.validate_options()

        else:
            self.gscale = gscale
            self.ignore_indels = ignore_indels
            self.min_length = min_length
            self.min_poly = min_poly
            self.min_score = min_score
            self.max_overlap = max_overlap

        self.raw_results = []
        self.results = []
        self.logger = logging.getLogger('recombination_detection')

    def set_options_from_config(self, settings):
        """
        Set the parameters of GENECONV from the config file
        :param settings: a dictionary of settings
        """
        self.gscale = int(settings['mismatch_penalty'])

        if settings['indels_as_polymorphisms'] == 'True':
            self.ignore_indels = False
        elif settings['indels_as_polymorphisms'] == 'False':
            self.ignore_indels = True
        else:
            self.ignore_indels = None

        self.min_length = int(settings['min_len'])
        self.min_poly = int(settings['min_poly'])
        self.min_score = int(settings['min_score'])
        self.max_overlap = int(settings['max_num'])

    def validate_options(self):
        """
        Check if the options from the config file are valid
        If the options are invalid, the default value will be used instead
        """
        if not isinstance(self.ignore_indels, bool):
            self.logger.warning("Invalid option for 'indels_as_polymorphisms'. Using default value (False) instead.")
            self.ignore_indels = False

        if self.min_length <= 0:
            self.logger.warning("Invalid option for 'min_len'. Using default value (1) instead.")
            self.min_length = 1

        if self.min_poly < 1:
            self.logger.warning("Invalid option for 'min_score'. Using default value (2) instead.")
            self.min_poly = 2

        if self.max_overlap < 0:
            self.logger.warning("Invalid option for 'max_num'. Using default value (1) instead.")
            self.max_overlap = 1

    def execute(self, in_path, output_dir=None):
        """
        Execute the GENECONV algorithm
            S. A. Sawyer (1999)
            GENECONV: A computer package for the statistical detection of gene conversion.
            Distributed by the author, Department of Mathematics, Washington University in St. Louis,
            Available at http://www.math.wustl.edu/~sawyer.
        :param in_path: Path to the input alignment file
        :param output_dir: Directory to save output files (default: directory of input file)
        :return: A list of results
        """
        # Set output directory 
        if output_dir is None:
            output_dir = os.path.dirname(os.path.abspath(in_path))
        else:
            # Ensure output directory exists
            os.makedirs(output_dir, exist_ok=True)
            
        self.logger.info(f"GENECONV output directory: {output_dir}")
        
        # Create temporary directory for GENECONV processing
        temp_dir = tempfile.mkdtemp(dir=output_dir)
        self.logger.info(f"Created temporary directory for GENECONV: {temp_dir}")
        
        try:
            # Copy input file to temp directory
            in_path_abs = os.path.abspath(in_path)
            in_name = os.path.basename(in_path)
            temp_input = os.path.join(temp_dir, in_name)
            shutil.copy2(in_path_abs, temp_input)
            self.logger.info(f"Copied input file to: {temp_input}")
            
            # Clear output files in temp directory
            out_files = glob.glob(os.path.join(temp_dir, '*.frags')) + glob.glob(os.path.join(temp_dir, '*.sum'))
            for f in out_files:
                try:
                    os.remove(f)
                    self.logger.info(f"Removed previous GENECONV output file: {f}")
                except OSError as e:
                    self.logger.warning(f"Could not remove file {f}: {e}")

            # Create config file in temp directory
            config_path = os.path.join(temp_dir, "geneconv.cfg")
            with open(config_path, 'w+') as cfg_handle:
                cfg_handle.write('#GCONV_CONFIG\n')
                cfg_handle.write('  -inputpath={}\n'.format(os.path.realpath(temp_input)))

                if not self.ignore_indels:
                    cfg_handle.write('  -Indel_blocs\n')

                cfg_handle.write('  -Gscale={}\n'.format(self.gscale))
                cfg_handle.write('  -Minlength={}\n'.format(self.min_length))
                cfg_handle.write('  -Minnpoly={}\n'.format(self.min_poly))
                cfg_handle.write('  -Minscore={}\n'.format(self.min_score))
                cfg_handle.write('  -Maxoverlap={}\n'.format(self.max_overlap))
            
            self.logger.info(f"Created GENECONV config file: {config_path}")

            # Path to GENECONV executables
            bin_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'bin')
            if sys.platform.startswith("win"):
                bin_path = os.path.join(bin_dir, 'GENECONV', 'windows_geneconv.exe')
            elif sys.platform == 'darwin':
                bin_path = os.path.join(bin_dir, 'GENECONV', 'geneconv.macOS')
            else:
                bin_path = os.path.join(bin_dir, 'GENECONV', 'geneconv.Unix')

            # Check if custom GENECONV executable exists in PATH
            custom_geneconv = shutil.which('geneconv')
            if custom_geneconv:
                bin_path = custom_geneconv
                self.logger.info(f"Using custom GENECONV executable from PATH: {bin_path}")

            if not os.path.isfile(bin_path):
                self.logger.error(f"No GENECONV executable file exists at {bin_path}")
                return []
                
            self.logger.info(f"Using GENECONV executable: {bin_path}")

            # Run GENECONV from the temp directory
            original_dir = os.getcwd()
            try:
                os.chdir(temp_dir)
                
                # Build command
                cmd = [bin_path, "-Seqfile={}".format(in_name),
                       "-Config={}".format(os.path.basename(config_path)), "-nolog"]
                
                self.logger.info(f"Running GENECONV command: {' '.join(cmd)}")
                
                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE,
                    text=True
                )
                
                stdout, stderr = process.communicate()
                
                # Log output
                if stdout:
                    self.logger.info(f"GENECONV stdout: {stdout}")
                if stderr:
                    self.logger.warning(f"GENECONV stderr: {stderr}")
                
                if process.returncode != 0:
                    self.logger.error(f"GENECONV exited with code {process.returncode}")
                    if "GLIBC" in stderr:
                        self.logger.error("GLIBC compatibility issue detected.")
                    return []
                
                # Get output file name
                output_base = os.path.splitext(in_name)[0]  # Remove the extension (.fa/.fasta)
                output_name = output_base + '.frags'
                
                # Check if output file exists
                if not os.path.exists(output_name):
                    self.logger.error(f"Expected GENECONV output file not found: {output_name}")
                    # List directory contents
                    self.logger.info(f"Files in directory: {os.listdir('.')}")
                    return []
                
                # Parse GENECONV output
                gc_results = self.parse_output(output_name)
                
                # Copy output files to the specified output directory
                for ext in ['.frags', '.sum']:
                    output_file = in_name + ext
                    if os.path.exists(output_file):
                        dest = os.path.join(output_dir, output_file)
                        shutil.copy2(output_file, dest)
                        self.logger.info(f"Copied output file to: {dest}")
                
                return gc_results
                
            finally:
                # Restore original directory
                os.chdir(original_dir)

        except Exception as e:
            self.logger.error(f"Error running GENECONV: {e}", exc_info=True)
            return []
            
        finally:
            # Clean up temporary directory
            try:
                shutil.rmtree(temp_dir)
                self.logger.info(f"Cleaned up temporary directory: {temp_dir}")
            except Exception as e:
                self.logger.warning(f"Could not clean up temporary directory {temp_dir}: {e}")
                
        return []

    def parse_output(self, out_path):
        """
        Parse the output of the GENECONV analysis
        :param out_path: Path to the output file
        :return: List of results
        """
        self.logger.info(f"Parsing GENECONV output: {out_path}")
        self.raw_results = []

        # Check that the out file exists
        try:
            with open(out_path) as out_handle:
                line_count = 0
                for line in out_handle:
                    line_count += 1
                    if not line.startswith('#'):
                        line = line.strip()
                        line = line.split()
                        
                        if len(line) < 6:
                            self.logger.warning(f"Invalid line format in GENECONV output: {line}")
                            continue
                            
                        seqs = line[1].split(';')  # Sequences, potential recombinant is first
                        rec_name = seqs[0]
                        if len(seqs) == 2:
                            parents = [seqs[1], '-']
                        else:
                            parents = ['-', '-']
                        corr_p_value = line[3]  # Bonferroni Corrected - Karlin-Altschul
                        
                        try:
                            locations = (int(line[4]), int(line[5]))  # Locations in alignment
                            self.raw_results.append((rec_name, parents, locations, corr_p_value))
                            self.logger.info(f"Found recombination event: {rec_name}, {parents}, {locations}, {corr_p_value}")
                        except (ValueError, IndexError) as e:
                            self.logger.warning(f"Error parsing locations: {e}, {line}")
                
                self.logger.info(f"Parsed {line_count} lines from GENECONV output")

        except FileNotFoundError as e:
            self.logger.error(f"GENECONV output file not found: {e}")

        return self.raw_results