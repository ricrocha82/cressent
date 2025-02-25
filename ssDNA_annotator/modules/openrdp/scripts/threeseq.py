import glob
import logging
import os
import subprocess
import sys
import shutil
import tempfile


class ThreeSeq:
    def __init__(self, in_path, output_dir=None):
        """
        Initialize ThreeSeq object
        
        Parameters:
        -----------
        in_path : str
            Path to input alignment file
        output_dir : str, optional
            Directory where output files should be saved (default: current directory)
        """
        self.in_path = os.path.abspath(in_path)
        self.in_name = os.path.basename(in_path)
        self.output_dir = os.path.abspath(output_dir) if output_dir else os.getcwd()
        self.raw_results = []
        self.results = []
        self.logger = logging.getLogger('recombination_detection')
        
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)

    def execute(self):
        """
        Execute the 3Seq algorithm.
            Lam HM, Ratmann O, Boni MF.
            Improved algorithmic complexity for the 3SEQ recombination detection algorithm.
            Mol Biol Evol, 35(1):247-251, 2018.
        :return: A list containing the results of the 3Seq analysis
                    Format: [triplets, uncorrected p-value, corrected p-value, breakpoint locations]
        """
        # Use the specified output directory
        self.logger.info(f"3Seq output directory: {self.output_dir}")
        
        # Clear output files in output directory
        out_files = glob.glob(os.path.join(self.output_dir, '*.3s.*'))
        for f in out_files:
            try:
                os.remove(f)
                self.logger.info(f"Removed previous 3Seq output file: {f}")
            except OSError as e:
                self.logger.warning(f"Could not remove file {f}: {e}")
        
        # Create temporary directory within the output directory
        temp_dir = tempfile.mkdtemp(dir=self.output_dir)
        self.logger.info(f"Created temporary directory for 3Seq: {temp_dir}")
        
        try:
            # Copy input file to temp directory
            temp_input = os.path.join(temp_dir, self.in_name)
            shutil.copy2(self.in_path, temp_input)
            self.logger.info(f"Copied input file to: {temp_input}")
            
            # Set paths to 3Seq executables
            bin_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'bin')
            
            if sys.platform.startswith("win"):
                bin_path = os.path.join(bin_dir, '3Seq', 'windows_3seq.exe')
            elif sys.platform == 'darwin':
                bin_path = os.path.join(bin_dir, '3Seq', '3seq.macOS')
            else:
                bin_path = os.path.join(bin_dir, '3Seq', '3seq.Unix')
            
            # Check if custom 3Seq executable exists in PATH
            custom_3seq = shutil.which('3seq')
            if custom_3seq:
                bin_path = custom_3seq
                self.logger.info(f"Using custom 3Seq executable from PATH: {bin_path}")
            
            if not os.path.isfile(bin_path):
                self.logger.error(f"No 3Seq executable found at {bin_path}")
                return []

            self.logger.info(f"Using 3Seq executable: {bin_path}")
            
            # Get output basename without extension
            output_base = os.path.splitext(self.in_name)[0]
            
            # Run 3Seq without p-value table since it's causing issues
            cmd = [bin_path, "-f", temp_input, "-d", "-id", output_base]
            
            self.logger.info(f"Running 3Seq command: {' '.join(cmd)}")
            self.logger.info(f"Working directory: {temp_dir}")
            
            try:
                # Change to temp directory to run 3Seq
                original_dir = os.getcwd()
                os.chdir(temp_dir)
                
                # Run 3Seq process
                process = subprocess.Popen(
                    cmd, 
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE,
                    text=True
                )
                
                # Send 'Y' to any prompts
                stdout, stderr = process.communicate(input="Y\n")
                
                # Log output
                if stdout:
                    self.logger.info(f"3Seq stdout: {stdout}")
                if stderr:
                    self.logger.warning(f"3Seq stderr: {stderr}")
                
                if process.returncode != 0:
                    self.logger.error(f"3Seq exited with code {process.returncode}")
                    if "GLIBC" in stderr:
                        self.logger.error("GLIBC compatibility issue detected. Try compiling 3Seq on this system.")
                    return []
                
                # Look for output files in temp directory
                # First check for the traditional filename
                expected_output = f"{output_base}.3s.rec"
                
                # Then check for the .csv version which appears to be used by newer 3Seq versions
                expected_output_csv = f"{output_base}.3s.rec.csv"
                
                self.logger.info(f"Looking for 3Seq output: {expected_output} or {expected_output_csv}")
                
                # List directory contents
                dir_contents = os.listdir('.')
                self.logger.info(f"Files in directory: {dir_contents}")
                
                # Check for the CSV file format first (newer versions)
                if os.path.exists(expected_output_csv):
                    self.logger.info(f"Found 3Seq output file: {expected_output_csv}")
                    ts_results = self.parse_output_csv(expected_output_csv)
                # If not found, check for the original format
                elif os.path.exists(expected_output):
                    self.logger.info(f"Found 3Seq output file: {expected_output}")
                    ts_results = self.parse_output(expected_output)
                # If neither is found, look for any .3s files
                else:
                    self.logger.warning(f"Expected output files not found")
                    rec_files = glob.glob("*.3s.rec*")
                    if rec_files:
                        found_file = rec_files[0]
                        self.logger.info(f"Found alternative output file: {found_file}")
                        if found_file.endswith('.csv'):
                            ts_results = self.parse_output_csv(found_file)
                        else:
                            ts_results = self.parse_output(found_file)
                    else:
                        self.logger.error("No .3s.rec* files found")
                        return []
                
                # Copy all output files to the specified output directory
                for f in glob.glob(f"*.3s.*"):
                    dest = os.path.join(self.output_dir, f)
                    shutil.copy2(f, dest)
                    self.logger.info(f"Copied output file to: {dest}")
                
                return ts_results
                
            finally:
                # Restore original directory
                os.chdir(original_dir)
                
        except Exception as e:
            self.logger.error(f"Error running 3Seq: {e}", exc_info=True)
            return []
            
        finally:
            # Clean up temporary directory
            try:
                shutil.rmtree(temp_dir)
                self.logger.info(f"Cleaned up temporary directory: {temp_dir}")
            except Exception as e:
                self.logger.warning(f"Could not clean up temporary directory {temp_dir}: {e}")

    def parse_output(self, out_path):
        """
        Parse the output of the 3Seq analysis in original format
        :param out_path: Path to the output file containing information about recombinant sequences
        :return: List of triplets, corrected and uncorrected p-values, and breakpoint locations
        """
        self.logger.info(f"Parsing 3Seq output (original format): {out_path}")
        
        # Check that the out file exists
        if not os.path.exists(out_path):
            self.logger.error(f"3Seq output file not found: {out_path}")
            return []
            
        try:
            with open(out_path) as out_handle:
                # Read header line
                out_handle.readline()  # Read first line

                for line in out_handle:
                    line = line.split('\t')
                    line = [l.strip() for l in line]
                    rec = line[0]
                    ps = [line[1], line[2]]
                    corr_p_value = line[10]  # Dunn-Sidak corrected p-value
                        

                    loc_line = line[12:]    # Breakpoint locations
                    for loc in loc_line:
                        parts = loc.split(' & ')
                        # Take the widest interval 3Seq returns
                        start_pos = parts[0].split('-')
                        end_pos = parts[1].split('-')
                        self.raw_results.append((rec, ps, start_pos[0], end_pos[-1], corr_p_value))
                        
                        self.logger.info(f"Found recombination event: {rec}, {ps}, {start_pos[0]}, {end_pos[-1]}, {corr_p_value}")

        except Exception as e:
            self.logger.error(f"Error parsing 3Seq output: {e}", exc_info=True)
            return []

        self.results = self.merge_breakpoints()
        return self.results

    def parse_output_csv(self, out_path):
        """
        Parse the output of the 3Seq analysis in CSV format (newer versions)
        :param out_path: Path to the CSV output file containing information about recombinant sequences
        :return: List of triplets, corrected and uncorrected p-values, and breakpoint locations
        """
        self.logger.info(f"Parsing 3Seq output (CSV format): {out_path}")
        
        # Check that the out file exists
        if not os.path.exists(out_path):
            self.logger.error(f"3Seq output CSV file not found: {out_path}")
            return []
            
        try:
            with open(out_path) as out_handle:
                # Read header line
                header = out_handle.readline()
                self.logger.info(f"3Seq CSV header: {header.strip()}")

                for line in out_handle:
                    line = line.strip().split(',')
                    
                    rec = line[0]
                    ps = [line[1], line[2]] 
                    
                    # Get p-value - use Dunn-Sidak corrected p-value (column index 10)
                    corr_p_value = line[10]
                    
                    # Get breakpoints (column index 12)
                    loc_line = line[12:] 
                    for loc in loc_line:
                        parts = loc.split(' & ')
                        # Take the widest interval 3Seq returns
                        start_pos = parts[0].split('-')
                        end_pos = parts[1].split('-')
                        self.raw_results.append((rec, ps, start_pos[0], end_pos[-1], corr_p_value))
                    
        except Exception as e:
            self.logger.warning(f"Error parsing breakpoint region: {e}, {region}")
            self.logger.error(f"Error parsing 3Seq CSV output: {e}", exc_info=True)
            return []


        self.results = self.merge_breakpoints()
        return self.results

    def merge_breakpoints(self):
        """
        Merge overlapping breakpoint locations
        :return: list of breakpoint locations where overlapping intervals are merged
        """
        if not self.raw_results:
            self.logger.warning("No raw results to merge")
            return []
            
        results_dict = {}
        results = []

        # Gather all regions with the same recombinant
        for i, bp in enumerate(self.raw_results):
            rec_name = self.raw_results[i][0]
            parents = tuple(sorted(self.raw_results[i][1]))
            key = (rec_name, parents)
            if key not in results_dict:
                results_dict[key] = []
            results_dict[key].append(self.raw_results[i][2:])

        # Merge any locations that overlap - eg [1, 5] and [3, 7] would become [1, 7]
        for key in results_dict:
            merged_regions = []
            for region in results_dict[key]:
                region = list(region)
                old_regions = list(results_dict[key])
                for region2 in old_regions:
                    try:
                        start = int(region[0])
                        end = int(region[1])
                        start2 = int(region2[0])
                        end2 = int(region2[1])
                        
                        if start <= start2 <= end or start <= end2 <= end:
                            region[0] = str(min(start, start2))
                            region[1] = str(max(end, end2))
                            results_dict[key].remove(region2)
                    except (ValueError, IndexError) as e:
                        self.logger.warning(f"Error merging regions: {e}, {region}, {region2}")
                        continue
                        
                merged_regions.append(region)

            # Output the results
            for region in merged_regions:
                rec_name = key[0]
                parents = key[1]
                
                try:
                    start = region[0]
                    end = region[1]
                    p_value = region[2]
                    results.append((rec_name, parents, start, end, p_value))
                    self.logger.info(f"Merged region: {rec_name}, {parents}, {start}, {end}, {p_value}")
                except IndexError as e:
                    self.logger.warning(f"Invalid region format: {e}, {region}")
                    continue

        return results